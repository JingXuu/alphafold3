#!/bin/bash
# 批量：fold_try 结构 + 预先缓存的 A_data.json / C_data.json + binder FASTA
#      -> 为每个 binder MPNN 序列构造 (A,B,C) 三链 *_data.json，调用 AF3 fold（不跑MSA）
#
# 用法：
#   bash complex/code/run_af3_fold_add_binder.sh \
#       <fold_root_dir> \
#       <af3_dp_root> \
#       <af3_model_dir> \
#       <af3_output_root>
#
# 示例：
# CUDA_VISIBLE_DEVICES=1 bash complex/code/run_af3_fold_add_binder_split.sh \
#     ./complex/fold_try \
#     ./data/split_msas/af3_datapipeline_output/af3-dp-20251221-021152-272270354 \
#     ~/public_databases/models \
#     ./data/af3_fold_results_try
#
# 说明：
#   - datapipeline 已被拆分为 A-only / C-only：各自有独立 *_data.json（含 MSA + 模板）
#   - 每个 binder FASTA 在自己的 job 目录中生成一份三链 *_data.json（插入 B 链，MSA-free）
#   - AF3 用 job 目录作为 input_dir + --norun_data_pipeline 做推理
#   - 输出目录结构：
#       OUT_ROOT/PDB_ID/binder_root/job_name/
#
#   - OUT_ROOT/folded_jobs.txt 记录已完成的 job_name，便于断点续跑

set -euo pipefail

if [ "$#" -ne 4 ]; then
  echo "Usage: $(basename "$0") <fold_root_dir> <af3_dp_root> <model_dir> <output_root>"
  exit 1
fi

FOLD_ROOT="$1"
AF3_DP_ROOT="$2"
MODEL_DIR="$3"
OUTPUT_ROOT="$4"

# 规范为绝对路径
FOLD_ROOT="$(readlink -f "${FOLD_ROOT}")"
AF3_DP_ROOT="$(readlink -f "${AF3_DP_ROOT}")"
MODEL_DIR="$(readlink -f "${MODEL_DIR}")"
OUTPUT_ROOT="$(readlink -f "${OUTPUT_ROOT}")"

# AF3 根目录（假设此脚本放在 complex/code 下）
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
AF3_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

echo "[INFO] FOLD_ROOT   = ${FOLD_ROOT}"
echo "[INFO] AF3_DP_ROOT = ${AF3_DP_ROOT}"
echo "[INFO] MODEL_DIR   = ${MODEL_DIR}"
echo "[INFO] OUTPUT_ROOT = ${OUTPUT_ROOT}"
echo "[INFO] AF3_ROOT    = ${AF3_ROOT}"

mkdir -p "${OUTPUT_ROOT}"

# 维护一个全局 folded_list.txt，记录已经跑过的 job_name（比如 1NSG_l83_s113551_mpnn5）
FOLDED_LIST="${OUTPUT_ROOT}/folded_jobs.txt"
touch "${FOLDED_LIST}"
echo "[INFO] FOLDED_LIST = ${FOLDED_LIST}"

# 切换到 AF3 根目录，便于直接调用 run_alphafold.py
cd "${AF3_ROOT}"

########################################
# Helper: 在 split datapipeline 输出里定位 A/C 的 *_data.json
########################################
find_cached_data_json() {
  # $1: pdb_id_upper (e.g., 1A7X)
  # $2: pdb_id_lower (e.g., 1a7x)
  # $3: chain_id     (A or C)  -- 注意：我们内部会同时尝试 A/a
  local pdbU="$1"
  local pdbL="$2"
  local chain="$3"

  local chainU="${chain}"
  local chainL
  chainL="$(echo "${chain}" | tr '[:upper:]' '[:lower:]')"  # A->a, C->c

  # 你现在的真实输出形如：msa/1a7x_a/1a7x_a_data.json
  local cand=(
    # 1) 最精确：小写 pdb + 小写链
    "${AF3_DP_ROOT}/msa/${pdbL}_${chainL}/${pdbL}_${chainL}_data.json"
    # 2) 小写 pdb + 大写链（兼容）
    "${AF3_DP_ROOT}/msa/${pdbL}_${chainU}/${pdbL}_${chainU}_data.json"
    # 3) 大写 pdb + 小写链（兼容）
    "${AF3_DP_ROOT}/msa/${pdbU}_${chainL}/${pdbU}_${chainL}_data.json"
    # 4) 大写 pdb + 大写链（兼容）
    "${AF3_DP_ROOT}/msa/${pdbU}_${chainU}/${pdbU}_${chainU}_data.json"

    # 5) 目录存在但文件名不完全一致时的通配兜底
    "${AF3_DP_ROOT}/msa/${pdbL}_${chainL}"/*_data.json
    "${AF3_DP_ROOT}/msa/${pdbL}_${chainU}"/*_data.json
    "${AF3_DP_ROOT}/msa/${pdbU}_${chainL}"/*_data.json
    "${AF3_DP_ROOT}/msa/${pdbU}_${chainU}"/*_data.json

    # 6) 更宽松兜底：目录名包含 <pdb>_<chain> 即可（应对 job_name 不止 4 位）
    "${AF3_DP_ROOT}/msa/"*"${pdbL}_${chainL}"*/*_data.json
    "${AF3_DP_ROOT}/msa/"*"${pdbL}_${chainU}"*/*_data.json
    "${AF3_DP_ROOT}/msa/"*"${pdbU}_${chainL}"*/*_data.json
    "${AF3_DP_ROOT}/msa/"*"${pdbU}_${chainU}"*/*_data.json
  )

  shopt -s nullglob
  for p in "${cand[@]}"; do
    for f in $p; do
      if [ -f "$f" ]; then
        echo "$f"
        shopt -u nullglob
        return 0
      fi
    done
  done
  shopt -u nullglob
  return 1
}


########################################
# Helper: 由 A_data_json + C_data_json + binder_fasta 合成三链 data json
########################################
merge_ac_add_binder() {
  # $1 A_data_json
  # $2 C_data_json
  # $3 binder_fasta
  # $4 job_name
  # $5 out_json
  # $6 binder_chain_id (default B)
  local a_json="$1"
  local c_json="$2"
  local fasta="$3"
  local job_name="$4"
  local out_json="$5"
  local binder_chain="${6:-B}"

  python - "${a_json}" "${c_json}" "${fasta}" "${job_name}" "${out_json}" "${binder_chain}" <<'PY'
import json, sys, re

a_json, c_json, fasta, job_name, out_json, binder_chain = sys.argv[1:]

def read_fasta_one(path: str) -> str:
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq.append(line)
    s = "".join(seq).replace(" ", "").upper()
    if not s or re.search(r"[^ACDEFGHIKLMNPQRSTVWYBXZJUO]", s):
        # 这里不做过度严格限制；若包含非常规字符，AF3 可能会报错
        pass
    return s

def load_json(path: str):
    with open(path, "r") as f:
        return json.load(f)

def get_protein_entry(obj: dict, chain_id: str) -> dict:
    seqs = obj.get("sequences", [])
    for s in seqs:
        prot = s.get("protein")
        if prot and prot.get("id") == chain_id:
            return s
    raise RuntimeError(f"Cannot find protein chain id={chain_id} in {obj.get('name','<no-name>')}")

A = load_json(a_json)
C = load_json(c_json)

A_entry = get_protein_entry(A, "A")
C_entry = get_protein_entry(C, "C")

binder_seq = read_fasta_one(fasta)

B_entry = {
    "protein": {
        "id": binder_chain,
        "sequence": binder_seq,
        # 明确 MSA-free / template-free
        "unpairedMsa": "",
        "pairedMsa": "",
        "templates": []
    }
}

# 以 A 的 data.json 为基底，尽量保留 datapipeline 生成的其它顶层字段（如果存在）
out = dict(A)

# name 改为 job_name（不依赖文件名）
out["name"] = job_name

# 关键：合并成三链 sequences（顺序 A,B,C）
out["sequences"] = [A_entry, B_entry, C_entry]

# 如果 C 的 data.json 有一些顶层键 A 没有，但你希望保留，可以按需合并（默认不强行覆盖 A）
# 这里做一个温和合并：缺失键补上
for k, v in C.items():
    if k not in out:
        out[k] = v

with open(out_json, "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

print(f"[merge_ac_add_binder] wrote {out_json}")
PY
}

########################################
# 遍历 fold_root 下的 target 目录
########################################

shopt -s nullglob

for target_dir in "${FOLD_ROOT}"/*; do
  [ -d "${target_dir}" ] || continue

  target_base="$(basename "${target_dir}")"
  # 例如 target_base = 1NSG_time_check_af3_1
  pdb_id="${target_base:0:4}"
  pdb_id_upper="$(echo "${pdb_id}" | tr '[:lower:]' '[:upper:]')"  # 1NSG
  pdb_id_lower="$(echo "${pdb_id}" | tr '[:upper:]' '[:lower:]')"  # 1nsg

  echo
  echo "==============================="
  echo "[INFO] Target dir : ${target_base}"
  echo "[INFO] PDB_ID     : ${pdb_id_upper}"
  echo "==============================="

  # 定位 A/C 的 cached data json（split 后独立）
  cached_data_json_A="$(find_cached_data_json "${pdb_id_upper}" "${pdb_id_lower}" "A" || true)"
  cached_data_json_C="$(find_cached_data_json "${pdb_id_upper}" "${pdb_id_lower}" "C" || true)"

  if [ -z "${cached_data_json_A}" ] || [ ! -f "${cached_data_json_A}" ]; then
    echo "[WARN] Cached A data json not found for ${pdb_id_upper}. Skip this target."
    echo "[WARN]   Tried under: ${AF3_DP_ROOT}/msa/*${pdb_id_upper}_A*/*_data.json"
    continue
  fi
  if [ -z "${cached_data_json_C}" ] || [ ! -f "${cached_data_json_C}" ]; then
    echo "[WARN] Cached C data json not found for ${pdb_id_upper}. Skip this target."
    echo "[WARN]   Tried under: ${AF3_DP_ROOT}/msa/*${pdb_id_upper}_C*/*_data.json"
    continue
  fi

  echo "[INFO] Cached A data json: ${cached_data_json_A}"
  echo "[INFO] Cached C data json: ${cached_data_json_C}"

  # 该 target 对应的 MPNN 序列目录
  mpnn_seq_dir="${target_dir}/MPNN/Sequences"
  if [ ! -d "${mpnn_seq_dir}" ]; then
    echo "[WARN] No MPNN/Sequences dir under ${target_dir}, skip."
    continue
  fi

  # 为该 target 单独建一个输出子目录
  target_out_root="${OUTPUT_ROOT}/${pdb_id_upper}"
  mkdir -p "${target_out_root}"

  # 遍历所有 binder FASTA
  for fasta in "${mpnn_seq_dir}"/*.fasta; do
    [ -f "${fasta}" ] || continue

    fasta_base="$(basename "${fasta}")"           # 例如 1NSG_l83_s113551_mpnn5.fasta
    job_name="${fasta_base%.fasta}"              # 1NSG_l83_s113551_mpnn5
    binder_root="${job_name%_mpnn*}"             # 1NSG_l83_s113551

    # 利用 folded_list.txt 去重
    if grep -qx "${job_name}" "${FOLDED_LIST}"; then
      echo "[INFO]   Skip already folded job: ${job_name}"
      continue
    fi

    echo "[INFO]   New job: ${job_name}"
    echo "[INFO]     FASTA       : ${fasta}"
    echo "[INFO]     binder_root : ${binder_root}"

    # 输出目录：OUT_ROOT/1NSG/1NSG_l83_s113551/1NSG_l83_s113551_mpnn5/
    binder_out_root="${target_out_root}/${binder_root}"
    job_out_dir="${binder_out_root}/${job_name}"
    mkdir -p "${job_out_dir}"

    ########################################
    # 1) 为该 job 构造专用的三链 *_data.json
    #    从 cached A/C data json 拷贝信息 + 插入 B（MSA-free）
    ########################################
    job_data_json="${job_out_dir}/${pdb_id_lower}_data.json"

    merge_ac_add_binder \
      "${cached_data_json_A}" \
      "${cached_data_json_C}" \
      "${fasta}" \
      "${job_name}" \
      "${job_data_json}" \
      "B"

    echo "[INFO]     Patched 3-chain data json: ${job_data_json}"

    ########################################
    # 2) 调用 AF3 ：input_dir 指向 job_out_dir，norun_data_pipeline
    ########################################
    echo "[INFO]     Calling AF3 fold for ${job_name}"

    python run_alphafold.py \
      --model_dir "${MODEL_DIR}" \
      --input_dir "${job_out_dir}" \
      --output_dir "${job_out_dir}" \
      --norun_data_pipeline

    # 若 run 成功，则把 job_name 记入 folded_list
    echo "${job_name}" >> "${FOLDED_LIST}"
    echo "[INFO]     Finished job: ${job_name}"
  done

done

echo
echo "[INFO] All done."
echo "[INFO] Folded jobs are listed in: ${FOLDED_LIST}"
echo "[INFO] Outputs are under: ${OUTPUT_ROOT}"
