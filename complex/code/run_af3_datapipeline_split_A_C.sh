#!/bin/bash
# 批量：PDB -> AF3 input JSON(AC) -> 本地拆分 A-only/C-only JSON -> 远程 run_pipeline.sh 跑 MSA+模板 -> 本地缓存 _data.json
#
# 用法：
#    nohup bash ./complex/code/run_af3_datapipeline_split_A_C.sh ./complex/pdb_no_mg_AC_50_256 ./data/split_msas > msa-20-256-split.log 2>&1 &
#
# 产物：
#   LOCAL_JSON_DIR/<TIMESTAMP>/
#     <job>_A.json
#     <job>_C.json
#     (可选) <job>_AC.json  # 中间文件，可配置删除
#
#   LOCAL_DP_OUT/msa/ ...   # 远程 datapipeline 输出同步回来

set -euo pipefail

########################################
# 参数检查
########################################

if [ "$#" -ne 2 ]; then
  echo "Usage: $(basename "$0") <pdb_dir> <local_work_dir>"
  exit 1
fi

PDB_DIR="$1"
LOCAL_WORK_DIR="$2"

# 规范为绝对路径，避免 rsync / ssh 时出错
PDB_DIR="$(readlink -f "${PDB_DIR}")"
LOCAL_WORK_DIR="$(readlink -f "${LOCAL_WORK_DIR}")"

########################################
# 路径与环境变量
########################################

# 当前脚本所在目录（complex/code）
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# AF3 根目录，假设 complex/ 在 AF3 根目录下
AF3_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# 本地 python（只用于运行 pdb_to_af3_json.py 与拆分 JSON）
PYTHON_BIN="python"

# 本次 run 的时间戳目录名（只作为目录名使用）
TIMESTAMP="af3-dp-$(date +%Y%m%d-%H%M%S-%N)"

# 本地 JSON 输出根目录（历史总目录）
LOCAL_JSON_DIR="${LOCAL_WORK_DIR}/af3_input_json"
# 本次 run 的 JSON 子目录（关键修复：避免混入旧 JSON）
LOCAL_JSON_RUN_DIR="${LOCAL_JSON_DIR}/${TIMESTAMP}"

# 本地 datapipeline 输出根目录
LOCAL_DP_ROOT="${LOCAL_WORK_DIR}/af3_datapipeline_output"
LOCAL_DP_OUT="${LOCAL_DP_ROOT}/${TIMESTAMP}"

# 远程 MSA 节点配置（可通过环境变量覆盖）
REMOTE_SRV="${REMOTE_SRV:-dbu@inode19}"
REMOTE_AF3_HOME="${REMOTE_AF3_HOME:-/public/dbu/alphafold3}"

# 远程本次 run 的根目录（本来就是 timestamp 独立目录）
REMOTE_RUN_ROOT="${REMOTE_AF3_HOME}/data/${TIMESTAMP}"
REMOTE_JSON_DIR="${REMOTE_RUN_ROOT}/fasta"
REMOTE_MSA_DIR="${REMOTE_RUN_ROOT}/msa"

# 是否保留中间的 AC.json（默认不保留）
KEEP_AC_JSON="${KEEP_AC_JSON:-0}"

echo "[INFO] PDB_DIR            = ${PDB_DIR}"
echo "[INFO] LOCAL_WORK_DIR     = ${LOCAL_WORK_DIR}"
echo "[INFO] TIMESTAMP          = ${TIMESTAMP}"
echo "[INFO] LOCAL_JSON_DIR     = ${LOCAL_JSON_DIR}"
echo "[INFO] LOCAL_JSON_RUN_DIR = ${LOCAL_JSON_RUN_DIR}"
echo "[INFO] LOCAL_DP_ROOT      = ${LOCAL_DP_ROOT}"
echo "[INFO] LOCAL_DP_OUT       = ${LOCAL_DP_OUT}"
echo "[INFO] AF3_ROOT           = ${AF3_ROOT}"
echo "[INFO] REMOTE_SRV         = ${REMOTE_SRV}"
echo "[INFO] REMOTE_AF3_HOME    = ${REMOTE_AF3_HOME}"
echo "[INFO] REMOTE_RUN_ROOT    = ${REMOTE_RUN_ROOT}"
echo "[INFO] KEEP_AC_JSON       = ${KEEP_AC_JSON}"

# 只建本次 run 的 JSON 目录，确保与本次 PDB_DIR 一致，不混入旧文件
mkdir -p "${LOCAL_JSON_RUN_DIR}"
mkdir -p "${LOCAL_DP_ROOT}"   # 只建根目录，不提前建 LOCAL_DP_OUT（避免失败留下空目录）

########################################
# 小工具：从 AC.json 拆出单链 A/C.json
########################################
split_chain_json() {
  local in_json="$1"
  local chain_id="$2"
  local job_name="$3"
  local out_json="$4"

  "${PYTHON_BIN}" - "${in_json}" "${chain_id}" "${job_name}" "${out_json}" <<'PY'
import json, sys

in_json   = sys.argv[1]
chain_id  = sys.argv[2]
job_name  = sys.argv[3]
out_json  = sys.argv[4]

with open(in_json, "r") as f:
    obj = json.load(f)

seqs = obj.get("sequences", [])
keep = []
for s in seqs:
    prot = s.get("protein", None)
    if not prot:
        continue
    if prot.get("id", None) == chain_id:
        keep.append(s)

if len(keep) != 1:
    raise RuntimeError(
        f"[split_chain_json] Expected exactly 1 chain with id={chain_id} in {in_json}, got {len(keep)}. "
        f"Found ids={[ (x.get('protein',{}) or {}).get('id',None) for x in seqs ]}"
    )

obj["name"] = f"{job_name}_{chain_id}"
obj["sequences"] = keep

with open(out_json, "w") as f:
    json.dump(obj, f, ensure_ascii=False, indent=2)

print(f"[split_chain_json] Wrote {out_json} (chain {chain_id})")
PY
}

########################################
# 第一步：PDB -> AF3 input JSON（先生成 AC，再拆分 A/C）
########################################

echo "[INFO] Step 1: Convert PDBs to AF3 input JSON (AC) and split into A-only/C-only ..."

shopt -s nullglob
pdb_files=( "${PDB_DIR}"/*.pdb "${PDB_DIR}"/*.PDB )
if [ "${#pdb_files[@]}" -eq 0 ]; then
  echo "[ERROR] No PDB files found in ${PDB_DIR}"
  exit 1
fi

for pdb in "${pdb_files[@]}"; do
  base=$(basename "${pdb}")
  job_name="${base%.*}"

  # 先生成 AC 合并版 json（中间产物）—— 注意：写到本次 run 的子目录
  out_json_ac="${LOCAL_JSON_RUN_DIR}/${job_name}_AC.json"
  out_json_a="${LOCAL_JSON_RUN_DIR}/${job_name}_A.json"
  out_json_c="${LOCAL_JSON_RUN_DIR}/${job_name}_C.json"

  echo "[INFO]   Processing PDB: ${base}"
  echo "[INFO]     1) AC json -> ${out_json_ac}"

  "${PYTHON_BIN}" "${SCRIPT_DIR}/pdb_to_af3_json.py" \
    --pdb "${pdb}" \
    --out_json "${out_json_ac}"

  echo "[INFO]     2) Split to A-only -> ${out_json_a}"
  split_chain_json "${out_json_ac}" "A" "${job_name}" "${out_json_a}"

  echo "[INFO]     3) Split to C-only -> ${out_json_c}"
  split_chain_json "${out_json_ac}" "C" "${job_name}" "${out_json_c}"

  if [ "${KEEP_AC_JSON}" -eq 0 ]; then
    rm -f "${out_json_ac}"
  else
    echo "[INFO]     (KEEP_AC_JSON=1) Keeping intermediate ${out_json_ac}"
  fi

  echo "[INFO]   Wrote: ${out_json_a} and ${out_json_c}"
done

########################################
# 第二步：同步 JSON 到远程 & 跑 AF3 data pipeline
########################################

echo "[INFO] Step 2: Sync JSON to remote & run AF3 data pipeline ..."

# 在远程创建本次 run 的目录
ssh "${REMOTE_SRV}" "mkdir -p '${REMOTE_JSON_DIR}' '${REMOTE_MSA_DIR}'"

# （可选保险）确保远程目录干净，避免你手动复用 TIMESTAMP 时混入旧文件
ssh "${REMOTE_SRV}" "rm -f '${REMOTE_JSON_DIR}'/*.json 2>/dev/null || true"

# 只同步本次 run 的 JSON 子目录（关键修复点）
rsync -avzP "${LOCAL_JSON_RUN_DIR}/" "${REMOTE_SRV}:${REMOTE_JSON_DIR}/"

# 在远程简单检查一下是否有 json 文件（仅用于 sanity check / 日志）
REMOTE_JSON_FILES=$(ssh "${REMOTE_SRV}" "ls '${REMOTE_JSON_DIR}'/*.json 2>/dev/null" || true)

if [ -z "${REMOTE_JSON_FILES}" ]; then
  echo "[ERROR] No JSON files found on remote in ${REMOTE_JSON_DIR}"
  exit 1
fi

echo "[INFO] Remote JSON files to process:"
echo "${REMOTE_JSON_FILES}"

echo "[INFO] Start remote data pipeline (MSA + templates) ..."
ssh "${REMOTE_SRV}" "cd '${REMOTE_AF3_HOME}' && \
  bash -x ./run_pipeline.sh \
    -o '${REMOTE_MSA_DIR}' \
    ${REMOTE_JSON_DIR}/*.json"
echo "[INFO] Remote data pipeline done."

########################################
# 第三步：同步 MSA+模板结果回本地，缓存 _data.json
########################################

echo "[INFO] Step 3: Sync MSA/template results back to local ..."

# 只有当远程 datapipeline 完成后，才创建本地 TIMESTAMP 目录，避免空目录
mkdir -p "${LOCAL_DP_OUT}"

# 同步 JSON（用于记录本次实际输入）和 MSA 结果
rsync -avzP "${REMOTE_SRV}:${REMOTE_JSON_DIR}/" "${LOCAL_DP_OUT}/fasta/"
rsync -avzP "${REMOTE_SRV}:${REMOTE_MSA_DIR}/"  "${LOCAL_DP_OUT}/msa/"

# 创建一个 latest 软链接，指向最近一次成功的 run
ln -sfn "${LOCAL_DP_OUT}" "${LOCAL_DP_ROOT}/latest"

echo "[INFO] Local datapipeline output root : ${LOCAL_DP_ROOT}"
echo "[INFO] This run output directory      : ${LOCAL_DP_OUT}"
echo "[INFO] Latest run symlink             : ${LOCAL_DP_ROOT}/latest"

echo "[INFO] Example _data.json paths (A/C are separate now):"
echo "[INFO]   ${LOCAL_DP_OUT}/msa/<job>_A/<job>_A_data.json"
echo "[INFO]   ${LOCAL_DP_OUT}/msa/<job>_C/<job>_C_data.json"

echo "[INFO] Done."
