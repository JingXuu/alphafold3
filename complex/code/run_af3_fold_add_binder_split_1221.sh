#!/bin/bash
# cal RMSD has problems.


# 批量：fold_try 结构 + 预先缓存的 A_data.json / C_data.json + binder FASTA
#      -> 为每个 binder MPNN 序列构造 (A,B,C) 三链 *_data.json，调用 AF3 fold（不跑MSA）
#
# 用法：
#   bash complex/code/run_af3_fold_add_binder_split.sh \
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
#     ./data/af3_fold_results
#
# CUDA_VISIBLE_DEVICES=1 nohup bash complex/code/run_af3_fold_add_binder_split.sh \
#     ./data/test \
#     ./data/split_msas/af3_datapipeline_output/af3-dp-20251221-021152-272270354 \
#     ~/public_databases/models \
#     ./data/test_results \
#     > ./data/logs/af3_$(date +%Y%m%d-%H%M%S).log 2>&1 &
#
# 说明：
#   - datapipeline 已被拆分为 A-only / C-only：各自有独立 *_data.json（含 MSA + 模板）
#   - 每个 binder FASTA 在自己的 job 目录中生成一份三链 *_data.json（插入 B 链，MSA-free）
#   - AF3 用 job 目录作为 input_dir + --norun_data_pipeline 做推理
#   - 输出目录结构：
#       OUT_ROOT/PDB_ID/binder_root/job_name/
#         job_name_data.json        # 我们拼的三链数据
#         job_name/                 # AF3 真正的输出目录
#           job_name_summary_confidences.json
#           job_name_model.cif
#           ...
#   - OUT_ROOT/metrics.csv 记录每个 job 的 iptm / scRMSD 等指标
#   - 每次 fold 前，只根据 metrics.csv 中是否已有该 job_name 的记录来决定是否跳过。

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

# 全局指标文件：首行表头
METRICS_CSV="${OUTPUT_ROOT}/metrics.csv"
if [ ! -f "${METRICS_CSV}" ]; then
  echo "job_name,pdb_id,binder_root,iptm,sc_rmsd,model_rank,pred_path,ref_path" > "${METRICS_CSV}"
  echo "[INFO] Created metrics CSV: ${METRICS_CSV}"
else
  echo "[INFO] Reusing existing metrics CSV: ${METRICS_CSV}"
fi

# 切换到 AF3 根目录，便于直接调用 run_alphafold.py
cd "${AF3_ROOT}"

########################################
# Helper: 在 split datapipeline 输出里定位 A/C 的 *_data.json
########################################
find_cached_data_json() {
  # $1: pdb_id_upper (e.g., 1A7X)
  # $2: pdb_id_lower (e.g., 1a7x)
  # $3: chain_id     (A or C)
  local pdbU="$1"
  local pdbL="$2"
  local chain="$3"

  local chainU="${chain}"
  local chainL
  chainL="$(echo "${chain}" | tr '[:upper:]' '[:lower:]')"  # A->a, C->c

  local cand=(
    "${AF3_DP_ROOT}/msa/${pdbL}_${chainL}/${pdbL}_${chainL}_data.json"
    "${AF3_DP_ROOT}/msa/${pdbL}_${chainU}/${pdbL}_${chainU}_data.json"
    "${AF3_DP_ROOT}/msa/${pdbU}_${chainL}/${pdbU}_${chainL}_data.json"
    "${AF3_DP_ROOT}/msa/${pdbU}_${chainU}/${pdbU}_${chainU}_data.json"

    "${AF3_DP_ROOT}/msa/${pdbL}_${chainL}"/*_data.json
    "${AF3_DP_ROOT}/msa/${pdbL}_${chainU}"/*_data.json
    "${AF3_DP_ROOT}/msa/${pdbU}_${chainL}"/*_data.json
    "${AF3_DP_ROOT}/msa/${pdbU}_${chainU}"/*_data.json

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
    raise RuntimeError(
        f"Cannot find protein chain id={chain_id} in {obj.get('name','<no-name>')}"
    )

A = load_json(a_json)
C = load_json(c_json)

A_entry = get_protein_entry(A, "A")
C_entry = get_protein_entry(C, "C")

binder_seq = read_fasta_one(fasta)

B_entry = {
    "protein": {
        "id": binder_chain,
        "sequence": binder_seq,
        "unpairedMsa": "",
        "pairedMsa": "",
        "templates": []
    }
}

out = dict(A)
out["name"] = job_name
out["sequences"] = [A_entry, B_entry, C_entry]

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
  pdb_id="${target_base:0:4}"
  pdb_id_upper="$(echo "${pdb_id}" | tr '[:lower:]' '[:upper:]')"  # 1NSG
  pdb_id_lower="$(echo "${pdb_id}" | tr '[:upper:]' '[:lower:]')"  # 1nsg

  echo
  echo "==============================="
  echo "[INFO] Target dir : ${target_base}"
  echo "[INFO] PDB_ID     : ${pdb_id_upper}"
  echo "==============================="

  cached_data_json_A="$(find_cached_data_json "${pdb_id_upper}" "${pdb_id_lower}" "A" || true)"
  cached_data_json_C="$(find_cached_data_json "${pdb_id_upper}" "${pdb_id_lower}" "C" || true)"

  if [ -z "${cached_data_json_A}" ] || [ ! -f "${cached_data_json_A}" ]; then
    echo "[WARN] Cached A data json not found for ${pdb_id_upper}. Skip this target."
    continue
  fi
  if [ -z "${cached_data_json_C}" ] || [ ! -f "${cached_data_json_C}" ]; then
    echo "[WARN] Cached C data json not found for ${pdb_id_upper}. Skip this target."
    continue
  fi

  echo "[INFO] Cached A data json: ${cached_data_json_A}"
  echo "[INFO] Cached C data json: ${cached_data_json_C}"

  mpnn_seq_dir="${target_dir}/MPNN/Sequences"
  if [ ! -d "${mpnn_seq_dir}" ]; then
    echo "[WARN] No MPNN/Sequences dir under ${target_dir}, skip."
    continue
  fi

  target_out_root="${OUTPUT_ROOT}/${pdb_id_upper}"
  mkdir -p "${target_out_root}"

  for fasta in "${mpnn_seq_dir}"/*.fasta; do
    [ -f "${fasta}" ] || continue

    fasta_base="$(basename "${fasta}")"           # 1NSG_l100_s313180_mpnn1.fasta
    job_name="${fasta_base%.fasta}"              # 1NSG_l100_s313180_mpnn1
    binder_root="${job_name%_mpnn*}"             # 1NSG_l100_s313180

    # 只根据 metrics.csv 去重
    if grep -q "^${job_name}," "${METRICS_CSV}"; then
      echo "[INFO]   Skip already evaluated job (metrics.csv): ${job_name}"
      continue
    fi

    echo "[INFO]   New job: ${job_name}"
    echo "[INFO]     FASTA       : ${fasta}"
    echo "[INFO]     binder_root : ${binder_root}"

    traj_dir="${target_dir}/Trajectory"
    if [ ! -d "${traj_dir}" ]; then
      echo "[WARN]   No Trajectory dir under ${target_dir}, skip job ${job_name}."
      continue
    fi

    ref_pdb_for_job=""
    # 文件名精确匹配 binder_root.pdb，例如 1NSG_l44_s890019.pdb
    ref_pdb_for_job=$(find "${traj_dir}" -maxdepth 1 -type f -name "${binder_root}.pdb" | head -n 1 || true)

    if [ -z "${ref_pdb_for_job}" ]; then
      echo "[WARN]   No reference PDB ${binder_root}.pdb under ${traj_dir}, skip job ${job_name}."
      continue
    fi
    echo "[INFO]     Found reference PDB for this job: ${ref_pdb_for_job}"


    binder_out_root="${target_out_root}/${binder_root}"
    job_out_dir="${binder_out_root}/${job_name}"
    mkdir -p "${job_out_dir}"

    job_data_json="${job_out_dir}/${pdb_id_lower}_data.json"

    merge_ac_add_binder \
      "${cached_data_json_A}" \
      "${cached_data_json_C}" \
      "${fasta}" \
      "${job_name}" \
      "${job_data_json}" \
      "B"

    echo "[INFO]     Patched 3-chain data json: ${job_data_json}"

    echo "[INFO]     Calling AF3 fold for ${job_name}"

    python run_alphafold.py \
      --model_dir "${MODEL_DIR}" \
      --input_dir "${job_out_dir}" \
      --output_dir "${job_out_dir}" \
      --norun_data_pipeline

    ########################################
    # 3) 成功后计算 iptm + scRMSD，并写入 metrics.csv
    ########################################
    python - "${job_out_dir}" "${job_name}" "${pdb_id_upper}" "${target_dir}" "${METRICS_CSV}" "${binder_root}" <<'PY'
import os
import sys
import json
import csv
from typing import List, Tuple

import numpy as np
from biotite.structure.io import load_structure


job_out_dir, job_name, pdb_id, target_dir, metrics_csv, binder_root = sys.argv[1:]


##########################
# Kabsch + CA RMSD 工具
##########################
def kabsch_algorithm(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    C_P = np.mean(P, axis=0)
    C_Q = np.mean(Q, axis=0)
    P_centered = P - C_P
    Q_centered = Q - C_Q
    H = np.dot(P_centered.T, Q_centered)

    try:
        U, S, Vt = np.linalg.svd(H)
        R = np.dot(Vt.T, U.T)
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = np.dot(Vt.T, U.T)
    except np.linalg.LinAlgError:
        print("[WARN metrics] SVD did not converge. Returning identity rotation.")
        R = np.eye(3)
    return R, C_P, C_Q


def calculate_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    diff = P - Q
    return float(np.sqrt(np.sum(diff * diff) / len(P)))


def get_ca_coords_from_path(pdb_path: str, chain_list: List[str]) -> np.ndarray:
    atom_array = load_structure(pdb_path)
    mask = np.isin(atom_array.chain_id, chain_list) & (atom_array.atom_name == "CA")
    ca_atoms = atom_array[mask]
    if ca_atoms.array_length() == 0:
        raise RuntimeError(f"No CA atoms found in chains {chain_list} of {pdb_path}")
    return ca_atoms.coord


def complex_ca_rmsd_from(
    path_1: str,
    path_2: str,
    chain_1_list: List[str],
    chain_2_list: List[str],
) -> float:
    if set(chain_1_list) != set(chain_2_list):
        raise ValueError(f"chain lists must match: {chain_1_list} vs {chain_2_list}")

    coords1 = get_ca_coords_from_path(path_1, chain_1_list)
    coords2 = get_ca_coords_from_path(path_2, chain_2_list)

    if coords1.shape != coords2.shape:
        raise ValueError(
            f"CA coord shape mismatch: {coords1.shape} vs {coords2.shape} "
            f"({path_1} vs {path_2})"
        )

    R, C_P, C_Q = kabsch_algorithm(coords1, coords2)
    coords2_aligned = np.dot(coords2 - C_Q, R) + C_P
    rmsd = calculate_rmsd(coords1, coords2_aligned)

    return round(rmsd, 3)


##########################
# 1) summary_confidences 路径（AF3 输出在子目录 job_name/ 下）
##########################
job_dir = os.path.join(job_out_dir, job_name)
summary_path = os.path.join(job_dir, f"{job_name}_summary_confidences.json")

if not os.path.isfile(summary_path):
    print(f"[WARN metrics] summary_confidences not found: {summary_path}, skip metrics.")
    sys.exit(0)

try:
    with open(summary_path, "r") as f:
        data = json.load(f)
except Exception as e:
    print(f"[WARN metrics] Failed to load {summary_path}: {e}")
    sys.exit(0)

iptm = data.get("iptm")
if iptm is None:
    mc = data.get("model_confidence", {})
    iptm = mc.get("iptm")

if iptm is None:
    print(f"[WARN metrics] iptm not found in {summary_path}, skip metrics.")
    sys.exit(0)

best_iptm = iptm
best_rank = 0  # 这里只有一个模型，rank 固定为 0

##########################
# 2) 预测结构路径（job_name_model.cif / .pdb）
##########################
pred_path = os.path.join(job_dir, f"{job_name}_model.cif")
if not os.path.isfile(pred_path):
    alt = os.path.join(job_dir, f"{job_name}_model.pdb")
    if os.path.isfile(alt):
        pred_path = alt
    else:
        print(f"[WARN metrics] No predicted model file found for {job_name}, skip scRMSD.")
        pred_path = None

##########################
# 3) 找到参考结构：
#    从 target_dir/Trajectory/ 下按 binder_root 匹配 .pdb 文件
##########################
ref_path = None

traj_dir = os.path.join(target_dir, "Trajectory")
ref_path = None

if os.path.isdir(traj_dir):
    cand = os.path.join(traj_dir, f"{binder_root}.pdb")
    if os.path.isfile(cand):
        ref_path = cand
    else:
        print(
            f"[WARN metrics] Expected ref PDB {binder_root}.pdb not found under {traj_dir}, "
            "skip scRMSD."
        )
else:
    print(f"[WARN metrics] Trajectory dir not found under {target_dir}, skip scRMSD.")

sc_rmsd = ""

if pred_path is not None and ref_path is not None:
    # A, C, B 三条链整体计算 CA RMSD
    chains = ["A", "C", "B"]
    try:
        sc_rmsd = complex_ca_rmsd_from(
            ref_path,
            pred_path,
            chain_1_list=chains,
            chain_2_list=chains,
        )
    except Exception as e:
        print(f"[WARN metrics] Failed to compute scRMSD for {job_name}: {e}")
        sc_rmsd = ""
else:
    if pred_path is None:
        print(f"[WARN metrics] pred_path missing for {job_name}, only iptm recorded.")
    if ref_path is None:
        print(
            f"[WARN metrics] No reference PDB (.pdb) found under {traj_dir} for binder_root={binder_root}, "
            "only iptm recorded."
        )

##########################
# 4) 将指标写入 metrics.csv
##########################
row = [
    job_name,
    pdb_id,
    binder_root,
    best_iptm,
    sc_rmsd,
    best_rank,
    pred_path or "",
    ref_path or "",
]

try:
    with open(metrics_csv, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)
    print(
        f"[metrics] job={job_name}, iptm={best_iptm}, sc_rmsd={sc_rmsd}, "
        f"model_rank={best_rank}, pred={pred_path}, ref={ref_path}"
    )
except Exception as e:
    print(f"[WARN metrics] Failed to append to {metrics_csv}: {e}")

PY

    echo "[INFO]     Finished job: ${job_name}"
  done

done

echo
echo "[INFO] All done."
echo "[INFO] Metrics CSV: ${METRICS_CSV}"
echo "[INFO] Outputs are under: ${OUTPUT_ROOT}"
