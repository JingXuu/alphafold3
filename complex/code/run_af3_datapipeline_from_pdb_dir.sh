#!/bin/bash
# 批量：PDB -> AF3 input JSON -> 远程 run_pipeline.sh 跑 MSA+模板 -> 本地缓存 _data.json
#
# 用法：
#   bash complex/code/run_af3_datapipeline_from_pdb_dir.sh \
#       ./complex/PDBs_try \
#       ./data
#
# 其中：
#   PDB_DIR        : 存放两条 target 链的 PDB 文件（每个 PDB = 一对 target）
#   LOCAL_WORK_DIR : 本地工作目录（建议 ./data）

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

# 本地 python（只用于运行 pdb_to_af3_json.py）
PYTHON_BIN="python"

# 本地 JSON 输出目录
LOCAL_JSON_DIR="${LOCAL_WORK_DIR}/af3_input_json"

# 本地 datapipeline 输出根目录
LOCAL_DP_ROOT="${LOCAL_WORK_DIR}/af3_datapipeline_output"

# 本次 run 的时间戳目录名（只作为目录名使用）
TIMESTAMP="af3-dp-$(date +%Y%m%d-%H%M%S-%N)"
LOCAL_DP_OUT="${LOCAL_DP_ROOT}/${TIMESTAMP}"

# 远程 MSA 节点配置（可通过环境变量覆盖）
# 例如：
#   REMOTE_SRV="dbu@inode20" REMOTE_AF3_HOME="/public/dbu/alphafold3" bash ...
REMOTE_SRV="${REMOTE_SRV:-dbu@inode19}"
REMOTE_AF3_HOME="${REMOTE_AF3_HOME:-/public/dbu/alphafold3}"

# 远程本次 run 的根目录
REMOTE_RUN_ROOT="${REMOTE_AF3_HOME}/data/${TIMESTAMP}"
REMOTE_JSON_DIR="${REMOTE_RUN_ROOT}/fasta"
REMOTE_MSA_DIR="${REMOTE_RUN_ROOT}/msa"

echo "[INFO] PDB_DIR         = ${PDB_DIR}"
echo "[INFO] LOCAL_WORK_DIR  = ${LOCAL_WORK_DIR}"
echo "[INFO] LOCAL_JSON_DIR  = ${LOCAL_JSON_DIR}"
echo "[INFO] LOCAL_DP_ROOT   = ${LOCAL_DP_ROOT}"
echo "[INFO] TIMESTAMP       = ${TIMESTAMP}"
echo "[INFO] LOCAL_DP_OUT    = ${LOCAL_DP_OUT}"
echo "[INFO] AF3_ROOT        = ${AF3_ROOT}"
echo "[INFO] REMOTE_SRV      = ${REMOTE_SRV}"
echo "[INFO] REMOTE_AF3_HOME = ${REMOTE_AF3_HOME}"
echo "[INFO] REMOTE_RUN_ROOT = ${REMOTE_RUN_ROOT}"

mkdir -p "${LOCAL_JSON_DIR}"
mkdir -p "${LOCAL_DP_ROOT}"   # 只建根目录，不提前建 LOCAL_DP_OUT（避免失败留下空目录）

########################################
# 第一步：PDB -> AF3 input JSON（两条链）
########################################

echo "[INFO] Step 1: Convert PDBs to AF3 input JSON ..."

shopt -s nullglob
pdb_files=( "${PDB_DIR}"/*.pdb "${PDB_DIR}"/*.PDB )
if [ "${#pdb_files[@]}" -eq 0 ]; then
  echo "[ERROR] No PDB files found in ${PDB_DIR}"
  exit 1
fi

for pdb in "${pdb_files[@]}"; do
  base=$(basename "${pdb}")
  job_name="${base%.*}"
  out_json="${LOCAL_JSON_DIR}/${job_name}.json"

  echo "[INFO]   Processing PDB: ${base} -> ${out_json}"

  # 这里沿用你之前的 pdb_to_af3_json.py 接口：
  #   --pdb <path> --out_json <path>
  # 如需指定链 ID，可以在脚本内部固定 A/C，或者在这里加参数（按你的脚本为准）
  "${PYTHON_BIN}" "${SCRIPT_DIR}/pdb_to_af3_json.py" \
    --pdb "${pdb}" \
    --out_json "${out_json}"

  echo "[INFO]   Wrote AF3 input JSON for ${job_name} -> ${out_json}"
done

########################################
# 第二步：同步 JSON 到远程 & 跑 AF3 data pipeline
########################################

echo "[INFO] Step 2: Sync JSON to remote & run AF3 data pipeline ..."

# 在远程创建本次 run 的目录
ssh "${REMOTE_SRV}" "mkdir -p '${REMOTE_JSON_DIR}' '${REMOTE_MSA_DIR}'"

# 同步本地 JSON 到远程
rsync -avzP "${LOCAL_JSON_DIR}/" "${REMOTE_SRV}:${REMOTE_JSON_DIR}/"

# 在远程简单检查一下是否有 json 文件（仅用于 sanity check / 日志）
REMOTE_JSON_FILES=$(ssh "${REMOTE_SRV}" "ls '${REMOTE_JSON_DIR}'/*.json 2>/dev/null" || true)

if [ -z "${REMOTE_JSON_FILES}" ]; then
  echo "[ERROR] No JSON files found on remote in ${REMOTE_JSON_DIR}"
  exit 1
fi

echo "[INFO] Remote JSON files to process:"
echo "${REMOTE_JSON_FILES}"

echo "[INFO] Start remote data pipeline (MSA + templates) ..."
# 注意这里让 *.json 在远程展开，避免多行变量插入导致的 “Permission denied”
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

# 同步 JSON（可选，用于记录输入）和 MSA 结果
rsync -avzP "${REMOTE_SRV}:${REMOTE_JSON_DIR}/" "${LOCAL_DP_OUT}/fasta/"
rsync -avzP "${REMOTE_SRV}:${REMOTE_MSA_DIR}/"  "${LOCAL_DP_OUT}/msa/"

# 创建一个 latest 软链接，指向最近一次成功的 run
ln -sfn "${LOCAL_DP_OUT}" "${LOCAL_DP_ROOT}/latest"

echo "[INFO] Local datapipeline output root : ${LOCAL_DP_ROOT}"
echo "[INFO] This run output directory      : ${LOCAL_DP_OUT}"
echo "[INFO] Latest run symlink             : ${LOCAL_DP_ROOT}/latest"
echo "[INFO] Example _data.json path        : ${LOCAL_DP_OUT}/msa/<pdbid>/<pdbid>_data.json"
echo "[INFO] Done."

# 用法示例：
#   chmod +x complex/code/run_af3_datapipeline_from_pdb_dir.sh
#   nohup  bash ./complex/code/run_af3_datapipeline_from_pdb_dir.sh ./complex/pdb_no_mg_AC_50_256 ./data > msa-20-256.log 2>&1 &
