#!/usr/bin/env python
import os
import sys
import json
import csv
import argparse

import numpy as np
from biotite.structure.io import load_structure

# 指标子函数：iptm + 预测结构路径
def compute_iptm_and_pred_path(job_out_dir: str, job_name: str):
    """
    返回 (iptm, model_rank, pred_path)
    - job_out_dir/job_name 下：
      * {job_name}_summary_confidences.json
      * {job_name}_model.cif / .pdb
    """
    job_dir = os.path.join(job_out_dir, job_name)
    summary_path = os.path.join(job_dir, f"{job_name}_summary_confidences.json")

    if not os.path.isfile(summary_path):
        print(f"[WARN metrics] summary_confidences not found: {summary_path}, skip iptm.")
        return None, None, ""

    try:
        with open(summary_path, "r") as f:
            data = json.load(f)
    except Exception as e:
        print(f"[WARN metrics] Failed to load {summary_path}: {e}")
        return None, None, ""

    iptm = data.get("iptm")
    if iptm is None:
        mc = data.get("model_confidence", {})
        iptm = mc.get("iptm")

    if iptm is None:
        print(f"[WARN metrics] iptm not found in {summary_path}, skip iptm.")
        return None, None, ""

    best_iptm = iptm
    best_rank = 0  # 当前每个 job 只有一个模型

    # 预测结构路径
    pred_path = os.path.join(job_dir, f"{job_name}_model.cif")
    if not os.path.isfile(pred_path):
        alt = os.path.join(job_dir, f"{job_name}_model.pdb")
        if os.path.isfile(alt):
            pred_path = alt
        else:
            print(f"[WARN metrics] No predicted model file found for {job_name}, pred_path empty.")
            pred_path = ""

    return best_iptm, best_rank, pred_path


##########################
# 主逻辑
##########################
def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compute iptm for an AF3 job and write/update metrics.csv.\n"
            "mode=iptm-only: compute iptm only, append one row (sc_rmsd empty)"
        )
    )
    parser.add_argument("--job_out_dir", required=True, help="job output dir (contains job_name subdir)")
    parser.add_argument("--job_name", required=True, help="job name, e.g. 1NSG_l100_s313180_mpnn1")
    parser.add_argument("--pdb_id", required=True, help="PDB id (upper case), e.g. 1NSG")
    parser.add_argument("--metrics_csv", required=True, help="global metrics csv path")

    args = parser.parse_args()

    job_out_dir = args.job_out_dir
    job_name = args.job_name
    pdb_id = args.pdb_id
    metrics_csv = args.metrics_csv

    # 确保 metrics_csv 存在（由 bash 预先创建）
    if not os.path.isfile(metrics_csv):
        with open(metrics_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                ["job_name", "pdb_id", "iptm", "model_rank", "pred_path"]
            )
        print(f"[INFO metrics] Created metrics CSV: {metrics_csv}")

    # 计算 iptm
    iptm, model_rank, pred_path = compute_iptm_and_pred_path(job_out_dir, job_name)
    if iptm is None:
        return

    row = [
        job_name,
        pdb_id,
        iptm,
        model_rank if model_rank is not None else "",
        pred_path,
    ]
    try:
        with open(metrics_csv, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(row)
        print(
            f"[metrics iptm-only] job={job_name}, iptm={iptm}, "
            f"model_rank={model_rank}, pred={pred_path}"
        )
    except Exception as e:
        print(f"[WARN metrics] Failed to append to {metrics_csv}: {e}")


if __name__ == "__main__":
    main()
