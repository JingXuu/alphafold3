#!/usr/bin/env python
import os
import sys
import json
import csv
import argparse
from typing import List, Tuple

import numpy as np
from biotite.structure.io import load_structure


##### copy parted from kabsch_algorithm ######


def kabsch_algorithm(P, Q):
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
        print("Warning: SVD did not converge. Returning identity rotation.")
        R = np.eye(3)
    return R, C_P, C_Q


def calculate_rmsd(P, Q):
    diff = P - Q
    return np.sqrt(np.sum(diff * diff) / len(P))


def get_ca_coords_from_path(pdb_path: str, chain_list: List[str]):
    atom_array = load_structure(pdb_path)
    ca_atoms = atom_array[
        np.isin(atom_array.chain_id, chain_list) & (atom_array.atom_name == "CA")
    ]
    assert len(ca_atoms) > 0, f"No CA atoms found in chain {chain_list} of {pdb_path}"
    return ca_atoms.coord


def complex_ca_rmsd_from(
    path_1: str,
    path_2: str,
    chain_1_list: List[str] = ["A", "B"],
    chain_2_list: List[str] = ["A", "B"],
):
    """
    Input:
        path_1: str, path to the first pdb file
        path_2: str, path to the second pdb file
        chain_1_list: list of str, chain ids in the first pdb file
        chain_2_list: list of str, chain ids in the second pdb file
    Output:
        rmsd: float, CA RMSD between the two complexes
    """

    assert set(chain_1_list) == set(chain_2_list)

    coords1 = get_ca_coords_from_path(path_1, chain_1_list)
    coords2 = get_ca_coords_from_path(path_2, chain_2_list)

    assert len(coords1) == len(coords2)

    R, C_P, C_Q = kabsch_algorithm(coords1, coords2)
    coords2_aligned = np.dot(coords2 - C_Q, R) + C_P
    rmsd = calculate_rmsd(coords1, coords2_aligned)

    return round(rmsd, 3)


##########################
# 指标子函数：iptm + 预测结构路径
##########################
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
# 指标子函数：scRMSD
##########################
def compute_scrmsd(pred_path: str, target_dir: str, binder_root: str) -> (str, str):
    """
    根据 pred_path 和 target_dir/Trajectory/binder_root.pdb 计算 scRMSD。
    返回 (sc_rmsd_str, ref_path)
    """
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

    if not pred_path or not os.path.isfile(pred_path) or ref_path is None:
        if not pred_path:
            print("[WARN metrics] pred_path missing, cannot compute scRMSD.")
        elif not os.path.isfile(pred_path):
            print(f"[WARN metrics] pred_path not a file: {pred_path}, cannot compute scRMSD.")
        if ref_path is None:
            print(
                f"[WARN metrics] ref_path missing for binder_root={binder_root} under {target_dir}, "
                "cannot compute scRMSD."
            )
        return "", ref_path or ""

    # 这里按 A, C, B 三链计算；内部的 complex_ca_rmsd_from 逻辑完全照你给的版本
    chains = ["A", "C", "B"]
    try:
        scrmsd_val = complex_ca_rmsd_from(
            ref_path,
            pred_path,
            chain_1_list=chains,
            chain_2_list=chains,
        )
        return str(scrmsd_val), ref_path
    except AssertionError as e:
        # 这里不改你原函数逻辑，只在外面打个 warning
        print(f"[WARN metrics] Assertion failed in complex_ca_rmsd_from: {e}")
        return "", ref_path
    except Exception as e:
        print(f"[WARN metrics] Failed to compute scRMSD: {e}")
        return "", ref_path


##########################
# 主逻辑
##########################
def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compute iptm + scRMSD for an AF3 job and write/update metrics.csv.\n"
            "mode=all:       compute iptm+scRMSD and append one row\n"
            "mode=iptm-only: compute iptm only, append one row (sc_rmsd empty)\n"
            "mode=scrmsd-only: recompute scRMSD and update existing row"
        )
    )
    parser.add_argument("--job_out_dir", required=True, help="job output dir (contains job_name subdir)")
    parser.add_argument("--job_name", required=True, help="job name, e.g. 1NSG_l100_s313180_mpnn1")
    parser.add_argument("--pdb_id", required=True, help="PDB id (upper case), e.g. 1NSG")
    parser.add_argument("--target_dir", required=True, help="fold root target dir (contains Trajectory/)")
    parser.add_argument("--metrics_csv", required=True, help="global metrics csv path")
    parser.add_argument("--binder_root", required=True, help="binder root, e.g. 1NSG_l100_s313180")
    parser.add_argument(
        "--mode",
        default="all",
        choices=["all", "iptm-only", "scrmsd-only"],
        help="metrics mode (default: all)",
    )

    args = parser.parse_args()

    job_out_dir = args.job_out_dir
    job_name = args.job_name
    pdb_id = args.pdb_id
    target_dir = args.target_dir
    metrics_csv = args.metrics_csv
    binder_root = args.binder_root
    mode = args.mode

    # 确保 metrics_csv 存在（all / iptm-only 模式通常由 bash 预先创建）
    if not os.path.isfile(metrics_csv):
        with open(metrics_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                ["job_name", "pdb_id", "binder_root", "iptm", "sc_rmsd", "model_rank", "pred_path", "ref_path"]
            )
        print(f"[INFO metrics] Created metrics CSV: {metrics_csv}")

    # ---------------------
    # mode = iptm-only / all : 需要 (iptm, rank, pred_path)
    # ---------------------
    if mode in ("all", "iptm-only"):
        iptm, model_rank, pred_path = compute_iptm_and_pred_path(job_out_dir, job_name)
        if iptm is None:
            return

        if mode == "iptm-only":
            row = [
                job_name,
                pdb_id,
                binder_root,
                iptm,
                "",           # sc_rmsd
                model_rank if model_rank is not None else "",
                pred_path,
                "",           # ref_path
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
            return

        # mode == "all": 继续算 scRMSD
        sc_rmsd, ref_path = compute_scrmsd(pred_path, target_dir, binder_root)
        row = [
            job_name,
            pdb_id,
            binder_root,
            iptm,
            sc_rmsd,
            model_rank if model_rank is not None else "",
            pred_path,
            ref_path,
        ]
        try:
            with open(metrics_csv, "a", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(row)
            print(
                f"[metrics all] job={job_name}, iptm={iptm}, sc_rmsd={sc_rmsd}, "
                f"model_rank={model_rank}, pred={pred_path}, ref={ref_path}"
            )
        except Exception as e:
            print(f"[WARN metrics] Failed to append to {metrics_csv}: {e}")
        return

    # ---------------------
    # mode = scrmsd-only : 只重算 scRMSD，更新已有行
    # ---------------------
    if mode == "scrmsd-only":
        with open(metrics_csv, "r", newline="") as f:
            reader = csv.reader(f)
            rows = list(reader)

        if not rows:
            print(f"[WARN metrics] {metrics_csv} is empty, cannot update scRMSD.")
            return

        header = rows[0]
        data_rows = rows[1:]

        job_idx = None
        for i, row in enumerate(data_rows, start=1):
            if row and row[0] == job_name:
                job_idx = i
                break

        if job_idx is None:
            print(f"[WARN metrics] No row for job_name={job_name} in {metrics_csv}, skip scrmsd-only.")
            return

        row = rows[job_idx]
        pred_path = row[6]
        if not pred_path:
            job_dir = os.path.join(job_out_dir, job_name)
            cand = os.path.join(job_dir, f"{job_name}_model.cif")
            if os.path.isfile(cand):
                pred_path = cand
            else:
                alt = os.path.join(job_dir, f"{job_name}_model.pdb")
                if os.path.isfile(alt):
                    pred_path = alt

        sc_rmsd, ref_path = compute_scrmsd(pred_path, target_dir, binder_root)

        if sc_rmsd:
            row[4] = sc_rmsd
        if pred_path:
            row[6] = pred_path
        if ref_path:
            row[7] = ref_path

        rows[job_idx] = row

        try:
            with open(metrics_csv, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerows(rows)
            print(
                f"[metrics scrmsd-only] job={job_name}, sc_rmsd={row[4]}, "
                f"pred={row[6]}, ref={row[7]}"
            )
        except Exception as e:
            print(f"[WARN metrics] Failed to rewrite {metrics_csv}: {e}")


if __name__ == "__main__":
    main()
