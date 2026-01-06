#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
从一个 PDB 文件中提取两条链的序列，并生成 AlphaFold3 的 input JSON（仅包含两个 protein 链）。

用法示例：
  python pdb_to_af3_json.py \
    --pdb /path/to/targets/7S2S_7S2R.pdb \
    --out_json /path/to/af3_input_json/7S2S_7S2R.json

可选指定链 ID（默认取文件中前两个链）：
  python pdb_to_af3_json.py \
    --pdb targets/7S2S_7S2R.pdb \
    --chains A,C \
    --out_json af3_input_json/7S2S_7S2R.json
"""

import argparse
import json
import os
from collections import defaultdict

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    # 常见特殊氨基酸
    "SEC": "U",  # Selenocysteine
    "PYL": "O",  # Pyrrolysine
}


def parse_pdb_two_chain_seqs(pdb_path, chain_ids=None):
    """
    从 PDB 中解析出每条链的氨基酸序列（只使用 ATOM 记录）。
    返回: {chain_id: sequence}
    """
    residues = {}  # key: (chain, resseq, icode) -> resname
    chains_in_file = []

    with open(pdb_path, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if len(line) < 54:
                continue
            chain_id = line[21].strip()  # PDB 格式 chain ID 在 column 22
            if not chain_id:
                continue
            if chain_id not in chains_in_file:
                chains_in_file.append(chain_id)

            resseq = line[22:26].strip()
            icode = line[26].strip()
            resname = line[17:20].strip()

            key = (chain_id, resseq, icode)
            if key not in residues:
                residues[key] = resname

    if not chains_in_file:
        raise ValueError(f"No chains found in PDB: {pdb_path}")

    if chain_ids is None:
        # 默认取前两个链
        if len(chains_in_file) < 2:
            raise ValueError(
                f"PDB {pdb_path} has only {len(chains_in_file)} chain(s), "
                f"need at least 2."
            )
        chain_ids = chains_in_file[:2]
    else:
        chain_ids = [c.strip() for c in chain_ids.split(",") if c.strip()]
        for c in chain_ids:
            if c not in chains_in_file:
                raise ValueError(
                    f"Specified chain '{c}' not found in PDB {pdb_path}. "
                    f"Available chains: {chains_in_file}"
                )
        if len(chain_ids) != 2:
            raise ValueError(
                f"Expect exactly 2 chains, got {chain_ids}. "
                f"Please pass like: --chains A,C"
            )

    seqs = {}
    # 分链组装序列
    residues_by_chain = defaultdict(list)
    for (ch, resseq, icode), resname in residues.items():
        residues_by_chain[ch].append((int(resseq), icode, resname))

    for ch in chain_ids:
        if ch not in residues_by_chain:
            raise ValueError(f"Chain '{ch}' has no residues in {pdb_path}")

        # 按 (resseq, icode) 排序
        res_list = sorted(residues_by_chain[ch], key=lambda x: (x[0], x[1]))
        one_letter_seq = []
        for _, _, resname in res_list:
            resname = resname.upper()
            aa = THREE_TO_ONE.get(resname, "X")
            one_letter_seq.append(aa)

        seqs[ch] = "".join(one_letter_seq)

    return seqs, chain_ids


def make_af3_input_json(pdb_path, out_json, chain_ids=None, seed=1):
    """
    为给定 PDB 构造 AF3 的 input JSON（仅两条 protein 链）。
    """
    basename = os.path.basename(pdb_path)
    job_name = os.path.splitext(basename)[0]

    seqs, used_chain_ids = parse_pdb_two_chain_seqs(pdb_path, chain_ids)

    sequences_block = []
    for ch in used_chain_ids:
        sequences_block.append(
            {
                "protein": {
                    "id": ch,
                    "sequence": seqs[ch],
                    # 注意：这里只提供序列。MSA/模板会在 data pipeline 阶段自动补充。
                }
            }
        )

    data = {
        "name": job_name,
        "modelSeeds": [seed],
        "sequences": sequences_block,
        "dialect": "alphafold3",
        # 对于 AF3 v3.0.0+，version 建议 >= 1；如使用外部 MSA/模板时需 >=2
        "version": 1,
    }

    os.makedirs(os.path.dirname(out_json), exist_ok=True)
    with open(out_json, "w") as f:
        json.dump(data, f, indent=2)

    print(f"[INFO] Wrote AF3 input JSON for {job_name} -> {out_json}")
    print(f"[INFO] Chains used: {used_chain_ids[0]}, {used_chain_ids[1]}")
    for ch in used_chain_ids:
        print(f"  Chain {ch}: length {len(seqs[ch])}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract two-chain sequences from PDB and write AF3 input JSON."
    )
    parser.add_argument("--pdb", required=True, help="Input PDB path.")
    parser.add_argument("--out_json", required=True, help="Output JSON path.")
    parser.add_argument(
        "--chains",
        default=None,
        help="Comma-separated chain IDs, e.g., 'A,C'. "
             "If not set, will use the first two chains in the PDB.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Model seed for AF3 (default: 1).",
    )

    args = parser.parse_args()
    make_af3_input_json(args.pdb, args.out_json, chain_ids=args.chains, seed=args.seed)


if __name__ == "__main__":
    main()
