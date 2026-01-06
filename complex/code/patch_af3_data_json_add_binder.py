#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
在 AF3 datapipeline 输出的 *_data.json（含 A/C 两条链）基础上，
插入一条 binder 链 B，并修改 name 为 job_name。

假设 cached *_data.json 顶层大致包含：
{
  "name": "...",
  "sequences": [
    { "protein": { "id": "A", "sequence": "...", ... } },
    { "protein": { "id": "C", "sequence": "...", ... } },
    ...  # 可能还有其他条目
  ],
  ...
}

本脚本：
- 保持 A/C 及其它字段完整；
- 在 sequences 中补一条 B 链：
  {
    "protein": {
      "id": "B",
      "sequence": "<binder_seq>",
      "modifications": [],
      "unpairedMsa": "",
      "pairedMsa": "",
      "templates": []
    }
  }
- 顶层 "name" 改为 job_name。
"""

import argparse
import json
from typing import Any, Dict, List


def parse_args():
    p = argparse.ArgumentParser(
        description="在 AF3 *_data.json 中插入 binder 链 B（复用 A/C MSA & 模板）"
    )
    p.add_argument(
        "--cached_data_json",
        required=True,
        help="datapipeline 输出的 *_data.json 路径，如 msa/1nsg/1nsg_data.json",
    )
    p.add_argument(
        "--binder_fasta",
        required=True,
        help="binder 的 FASTA 文件路径（单序列）",
    )
    p.add_argument(
        "--job_name",
        required=True,
        help="AF3 job 名称，例如 1NSG_l83_s113551_mpnn5",
    )
    p.add_argument(
        "--binder_chain_id",
        default="B",
        help="binder 链 ID，默认 B",
    )
    p.add_argument(
        "--out_json",
        required=True,
        help="输出 JSON 路径（写到 job 目录，不覆盖原始 *_data.json）",
    )
    return p.parse_args()


def read_fasta_first_seq(path: str) -> str:
    """读取 FASTA 文件中的首条序列，返回大写氨基酸序列。"""
    seq_lines: List[str] = []
    in_seq = False
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if in_seq:
                    break
                in_seq = True
                continue
            if in_seq:
                seq_lines.append(line)
    if not seq_lines:
        raise ValueError(f"No sequence found in FASTA: {path}")
    return "".join(seq_lines).upper()


def get_protein_id(protein_obj: Dict[str, Any]) -> str:
    """
    兼容两种写法：
    - "id": "A"
    - "id": ["A"]
    """
    pid = protein_obj.get("id")
    if isinstance(pid, list):
        if not pid:
            raise ValueError("protein.id is an empty list")
        return str(pid[0])
    return str(pid)


def main():
    args = parse_args()

    # 1) 读取 cached *_data.json（只读）
    with open(args.cached_data_json, "r") as f:
        data = json.load(f)

    if "sequences" not in data or not isinstance(data["sequences"], list):
        raise ValueError(
            f"{args.cached_data_json} 中缺少顶层 'sequences' 列表，请先查看实际结构再调整脚本。"
        )

    sequences: List[Dict[str, Any]] = data["sequences"]

    # 2) 拿出 A/C 对应的 entry，保留其他 entry 备用
    a_entry = None
    c_entry = None
    other_entries: List[Dict[str, Any]] = []

    for entry in sequences:
        protein = entry.get("protein")
        if not isinstance(protein, dict):
            other_entries.append(entry)
            continue

        pid = get_protein_id(protein)
        if pid == "A":
            a_entry = entry
        elif pid == "C":
            c_entry = entry
        else:
            other_entries.append(entry)

    if a_entry is None or c_entry is None:
        raise ValueError(
            f"{args.cached_data_json} 中未找到 protein.id 为 'A' 或 'C' 的条目，"
            "请检查文件内容。"
        )

    # 3) 读取 binder 序列
    binder_seq = read_fasta_first_seq(args.binder_fasta)

    # 4) 构造 B 链的 protein 条目
    binder_protein = {
        "id": args.binder_chain_id,  # 默认 "B"
        "sequence": binder_seq,
        "modifications": [],
        "unpairedMsa": "",
        "pairedMsa": "",
        "templates": [],
    }
    binder_entry = {"protein": binder_protein}

    # 5) 新 sequences：A, B, C（如需保留其他条目，可附加 other_entries）
    new_sequences: List[Dict[str, Any]] = [a_entry, binder_entry, c_entry]
    # 若你确认还需要保留其他链（如配体），可以取消下面一行注释：
    # new_sequences.extend(other_entries)

    data["sequences"] = new_sequences
    data["name"] = args.job_name

    # 6) 写出新的 JSON（写到 job 目录，不覆盖原始 cached_data_json）
    with open(args.out_json, "w") as f:
        json.dump(data, f, indent=2)

    print(
        f"[INFO] Patched {args.cached_data_json} -> {args.out_json} "
        f"(job_name={args.job_name}, binder_chain_id={args.binder_chain_id})"
    )


if __name__ == "__main__":
    main()
