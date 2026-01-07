#!/usr/bin/env python
import json
import argparse


def read_fasta_one(path: str) -> str:
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
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


def main():
    parser = argparse.ArgumentParser(
        description="Merge A/C datapipeline json with binder FASTA to build 3-chain data.json"
    )
    parser.add_argument("--a_json", required=True, help="A-chain datapipeline json")
    parser.add_argument("--c_json", required=True, help="C-chain datapipeline json")
    parser.add_argument("--fasta", required=True, help="binder FASTA file (B chain)")
    parser.add_argument("--job_name", required=True, help="job name used as data.json name")
    parser.add_argument("--out_json", required=True, help="output 3-chain data.json path")
    parser.add_argument("--binder_chain", default="B", help="binder chain id, default B")

    args = parser.parse_args()

    A = load_json(args.a_json)
    C = load_json(args.c_json)

    A_entry = get_protein_entry(A, "A")
    C_entry = get_protein_entry(C, "C")

    binder_seq = read_fasta_one(args.fasta)

    B_entry = {
        "protein": {
            "id": args.binder_chain,
            "sequence": binder_seq,
            "unpairedMsa": "",
            "pairedMsa": "",
            "templates": [],
        }
    }

    out = dict(A)
    out["name"] = args.job_name
    out["sequences"] = [A_entry, B_entry, C_entry]

    for k, v in C.items():
        if k not in out:
            out[k] = v

    with open(args.out_json, "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    print(f"[merge_ac_add_binder] wrote {args.out_json}")


if __name__ == "__main__":
    main()
