from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a normalized domain mapping CSV from SCOPe, CATH, and ECOD."
    )
    parser.add_argument(
        "--classifications-dir",
        type=Path,
        default=Path("data/raw/classifications"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/domain_mapping.csv"),
    )
    parser.add_argument(
        "--min-length", type=int, default=60,
    )
    parser.add_argument(
        "--max-length", type=int, default=180,
    )
    return parser.parse_args()


def parse_scope_cla(path: Path) -> dict[str, dict]:
    domains = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            sid = parts[0]
            pdb_id = parts[1]
            chain_range = parts[2].strip()
            sccs = parts[3]
            sunid = parts[4]
            hierarchy = parts[5]

            chain_id = ""
            res_start = ""
            res_end = ""
            if chain_range and chain_range != "-":
                m = re.match(r"([A-Za-z0-9]):(-?\d+)-(-?\d+)", chain_range)
                if m:
                    chain_id = m.group(1)
                    res_start = m.group(2)
                    res_end = m.group(3)
                else:
                    m2 = re.match(r"([A-Za-z0-9]):", chain_range)
                    if m2:
                        chain_id = m2.group(1)

            domains[sid] = {
                "domain_id": sid,
                "structure_id": pdb_id.lower(),
                "chain_id": chain_id,
                "scope_id": sccs,
                "res_start": res_start,
                "res_end": res_end,
                "hierarchy": hierarchy,
            }
    return domains


def parse_scope_des(path: Path) -> dict[str, dict]:
    entries = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            sunid = parts[0]
            level = parts[1]
            sccs = parts[2]
            domain_or_dash = parts[3]
            desc = parts[4] if len(parts) > 4 else ""
            entries[sunid] = {
                "sunid": sunid,
                "level": level,
                "sccs": sccs,
                "desc": desc,
            }
    return entries


def parse_astral_fasta(path: Path) -> dict[str, int]:
    lengths = {}
    current_id = None
    current_len = 0
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                current_id = line[1:].split()[0].strip()
                current_len = 0
            else:
                current_len += len(line.strip())
    if current_id is not None:
        lengths[current_id] = current_len
    return lengths


def parse_cath_json(path: Path) -> dict[str, list[dict]]:
    with path.open("r", encoding="utf-8") as f:
        raw = json.load(f)
    entries = raw.get("data", raw if isinstance(raw, list) else [])

    by_pdb_chain: dict[str, list[dict]] = {}
    for entry in entries:
        domain_id = entry.get("domain_id", "")
        pdb_code = entry.get("pdb_code", "").lower()
        if not domain_id or not pdb_code:
            continue
        chain_id = domain_id[4] if len(domain_id) > 4 else ""
        sf_id = entry.get("superfamily_id", "")
        atom_length = entry.get("atom_length", "")
        key = f"{pdb_code}_{chain_id}"
        by_pdb_chain.setdefault(key, []).append({
            "cath_domain_id": domain_id,
            "cath_id": sf_id,
            "atom_length": atom_length,
        })
    return by_pdb_chain


def parse_ecod(path: Path) -> dict[str, list[dict]]:
    by_pdb_chain: dict[str, list[dict]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header_cols = None
        for line in f:
            stripped = line.lstrip("#").strip()
            if header_cols is None:
                if stripped.startswith("uid"):
                    header_cols = stripped.split("\t")
                continue
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < len(header_cols):
                continue
            row = dict(zip(header_cols, parts))
            pdb_code = row.get("pdb", "").lower()
            chain_id = row.get("chain", "")
            ecod_id = row.get("ecod_domain_id", "")
            arch = row.get("architecture_name", "")
            x_name = row.get("x_name", "")
            t_name = row.get("t_name", "")
            f_id = row.get("f_id", "")
            pdb_range = row.get("pdb_range", "")

            if not pdb_code:
                continue
            key = f"{pdb_code}_{chain_id}"
            by_pdb_chain.setdefault(key, []).append({
                "ecod_domain_id": ecod_id,
                "ecod_id": f_id,
                "arch": arch,
                "x_name": x_name,
                "t_name": t_name,
                "pdb_range": pdb_range,
            })
    return by_pdb_chain


def extract_kingdom(desc: str) -> str:
    desc_lower = desc.lower()
    bacteria_indicators = [
        "escherichia", "bacillus", "staphylococcus", "streptococcus",
        "salmonella", "mycobacterium", "pseudomonas", "helicobacter",
        "clostridium", "vibrio", "neisseria", "campylobacter",
        "bacteria", "thermus", "thermotoga", "aquifex",
    ]
    archaea_indicators = [
        "archaea", "archaeon", "methanobacterium", "halobacterium",
        "sulfolobus", "pyrococcus", "thermococcus", "methanococcus",
    ]
    eukaryote_indicators = [
        "human", "mouse", "rat", "homo sapiens", "mus musculus",
        "saccharomyces", "drosophila", "caenorhabditis", "arabidopsis",
        "xenopus", "danio", "gallus", "bos", "sus", "oryctolagus",
        "eukaryot", "fungal", "fungi", "plant", "mammal",
    ]
    for kw in archaea_indicators:
        if kw in desc_lower:
            return "archaea"
    for kw in bacteria_indicators:
        if kw in desc_lower:
            return "bacteria"
    for kw in eukaryote_indicators:
        if kw in desc_lower:
            return "eukaryotes"
    return ""


def main() -> int:
    args = parse_args()
    cls_dir = args.classifications_dir

    print("Parsing SCOPe classification ...")
    scope_domains = parse_scope_cla(cls_dir / "scope_cla.txt")
    print(f"  {len(scope_domains)} domains in SCOPe classification")

    print("Parsing SCOPe descriptions ...")
    scope_des = parse_scope_des(cls_dir / "scope_des.txt")
    print(f"  {len(scope_des)} description entries")

    print("Parsing ASTRAL 40% representatives ...")
    astral_lengths = parse_astral_fasta(cls_dir / "scope_astral_40.fa")
    print(f"  {len(astral_lengths)} domain sequences")

    print("Parsing CATH domain data ...")
    cath_by_chain = parse_cath_json(cls_dir / "cath_domain_list.txt")
    cath_total = sum(len(v) for v in cath_by_chain.values())
    print(f"  {cath_total} CATH domains across {len(cath_by_chain)} PDB chains")

    print("Parsing ECOD domain data ...")
    ecod_by_chain = parse_ecod(cls_dir / "ecod_domains.txt")
    ecod_total = sum(len(v) for v in ecod_by_chain.values())
    print(f"  {ecod_total} ECOD domains across {len(ecod_by_chain)} PDB chains")

    sunid_to_species = {}
    for sunid, entry in scope_des.items():
        if entry["level"] == "sp":
            sunid_to_species[sunid] = entry["desc"]

    print(f"\nBuilding domain mapping (length window: {args.min_length}-{args.max_length}) ...")
    rows = []
    filtered_out = 0
    no_length = 0
    for sid, info in scope_domains.items():
        if sid not in astral_lengths:
            no_length += 1
            continue
        seq_len = astral_lengths[sid]
        if seq_len < args.min_length or seq_len > args.max_length:
            filtered_out += 1
            continue

        pdb_id = info["structure_id"]
        chain_id = info["chain_id"]
        chain_key = f"{pdb_id}_{chain_id}"

        cath_entries = cath_by_chain.get(chain_key, [])
        cath_id = cath_entries[0]["cath_id"] if cath_entries else ""

        ecod_entries = ecod_by_chain.get(chain_key, [])
        ecod_id = ecod_entries[0]["ecod_id"] if ecod_entries else ""

        hierarchy = info.get("hierarchy", "")
        kingdom = ""
        sp_match = re.search(r"sp=(\d+)", hierarchy)
        if sp_match:
            sp_sunid = sp_match.group(1)
            species_desc = sunid_to_species.get(sp_sunid, "")
            kingdom = extract_kingdom(species_desc)

        res_start = info["res_start"]
        res_end = info["res_end"]

        rows.append({
            "domain_id": sid,
            "structure_id": pdb_id,
            "chain_id": chain_id,
            "residue_start": res_start,
            "residue_end": res_end,
            "domain_length": str(seq_len),
            "kingdom": kingdom,
            "scope_id": info["scope_id"],
            "cath_id": cath_id,
            "ecod_id": ecod_id,
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "domain_id", "structure_id", "chain_id",
        "residue_start", "residue_end", "domain_length",
        "kingdom", "scope_id", "cath_id", "ecod_id",
    ]
    with args.output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    kingdoms = {}
    for r in rows:
        k = r["kingdom"] or "unknown"
        kingdoms[k] = kingdoms.get(k, 0) + 1
    has_cath = sum(1 for r in rows if r["cath_id"])
    has_ecod = sum(1 for r in rows if r["ecod_id"])

    print(f"\nResults:")
    print(f"  Total ASTRAL 40% domains: {len(astral_lengths)}")
    print(f"  Outside length window: {filtered_out}")
    print(f"  Not in ASTRAL 40%: {no_length}")
    print(f"  Final domain mapping rows: {len(rows)}")
    print(f"  With CATH cross-reference: {has_cath}")
    print(f"  With ECOD cross-reference: {has_ecod}")
    print(f"  Kingdom distribution: {kingdoms}")
    print(f"  Saved to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
