"""A6: extract AlphaMissense and ESM1b scores for the manifest proteins.

One streaming pass over AlphaMissense_aa_substitutions.tsv.gz (1.2 GB gzipped)
pulls every needed accession at once, so the file is never scanned twice.
ESM1b comes from per-accession CSVs inside ESM1b_ALL_hum_isoforms_LLR.zip by
random access.

Scope discipline. Confirmation-protein scores are written to a SEPARATE
directory, `features/confirmation_held/`, which no development or fitting
script reads. Extracting a predictor score is not looking at an outcome and is
not a model-versus-AM comparison, but keeping the two sets physically apart
makes an accidental peek take a deliberate act.

Orientation (locked decision 1): AM pathogenicity is already damage-oriented
(higher = more pathogenic). ESM1b LLR is fitness-oriented (alternative minus
wild type), so the damage score written here is NEGATED: damage = -LLR.

AlphaMissense data is (c) 2023 DeepMind Technologies Limited. The file header
declares CC BY-NC-SA 4.0 (note: the project prompt says CC BY 4.0; the file
itself is the more restrictive NC-SA). Not validated for, and not approved for,
any clinical use.

Creates new files only.
"""

from __future__ import annotations

import csv
import gzip
import json
import time
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INV = ROOT / "results" / "beat_am_20261016" / "inventory"
FEAT = ROOT / "results" / "beat_am_20261016" / "features"

AM_GZ = ROOT / "data" / "raw" / "AlphaMissense_aa_substitutions.tsv.gz"
ESM_ZIP = ROOT / "data" / "raw" / "ESM1b_ALL_hum_isoforms_LLR.zip"
ESM_PREFIX = "content/ALL_hum_isoforms_ESM1b_LLR/"

CONFIRM_ROLES = {"confirmation_membrane", "confirmation_allhuman"}


def load_manifest() -> list[dict]:
    with (INV / "group_manifest.csv").open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def extract_am(accessions: dict[str, Path]) -> dict[str, int]:
    """Single streaming pass. accessions maps accession -> output directory."""
    handles: dict[str, object] = {}
    counts: dict[str, int] = {a: 0 for a in accessions}
    t0 = time.time()
    lines = 0
    try:
        with gzip.open(AM_GZ, "rt") as fh:
            for line in fh:
                lines += 1
                if lines % 20_000_000 == 0:
                    print(f"  {lines/1e6:.0f}M lines, {time.time()-t0:.0f}s", flush=True)
                if line.startswith("#") or line.startswith("uniprot_id"):
                    continue
                tab = line.find("\t")
                acc = line[:tab]
                if acc not in counts:
                    continue
                h = handles.get(acc)
                if h is None:
                    out = accessions[acc] / f"{acc}.am.tsv"
                    h = out.open("w", encoding="utf-8")
                    h.write("protein_variant\tam_pathogenicity\n")
                    handles[acc] = h
                rest = line[tab + 1:].rstrip("\n").split("\t")
                h.write(f"{rest[0]}\t{rest[1]}\n")
                counts[acc] += 1
    finally:
        for h in handles.values():
            h.close()
    print(f"  AM pass done: {lines/1e6:.1f}M lines in {time.time()-t0:.0f}s")
    return counts


def extract_esm(accessions: dict[str, Path]) -> dict[str, int]:
    """Per-accession CSV out of the zip. Writes damage = -LLR."""
    counts: dict[str, int] = {}
    with zipfile.ZipFile(ESM_ZIP) as z:
        available = set(z.namelist())
        for acc, outdir in accessions.items():
            name = f"{ESM_PREFIX}{acc}_LLR.csv"
            if name not in available:
                counts[acc] = 0
                continue
            with z.open(name) as fh:
                rows = fh.read().decode("utf-8").splitlines()
            if not rows:
                counts[acc] = 0
                continue
            header = rows[0].split(",")
            # header[0] is empty; the rest are "<native_aa> <position>"
            sites = []
            for cell in header[1:]:
                parts = cell.strip().split()
                sites.append((parts[0], int(parts[1])) if len(parts) == 2 else (None, None))
            n = 0
            out = outdir / f"{acc}.esm1b.tsv"
            with out.open("w", encoding="utf-8") as oh:
                oh.write("protein_variant\tesm1b_damage\n")
                for row in rows[1:]:
                    cells = row.split(",")
                    alt = cells[0].strip()
                    if not alt:
                        continue
                    for (native, pos), value in zip(sites, cells[1:]):
                        if native is None or alt == native:
                            continue
                        try:
                            llr = float(value)
                        except ValueError:
                            continue
                        # fitness-oriented LLR -> damage score
                        oh.write(f"{native}{pos}{alt}\t{-llr:.4f}\n")
                        n += 1
            counts[acc] = n
    return counts


def main() -> None:
    manifest = load_manifest()
    fit_dir = FEAT / "fitting_pool"
    hold_dir = FEAT / "confirmation_held"
    fit_dir.mkdir(parents=True, exist_ok=True)
    hold_dir.mkdir(parents=True, exist_ok=True)

    routing = {
        r["accession"]: (hold_dir if r["role"] in CONFIRM_ROLES else fit_dir)
        for r in manifest
    }
    print(f"accessions to extract: {len(routing)} "
          f"({sum(1 for v in routing.values() if v == hold_dir)} held, "
          f"{sum(1 for v in routing.values() if v == fit_dir)} fitting)")

    print("ESM1b...")
    esm_counts = extract_esm(routing)
    print(f"  ESM1b: {sum(1 for v in esm_counts.values() if v)} accessions found, "
          f"{sum(esm_counts.values()):,} variant rows")

    print("AlphaMissense (single streaming pass)...")
    am_counts = extract_am(routing)
    print(f"  AM: {sum(1 for v in am_counts.values() if v)} accessions found, "
          f"{sum(am_counts.values()):,} variant rows")

    by_acc = {r["accession"]: r for r in manifest}
    rows = []
    for acc in sorted(routing):
        r = by_acc[acc]
        rows.append(
            {
                "accession": acc,
                "gene": r["gene"],
                "role": r["role"],
                "fitting_eligible": r["fitting_eligible"],
                "held_separately": routing[acc] == hold_dir,
                "am_rows": am_counts.get(acc, 0),
                "esm1b_rows": esm_counts.get(acc, 0),
                "am_found": am_counts.get(acc, 0) > 0,
                "esm1b_found": esm_counts.get(acc, 0) > 0,
            }
        )
    with (FEAT / "predictor_coverage.csv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    missing_am = [r["gene"] for r in rows if not r["am_found"]]
    missing_esm = [r["gene"] for r in rows if not r["esm1b_found"]]
    (FEAT / "extraction_manifest.json").write_text(
        json.dumps(
            {
                "orientation": {
                    "am_pathogenicity": "already damage-oriented, used as-is",
                    "esm1b": "source LLR is fitness-oriented (alt - wt); written "
                             "here as damage = -LLR per locked decision 1",
                },
                "scope_discipline": "confirmation-protein scores live in "
                                    "features/confirmation_held/ and are read by no "
                                    "fitting or development script",
                "attribution": "AlphaMissense data (c) 2023 DeepMind Technologies "
                               "Limited. The distributed file header declares "
                               "CC BY-NC-SA 4.0 (the project prompt states CC BY 4.0; "
                               "the file itself is the more restrictive NC-SA licence). "
                               "Not validated for, and not approved for, any clinical use.",
                "n_accessions": len(rows),
                "am_accessions_found": sum(1 for r in rows if r["am_found"]),
                "esm1b_accessions_found": sum(1 for r in rows if r["esm1b_found"]),
                "missing_am": missing_am,
                "missing_esm1b": missing_esm,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    print(f"\nmissing AM ({len(missing_am)}): {missing_am}")
    print(f"missing ESM1b ({len(missing_esm)}): {missing_esm}")


if __name__ == "__main__":
    main()
