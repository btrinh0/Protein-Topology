"""Controls and scale fixes for the within-site charge audit.

Five things the headline result needs before it can be frozen.

**Standardised units.** AlphaMissense scores run 0 to 1; ESM1b is an unbounded
log-likelihood ratio. "+1.11 versus +0.051" compares nothing. Every side
difference here is also divided by the standard deviation of that model's own
within-site penalties, and reported as a rank statistic P(cytosolic > outer),
which assumes no scale at all.

**Glycosylation control.** On the outer side, S/T to K/R destroys an N-X-S/T
sequon, and N to K/R destroys it too. Losing glycosylation is damaging for
reasons that have nothing to do with the positive-inside rule, and it only
happens on the outer side, so it biases the gain contrast toward the hypothesis.
Sites participating in a sequon are excluded here.

**Lysine versus arginine.** Lysine is a ubiquitination target, so losing one can
look benign and gaining one damaging for degradation reasons. That biases against
the hypothesis. Arginine-only is the cleaner test.

**Floor effects.** AlphaMissense is bounded, so a side difference could be
compression near an endpoint. Repeated on the logit scale and on ranks.

**Ion channels.** The census carries no pore or re-entrant loop annotation, so
pore loops cannot be excluded directly. Excluding ion channels entirely is the
available proxy, and it matters because Kir2.1 has a pore loop.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from run_gate0 import load_frame
from run_gate0_within_site import estimate, site_gain_penalties, site_penalties
from run_gate1 import read_topology

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TOPOLOGY = ROOT / "data" / "raw" / "topology_viewer.html"
DEFAULT_INTRAMEM = ROOT / "data" / "external" / "uniprot" / "human_intramem.tsv"
PORE_MARGIN = 10  # a pore span plus the flank window either side of it
PRIMARY_WINDOW = 10
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def sequon_positions(sequence: str) -> set[int]:
    """1-based positions taking part in an N-X-S/T sequon (X is anything but proline)."""
    positions: set[int] = set()
    for i in range(len(sequence) - 2):
        if sequence[i] == "N" and sequence[i + 1] != "P" and sequence[i + 2] in "ST":
            positions.add(i + 1)
            positions.add(i + 3)
    return positions


def load_pore_spans(path: Path) -> dict[str, list[tuple[int, int]]]:
    """UniProt INTRAMEM features: pore helices, selectivity filters, re-entrant loops.

    The census annotates only transmembrane segments and loops, so a pore loop looks
    like ordinary aqueous flank to every analysis here. These spans are the real
    exclusion the ion-channel proxy was standing in for."""
    if not path.exists():
        return {}
    spans: dict[str, list[tuple[int, int]]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="	"):
            found = [(int(a), int(b)) for a, b in
                     re.findall(r"INTRAMEM (\d+)\.\.(\d+)", row.get("Intramembrane") or "")]
            if found:
                spans[row["Entry"]] = found
    return spans


def build_annotations(proteins: dict, accessions: list[str],
                      pore_spans: dict) -> tuple[dict, dict, dict]:
    sequons, ion_channel, pore = {}, {}, {}
    for index, accession in enumerate(accessions):
        record = proteins.get(accession)
        if record is None:
            continue
        sequons[index] = sequon_positions(record["sequence"])
        ion_channel[index] = bool(record.get("isIonChannel"))
        blocked: set[int] = set()
        for start, end in pore_spans.get(accession, ()):
            blocked |= set(range(start - PORE_MARGIN, end + PORE_MARGIN + 1))
        if blocked:
            pore[index] = blocked
    return sequons, ion_channel, pore


def annotate(sites: pd.DataFrame, sequons: dict, ion_channel: dict,
             pore: dict) -> pd.DataFrame:
    """Annotate a per-site table. site_penalties aggregates away extra columns, so
    this has to run on its output rather than on the raw variant frame."""
    if sites.empty:
        return sites
    sites = sites.copy()
    keys = list(zip(sites["protein"].to_numpy(), sites["position"].to_numpy()))
    sites["in_sequon"] = [position in sequons.get(protein, ()) for protein, position in keys]
    sites["ion_channel"] = [bool(ion_channel.get(protein, False))
                            for protein in sites["protein"].to_numpy()]
    sites["near_pore"] = [position in pore.get(protein, ()) for protein, position in keys]
    return sites


def logit(values: np.ndarray) -> np.ndarray:
    clipped = np.clip(values, 1e-4, 1 - 1e-4)
    return np.log(clipped / (1 - clipped))


def standardised(sites: pd.DataFrame, result: dict) -> dict:
    spread = float(sites["penalty"].std(ddof=1)) if len(sites) > 1 else float("nan")
    result = dict(result)
    result["penalty_sd"] = spread
    difference = result.get("side_difference")
    result["standardised_side_difference"] = (
        difference / spread if spread and np.isfinite(spread) and spread > 0 else None)
    if result.get("ci95") and all(v is not None for v in result["ci95"]) and spread:
        result["standardised_ci95"] = [v / spread for v in result["ci95"]]
    return result


def block(sites: pd.DataFrame, n_proteins: int, draws: int, seed: int,
         label: str) -> dict:
    if sites.empty:
        return {"scope": label, "n_sites": 0}
    result = standardised(sites, estimate(sites, "penalty", n_proteins, draws, seed))
    result["scope"] = label
    return result


def run(args: argparse.Namespace) -> dict:
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame, _ = load_frame(Path(args.scores))
    blob = np.load(Path(args.scores), allow_pickle=False)
    accessions = list(blob["accessions"])
    proteins = read_topology(Path(args.topology))
    pore_spans = load_pore_spans(Path(args.intramem))
    sequons, ion_channel, pore = build_annotations(proteins, accessions, pore_spans)
    n_proteins = int(frame["protein"].max()) + 1
    bounded = bool(frame["score"].min() >= 0.0 and frame["score"].max() <= 1.0)

    results: dict = {"loss": {}, "gain": {}}
    code = {aa: i for i, aa in enumerate(ALPHABET)}

    def window(sites: pd.DataFrame) -> pd.DataFrame:
        sites = sites[sites["distance"] <= PRIMARY_WINDOW]
        return annotate(sites, sequons, ion_channel, pore)

    loss_all = window(site_penalties(frame, set("KR")))
    results["loss"]["all"] = block(loss_all, n_proteins, args.draws, args.seed, "all K/R sites")
    for residue in ("K", "R"):
        subset = loss_all[loss_all["ref_aa"] == code[residue]]
        results["loss"][f"{residue}_only"] = block(
            subset, n_proteins, args.draws, args.seed, f"{residue} sites only")
    results["loss"]["excluding_sequon_sites"] = block(
        loss_all[~loss_all["in_sequon"]], n_proteins, args.draws, args.seed,
        "K/R sites not in an N-X-S/T sequon")
    results["loss"]["excluding_ion_channels"] = block(
        loss_all[~loss_all["ion_channel"]], n_proteins, args.draws, args.seed,
        "K/R sites outside ion channels")
    results["loss"]["excluding_pore_regions"] = block(
        loss_all[~loss_all["near_pore"]], n_proteins, args.draws, args.seed,
        "K/R sites away from UniProt intramembrane/pore spans")

    gain_all = window(site_gain_penalties(frame, set("KR")))
    results["gain"]["all"] = block(gain_all, n_proteins, args.draws, args.seed,
                                   "all uncharged sites")
    results["gain"]["excluding_sequon_sites"] = block(
        gain_all[~gain_all["in_sequon"]], n_proteins, args.draws, args.seed,
        "uncharged sites not in an N-X-S/T sequon")
    results["gain"]["excluding_ion_channels"] = block(
        gain_all[~gain_all["ion_channel"]], n_proteins, args.draws, args.seed,
        "uncharged sites outside ion channels")
    results["gain"]["excluding_pore_regions"] = block(
        gain_all[~gain_all["near_pore"]], n_proteins, args.draws, args.seed,
        "uncharged sites away from UniProt intramembrane/pore spans")

    if bounded:
        logit_frame = frame.copy()
        logit_frame["score"] = logit(logit_frame["score"].to_numpy())
        results["loss"]["logit_scale"] = block(
            window(site_penalties(logit_frame, set("KR"))), n_proteins,
            args.draws, args.seed, "K/R loss, logit-transformed scores")
        results["gain"]["logit_scale"] = block(
            window(site_gain_penalties(logit_frame, set("KR"))), n_proteins,
            args.draws, args.seed, "K/R gain, logit-transformed scores")

    summary = {
        "analysis": "Within-site charge audit: controls, scale checks, standardised units",
        "clinical_labels_used": False,
        "source_scores": str(args.scores),
        "scores_bounded_0_1": bounded,
        "as_of_utc": datetime.now(timezone.utc).isoformat(),
        "primary_window_residues": PRIMARY_WINDOW,
        "results": results,
    }
    (out_dir / "within_site_controls.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    for direction, blocks in results.items():
        print(f"\n=== {direction} ===")
        for key, value in blocks.items():
            if not value.get("n_sites"):
                continue
            std = value.get("standardised_side_difference")
            print(f"  {key:<28} raw {value['side_difference']:+.4f}  "
                  f"std {std:+.4f}  " if std is not None else
                  f"  {key:<28} raw {value['side_difference']:+.4f}  std n/a  ", end="")
            print(f"P(cyto>outer) {value['p_cytosolic_higher']:.4f}  sites {value['n_sites']:>7}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--topology", default=str(DEFAULT_TOPOLOGY))
    parser.add_argument("--intramem", default=str(DEFAULT_INTRAMEM))
    parser.add_argument("--draws", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260927)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
