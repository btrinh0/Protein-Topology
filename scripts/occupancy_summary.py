from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Produce the first quantitative occupancy summary from the atlas table."
    )
    parser.add_argument(
        "--atlas",
        type=Path,
        default=Path("data/processed/domain_atlas_table.csv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/occupancy_summary.csv"),
    )
    return parser.parse_args()


def assign_occupancy_tier(count: int) -> str:
    if count >= 10:
        return "common"
    if count >= 3:
        return "rare"
    if count >= 1:
        return "ultra-rare"
    return "unoccupied"


def main() -> int:
    args = parse_args()
    df = pd.read_csv(args.atlas)

    print(f"Atlas table: {len(df)} rows")
    print(f"Unique SCOPe fold classes: {df['scope_id'].nunique()}")
    print()

    for col in ["ps_fraction", "pp_fraction", "x_fraction"]:
        vals = pd.to_numeric(df[col], errors="coerce").dropna()
        print(f"{col}: mean={vals.mean():.4f}, std={vals.std():.4f}, "
              f"min={vals.min():.4f}, max={vals.max():.4f}")
    print()

    print("Kingdom distribution:")
    kingdom_counts = df["kingdom"].fillna("unknown").value_counts()
    for kingdom, count in kingdom_counts.items():
        print(f"  {kingdom}: {count}")
    print()

    scope_fold_counts = df["scope_id"].value_counts()
    fold_tier_map = {fold: assign_occupancy_tier(count)
                     for fold, count in scope_fold_counts.items()}
    tier_dist = pd.Series(fold_tier_map).value_counts()

    print("SCOPe fold occupancy tiers (by unique fold class):")
    for tier in ["common", "rare", "ultra-rare"]:
        count = tier_dist.get(tier, 0)
        print(f"  {tier}: {count} fold classes")
    total_folds = len(scope_fold_counts)
    print(f"  Total occupied fold classes: {total_folds}")
    print()

    cath_coverage = (df["cath_id"].astype(str).str.len() > 0).sum()
    ecod_coverage = (df["ecod_id"].astype(str).str.len() > 0).sum()
    print(f"CATH cross-reference coverage: {cath_coverage}/{len(df)} ({100*cath_coverage/len(df):.1f}%)")
    print(f"ECOD cross-reference coverage: {ecod_coverage}/{len(df)} ({100*ecod_coverage/len(df):.1f}%)")
    print()

    scope_class = df["scope_id"].str[0]
    print("SCOPe major class distribution:")
    class_names = {
        "a": "all-alpha",
        "b": "all-beta",
        "c": "alpha/beta (a/b)",
        "d": "alpha+beta (a+b)",
        "e": "multi-domain",
        "f": "membrane/cell surface",
        "g": "small proteins",
    }
    for cls, name in class_names.items():
        count = (scope_class == cls).sum()
        if count > 0:
            folds_in_class = df[scope_class == cls]["scope_id"].nunique()
            print(f"  {cls} ({name}): {count} domains, {folds_in_class} folds")
    print()

    df["occupancy_tier"] = df["scope_id"].map(fold_tier_map)
    df.to_csv(args.atlas, index=False)

    summary_rows = []
    for fold, count in scope_fold_counts.items():
        tier = fold_tier_map[fold]
        subset = df[df["scope_id"] == fold]
        kingdoms = subset["kingdom"].fillna("unknown").value_counts().to_dict()
        summary_rows.append({
            "scope_id": fold,
            "domain_count": count,
            "occupancy_tier": tier,
            "kingdoms": str(kingdoms),
        })
    summary_df = pd.DataFrame(summary_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(args.output, index=False)
    print(f"Occupancy summary saved to {args.output}")

    print()
    print("--- HEADLINE NUMBERS ---")
    print(f"Under a bounded topology definition for {len(df)} single-domain proteins")
    print(f"(ASTRAL 40%, 60-180 residues, SCOPe 2.08),")
    print(f"natural structures occupy {total_folds} distinct fold classes:")
    common = tier_dist.get("common", 0)
    rare = tier_dist.get("rare", 0)
    ultra_rare = tier_dist.get("ultra-rare", 0)
    print(f"  {common} common, {rare} rare, {ultra_rare} ultra-rare.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
