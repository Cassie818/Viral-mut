#!/usr/bin/env python3
"""ESM-2 650M scoring for ClinMAVE missense variants.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from esm import pretrained


AA_COLS = list("ACDEFGHIKLMNPQRSTVWY")
DEFAULT_INPUT_GLOB = "Results/ClinMAVE/*/missense/*_LLR_results.csv"
DEFAULT_PROTEIN_DIR = Path("/Users/cassie/Desktop/Protein")
DEFAULT_OUT_DIR = Path("Results/Revision/ClinMAVE_ESM2_650M")


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip().upper()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith(">")
    )


def load_inputs(pattern: str) -> pd.DataFrame:
    paths = sorted(Path(".").glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No input files matched {pattern!r}")
    frames = []
    for path in paths:
        df = pd.read_csv(path)
        df["source_file"] = str(path)
        df["assay"] = path.parts[2] if len(path.parts) > 2 else "unknown"
        name = path.name.replace("_LLR_results.csv", "")
        df["effect_class"] = name.replace("_dms", "").replace("_cbge", "")
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    required = ["Identifier", "Gene", "Site", "Ref", "Mut", "LLR", "source_file"]
    missing = [col for col in required if col not in out.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    out["Site"] = pd.to_numeric(out["Site"], errors="coerce").astype("Int64")
    return out


def variant_key(row: pd.Series) -> str:
    return "|".join(
        [
            str(row["Gene"]),
            str(int(row["Site"])),
            str(row["Ref"]),
            str(row["Mut"]),
        ]
    )


def make_variant_table(df: pd.DataFrame, protein_dir: Path, max_len: int) -> pd.DataFrame:
    variants = (
        df.dropna(subset=["Gene", "Site", "Ref", "Mut"])
        .drop_duplicates(subset=["Gene", "Site", "Ref", "Mut"])
        .copy()
    )
    rows = []
    seq_cache: dict[str, str | None] = {}
    for _, row in variants.iterrows():
        gene = str(row["Gene"])
        if gene not in seq_cache:
            fasta = protein_dir / f"{gene}_protein.fasta"
            seq_cache[gene] = read_fasta(fasta) if fasta.exists() else None
        seq = seq_cache[gene]
        site = int(row["Site"])
        rows.append(
            {
                **row.to_dict(),
                "variant_key": variant_key(row),
                "protein_length": len(seq) if seq is not None else np.nan,
                "has_fasta": seq is not None,
                "length_compatible": seq is not None and len(seq) <= max_len,
                "site_in_range": seq is not None and 1 <= site <= len(seq),
                "ref_matches_fasta": seq is not None and 1 <= site <= len(seq) and seq[site - 1] == str(row["Ref"]),
            }
        )
    out = pd.DataFrame(rows)
    return out


def append_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def load_done(path: Path) -> set[str]:
    if not path.exists() or path.stat().st_size == 0:
        return set()
    return set(pd.read_csv(path, usecols=["variant_key"])["variant_key"].astype(str))


def score_gene(
    model: torch.nn.Module,
    alphabet,
    batch_converter,
    gene: str,
    group: pd.DataFrame,
    protein_dir: Path,
    device: torch.device,
) -> list[dict[str, object]]:
    sequence = read_fasta(protein_dir / f"{gene}_protein.fasta")
    _, _, tokens = batch_converter([(gene, sequence)])
    tokens = tokens.to(device)
    with torch.no_grad():
        logits = model(tokens, repr_layers=[], return_contacts=False)["logits"][0, 1 : len(sequence) + 1]
        log_probs = F.log_softmax(logits, dim=-1).detach().cpu().numpy()
    aa_to_idx = alphabet.tok_to_idx
    rows = []
    for _, row in group.iterrows():
        site = int(row["Site"])
        ref = str(row["Ref"])
        mut = str(row["Mut"])
        if ref not in aa_to_idx or mut not in aa_to_idx:
            continue
        llr = float(log_probs[site - 1, aa_to_idx[mut]] - log_probs[site - 1, aa_to_idx[ref]])
        rows.append(
            {
                "variant_key": row["variant_key"],
                "Gene": gene,
                "Site": site,
                "Ref": ref,
                "Mut": mut,
                "protein_length": int(row["protein_length"]),
                "esm2_650m_llr": llr,
                "esm2_650m_score": -llr,
            }
        )
    return rows


def choose_device(requested: str) -> torch.device:
    if requested != "auto":
        return torch.device(requested)
    if torch.backends.mps.is_available():
        return torch.device("mps")
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-glob", default=DEFAULT_INPUT_GLOB)
    parser.add_argument("--protein-dir", type=Path, default=DEFAULT_PROTEIN_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--max-len", type=int, default=1022)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--max-genes", type=int, default=None)
    parser.add_argument("--report-every", type=int, default=10)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    all_inputs = load_inputs(args.input_glob)
    variant_table = make_variant_table(all_inputs, args.protein_dir, args.max_len)
    variant_table.to_csv(args.out_dir / "clinmave_missense_variant_audit.csv", index=False)

    score_path = args.out_dir / "clinmave_missense_esm2_650m_variant_scores.csv"
    if args.force and score_path.exists():
        score_path.unlink()
    done = load_done(score_path)
    scorable = variant_table[
        variant_table["length_compatible"]
        & variant_table["site_in_range"]
        & variant_table["ref_matches_fasta"]
        & variant_table["Ref"].isin(AA_COLS)
        & variant_table["Mut"].isin(AA_COLS)
    ].copy()
    scorable = scorable[~scorable["variant_key"].astype(str).isin(done)]
    genes = sorted(scorable["Gene"].astype(str).unique())
    if args.max_genes is not None:
        genes = genes[: args.max_genes]

    audit = {
        "input_rows": len(all_inputs),
        "unique_variants": len(variant_table),
        "scorable_unique_variants": int(
            (
                variant_table["length_compatible"]
                & variant_table["site_in_range"]
                & variant_table["ref_matches_fasta"]
                & variant_table["Ref"].isin(AA_COLS)
                & variant_table["Mut"].isin(AA_COLS)
            ).sum()
        ),
        "already_scored_variants": len(done),
        "remaining_variants_this_run": len(scorable),
        "remaining_genes_this_run": len(genes),
        "missing_fasta_variants": int((~variant_table["has_fasta"]).sum()),
        "length_incompatible_variants": int((variant_table["has_fasta"] & ~variant_table["length_compatible"]).sum()),
        "ref_mismatch_variants": int((variant_table["length_compatible"] & variant_table["site_in_range"] & ~variant_table["ref_matches_fasta"]).sum()),
    }
    pd.DataFrame([audit]).to_csv(args.out_dir / "clinmave_missense_esm2_650m_run_audit.csv", index=False)
    print(pd.DataFrame([audit]).to_string(index=False))

    if genes:
        device = choose_device(args.device)
        print(f"Loading ESM-2 650M on {device}...")
        model, alphabet = pretrained.esm2_t33_650M_UR50D()
        model.eval().to(device)
        batch_converter = alphabet.get_batch_converter()
        for idx, gene in enumerate(genes, start=1):
            group = scorable[scorable["Gene"].astype(str) == gene]
            rows = score_gene(model, alphabet, batch_converter, gene, group, args.protein_dir, device)
            append_rows(score_path, rows)
            if idx % args.report_every == 0 or idx == len(genes):
                print(f"Scored {idx}/{len(genes)} genes; latest={gene}; variants_written={len(rows)}")

    scores = pd.read_csv(score_path) if score_path.exists() else pd.DataFrame()
    if not scores.empty:
        merged = all_inputs.copy()
        merged["variant_key"] = merged.apply(variant_key, axis=1)
        merged = merged.merge(scores, on="variant_key", how="left", suffixes=("", "_scored"))
        merged.to_csv(args.out_dir / "clinmave_missense_all_with_esm2_650m.csv", index=False)
        for source, group in merged.groupby("source_file", sort=True):
            source_path = Path(source)
            out_name = source_path.name.replace("_LLR_results.csv", "_LLR_ESM2_650M_results.csv")
            out_subdir = args.out_dir / source_path.parts[2] / source_path.parts[3]
            out_subdir.mkdir(parents=True, exist_ok=True)
            cols = [
                "Identifier",
                "Chr",
                "Gene",
                "Transcriptid",
                "Site",
                "Ref",
                "Mut",
                "LLR",
                "esm2_650m_llr",
                "esm2_650m_score",
                "protein_length",
            ]
            group[[col for col in cols if col in group.columns]].to_csv(out_subdir / out_name, index=False)
        print(f"Wrote {score_path}")
        print(f"Wrote {args.out_dir / 'clinmave_missense_all_with_esm2_650m.csv'}")


if __name__ == "__main__":
    main()
