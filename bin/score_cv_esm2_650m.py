#!/usr/bin/env python3
"""ClinVar missense scoring with ESM-2 650M.
"""

from __future__ import annotations

import argparse
import csv
import gc
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F


AA_ORDER = set("ACDEFGHIKLMNPQRSTVWY")
PATHOGENIC_LABELS = {"pathogenic", "likely_pathogenic"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clinvar-dir", default="Results/ClinVar/missense")
    parser.add_argument("--protein-dir", default="/Users/cassie/Desktop/Protein")
    parser.add_argument("--cache-dir", default="Results/Revision/model_cache")
    parser.add_argument(
        "--output",
        default="Results/Revision/esm2_650m_full/clinvar_missense_esm2_650m_scores.csv",
    )
    parser.add_argument(
        "--failed-output",
        default="Results/Revision/esm2_650m_full/failed_genes.csv",
    )
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda", "mps"])
    parser.add_argument(
        "--max-genes",
        type=int,
        default=0,
        help="Optional smoke-test limit. 0 means all remaining genes.",
    )
    parser.add_argument(
        "--report-every",
        type=int,
        default=25,
        help="Print progress after this many completed genes.",
    )
    parser.add_argument(
        "--min-pos",
        type=int,
        default=0,
        help="Require at least this many pathogenic/likely pathogenic variants per gene.",
    )
    parser.add_argument(
        "--min-neg",
        type=int,
        default=0,
        help="Require at least this many benign/likely benign variants per gene.",
    )
    parser.add_argument(
        "--sort-by-length",
        action="store_true",
        help="Run shorter proteins first after applying label-count filters.",
    )
    return parser.parse_args()


def choose_device(requested: str) -> torch.device:
    if requested == "cuda":
        return torch.device("cuda")
    if requested == "mps":
        return torch.device("mps")
    if requested == "cpu":
        return torch.device("cpu")
    if torch.cuda.is_available():
        return torch.device("cuda")
    if getattr(torch.backends, "mps", None) is not None and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip().upper()
        for line in path.read_text(errors="ignore").splitlines()
        if not line.startswith(">")
    )


def load_clinvar(clinvar_dir: Path) -> pd.DataFrame:
    files = [
        clinvar_dir / "missense_benign.csv",
        clinvar_dir / "missense_likely_benign.csv",
        clinvar_dir / "missense_likely_pathogenic.csv",
        clinvar_dir / "missense_pathogenic.csv",
    ]
    missing = [str(path) for path in files if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing ClinVar files: " + ", ".join(missing))
    df = pd.concat([pd.read_csv(path) for path in files], ignore_index=True)
    df = df.dropna(
        subset=["Gene_prot", "Site_prot", "Ref_prot", "Mut_prot", "LLR_prot", "LLR_gene"]
    ).copy()
    df["Site_prot"] = df["Site_prot"].astype(int)
    df["label"] = df["Label_prot"].isin(PATHOGENIC_LABELS).astype(int)
    df["variant_id"] = (
        df["Gene_prot"].astype(str)
        + ":"
        + df["Site_prot"].astype(str)
        + ":"
        + df["Ref_prot"].astype(str)
        + ">"
        + df["Mut_prot"].astype(str)
        + ":"
        + df["Ref_gene"].astype(str)
        + ">"
        + df["Mut_gene"].astype(str)
    )
    return df


def load_done_genes(output: Path) -> set[str]:
    if not output.exists() or output.stat().st_size == 0:
        return set()
    try:
        return set(pd.read_csv(output, usecols=["Gene_prot"])["Gene_prot"].astype(str).unique())
    except Exception:
        return set()


def append_rows(output: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    write_header = not output.exists() or output.stat().st_size == 0
    with output.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def append_failure(failed_output: Path, gene: str, n_variants: int, seq_len: int | None, error: str) -> None:
    failed_output.parent.mkdir(parents=True, exist_ok=True)
    write_header = not failed_output.exists() or failed_output.stat().st_size == 0
    with failed_output.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Gene_prot", "n_variants", "seq_len", "error"])
        if write_header:
            writer.writeheader()
        writer.writerow(
            {
                "Gene_prot": gene,
                "n_variants": n_variants,
                "seq_len": "" if seq_len is None else seq_len,
                "error": error,
            }
        )


def load_model(cache_dir: Path, device: torch.device):
    import esm

    torch.hub.set_dir(str(cache_dir / "torch_hub"))
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.eval().to(device)
    return model, alphabet, alphabet.get_batch_converter()


def score_gene(
    gene: str,
    group: pd.DataFrame,
    sequence: str,
    model,
    alphabet,
    batch_converter,
    device: torch.device,
) -> list[dict[str, object]]:
    aa_to_idx = {alphabet.get_tok(i): i for i in range(len(alphabet))}
    valid_indices = []
    for idx, row in group.iterrows():
        site = int(row["Site_prot"])
        ref = str(row["Ref_prot"])
        mut = str(row["Mut_prot"])
        if (
            1 <= site <= len(sequence)
            and ref in AA_ORDER
            and mut in AA_ORDER
            and sequence[site - 1] == ref
        ):
            valid_indices.append(idx)

    if not valid_indices:
        return []

    valid = group.loc[valid_indices]
    with torch.no_grad():
        _, _, tokens = batch_converter([(gene, sequence)])
        tokens = tokens.to(device)
        logits = model(tokens, repr_layers=[], return_contacts=False)["logits"][0, 1 : len(sequence) + 1]
        log_probs = F.log_softmax(logits, dim=-1).detach().cpu().numpy()

    rows = []
    for _, row in valid.iterrows():
        site_zero = int(row["Site_prot"]) - 1
        ref_idx = aa_to_idx[row["Ref_prot"]]
        mut_idx = aa_to_idx[row["Mut_prot"]]
        llr = float(log_probs[site_zero, mut_idx] - log_probs[site_zero, ref_idx])
        rows.append(
            {
                "Label_prot": row["Label_prot"],
                "Gene_prot": row["Gene_prot"],
                "Site_prot": int(row["Site_prot"]),
                "Ref_prot": row["Ref_prot"],
                "Mut_prot": row["Mut_prot"],
                "LLR_prot_150m": float(row["LLR_prot"]),
                "Label_gene": row["Label_gene"],
                "Gene_gene": row["Gene_gene"],
                "Site_gene": int(row["Site_gene"]),
                "Ref_gene": row["Ref_gene"],
                "Mut_gene": row["Mut_gene"],
                "LLR_gene": float(row["LLR_gene"]),
                "label": int(row["label"]),
                "variant_id": row["variant_id"],
                "protein_length": len(sequence),
                "esm2_150m_score": -float(row["LLR_prot"]),
                "calm_existing_score": -float(row["LLR_gene"]),
                "esm2_650m_llr": llr,
                "esm2_650m_score": -llr,
            }
        )
    return rows


def filter_genes_by_label_counts(df: pd.DataFrame, min_pos: int, min_neg: int) -> set[str]:
    counts = df.groupby("Gene_prot")["label"].agg(n="size", n_pos="sum")
    counts["n_neg"] = counts["n"] - counts["n_pos"]
    keep = counts[(counts["n_pos"] >= min_pos) & (counts["n_neg"] >= min_neg)]
    return set(keep.index.astype(str))


def protein_lengths(genes: list[str], protein_dir: Path) -> dict[str, int]:
    lengths = {}
    for gene in genes:
        fasta = protein_dir / f"{gene}_protein.fasta"
        if fasta.exists():
            try:
                lengths[gene] = len(read_fasta(fasta))
            except Exception:
                lengths[gene] = 10**9
        else:
            lengths[gene] = 10**9
    return lengths


def main() -> None:
    args = parse_args()
    output = Path(args.output)
    failed_output = Path(args.failed_output)
    device = choose_device(args.device)

    df = load_clinvar(Path(args.clinvar_dir))
    done = load_done_genes(output)
    eligible = filter_genes_by_label_counts(df, args.min_pos, args.min_neg)
    genes = sorted(gene for gene in df["Gene_prot"].astype(str).unique() if gene in eligible and gene not in done)
    if args.sort_by_length:
        lengths = protein_lengths(genes, Path(args.protein_dir))
        genes = sorted(genes, key=lambda gene: (lengths.get(gene, 10**9), gene))
    if args.max_genes > 0:
        genes = genes[: args.max_genes]

    print(
        f"Loaded {len(df)} variants across {df['Gene_prot'].nunique()} genes. "
        f"Filter min_pos={args.min_pos}, min_neg={args.min_neg}. "
        f"{len(done)} genes already scored in output; {len(genes)} genes to run on {device}.",
        flush=True,
    )

    model, alphabet, batch_converter = load_model(Path(args.cache_dir), device)
    fieldnames = [
        "Label_prot",
        "Gene_prot",
        "Site_prot",
        "Ref_prot",
        "Mut_prot",
        "LLR_prot_150m",
        "Label_gene",
        "Gene_gene",
        "Site_gene",
        "Ref_gene",
        "Mut_gene",
        "LLR_gene",
        "label",
        "variant_id",
        "protein_length",
        "esm2_150m_score",
        "calm_existing_score",
        "esm2_650m_llr",
        "esm2_650m_score",
    ]

    start = time.time()
    completed = 0
    written_variants = 0
    for gene in genes:
        group = df[df["Gene_prot"].astype(str) == gene]
        fasta = Path(args.protein_dir) / f"{gene}_protein.fasta"
        seq_len = None
        try:
            if not fasta.exists():
                raise FileNotFoundError(f"Missing FASTA: {fasta}")
            sequence = read_fasta(fasta)
            seq_len = len(sequence)
            rows = score_gene(gene, group, sequence, model, alphabet, batch_converter, device)
            if rows:
                append_rows(output, rows, fieldnames)
                written_variants += len(rows)
            else:
                append_failure(failed_output, gene, len(group), seq_len, "No valid reference-matching variants")
        except Exception as exc:
            append_failure(failed_output, gene, len(group), seq_len, repr(exc))
        finally:
            completed += 1
            gc.collect()

        if completed == 1 or completed % args.report_every == 0:
            elapsed = time.time() - start
            rate = completed / elapsed if elapsed else 0.0
            remaining = len(genes) - completed
            eta_hours = remaining / rate / 3600 if rate else float("nan")
            print(
                f"Progress: {completed}/{len(genes)} genes this run, "
                f"{written_variants} variants written, {rate*3600:.1f} genes/hour, "
                f"ETA {eta_hours:.2f} h",
                flush=True,
            )

    print(f"Done. Wrote {written_variants} variants to {output}", flush=True)
    print(f"Failures, if any, are in {failed_output}", flush=True)


if __name__ == "__main__":
    main()
