#!/usr/bin/env python3
"""Scrape ClinMAVE API for DMS/CBGE matched missense candidates."""

from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import requests


BASE_URL = "https://ngdc.cncb.ac.cn"
OUTDIR = Path("Results/Revision/ClinMAVE_api_cross_platform")
SESSION = requests.Session()


def get_json(path: str, params: dict, timeout: int = 60) -> dict:
    url = BASE_URL + path
    for attempt in range(4):
        try:
            response = SESSION.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response.json()
        except Exception:
            if attempt == 3:
                raise
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError("unreachable")


def fetch_genes() -> pd.DataFrame:
    first = get_json(
        "/clinmave/api/fetch/table/genesummary",
        {"page": 0, "size": 1, "sort": "id,asc"},
    )
    total = int(first["totalRows"])
    data = get_json(
        "/clinmave/api/fetch/table/genesummary",
        {"page": 0, "size": total, "sort": "id,asc"},
    )
    return pd.DataFrame(data["data"])


def fetch_gene_variants(gene: str) -> pd.DataFrame:
    page_size = 1000
    params = {
        "page": 0,
        "size": page_size,
        "sort": "id,asc",
        "geneName": gene,
        "molecularConsequence": "Missense",
    }
    first = get_json("/clinmave/api/fetch/table/variant", params)
    total = int(first.get("totalRows", 0))
    if total == 0:
        return pd.DataFrame()
    rows = first.get("data", [])
    for page in range(1, (total + page_size - 1) // page_size):
        params["page"] = page
        rows.extend(get_json("/clinmave/api/fetch/table/variant", params).get("data", []))
        time.sleep(0.15)
    return pd.DataFrame(rows)


def platform(value: str) -> str | None:
    if value == "Deep Mutational Scanning":
        return "DMS"
    if value == "CRISPR-Based Genome Editing":
        return "CBGE"
    return None


def label(value: str) -> int | None:
    if value == "Functionally normal":
        return 0
    if value == "Loss-of-function":
        return 1
    return None


def summarize_gene(df: pd.DataFrame) -> dict:
    gene = df["geneName"].iloc[0]
    df = df.copy()
    df["platform"] = df["maveTechnique"].map(platform)
    df["binary_label"] = df["consequenceClass"].map(label)
    df = df.dropna(subset=["platform", "binary_label"])
    df["binary_label"] = df["binary_label"].astype(int)

    grouped = (
        df.groupby(["identifier", "platform"])["binary_label"]
        .agg(lambda x: "|".join(map(str, sorted(set(x)))))
        .reset_index()
    )
    strict = grouped[grouped["binary_label"].isin(["0", "1"])].copy()
    strict["binary_label"] = strict["binary_label"].astype(int)
    wide = strict.pivot_table(index="identifier", columns="platform", values="binary_label", aggfunc="first")
    if "DMS" not in wide.columns or "CBGE" not in wide.columns:
        matched = pd.DataFrame()
    else:
        matched = wide.dropna(subset=["DMS", "CBGE"]).copy()
        matched[["DMS", "CBGE"]] = matched[["DMS", "CBGE"]].astype(int)

    row = {
        "Gene": gene,
        "raw_missense_records": len(df),
        "unique_missense_identifiers": df["identifier"].nunique(),
        "n_matched_strict": len(matched),
        "DMS_cases": int(matched["DMS"].sum()) if len(matched) else 0,
        "DMS_controls": int((1 - matched["DMS"]).sum()) if len(matched) else 0,
        "CBGE_cases": int(matched["CBGE"].sum()) if len(matched) else 0,
        "CBGE_controls": int((1 - matched["CBGE"]).sum()) if len(matched) else 0,
        "label_concordance": float((matched["DMS"] == matched["CBGE"]).mean()) if len(matched) else float("nan"),
        "n_within_platform_conflict": int((grouped["binary_label"].str.contains("\\|", regex=True)).sum()),
    }
    return row


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    genes = fetch_genes()
    genes.to_csv(OUTDIR / "clinmave_genesummary_api.csv", index=False)
    both = genes[
        genes["infoTechnique"].str.contains("Deep Mutational Scanning", na=False)
        & genes["infoTechnique"].str.contains("CRISPR-Based Genome Editing", na=False)
    ].copy()
    both.to_csv(OUTDIR / "clinmave_genes_with_dms_and_cbge.csv", index=False)

    rows = []
    failures = []
    raw_path = OUTDIR / "clinmave_dms_cbge_gene_variants.jsonl"
    with raw_path.open("w") as handle:
        for i, gene in enumerate(both["geneName"], start=1):
            print(f"[{i}/{len(both)}] {gene}", flush=True)
            try:
                variants = fetch_gene_variants(gene)
                if not variants.empty:
                    for record in variants.to_dict(orient="records"):
                        handle.write(json.dumps(record, ensure_ascii=False) + "\n")
                    rows.append(summarize_gene(variants))
            except Exception as exc:
                failures.append({"Gene": gene, "error": str(exc)})
            time.sleep(0.2)

    summary = pd.DataFrame(rows).sort_values(["n_matched_strict", "Gene"], ascending=[False, True])
    summary.to_csv(OUTDIR / "clinmave_cross_platform_candidate_summary.csv", index=False)
    pd.DataFrame(failures).to_csv(OUTDIR / "clinmave_cross_platform_scrape_failures.csv", index=False)

    eligible = summary[
        (summary["DMS_cases"] >= 5)
        & (summary["DMS_controls"] >= 5)
        & (summary["CBGE_cases"] >= 5)
        & (summary["CBGE_controls"] >= 5)
    ].copy()
    eligible.to_csv(OUTDIR / "clinmave_cross_platform_eligible_old_style.csv", index=False)

    print("\nGenes with both DMS and CBGE:", len(both))
    print("Eligible strict old-style matched genes:", len(eligible))
    print(eligible.to_string(index=False))
    print(f"\nWrote outputs to {OUTDIR}")


if __name__ == "__main__":
    main()
