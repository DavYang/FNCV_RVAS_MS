import marimo

__generated_with = "0.20.4"
app = marimo.App(
    width="full",
    app_title="Phase 2d: GWAS Catalog EUR MS Loci (HPC)",
)


@app.cell(hide_code=True)
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    # Phase 2d: GWAS Catalog EUR MS Loci

    **Environment:** HPC (`/gs/gsfs0/shared-lab/greally-lab/David/GA4K/6_AoU-AxA/FNCV_RVAS_MS`)

    This notebook runs **step 2d** of the Phase 2 pipeline on the HPC — parsing the raw GWAS
    Catalog TSV to extract EUR-only Multiple Sclerosis loci, lifting coordinates from GRCh37
    to GRCh38, and uploading the result to the AoU workspace GCS bucket for use in step 2e.

    | Step | Script | Environment |
    |---|---|---|
    | 2a | `02a_export_gwas_ma.sh` | AoU VM |
    | 2c | `02c_export_top_gwas_snps.sh` | AoU VM |
    | **2d (this notebook)** | `02d_parse_gwas_catalog.sh` | **HPC** |
    | 2e | `02e_define_loci_catalog.sh` | AoU VM |

    ## Filtering criteria
    - `DISEASE/TRAIT` is **exactly** `"Multiple sclerosis"` (case-insensitive exact match;
      excludes MTAG, OCB status, drug-induced, severity, and relapse sub-traits)
    - `P-VALUE` < 5e-8, autosomes only, non-null `CHR_POS`
    - **EUR-only:** `INITIAL SAMPLE SIZE` must contain `"European ancestry"` and must not
      mention any non-EUR or unknown ancestry group — multi-ancestry GWAS excluded entirely
    - Liftover `CHR_POS` GRCh37 → GRCh38 via `pyliftover` + UCSC `hg19ToHg38.over.chain.gz`
    - Deduplicated on `chrom + pos_hg38 + rsid` (lowest p-value kept)
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("## 0. Setup & Configuration")
    return


@app.cell
def _():
    import os
    import json
    import subprocess

    PROJECT_ROOT = "/gs/gsfs0/shared-lab/greally-lab/David/GA4K/6_AoU-AxA/FNCV_RVAS_MS"
    os.chdir(PROJECT_ROOT)

    CONFIG_PATH = "config/config.json"
    with open(CONFIG_PATH) as _f:
        config = json.load(_f)

    CATALOG_RAW   = config["inputs"]["gwas_catalog_raw"]
    CHAIN_FILE    = config["inputs"]["liftover_chain"]
    LOCUS_DEF_DIR = config["outputs"]["locus_def_dir"]
    GCS_DEST      = config["params"]["gwas_catalog_gcs_path"]
    P_THRESHOLD   = config["params"]["gwas_catalog_p_threshold"]
    MHC_INTERVAL  = config["params"]["mhc_interval"]
    OUTPUT_FILE   = os.path.join(LOCUS_DEF_DIR, "gwas_catalog_ms_eur_hg38.tsv")
    PLOT_PATH     = os.path.join(LOCUS_DEF_DIR, "manhattan_gwas_catalog_ms_eur.png")

    print(f"CWD           : {os.getcwd()}")
    print(f"Catalog input : {CATALOG_RAW}")
    print(f"Chain file    : {CHAIN_FILE}")
    print(f"Output        : {OUTPUT_FILE}")
    print(f"GCS dest      : {GCS_DEST}")
    print(f"P threshold   : {P_THRESHOLD}")
    print(f"MHC interval  : {MHC_INTERVAL}")
    return (
        CATALOG_RAW,
        CHAIN_FILE,
        CONFIG_PATH,
        GCS_DEST,
        LOCUS_DEF_DIR,
        MHC_INTERVAL,
        OUTPUT_FILE,
        P_THRESHOLD,
        PLOT_PATH,
        config,
        json,
        os,
        subprocess,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md("## 1. Pre-flight Checks")
    return


@app.cell
def _(CATALOG_RAW, CHAIN_FILE, OUTPUT_FILE, os):
    import shutil

    checks_passed = True
    results = []

    if os.path.exists(CATALOG_RAW):
        size_gb = os.path.getsize(CATALOG_RAW) / 1e9
        results.append(f"[OK]   Catalog TSV   : {CATALOG_RAW} ({size_gb:.2f} GB)")
    else:
        results.append(f"[FAIL] Catalog TSV not found: {CATALOG_RAW}")
        checks_passed = False

    py_script = "python/02d_parse_gwas_catalog.py"
    if os.path.exists(py_script):
        results.append(f"[OK]   Python script : {py_script}")
    else:
        results.append(f"[FAIL] Python script not found: {py_script}")
        checks_passed = False

    try:
        import pyliftover
        results.append(f"[OK]   pyliftover    : installed")
    except ImportError:
        results.append("[FAIL] pyliftover not installed -- pip install pyliftover")
        checks_passed = False

    try:
        import pandas as _pd
        results.append(f"[OK]   pandas        : {_pd.__version__}")
    except ImportError:
        results.append("[FAIL] pandas not installed")
        checks_passed = False

    if shutil.which("gsutil"):
        results.append("[OK]   gsutil        : found")
    else:
        results.append("[WARN] gsutil not found -- upload step must be run separately")

    if os.path.exists(CHAIN_FILE):
        results.append(f"[OK]   Chain file    : cached at {CHAIN_FILE}")
    else:
        results.append(f"[INFO] Chain file    : not cached, will download from UCSC at runtime")

    if os.path.exists(OUTPUT_FILE):
        results.append(f"[INFO] Output file   : already exists (use --force to overwrite)")
    else:
        results.append(f"[INFO] Output file   : not yet generated")

    print("\n".join(results))
    print()
    print("All checks passed -- ready to run." if checks_passed else "One or more checks FAILED.")
    return checks_passed, results, shutil


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## 2. Run GWAS Catalog Parser

    Runs `bash/02d_parse_gwas_catalog.sh` which:
    1. Validates inputs from `config/config.json`
    2. Calls `python/02d_parse_gwas_catalog.py` (chunked 50k-row reads, ~500 MB TSV)
    3. Filters to EUR-only MS associations at p < 5e-8 (exact trait match)
    4. Downloads chain file if not cached, lifts GRCh37 → GRCh38
    5. Writes `results/2-locus_definition/gwas_catalog_ms_eur_hg38.tsv`

    Runtime: typically **2–5 minutes**. Logs: `logs/02d_parse_gwas_catalog_<timestamp>.log`

    > To overwrite an existing output, set `FORCE = True` below.
    """)
    return


@app.cell
def _():
    FORCE = False
    return (FORCE,)


@app.cell
def _(FORCE, subprocess):
    _force_flag = "--force" if FORCE else ""
    _cmd = f"bash bash/02d_parse_gwas_catalog.sh {_force_flag}".strip()
    print(f"Running: {_cmd}")
    _result = subprocess.run(_cmd, shell=True, text=True)
    if _result.returncode != 0:
        raise RuntimeError(f"Script exited with code {_result.returncode}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("## 3. Inspect Output")
    return


@app.cell
def _(OUTPUT_FILE, os):
    import pandas as pd

    if not os.path.exists(OUTPUT_FILE):
        raise FileNotFoundError(
            f"Output not found: {OUTPUT_FILE} -- did step 2 complete successfully?"
        )

    df = pd.read_csv(OUTPUT_FILE, sep="\t")
    print(f"Rows    : {len(df):,}")
    print(f"Columns : {list(df.columns)}")
    df.head(10)
    return df, pd


@app.cell
def _(df):
    print("Loci per chromosome:")
    print(df.groupby("chrom").size().to_string())
    print()
    print(f"Unique PubMed IDs (studies) : {df['pubmed_id'].nunique()}")
    print(f"P-value range               : {df['p_value'].min():.2e} -- {df['p_value'].max():.2e}")
    print(f"Chromosomes represented     : {sorted(df['chrom'].unique())}")
    return


@app.cell
def _(df):
    _non_eur_kw = [
        "african", "asian", "japanese", "hispanic", "admixed",
        "latino", "latina", "south asian", "east asian", "multi",
        "chinese", "korean", "indigenous", "amerindian",
        "native american", "unknown ancestry",
    ]
    _contaminated = df[
        df["sample_size"].str.lower().str.contains("|".join(_non_eur_kw), na=False)
    ]
    if len(_contaminated) == 0:
        print("[OK] No multi-ancestry rows found in output.")
    else:
        print(f"[WARN] {len(_contaminated)} rows with non-EUR ancestry terms survived filtering:")
        print(_contaminated[["chrom", "pos_hg38", "rsid", "sample_size"]].to_string())
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("## 4. Manhattan Plot")
    return


@app.cell
def _(MHC_INTERVAL, PLOT_PATH, df, os, pd):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _mhc_chrom, _mhc_coords = MHC_INTERVAL.split(":")
    _mhc_start, _mhc_end = [int(x) for x in _mhc_coords.split("-")]

    CHROM_ORDER = [f"chr{i}" for i in range(1, 23)]

    _plot_df = df[df["chrom"].isin(CHROM_ORDER)].copy()
    _mhc_mask = (
        (_plot_df["chrom"] == _mhc_chrom)
        & (_plot_df["pos_hg38"] >= _mhc_start)
        & (_plot_df["pos_hg38"] <= _mhc_end)
    )
    _plot_df = _plot_df[~_mhc_mask].copy()
    print(f"MHC loci excluded ({MHC_INTERVAL}): {_mhc_mask.sum():,}")
    print(f"Loci retained for plot             : {len(_plot_df):,}")

    _plot_df["chrom"] = pd.Categorical(
        _plot_df["chrom"], categories=CHROM_ORDER, ordered=True
    )
    _plot_df = _plot_df.sort_values(["chrom", "pos_hg38"]).reset_index(drop=True)

    CHROM_SIZES = {
        "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
        "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
        "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
        "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
        "chr21": 46709983, "chr22": 50818468,
    }
    CHROM_GAP = 5_000_000

    _offsets = {}
    _offset = 0
    for _c in CHROM_ORDER:
        _offsets[_c] = _offset
        _offset += CHROM_SIZES.get(_c, 0) + CHROM_GAP

    _plot_df["x"] = _plot_df.apply(
        lambda r: _offsets[r["chrom"]] + r["pos_hg38"], axis=1
    )
    _plot_df["neg_log10_p"] = -np.log10(_plot_df["p_value"].clip(lower=1e-300))

    _xtick_pos = [_offsets[c] + CHROM_SIZES[c] / 2 for c in CHROM_ORDER]
    _xtick_labels = [c.replace("chr", "") for c in CHROM_ORDER]

    PALETTE = ["#1f4e79", "#2e86ab"]
    _colour_map = {c: PALETTE[i % 2] for i, c in enumerate(CHROM_ORDER)}

    _fig, _ax = plt.subplots(figsize=(18, 5))

    for _chrom in CHROM_ORDER:
        _sub = _plot_df[_plot_df["chrom"] == _chrom]
        _ax.scatter(
            _sub["x"], _sub["neg_log10_p"],
            c=_colour_map[_chrom], s=12, alpha=0.85, linewidths=0, rasterized=True
        )

    _sig_y = -np.log10(5e-8)
    _ax.axhline(_sig_y, color="crimson", linewidth=0.8, linestyle="--", zorder=3)
    _ax.text(
        _plot_df["x"].max() * 1.002, _sig_y, "p=5e-8",
        va="center", ha="left", fontsize=7, color="crimson"
    )

    _ax.set_xticks(_xtick_pos)
    _ax.set_xticklabels(_xtick_labels, fontsize=7)
    _ax.set_xlabel("Chromosome", fontsize=10)
    _ax.set_ylabel("$-\\log_{10}(p)$", fontsize=10)
    _ax.set_title(
        f"EUR-only MS GWAS Catalog Loci (GRCh38, MHC excluded)  |  n={len(_plot_df):,} loci",
        fontsize=11,
    )
    _ax.set_xlim(0, _plot_df["x"].max() * 1.01)
    _ax.tick_params(axis="x", length=3)
    _ax.spines["top"].set_visible(False)
    _ax.spines["right"].set_visible(False)

    plt.tight_layout()

    os.makedirs(os.path.dirname(PLOT_PATH), exist_ok=True)
    _fig.savefig(PLOT_PATH, dpi=200, bbox_inches="tight", format="png")
    print(f"Saved: {PLOT_PATH}")

    _fig
    return CHROM_ORDER, CHROM_SIZES, PALETTE, matplotlib, np, plt


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## 5. Upload to AoU GCS Workspace Bucket

    Uploads `gwas_catalog_ms_eur_hg38.tsv` to the AoU workspace bucket so that
    `02e_define_loci_catalog.sh` on the AoU workbench can download and use it.
    """)
    return


@app.cell
def _(GCS_DEST, OUTPUT_FILE, os, subprocess):
    if not os.path.exists(OUTPUT_FILE):
        raise FileNotFoundError(f"Output not found: {OUTPUT_FILE}")

    print(f"Uploading : {OUTPUT_FILE}")
    print(f"        to: {GCS_DEST}")
    print()

    _upload = subprocess.run(
        ["gsutil", "cp", OUTPUT_FILE, GCS_DEST],
        text=True, capture_output=True
    )
    if _upload.returncode == 0:
        print("[OK] Upload complete.")
    else:
        print(f"[FAIL] gsutil exited {_upload.returncode}:")
        print(_upload.stderr)
    return


@app.cell
def _(GCS_DEST, subprocess):
    _verify = subprocess.run(
        ["gsutil", "ls", "-lh", GCS_DEST],
        text=True, capture_output=True
    )
    if _verify.returncode == 0:
        print("[OK] Upload verified:")
        print(_verify.stdout.strip())
    else:
        print("[FAIL] Could not verify upload:")
        print(_verify.stderr.strip())
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## 6. Next Steps

    Once the upload is confirmed, switch to the **AoU Researcher Workbench** and run:

    ```bash
    # Step 2e: cross-reference AoU top SNPs vs catalog, define ±250 kb windows
    nohup bash bash/02e_define_loci_catalog.sh > /dev/null 2>&1 &
    tail -f logs/02e_define_loci_catalog_*.log
    ```

    | Output file | Description |
    |---|---|
    | `target_loci.bed` | Final RVAS loci (BED, GRCh38, 0-based) |
    | `gwas_catalog_validated_snps.tsv` | AoU top SNPs matched to catalog |
    | `novel_snps.tsv` | AoU top SNPs with no catalog support within ±250 kb |
    """)
    return


if __name__ == "__main__":
    app.run()
