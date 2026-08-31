"""Design of the tiling ABE/CBE sgRNA library the MAPK screens use.

This is the step upstream of everything else in this repository: the guides
counted in ``TableS1-ScreenData.xlsx`` are the guides this module produces.

Seven steps, run in order by ``MAPK_Library_Design.ipynb``:

1. Map each target gene symbol to its MANE Select transcript.
2. Design ABE and CBE sgRNAs against those transcripts.
3. Merge the two design tables and filter the result.
4. Summarize what the library can install, per gene.
5. Predict editing efficiency with BEhive.
6. Score off-target specificity with FlashFry.
7. Join the three together into the final library.

Inputs
------
Everything a run reads is in ``TableS6-LibraryDesign.xlsx``: the base-editor
parameters, the target-gene transcript mapping, the two design tables, the
precomputed BEhive and FlashFry scores, and the exon sequences. The raw files
it was built from stay under ``library_design_JW/inputs/`` as provenance, and
``build_inputs`` is the only module that understands their formats.

Heavy steps
-----------
Steps 2, 5 and 6 need resources the rest of the pipeline does not: a live
Ensembl REST endpoint and a 3.9 GB ClinVar dump for step 2, a pinned Python 3.7
environment for step 5, a multi-gigabyte hg38 off-target database for step 6.
Results for all three are in the workbook, so the default run reproduces the
published library in seconds.

Each heavy step is therefore reached through a ``resolve_*`` function that
reports what it can and cannot do rather than raising: it returns a record
naming the source of the results it used and every reason it could not
recompute them. The notebook prints that record. A run always completes.

Nothing here prints; every function returns what the notebook reports.
"""

from __future__ import annotations

import datetime
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import config as cfg

# ---------------------------------------------------------------------------
# STEP 1 - GENE SELECTION AND TRANSCRIPT MAPPING
# ---------------------------------------------------------------------------


def load_transcripts() -> pd.DataFrame:
    """The target genes and their MANE Select transcripts, from the workbook.

    The MANE v1.0 summary is restricted to ``MANE Select`` and joined to the
    target gene list by ``build_inputs.read_mane_transcripts``, which is where
    that file's format is understood.
    """
    return pd.read_excel(cfg.LIBRARY_DESIGN_XLSX,
                         sheet_name=cfg.LIB_TRANSCRIPT_SHEET)


def map_genes_to_transcripts(out_dir: Path) -> tuple[pd.DataFrame, Path]:
    """Write the gene table and the design script's input file.

    Returns ``(transcripts, design_input_path)``: the target genes with their
    transcripts, and the tab-separated transcript/symbol file step 2 reads.
    """
    transcripts = load_transcripts()
    assert transcripts["Ensembl_nuc"].notna().all(), \
        "some genes did not map to a MANE transcript"

    transcripts.to_csv(out_dir / "01_gene_list_with_transcripts.csv", index=False)

    design_input_path = out_dir / "01_genes_for_lib_design.txt"
    (transcripts[["Ensembl_nuc", "Gene"]]
        .rename(columns={"Ensembl_nuc": "Transcript ID", "Gene": "Gene Symbol"})
        .to_csv(design_input_path, sep="\t", index=False))

    return transcripts, design_input_path


# ---------------------------------------------------------------------------
# STEP 2 - SGRNA DESIGN
# ---------------------------------------------------------------------------
#
# base_editing_guide_designs.py tiles every sgRNA whose editing window falls
# within a coding exon of the target transcript, plus an intron buffer, and
# annotates the resulting amino-acid changes. It is run once per base editor,
# with PAM, window, length, edit type and buffer read from BEs_to_use.txt.


def load_be_parameters() -> pd.DataFrame:
    """The base-editor parameter table the design script is driven by."""
    return pd.read_excel(cfg.LIBRARY_DESIGN_XLSX,
                         sheet_name=cfg.LIB_BASE_EDITOR_SHEET)


def supplied_designs() -> dict[str, pd.DataFrame]:
    """The designs used for the published library, one sheet per base editor."""
    return {be: pd.read_excel(cfg.LIBRARY_DESIGN_XLSX, sheet_name=be)
            for be in (cfg.LIB_ABE_NAME, cfg.LIB_CBE_NAME)}


def read_design_file(path: Path) -> pd.DataFrame:
    """A design table straight from the design script's tab-separated output.

    Regenerated designs never reach the workbook, so this is the one place a
    run reads the raw format that ``build_inputs`` otherwise owns.
    """
    return pd.read_table(path)


def find_existing_designs(out_dir: Path) -> dict[str, Path] | None:
    """Designs left in ``out_dir`` by a previous step-2 run, or None.

    Both editors have to be present: a half-regenerated pair would mix designs
    from two Ensembl releases into one library.
    """
    found = {}
    for be in (cfg.LIB_ABE_NAME, cfg.LIB_CBE_NAME):
        hits = list(out_dir.glob(f"MAPK_{be}_*/sgrna_designs_MAPK_{be}.txt"))
        if hits:
            found[be] = max(hits, key=lambda p: p.stat().st_mtime)
    return found if len(found) == 2 else None


def design_preflight(transcript_id: str, server: str = None) -> tuple[list[str], list[str]]:
    """Check both things step 2 needs. Returns ``(blockers, notes)``.

    An empty ``blockers`` list means the design run can go ahead. ``notes``
    carries what was learned on the way - the Ensembl release answering, the
    size of the test sequence fetch - which is worth reporting either way.
    """
    server = server or cfg.LIB_ENSEMBL_SERVER
    problems, notes = [], [f"Ensembl server: {server}"]

    if not cfg.LIB_CLINVAR_SUMMARY.exists():
        problems.append(
            "ClinVar variant_summary.txt not found at %s\n"
            "      cd %s\n"
            "      curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/"
            "variant_summary.txt.gz && gunzip variant_summary.txt.gz"
            % (cfg.LIB_CLINVAR_SUMMARY, cfg.LIB_DESIGN_SCRIPTS_DIR))

    try:
        import requests

        # /info/ping is answered by the front end even when the databases are
        # down, so query /info/data instead: it reports the release and only
        # succeeds if the databases are actually attached.
        response = requests.get(
            f"{server}/info/data?content-type=application/json", timeout=30)
        response.raise_for_status()
        releases = response.json().get("releases", [])
        notes.append(f"Ensembl release: {', '.join(str(x) for x in releases) or 'unknown'}")

        # Exercise the endpoint the design script actually depends on.
        url = (f"{server}/sequence/id/{transcript_id}"
               "?content-type=text/plain&expand_5prime=40&expand_3prime=40")
        response = requests.get(url, headers={"Content-Type": "text/plain"}, timeout=60)
        response.raise_for_status()
        notes.append(f"sequence fetch OK: {len(response.text):,} nt for {transcript_id}")
    except ImportError:
        problems.append("the `requests` package is not installed")
    except Exception as exc:
        problems.append(
            f"Ensembl fetch failed ({type(exc).__name__}: {exc})\n"
            f"      Try another endpoint by changing ENSEMBL_SERVER, or check\n"
            f"      https://status.ensembl.org")

    return problems, notes


def run_design(be_row: pd.Series, design_input_path: Path, out_dir: Path,
               server: str = None) -> Path:
    """Run the design script for one base editor; return the design file.

    Called with ``cwd=out_dir`` and a bare ``--output-name``, because the
    script embeds ``output_name`` in a filename as well as in the directory
    name it creates.
    """
    server = server or cfg.LIB_ENSEMBL_SERVER
    output_name = f"MAPK_{be_row['BEs']}"
    cmd = [
        sys.executable, str(cfg.LIB_DESIGN_SCRIPT),
        "--input-file", str(design_input_path),
        "--input-type", "tid",
        "--variant-file", str(cfg.LIB_CLINVAR_SUMMARY),
        "--pam", str(be_row["PAM"]),
        "--edit-window", str(be_row["Window"]),
        "--sg-len", str(int(be_row["Sgrna length"])),
        "--edit", str(be_row["Edit"]),
        "--intron-buffer", str(int(be_row["Intron buffer"])),
        "--filter-gc", str(be_row["Filter GC"]),
        "--output-name", output_name,
    ]
    # The design script reads its REST endpoint from the environment.
    env = dict(os.environ, ENSEMBL_REST_SERVER=server)
    subprocess.run(cmd, check=True, cwd=str(out_dir), env=env)

    # The script appends a run timestamp to the directory name.
    made = list(out_dir.glob(f"{output_name}_*/sgrna_designs_{output_name}.txt"))
    if not made:
        raise RuntimeError(f"{be_row['BEs']}: design script produced no sgrna_designs file")
    return max(made, key=lambda p: p.stat().st_mtime)


def resolve_designs(transcripts: pd.DataFrame, design_input_path: Path,
                    out_dir: Path, run_heavy: bool = False,
                    reuse_existing: bool = True, server: str = None) -> dict:
    """Decide which sgRNA designs this run uses, regenerating them if asked.

    Returns a record with ``source`` ('supplied', 'reused' or 'regenerated'),
    ``designs`` (base-editor name to design table), ``origins`` (where each
    table was read from, for reporting), ``regenerated`` (whether the designs
    came from a live Ensembl rather than the workbook, which is what tells the
    verification step that differences are expected), ``blockers`` and
    ``notes``.

    A failure inside the hour-long design run is caught and reported rather
    than raised, so it does not cost the rest of the pipeline.
    """
    server = server or cfg.LIB_ENSEMBL_SERVER
    workbook = cfg.LIBRARY_DESIGN_XLSX.name
    record = {"source": "supplied", "designs": None, "regenerated": False,
              "origins": {}, "blockers": [], "notes": []}

    def use_supplied():
        record["designs"] = supplied_designs()
        record["origins"] = {be: f"{workbook} [{be}]" for be in record["designs"]}
        return record

    if not run_heavy:
        record["notes"].append(
            f"RUN_HEAVY_STEPS is False - using the designs in {workbook}")
        return use_supplied()

    existing = find_existing_designs(out_dir) if reuse_existing else None
    if existing:
        record.update(source="reused", regenerated=True,
                      designs={be: read_design_file(path)
                               for be, path in existing.items()},
                      origins={be: str(path) for be, path in existing.items()})
        record["notes"].append(
            "Reusing regenerated designs already in the output directory. Set "
            "REUSE_EXISTING_DESIGNS = False (or delete those directories) to "
            "force a fresh run.")
        return record

    blockers, notes = design_preflight(transcripts["Ensembl_nuc"].iloc[0], server)
    record["notes"] += notes
    if blockers:
        record["blockers"] = blockers
        record["notes"].append(f"Using the designs in {workbook}")
        return use_supplied()

    try:
        produced = {row["BEs"]: run_design(row, design_input_path, out_dir, server)
                    for _, row in load_be_parameters().iterrows()}
        releases = _ensembl_releases(server)
        (out_dir / "02_design_provenance.json").write_text(json.dumps({
            "ensembl_server": server,
            "ensembl_release": releases,
            "generated": datetime.datetime.now().isoformat(timespec="seconds"),
            "designs": {be: str(p) for be, p in produced.items()},
        }, indent=2))
        record.update(source="regenerated", regenerated=True,
                      designs={be: read_design_file(path)
                               for be, path in produced.items()},
                      origins={be: str(path) for be, path in produced.items()})
        record["notes"].append(
            f"Ensembl release {releases} recorded in 02_design_provenance.json")
        return record
    except Exception as exc:
        record["blockers"].append(
            f"design regeneration FAILED: {type(exc).__name__}: {exc}")
        record["notes"].append(f"Using the designs in {workbook}. "
                               "Re-run this cell to retry.")
        return use_supplied()


def _ensembl_releases(server: str) -> list:
    """Releases served by an Ensembl endpoint, for the provenance record."""
    import requests
    return requests.get(f"{server}/info/data?content-type=application/json",
                        timeout=30).json().get("releases", [])


# ---------------------------------------------------------------------------
# STEP 3 - MERGE AND FILTER
# ---------------------------------------------------------------------------
#
# Both editors are designed against the same NG PAM guide space, so most
# sgRNAs carry both an A->G and a C->T annotation. Merging on sequence gives
# one row per guide with both.

# Columns taken from each design table, and what they become after the merge.
_ABE_COLUMNS = {
    "sgRNA sequence": "sgRNA sequence",
    "Gene Symbol": "Gene Symbol",
    "Ensembl Gene ID": "Ensembl Gene ID",
    "Ensembl transcript ID": "Ensembl transcript ID",
    "sgRNA Strand": "sgRNA Strand",
    "PAM": "PAM",
    "# edits": "ABE # edits",
    "Amino acid edits": "ABE Amino acid edits",
    "Mutation category": "ABE Mutation category",
    "Clinical significance": "ABE Clinical significance",
    "4T flag": "4T flag",
}
_CBE_COLUMNS = {
    "sgRNA sequence": "sgRNA sequence",
    "# edits": "CBE # edits",
    "Amino acid edits": "CBE Amino acid edits",
    "Mutation category": "CBE Mutation category",
    "Clinical significance": "CBE Clinical significance",
}


def _select(designs: pd.DataFrame, columns: dict[str, str]) -> pd.DataFrame:
    """One design table, restricted and renamed to its merged column names."""
    frame = designs[list(columns)].copy()
    frame.columns = list(columns.values())
    return frame


def merge_designs(abe_designs: pd.DataFrame, cbe_designs: pd.DataFrame,
                  out_dir: Path) -> tuple[pd.DataFrame, list[str], Path]:
    """Merge the ABE and CBE design tables into one row per guide.

    Returns the merged table, the gene order taken from it, and the path the
    merge was written to. That order is captured here, before filtering,
    because guide numbering follows it.

    The unfiltered merge is written out as well, and step 5 reads it back:
    BEhive scores every designed guide, not only the ones that survive step 3.
    """
    merged = _select(abe_designs, _ABE_COLUMNS).merge(
        _select(cbe_designs, _CBE_COLUMNS), "left", on="sgRNA sequence")

    gene_list = list(merged["Gene Symbol"].drop_duplicates())

    merged["ABE Amino acid edits"] = merged["ABE Amino acid edits"].astype(str)
    merged["CBE Amino acid edits"] = merged["CBE Amino acid edits"].astype(str)

    unfiltered_path = out_dir / "03_MAPK_lib_unfiltered.csv"
    merged.to_csv(unfiltered_path)
    return merged, gene_list, unfiltered_path


def filter_library(merged: pd.DataFrame,
                   out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drop the guides the library cannot use. Returns ``(kept, counts)``.

    Removed are guides that edit the start codon, create or remove a stop
    codon, edit only intronic or only UTR sequence, or occur more than once in
    the library. Every copy of a non-unique sequence goes: such a guide cannot
    be attributed to a single gene. The duplicates are written out first.
    """
    kept = merged
    rows = [{"filter": "designed guides", "dropped": 0, "remaining": len(kept)}]

    for pattern, reason in cfg.LIB_FILTERS:
        keep_mask = ~(kept["ABE Amino acid edits"].str.contains(pattern)
                      | kept["CBE Amino acid edits"].str.contains(pattern))
        rows.append({"filter": reason, "dropped": int((~keep_mask).sum()),
                     "remaining": int(keep_mask.sum())})
        kept = kept[keep_mask]

    kept[kept.duplicated("sgRNA sequence")].to_csv(out_dir / "03_duplicate_sgRNAs.csv")
    before = len(kept)
    kept = kept.drop_duplicates(subset="sgRNA sequence", keep=False).copy()
    rows.append({"filter": "non-unique sequences", "dropped": before - len(kept),
                 "remaining": len(kept)})

    return kept, pd.DataFrame(rows)


def number_guides(kept: pd.DataFrame, gene_list: list[str],
                  out_dir: Path) -> tuple[pd.DataFrame, Path]:
    """Add ``<GENE>_<n>`` guide IDs and write the filtered library.

    Guides are contiguous by gene, so numbering follows the gene order
    captured before filtering.
    """
    guide_ids = []
    for gene in gene_list:
        gene_guides = kept.loc[kept["Gene Symbol"] == gene, "sgRNA sequence"]
        guide_ids += [f"{gene}_{i + 1}" for i in range(len(gene_guides))]

    assert len(guide_ids) == len(kept)
    numbered = kept.copy()
    numbered.insert(0, "guide_id", guide_ids)

    assert (numbered["guide_id"].str.rsplit("_", n=1).str[0]
            == numbered["Gene Symbol"]).all(), "guide_id / gene mismatch"

    library_path = out_dir / "03_MAPK_lib_filtered.csv"
    numbered.to_csv(library_path)
    return numbered, library_path


# ---------------------------------------------------------------------------
# STEP 4 - PER-GENE STATISTICS
# ---------------------------------------------------------------------------


def _unique_missense(edit_series) -> set[str]:
    """Unique missense mutations in a set of ';'-separated edit annotations."""
    muts = []
    for entry in edit_series:
        for mut in str(entry).split(";"):
            mut = mut.strip()
            if not mut or mut == "flankseq":
                continue
            if mut[:3] == mut[-3:]:      # e.g. Ala12Ala -> silent
                continue
            muts.append(mut)
    return set(muts)


def _residues(muts: set[str]) -> set[int]:
    """Residue numbers touched by a set of mutations (Ala12Val -> 12)."""
    return {int(m[3:-3]) for m in muts}


def gene_statistics(library: pd.DataFrame, gene_list: list[str],
                    transcripts: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    """Unique missense mutations and residues reachable, per gene.

    Reported for the CBE and the ABE separately and for the two combined.
    Silent substitutions and non-coding annotations are excluded.
    """
    rows = []
    for gene in gene_list:
        gene_df = library[library["Gene Symbol"] == gene]
        cbe = _unique_missense(gene_df["CBE Amino acid edits"])
        abe = _unique_missense(gene_df["ABE Amino acid edits"])
        cbe_res, abe_res = _residues(cbe), _residues(abe)
        rows.append([gene, len(gene_df),
                     len(cbe), len(cbe_res),
                     len(abe), len(abe_res),
                     len(cbe | abe), len(cbe_res | abe_res)])

    stats = pd.DataFrame(rows, columns=[
        "Gene", "Total guides",
        "CBE: unique missense mutations", "CBE: residues mutated",
        "ABE: unique missense mutations", "ABE: residues mutated",
        "Combined CBE/ABE: unique missense mutations",
        "Combined CBE/ABE: residues mutated"])

    stats = (transcripts[["Gene", "name", "Ensembl_nuc"]]
             .merge(stats, "right", on="Gene")
             .sort_values("Total guides", ascending=False))

    stats.to_csv(out_dir / "04_MAPK_lib_stats.csv", index=False)
    return stats


# ---------------------------------------------------------------------------
# STEP 5 - BEHIVE PREDICTED EDITING EFFICIENCY
# ---------------------------------------------------------------------------
#
# BEhive (https://github.com/maxwshen/be_predict_efficiency) predicts editing
# efficiency from the 50 bp of genomic context centred on the protospacer.
# Scores are computed for every designed guide, before filtering, so no guide
# reaches step 7 without one.

# Where setup_behive_env.sh clones it, then the original location.
BEHIVE_SEARCH_DIRS = [cfg.LIB_BEHIVE_DIR, Path.home()]


def pull_out_50bp_region(exon_list, guide_seq: str) -> str:
    """The 50 bp context around ``guide_seq``, searched across ``exon_list``.

    That is 20 nt upstream, the 20 nt guide, and 10 nt downstream. A hit too
    close to an exon edge to yield the full window is skipped and the search
    continues, since the same guide may appear again with more flanking
    sequence.
    """
    guide_seq = guide_seq.upper()
    upstream, downstream = cfg.LIB_BEHIVE_CONTEXT

    for exon_seq in exon_list:
        exon_seq = exon_seq.upper()
        search_pos = 0
        while True:
            guide_start = exon_seq.find(guide_seq, search_pos)
            if guide_start == -1:
                break                                   # no more hits in this exon
            context_start = guide_start - upstream
            context_end = guide_start + downstream
            if context_start >= 0 and context_end <= len(exon_seq):
                return exon_seq[context_start:context_end]
            search_pos = guide_start + 1                # hit too close to the edge

    return "None"


def load_exon_records() -> pd.DataFrame:
    """The exon sequences BEhive reads a guide's genomic context from."""
    return pd.read_excel(cfg.LIBRARY_DESIGN_XLSX, sheet_name=cfg.LIB_EXON_SHEET)


def load_exons_by_gene(exon_records: pd.DataFrame,
                       refseq_to_symbol: dict) -> dict:
    """Group exon records into ``{gene_symbol: [exon, exon_revcomp, ...]}``.

    Record IDs look like ``hg38_ncbiRefSeq_NM_001654.5_3``; the RefSeq
    accession is recovered from that and translated through
    ``refseq_to_symbol``. Both strands are kept, because a guide may match
    either.
    """
    pattern = re.compile(r"^hg38_ncbiRefSeq_(.+?)_[^_]+$")

    exons_by_refseq: dict[str, list[str]] = {}
    for record_id, sequence in zip(exon_records["Record"],
                                   exon_records["Sequence"]):
        match = pattern.match(str(record_id))
        if not match:
            continue
        sequence = Seq(str(sequence))
        exons_by_refseq.setdefault(match.group(1), []).extend(
            [str(sequence), str(sequence.reverse_complement())])

    exons_by_gene = {refseq_to_symbol[r]: exons
                     for r, exons in exons_by_refseq.items() if r in refseq_to_symbol}

    # Fail loudly rather than silently returning empty contexts for everything.
    if not exons_by_gene:
        raise KeyError(
            "No FASTA record matched the RefSeq->symbol mapping.\n"
            "  FASTA keys look like: %s\n"
            "  mapping keys look like: %s\n"
            "Version suffixes must match on both sides."
            % (sorted(exons_by_refseq)[:3], sorted(refseq_to_symbol)[:3]))

    return exons_by_gene


def find_behive(behive_parent=None) -> Path | None:
    """The directory *containing* ``be_predict_efficiency/``, or None."""
    candidates = [Path(behive_parent)] if behive_parent else BEHIVE_SEARCH_DIRS
    for parent in candidates:
        if (parent / "be_predict_efficiency").is_dir():
            return parent
    return None


def behive_preflight(behive_parent=None) -> list[str]:
    """Reasons BEhive cannot run in this kernel. Empty list means it can."""
    problems = []

    try:
        import sklearn
        if sklearn.__version__ != cfg.LIB_BEHIVE_SKLEARN_VERSION:
            problems.append(
                "scikit-learn is %s, but BEhive's pickled models need %s "
                "(they reference sklearn.ensemble.gradient_boosting, removed in 0.22)"
                % (sklearn.__version__, cfg.LIB_BEHIVE_SKLEARN_VERSION))
    except ImportError:
        problems.append("scikit-learn is not installed in this kernel")

    if find_behive(behive_parent) is None:
        searched = [Path(behive_parent)] if behive_parent else BEHIVE_SEARCH_DIRS
        problems.append("no be_predict_efficiency/ found in: %s"
                        % ", ".join(str(c) for c in searched))

    return problems


def import_behive(behive_parent=None):
    """Import BEhive, adding its parent directory to ``sys.path``.

    The package is imported by name, so the directory placed on the path is
    the one *containing* ``be_predict_efficiency/``, not the clone itself.
    """
    parent = find_behive(behive_parent)
    if parent is None:
        raise ImportError("be_predict_efficiency not found; run "
                          "inputs/behive/setup_behive_env.sh")
    if str(parent) not in sys.path:
        sys.path.insert(0, str(parent))
    from be_predict_efficiency import predict as be_efficiency_model
    return be_efficiency_model


def calc_behive_scores(guides: pd.DataFrame, exon_records: pd.DataFrame,
                       refseq_to_symbol: dict,
                       behive_parent=None) -> tuple[pd.DataFrame, list[str]]:
    """Score every guide with BEhive. Returns ``(scored, warnings)``.

    Adds ``context_50bp`` and one column per model. Guides whose context
    cannot be recovered are scored NaN rather than aborting the run: BEhive
    asserts on any input that is not exactly 50 nt.
    """
    be_efficiency_model = import_behive(behive_parent)
    exons_by_gene = load_exons_by_gene(exon_records, refseq_to_symbol)
    warnings = []

    missing_genes = sorted(set(guides["Gene Symbol"]) - set(exons_by_gene))
    if missing_genes:
        warnings.append(f"no exons in the FASTA for: {', '.join(missing_genes)}")

    scored = guides.copy()
    scored["context_50bp"] = [
        pull_out_50bp_region(exons_by_gene.get(gene, []), sgrna)
        for sgrna, gene in zip(scored["sgRNA sequence"], scored["Gene Symbol"])]

    scorable = scored["context_50bp"].str.len() == 50
    if (~scorable).any():
        warnings.append(f"{(~scorable).sum()} guides have no usable 50bp context "
                        f"and are scored NaN")

    for model, column in cfg.LIB_BEHIVE_MODELS:
        be_efficiency_model.init_model(base_editor=model,
                                       celltype=cfg.LIB_BEHIVE_CELLTYPE)
        scored[column] = [
            be_efficiency_model.predict(x)["Predicted logit score"] if ok else np.nan
            for x, ok in zip(scored["context_50bp"], scorable)]

    return scored, warnings


def cross_check_behive(scored: pd.DataFrame) -> pd.DataFrame:
    """Compare freshly computed BEhive scores against the published ones.

    Joined on sequence *and* gene, not sequence alone: 93 guides occur in more
    than one gene with a different genomic context - and therefore a different
    score - in each. Deduplicating by sequence picks one arbitrarily and makes
    the comparison report spurious differences.
    """
    reference = load_behive_sheet()
    key = ["sgRNA sequence", "Gene Symbol"]
    score_columns = [column for _, column in cfg.LIB_BEHIVE_MODELS]

    check = scored.merge(reference[key + ["context_50bp"] + score_columns],
                         "inner", on=key, suffixes=("", "_published"))
    if check.empty:
        return pd.DataFrame(columns=["check", "n_rows", "max_abs_difference", "agrees"])

    rows = [{"check": "50bp context", "n_rows": len(check), "max_abs_difference": np.nan,
             "agrees": bool((check["context_50bp"]
                             == check["context_50bp_published"]).all())}]
    for column in score_columns:
        delta = (check[column] - check[f"{column}_published"]).abs().max()
        rows.append({"check": column, "n_rows": len(check),
                     "max_abs_difference": delta, "agrees": bool(delta < 1e-6)})
    return pd.DataFrame(rows)


def load_behive_sheet() -> pd.DataFrame:
    """The precomputed BEhive scores, one row per designed guide."""
    return pd.read_excel(cfg.LIBRARY_DESIGN_XLSX, sheet_name=cfg.LIB_BEHIVE_SHEET)


def resolve_behive_scores(unfiltered_path: Path, transcripts: pd.DataFrame,
                          out_dir: Path, run_heavy: bool = False,
                          behive_parent=None) -> dict:
    """Decide which BEhive scores this run uses, recomputing them if asked.

    Returns a record with ``source`` ('supplied' or 'computed'), ``scores``
    (one row per unique sgRNA, for the step-7 join), ``blockers``, ``notes``
    and, when recomputed, a ``cross_check`` frame against the published
    scores.
    """
    record = {"source": "supplied", "blockers": [], "notes": [],
              "cross_check": None}

    blockers = behive_preflight(behive_parent) if run_heavy else []
    if not run_heavy:
        record["notes"].append("RUN_HEAVY_STEPS is False - using the BEhive "
                               f"scores in {cfg.LIBRARY_DESIGN_XLSX.name}")
    elif blockers:
        record["blockers"] = blockers
        record["notes"] += [
            f"current kernel: {sys.executable}",
            f"create the environment: bash {cfg.LIB_BEHIVE_DIR / 'setup_behive_env.sh'}",
            "then in Jupyter: Kernel > Change kernel > Python (BEhive_env)",
            f"Using the BEhive scores in {cfg.LIBRARY_DESIGN_XLSX.name}."]

    if run_heavy and not blockers:
        # Keep the RefSeq version suffix. UCSC record IDs are of the form
        # hg38_ncbiRefSeq_NM_001654.5_3, and load_exons_by_gene() parses out
        # "NM_001654.5" - so the dict keys have to carry the version too.
        refseq_to_symbol = dict(zip(transcripts["RefSeq_nuc"], transcripts["Gene"]))
        scored, warnings = calc_behive_scores(
            pd.read_csv(unfiltered_path, index_col=0), load_exon_records(),
            refseq_to_symbol, behive_parent)
        scored.to_csv(out_dir / "05_MAPK_lib_BEhive_scores.csv")
        record.update(source="computed", cross_check=cross_check_behive(scored))
        record["notes"] += warnings
    else:
        scored = load_behive_sheet()

    score_columns = [column for _, column in cfg.LIB_BEHIVE_MODELS]
    scores = scored[["sgRNA sequence"] + score_columns].drop_duplicates(
        subset="sgRNA sequence")
    scores.to_csv(out_dir / "05_BEhive_scores_by_sgRNA.csv", index=False)

    record["scored"] = scored
    record["scores"] = scores
    return record


# ---------------------------------------------------------------------------
# STEP 6 - FLASHFRY OFF-TARGET SPECIFICITY
# ---------------------------------------------------------------------------
#
# FlashFry accepts only NGG guides and discards anything else, which would
# exclude most of this library, since SpG guides require only NG. A G is
# therefore appended to each guide's NG PAM to form a 23 nt sequence FlashFry
# will accept. Only the 20 nt protospacer determines the off-target search;
# the appended base is stripped again when the scores are merged back.


def df_to_fasta(frame: pd.DataFrame, id_col: str, seq_col: str,
                output_file: Path) -> None:
    """Write two columns of a dataframe out as a FASTA file."""
    records = [SeqRecord(Seq(row[seq_col]), id=str(row[id_col]), description="_")
               for _, row in frame.iterrows()]
    with open(output_file, "w") as handle:
        SeqIO.write(records, handle, "fasta")


def write_flashfry_fasta(library_path: Path, out_dir: Path) -> tuple[pd.DataFrame, Path]:
    """Write the padded 23 nt query FASTA FlashFry is run against."""
    library = pd.read_csv(library_path, index_col=0)
    library["sgRNA sequence w NGG"] = library["sgRNA sequence"] + library["PAM"] + "G"

    fasta_path = out_dir / "06_MAPK_sgRNAs_w_NGGs.fasta"
    df_to_fasta(library, id_col="guide_id", seq_col="sgRNA sequence w NGG",
                output_file=fasta_path)
    return library, fasta_path


def load_flashfry_sheet() -> pd.DataFrame:
    """The precomputed FlashFry off-target scores, one row per scored guide."""
    return pd.read_excel(cfg.LIBRARY_DESIGN_XLSX, sheet_name=cfg.LIB_FLASHFRY_SHEET)


def resolve_specificity_scores(out_dir: Path, scored_path: Path = None,
                               run_heavy: bool = False) -> dict:
    """Load FlashFry specificity scores. Returns a record like the others.

    FlashFry needs a multi-gigabyte hg38 NGG database, so this step is never
    run inline and ``RUN_HEAVY_STEPS`` does not change its behaviour. Pass
    ``scored_path`` to use your own ``*_gRNA_scored.txt``; otherwise the scores
    come from the workbook.
    """
    workbook = cfg.LIBRARY_DESIGN_XLSX.name
    record = {"source": "provided" if scored_path else "supplied",
              "notes": [], "blockers": []}

    if run_heavy:
        record["blockers"].append(
            "FlashFry cannot run inline; it needs an hg38 NGG database. To "
            "regenerate, upload the query FASTA to the cluster, run %s, then "
            "pass the resulting *_gRNA_scored.txt as FLASHFRY_SCORED."
            % " and ".join(str(cfg.LIB_FLASHFRY_DIR / s)
                           for s in cfg.LIB_FLASHFRY_SCRIPTS))
        record["notes"].append(f"Using the scores in {workbook}.")
    elif not scored_path:
        record["notes"].append("RUN_HEAVY_STEPS is False - using the FlashFry "
                               f"scores in {workbook}")

    scored = pd.read_table(scored_path) if scored_path else load_flashfry_sheet()
    # Strip the 3 nt PAM back off to recover the original 20 nt guide.
    scored["sgRNA sequence"] = [x[:-3] for x in scored["target"]]

    scores = scored[["sgRNA sequence", cfg.LIB_SPECIFICITY_COLUMN]].drop_duplicates(
        subset="sgRNA sequence")
    scores.to_csv(out_dir / "06_specificity_scores_by_sgRNA.csv", index=False)

    record["scored"] = scored
    record["scores"] = scores
    return record


# ---------------------------------------------------------------------------
# STEP 7 - FINAL LIBRARY
# ---------------------------------------------------------------------------


def assemble_library(library_path: Path, behive_scores: pd.DataFrame,
                     specificity_scores: pd.DataFrame,
                     out_dir: Path) -> tuple[pd.DataFrame, Path]:
    """The filtered library joined to its efficiency and specificity scores.

    ``07_MAPK_lib_final.csv`` is the deliverable: the library that was ordered
    and screened.
    """
    final = (pd.read_csv(library_path, index_col=0)
             .merge(behive_scores, "left", on="sgRNA sequence")
             .merge(specificity_scores, "left", on="sgRNA sequence"))

    final_path = out_dir / "07_MAPK_lib_final.csv"
    final.to_csv(final_path, index=False)
    return final, final_path


# ---------------------------------------------------------------------------
# VERIFICATION
# ---------------------------------------------------------------------------


def structural_checks(final: pd.DataFrame) -> pd.DataFrame:
    """Invariants that must hold no matter where the designs came from."""
    checks = [
        ("unique sgRNA sequences", not final["sgRNA sequence"].duplicated().any()),
        ("unique guide_ids", not final["guide_id"].duplicated().any()),
        ("every guide is 20 nt", bool(final["sgRNA sequence"].str.len().eq(20).all())),
        ("guide_id agrees with gene",
         bool((final["guide_id"].str.rsplit("_", n=1).str[0]
               == final["Gene Symbol"]).all())),
    ]
    return pd.DataFrame(checks, columns=["check", "passed"])


def published_comparison(final: pd.DataFrame, stats: pd.DataFrame) -> pd.DataFrame:
    """This run's counts beside the published library's."""
    score_columns = [column for _, column in cfg.LIB_BEHIVE_MODELS]
    observed = {
        "guides in final library": len(final),
        "genes targeted": final["Gene Symbol"].nunique(),
        "rows in per-gene stats": len(stats),
        "guides with a BEhive CBE score": int(final[score_columns[0]].notna().sum()),
        "guides with a BEhive ABE score": int(final[score_columns[1]].notna().sum()),
        "guides with a Hsu2013 score":
            int(final[cfg.LIB_SPECIFICITY_COLUMN].notna().sum()),
    }
    return pd.DataFrame(
        [{"check": label, "this_run": got,
          "published": cfg.LIB_PUBLISHED_COUNTS[label],
          "matches": got == cfg.LIB_PUBLISHED_COUNTS[label]}
         for label, got in observed.items()])
