#!/usr/bin/env python3
"""
blast_remote.py
Run remote blastn / megablast, save TSV, then convert to XLSX.

(python >=3.8, pandas & openpyxl required for the Excel step)
Author: Kimhun Tuntikawinwong
"""
import argparse, os, random, subprocess, time
from pathlib import Path
import datetime
import pandas as pd     

# ---------------------------------------------------------------------
def parse_args():
    """
    Parse the command-line interface for *blast_remote.py*.

    Returns
    -------
    argparse.Namespace
        The populated namespace whose attributes correspond to the
        command-line options (e.g. ``namespace.query``, ``namespace.db`` …).
    """
    p = argparse.ArgumentParser(
        description="Remote BLAST → TSV → XLSX (short-sequence wrapper)")
    p.add_argument("query", help="FASTA file containing ONE query sequence")
    p.add_argument("-d", "--db", default="core_nt",
                   choices=["nt", "core_nt"],
                   help="BLAST database (default: core_nt)")
    p.add_argument("-o", "--outdir", default=".",
                   help="output folder (default: current dir)")
    p.add_argument("-u", "--outfile", default="", metavar="NAME.tsv",
                   help="TSV file name; default = auto-generated")
    p.add_argument("-t", "--task", default="blastn",
                   choices=["blastn", "megablast", "dc-megablast", "blastn-short"],
                   help="BLAST task (default: blastn)")
    p.add_argument("-f", "--filter", default=None,
                   help="filter (default: none)")
    p.add_argument("-n", "--filter_name", default=None,
                   help="filter name(default: none)")
    p.add_argument("-s","--sleep", nargs=2, type=int, default=[11, 15],
                   metavar=("MIN", "MAX"),
                   help="random pause after job (default: 11-15 s)")
    return p.parse_args()

# ---------------------------------------------------------------------
def auto_name(query: Path, outdir: Path, filter_name: str) -> Path:
    """
    Build a timestamped default output filename.

    The name is constructed as  
    ``<query stem>_vs_<filter_name>_<YYYYMMDD_HHMMSS>.tsv`` when
    *filter_name* is supplied, otherwise  
    ``<query stem>_<YYYYMMDD_HHMMSS>.tsv``.

    Parameters
    ----------
    query : pathlib.Path
        The original FASTA file provided by the user.
    outdir : pathlib.Path
        Destination directory in which the file will be created.
    filter_name : str
        Readable label for the Entrez filter (e.g. ``"Mycobacteriaceae"``).

    Returns
    -------
    pathlib.Path
        Full path of the output *.tsv* file.
    """
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    base = query.stem.split(".")[0]
    if filter_name:
        return outdir / f"{base}_vs_{filter_name}_{ts}.tsv"
    else:
        return outdir / f"{base}_{ts}.tsv"
# ---------------------------------------------------------------------
def convert_tsv_to_xlsx(tsv: Path):
    """
    Load a BLAST *outfmt 6* TSV, add column headers, and save as Excel.

    Parameters
    ----------
    tsv : pathlib.Path
        Path to the BLAST results in *outfmt 6* (tab-delimited).

    Returns
    -------
    pathlib.Path
        Path to the newly-created *.xlsx* file (same stem as *tsv*).

    Notes
    -----
    The column list (`cols`) must match the fields requested in the
    ``-outfmt`` string, **in the same order**; otherwise Pandas will
    mis-align the data.
    """
    cols = [
        "query_seq_id","subject_seq_id","query_coverage","percent_identity","length","mismatch","gapopen",
        "query_start","query_end","subject_start","send","evalue","bitscore","score",
        "qcovhsp","query_length","subject_length(Acc. Len)","staxids","sscinames","sskingdoms"
    ]
    df = pd.read_csv(tsv, sep="\t", header=None, names=cols)
    xlsx = tsv.with_suffix(".xlsx")
    df.to_excel(xlsx, index=False)
    print(f"Create excel: {xlsx.name}")
    return xlsx
# -------------------------------------------------------
def warn_if_short(query_fa: Path,
                  requested_task: str,
                  threshold: int = 50) -> None:
    """
    Emit a polite hint when a very short query is run with the wrong task.

    If the first sequence in *query_fa* is ≤ *threshold* bp and the user
    did **not** request ``-task blastn-short``, a warning message is printed
    to stderr suggesting that task.

    Parameters
    ----------
    query_fa : pathlib.Path
        FASTA file provided on the command line.
    requested_task : str
        Value of the ``--task`` option the user supplied.
    threshold : int, optional
        Length (in bp) below which the message should be shown
        (default = 50 bp).

    Returns
    -------
    None
        The function only prints a message; it does not modify state.

    Examples
    --------
    >>> warn_if_short(Path("probe.fa"), "blastn")
    >  Query length = 35 bp; for sequences ≤ 50 bp [...]
    """
    if requested_task == "blastn-short":
        return                                     # user already chose it

    # count only the first sequence (fast)
    length = 0
    with query_fa.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if length:                         # we’ve finished 1st record
                    break
                continue
            length += len(line.strip())
            if length > threshold:                 # long enough, stop early
                break

    if length <= threshold:
        print(
            f">  Query length = {length} bp; "
            "for sequences ≤ {threshold} bp NCBI recommends `-task blastn-short`.\n"
            f">   You requested `-task {requested_task}`. "
            "If you want optimal sensitivity for very short probes/primers, "
            "consider rerunning with `-t blastn-short`."
        )
# -------------------------------------------------------
def main():
    """
    Orchestrate the batch run:

    1. Parse CLI arguments.  
    2. Discover FASTA files (with optional keyword whitelist).  
    3. Split any multi-record FASTA.  
    4. Invoke *blast_remote.py* for each sequence not yet processed.  
    5. Log progress and sleep a random interval between jobs.

    Notes
    -----
    * Logging goes both to the console **and** to a timestamped text file
      in *out_dir*.
    * The wrapper passes ``--sleep 0 0`` to *blast_remote.py* because the
      outer script already handles the polite back-off.
    """
    a = parse_args()

    query  = Path(a.query).expanduser().resolve()
    if not query.is_file():
        raise FileNotFoundError(query)
    warn_if_short(query, a.task)
    outdir = Path(a.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    tsv_path = outdir / a.outfile if a.outfile else auto_name(query, outdir, a.filter_name)

    # Where BLAST looks for taxdb.* if you have a local copy (not nessecary for remote BLAST)
    os.environ.setdefault("BLASTDB", str(Path.home() / "blastdb"))

    cmd = [
        "blastn",
        "-task", a.task,
        "-query", str(query),
        "-db", a.db,
        "-remote",
        "-max_target_seqs", "200",
        "-outfmt", 
        (
        "6 qseqid saccver "      # query / subject IDs
        "qcovs "                # Query Cover (sum of all HSPs, %)
        "pident length mismatch gapopen "
        "qstart qend sstart send "
        "evalue "
        "bitscore "             # Max Score  (bit-score of the best HSP)
        "score "                # Total Score (raw score of the best HSP)
        "qcovhsp "              # coverage of this *HSP* (%), optional
        "qlen slen "            # Acc. Len  ≈ subject length (and query length)
        "staxids sscinames sskingdoms"),
        "-out", str(tsv_path),
    ]

    if a.filter or a.filter == "":                      
        cmd.extend(["-entrez_query", a.filter])

    print("> Running BLAST ⇒\n  " + "   ".join(cmd) + "\n")
    t0 = time.perf_counter()
    print(f"> Waiting for BLAST result ...")
    res = subprocess.run(cmd, capture_output=True, text=True)
    t1 = time.perf_counter()                    
    elapsed_sec = t1 - t0                        # float, seconds
    elapsed_hms = str(datetime.timedelta(seconds=int(elapsed_sec)))
    print(f"> Receive BLAST result! done in {elapsed_hms}")
    print(f"> Query time: {elapsed_hms}")


    if res.returncode != 0:
        print(res.stderr)
        raise SystemExit(f"BLAST failed (exit {res.returncode})")

    # ── sanity-check the output file ────────────────────────────────
    if not tsv_path.exists() or tsv_path.stat().st_size == 0:
        raise SystemExit(
            "BLAST finished but produced no results — "
            "check your query and/or filter string."
        )

    print(f"> TSV written → {tsv_path.name}  ({tsv_path.stat().st_size} bytes)")

    #convert to excel
    convert_tsv_to_xlsx(tsv_path)

    # polite back-off
    sleep_for = random.randint(*a.sleep)
    print(f"Sleeping {sleep_for}s …")
    time.sleep(sleep_for)

# ---------------------------------------------------------------------
if __name__ == "__main__":
    main()
