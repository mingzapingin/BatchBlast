#!/usr/bin/env python3
"""
run_batch_blast.py
Call blast_remote.py sequentially for many query FASTA files.

> skips files already processed (xlsx present)
> logs every step
> 11-15 s random pause between jobs

Author: Kimhun Tuntikawinwong
"""

from pathlib import Path
from textwrap import wrap
import subprocess, sys, time, datetime, random, logging, argparse

# ─── command-line options ──────────────────────────────────────────────
def get_args():
    """
    Parse all CLI options for *run_batch_blast.py*.

    Returns
    -------
    argparse.Namespace
        Populated namespace with attributes:

        * **input_dir** – folder containing FASTA files  
        * **out_dir** – destination for *.tsv / .xlsx*  
        * **blast_script** – path to *blast_remote.py*  
        * **include** – list[str] keywords for filename whitelist  
        * **filter** – Entrez query string (or *None*)  
        * **filter_name** – short human-readable label for the filter  
        * **sleep** – two-element list [min, max] pause in seconds
    """
    ap = argparse.ArgumentParser(
        description="Batch wrapper around blast_remote.py")
    ap.add_argument("input_dir",  help="folder with *.fa / *.fna / *.fasta")
    ap.add_argument("out_dir",    help="folder to write .tsv / .xlsx results")
    ap.add_argument("--blast_script", default="blast_remote.py",
                    help="path to blast_remote.py (default: same folder)")
    ap.add_argument("--include",nargs="*",             # 0, 1, or many strings
                    metavar="Keyword",help="whitelist of keyword; keep FASTA whose *filename* contains "
                    "any Keywords (case-insensitive). If omitted, every FASTA is kept."
                    "Examples: ['marinum'], or['marinum','fortuitum'] ")
    ap.add_argument("--filter", default=None,
                    help="Entrez filter string, Example: 'txid1762[Organism]' for Mycobacteriaceae")
    ap.add_argument("--filter_name", default=None,
                    help="Entrez filter string, Example: 'Mycobacteriaceae'")
    ap.add_argument("--sleep", nargs=2, type=int, default=[11,15],
                    metavar=("MIN","MAX"),
                    help="random delay (s) between jobs")
    return ap.parse_args()

# ─── helpers ───────────────────────────────────────────────────────────
def human(seconds: float) -> str:
    """
    Convert a duration in seconds into *HH:MM:SS*.

    Parameters
    ----------
    seconds : float
        Number of seconds.

    Returns
    -------
    str
        Human-readable time string, e.g. ``'02:17:45'``.
    """
    return str(datetime.timedelta(seconds=int(seconds)))

def done_already(fasta: Path, out_dir: Path) -> bool:
    """
    Check whether this FASTA has been processed before.

    The rule is: if any *.xlsx* in *out_dir* starts with the same **stem**
    (basename without extension) as *fasta*, we consider the job done.

    Parameters
    ----------
    fasta : pathlib.Path
        Query FASTA file.
    out_dir : pathlib.Path
        Output folder that may already contain results.

    Returns
    -------
    bool
        ``True`` if a matching Excel file exists, else ``False``.
    """
    stem = fasta.stem
    return any(p.name.startswith(stem) and p.suffix == ".xlsx"
               for p in out_dir.iterdir())

def split_multifasta(fasta: Path, out_folder: Path, line_width: int = 70,):
    """
    Split a multi-record FASTA into individual one-record files.

    If *fasta* contains only a single header, the same path is returned
    unchanged. Otherwise each record is written to  
    ``<out_folder>/<stem>_<idx>.fasta`` and a list of those paths is
    returned.

    Parameters
    ----------
    fasta : pathlib.Path
        Input FASTA (may have 1 or many sequences).
    out_folder : pathlib.Path
        Folder to hold the generated single-record FASTA files.
    line_width : int, optional
        Column width used when re-wrapping sequence lines (default 70).

    Returns
    -------
    list[Path]
        List of FASTA paths that the caller should submit to BLAST.
    """

    # --- quick scan to find deflines ------------------------------
    with fasta.open() as fh:
        headers = [i for i, line in enumerate(fh) if line.startswith(">")]

    if len(headers) <= 1:                       # single-record FASTA
        return [fasta]   
    log = logging.getLogger("batch")
    log.info("Detected  multiple headers in one FASTA files, splitting into %d sequences", len(headers))

    out_folder.mkdir(parents=True, exist_ok=True) 
    
    # --- multi-record: split --------------------------------------
    split_files = []
    record_idx = 0
    seq_lines  = []
    header     = None

    def flush():
        nonlocal record_idx, seq_lines, header
        if header is None:
            return
        record_idx += 1
        out_path = out_folder / f"{fasta.stem}_{record_idx}.fasta"
        with out_path.open("w") as out:
            out.write(header)
            for chunk in wrap("".join(seq_lines).replace("\n", ""), line_width):
                out.write(chunk + "\n")
        split_files.append(out_path)
        seq_lines, header = [], None

    with fasta.open() as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                header = line
            else:
                seq_lines.append(line.strip())
        flush()                                 # last record

    return split_files

# ─── main ──────────────────────────────────────────────────────────────
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
    a = get_args()

    blast_script = Path(a.blast_script).expanduser().resolve()
    if not blast_script.is_file():
        sys.exit(f"> blast_remote.py not found: {blast_script}")

    in_dir  = Path(a.input_dir).expanduser().resolve()
    out_dir = Path(a.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # logging to screen + file
    log_path = out_dir / f"batch_log_{datetime.datetime.now():%Y%m%d_%H%M%S}.txt"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s",
        datefmt="%H:%M:%S",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_path)]
    )
    log = logging.getLogger("batch")

    # collect candidate FASTA files
    fasta_ext = {".fa", ".fna", ".fasta", ".fas"}
    fastas = sorted(
        f for f in in_dir.rglob("*")
        if f.suffix.lower() in fasta_ext
        and (a.include is None or any(kw.lower() in f.name.lower() for kw in a.include)))

    if not fastas:
        sys.exit("> No matching FASTA files found.")

    log.info("Found %d FASTA files", len(fastas))

    split_dir  = out_dir / "split_seqs"
    split_map  = {}                              # {fasta: [Path, Path, …]}
    for fa in fastas:
        split_map[fa] = split_multifasta(fa, split_dir)

    seq_total = sum(len(v) for v in split_map.values())
    log.info("Total sequences to BLAST: %d", seq_total)

    t_batch0 = time.perf_counter()
    seq_done = 0

    for fa in fastas:
        for q in split_multifasta(fa, out_dir / "split_seqs"):
            if done_already(q, out_dir):
                seq_done += 1
                log.info("> (%d/%d) %-35s  (already done)",
                         seq_done, seq_total, q.name)
                continue

            seq_done += 1
            log.info("> (%d/%d) BLAST %s", seq_done, seq_total, q.name)

            cmd = [
                sys.executable, str(blast_script), str(q),         
                "-d", "core_nt",
                "-o", str(out_dir),
                "-t", "megablast",
                "--sleep", "0", "0"
            ]

            if a.filter:                                 
                cmd.extend(["-f", a.filter])
            
            if a.filter_name:
                cmd.extend(["-n", a.filter_name])

            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                log.error(">  %s failed (exit %d)", q.name, e.returncode)
                continue
            log.info(">  %s done", q.name)

            nap = random.randint(*a.sleep)
            log.info("> sleeping %d s …", nap)
            time.sleep(nap)

    log.info("> All jobs finished in %s", human(time.perf_counter()-t_batch0))
    log.info("Log saved to %s", log_path)

# ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    main()
