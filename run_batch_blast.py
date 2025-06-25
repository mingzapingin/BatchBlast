#!/usr/bin/env python3
"""
run_batch_blast.py
Call blast_remote.py sequentially for many query FASTA files.

✓ skips files already processed (xlsx present)
✓ logs every step
✓ 11-15 s random pause between jobs

Author: Kimhun Tuntikawinwong
"""

from pathlib import Path
from textwrap import wrap
import subprocess, sys, time, datetime, random, logging, argparse

# ─── command-line options ──────────────────────────────────────────────
def get_args():
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
                    help="Entrez filter string (default: Mycobacteriaceae)")
    ap.add_argument("--sleep", nargs=2, type=int, default=[11,15],
                    metavar=("MIN","MAX"),
                    help="random delay (s) between jobs")
    return ap.parse_args()

# ─── helpers ───────────────────────────────────────────────────────────
def human(seconds: float) -> str:
    return str(datetime.timedelta(seconds=int(seconds)))

def done_already(fasta: Path, out_dir: Path) -> bool:
    """
    True if an .xlsx whose stem starts with the FASTA stem already exists.
    """
    stem = fasta.stem
    return any(p.name.startswith(stem) and p.suffix == ".xlsx"
               for p in out_dir.iterdir())

def split_multifasta(fasta: Path, out_folder: Path, line_width: int = 70,):
    """
    If *fasta* contains > 1 record, write each record to its own
    FASTA file in *out_folder* and return the list of new paths.
    If the file already has a single record, return [fasta].

    Parameters
    ----------
    fasta : Path
        Path to the input FASTA.
    out_folder : Path
        Where the per-sequence files will be written.
    line_width : int
        Wrap sequence lines to this many characters (default 70).

    Returns
    -------
    list[Path]
        Paths of FASTA files that should be BLASTed.
    """
    out_folder.mkdir(parents=True, exist_ok=True)

    # --- quick scan to find deflines ------------------------------
    with fasta.open() as fh:
        headers = [i for i, line in enumerate(fh) if line.startswith(">")]

    if len(headers) <= 1:                       # single-record FASTA
        return [fasta]
    log = logging.getLogger("batch")
    log.info("Detected  multiple headers in one FASTA files, splitting into %d sequences", len(headers))

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
    a = get_args()

    blast_script = Path(a.blast_script).expanduser().resolve()
    if not blast_script.is_file():
        sys.exit(f"❌ blast_remote.py not found: {blast_script}")

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
        sys.exit("❌ No matching FASTA files found.")

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
                log.info("⏩ (%d/%d) %-35s  (already done)",
                         seq_done, seq_total, q.name)
                continue

            seq_done += 1
            log.info("▶︎ (%d/%d) BLAST %s", seq_done, seq_total, q.name)

            cmd = [
                sys.executable, str(blast_script), str(q),         
                "-d", "core_nt",
                "-o", str(out_dir),
                "-t", "megablast",
                "-n", "Mycobacteriaceae",
                "--sleep", "0", "0"
            ]

            if a.filter:                                 
                cmd.extend(["-f", a.filter])
            
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                log.error("✖  %s failed (exit %d)", q.name, e.returncode)
                continue
            log.info("✔  %s done", q.name)

            nap = random.randint(*a.sleep)
            log.info("🕒 sleeping %d s …", nap)
            time.sleep(nap)

    log.info("🎉 All jobs finished in %s", human(time.perf_counter()-t_batch0))
    log.info("Log saved to %s", log_path)

# ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    main()
