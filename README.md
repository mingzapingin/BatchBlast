# BLAST batch wrapper

`blast_remote.py` – run a single remote BLAST (`blastn`, `megablast`, etc.),
save results in TSV **and** Excel.

`run_batch_blast.py` – recurse through a folder, split multi-record FASTA
files, skip finished jobs, and call `blast_remote.py` for every sequence.

## Quick start

```bash
# 1. create env (conda or venv)
conda create -n blastwrap python=3.10
conda activate blastwrap

# 2. install Python deps
pip install -r requirements.txt          # installs pandas & openpyxl

# 3. install NCBI BLAST+
conda install -c bioconda blast          # (or grab binaries from NCBI)

# 4. run a single query
python blast_remote.py query.fasta -t megablast

# 5. or batch-process a folder
python run_batch_blast.py ./input_fasta ./blast_results
