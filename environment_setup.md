# 1. Install miniforge
brew install miniforge

# 2. Initialize conda (for zsh - default on modern macOS)
# Initialize conda for zsh (required so `conda activate` works in new terminals)

conda init zsh

# 3. Restart terminal 

# 4. Create environment without quast first
mamba create -n wgs_pipeline \
  -c conda-forge -c bioconda \
  python=3.9 \
  fastqc \
  trimmomatic \
  openjdk \
  spades \
  kraken2 \
  setuptools

# 5. Activate it if any error run step 2 and step 3 again and then activate conda below
conda activate wgs_pipeline

# 6. Install Docker for Prokka(container in the script)
# Prokka is run via Docker to avoid dependency conflicts.
# Ensure Docker Desktop is installed and running:

open /Applications/Docker.app 

# 7. Install QUAST via pip (more reliable on macOS than conda)
pip install quast

# Additional utilities used by downstream steps
conda install conda-forge::ncbi-datasets-cli
conda install conda-forge::matplotlib
mamba install -c bioconda minimap2

# 8. # Flye is installed via Homebrew for long-read assembly support
brew install flye

# 9. Verify everything works
fastqc --version
trimmomatic --version
spades.py --version
kraken2 --version
quast --version
ncbi-datasets-cli --version
minimap2 --version
flye --version
