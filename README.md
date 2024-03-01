# levantovsky-2024
Published code and computational methods for Levantovsky et. al. 2024

This codebase include submodules, so the `--recursive` flag is required.

[![DOI](https://zenodo.org/badge/624900220.svg)](https://zenodo.org/doi/10.5281/zenodo.10729881)

```sh
git clone --recursive https://github.com/ChoBioLab/levantovsky-2024.git
```

## Tools
- [CellPhoneDB](https://github.com/Teichlab/cellphonedb) - cell-ligand analysis
- [coreSC](https://github.com/ChoBioLab/coreSC) - scRNAseq and scATAC-multiome preprocessing (Seurat & Signac wrapper)
- [DESeq2](https://github.com/mikelove/DESeq2) - DEG analysis
- [Scenic](https://github.com/aertslab/SCENIC) - gene regulatory network analysis
- [Scenicplus](https://github.com/aertslab/scenicplus) - enhancer driven gene regulatory network analysis
- [scVelo](https://github.com/theislab/scvelo) - scRNAseq velocity analysis
- [Seurat](https://github.com/satijalab/seurat) - scRNAseq preprocessing
- [Signac](https://github.com/stuart-lab/signac) - scATAC-multiome preprocessing
- [Tobias](https://github.com/loosolab/TOBIAS) - Chromatin availability analysis

## Dependencies
Package management in this project was handled entirely by some form of virtualisation framework. The primary dependencies given below are necessary to orchestrate those systems. All other package requirements are generally accounted for.

- [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) (coreSC)
- [Conda](https://docs.conda.io/en/latest/miniconda.html) (scenic, scenicplus, scVelo, tobias, cellphoneDB)
- [Docker](https://docs.docker.com/engine/install/) (R packages image)
- [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) (data)
- A development environment that can connect to a container such as [VSCode](https://code.visualstudio.com/download)

## Environments
All environments outlined in this codebase were originally assembled in Linux. It is generally advised that any use of these workflows be carried out in a comparable architecture.

**R Packages**

(Seurat, Scenic, DESeq2)

Every method using Seurat, Signac, and DESeq2 (scRNAseq, ATAC-multiome, etc) was executed in a containerized environment. The setup for this environment can be instantiated with a single script as shown below. It's generally recommended to use the `-v` parameter to mount a drive from the host system to pass content into the container. 

```sh
# running the container
# the first-time setup will take a while as the image is quite large (~7.5GB)
## directly from the cli

docker run \
    -it \
    -v `pwd`:/data \
    public.ecr.aws/chobiolab/r-dev:levantovsky-2024 \
    /bin/bash

## in the background to attach to an IDE

docker run \
    -it -d \
    -v `pwd`:/data \
    public.ecr.aws/chobiolab/r-dev:levantovsky-2024

# how to connect to a container running in the background https://code.visualstudio.com/docs/devcontainers/attach-container
```

**conda** 

(scenic, scVelo, tobias, cellphoneDB)

Several methods were managed with Conda environments. Deployment of these is typically tool-specific and follows the standard Conda mechanisms. 

```sh
conda create --name <env> --file <this file>

# example
conda create --name scenic --file levantovsky-2024/env/conda/scenic-reqs.txt
```

**scenicplus**

At the time of use, the current release of scenic+ contained dependencies that could not be cleanly captured by conda. The most effective way to recreate the published method is to follow the scenic+ conda install method from the commit at the time of publication (e4bdd9f).

We've also included an optional conda env file for posterity, but this will likely suffer from dependency issues with `pycisTopic` and `pycistarget` that would need to be resolved manually. The env file is located at `levantovsky-2024/env/conda/scenicplus_env.yml`

```sh
# fresh scenicplus conda install

conda create --name scenicplus python=3.8
conda activate scenicplus
git clone https://github.com/aertslab/scenicplus
cd scenicplus
git checkout e4bdd9f
pip install -e .
```
