# levantovsky-2023
Published code and computational methods for Levantovsky et. al. 2023

## Tools
- [CellPhoneDB](https://github.com/Teichlab/cellphonedb)
- [coreSC](https://github.com/ChoBioLab/coreSC)
- [DESeq2](https://github.com/mikelove/DESeq2)
- [Scenic](https://github.com/aertslab/SCENIC)
- [scVelo](https://github.com/theislab/scvelo)
- [Seurat](https://github.com/satijalab/seurat)
- [Signac](https://github.com/stuart-lab/signac)
- [Tobias](https://github.com/loosolab/TOBIAS)

## Dependencies
Package management in this project was handled entirely by some form of virtualisation framework. The primary dependencies given below are necessary to orchestrate those systems.

- [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) (coreSC)
- [Conda](https://docs.conda.io/en/latest/miniconda.html) (scenic, scVelo, tobias, cellphoneDB)
- [Docker](https://docs.docker.com/engine/install/) (seurat-deseq)

## Environments
All environments outlined in this codebase were originally assembled in a Unix environment (linux / mac os). It is generally advised that any use of these workflows be carried out in a comparable architecture.

**seurat-deseq**

**conda** (scenic, scVelo, tobias, cellphoneDB)

```sh
conda create --name <env> --file <this file>
conda create --name scenic --file levantovsky-2023/env/conda/scenic-reqs.txt
```

## Data
All associated data inputs can be found in the AWS bucket at the link below.

- https://levantovsky-2023.s3.amazonaws.com/data/

