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
Package management in this project was handled entirely by some form of virtualisation framework. The primary dependencies given below are necessary to orchestrate those systems. All other package requirements are accounted for.

- [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) (coreSC)
- [Conda](https://docs.conda.io/en/latest/miniconda.html) (scenic, scVelo, tobias, cellphoneDB)
- [Docker](https://docs.docker.com/engine/install/) (seurat-deseq)
- [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)
- A development environment that an connect to a container such as [VSCode](https://code.visualstudio.com/download) (seurat-deseq)

## Environments
All environments outlined in this codebase were originally assembled in Unix (linux / mac os). It is generally advised that any use of these workflows be carried out in a comparable architecture.

**seurat-deseq**
Every element within the Seurat framework was executed inside a container. The setup for this environment can be executed with a single script as shown below.

```sh
# docker must be installed!

cd levantovsky-2023/env/seurat-deseq    # navigate to the env files
vim .env file                           # edit the .env file
sudo ./rebuild                          # execute the container build - admin privs required

# the first-time setup will take a while
# following completion, a running container should be available to the host system for connection
# how to connect to a container https://code.visualstudio.com/docs/devcontainers/attach-container
```

**conda** (scenic, scVelo, tobias, cellphoneDB)

```sh
conda create --name <env> --file <this file>

# example
conda create --name scenic --file levantovsky-2023/env/conda/scenic-reqs.txt
```

## Data
All associated data inputs can be accessed from the public AWS bucket given below. Total assets are moderately large (~350GB).

```sh
# all data
aws s3 cp s3://levantovsky-2023/ . --recurisve --no-sign-request
```

Smaller data subsets can be accessed by substituting the s3 path with one of the following.

```sh
s3://levantovsky-2023/data/cpdb/
s3://levantovsky-2023/data/scenic/
s3://levantovsky-2023/data/scVelo/
s3://levantovsky-2023/data/seurat-deseq/
s3://levantovsky-2023/data/tobias/
```
