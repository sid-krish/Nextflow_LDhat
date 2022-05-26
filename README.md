<h1 align="center">LDhat Nextflow Pipeline</h1>
  <p align="center">
    LDhat pipeline configured to work with gene-conversion type recombination

- [About](#about)
- [Getting Started](#getting-started)
  - [Set up using conda](#set-up-using-conda)
  - [Set up using docker](#set-up-using-docker)

## About
Nextflow pipeline for LDhat, specifically for the LDhat pairwise module, configured to work with gene-conversion type recombination.

<!-- GETTING STARTED -->
## Getting Started

This LDhat nextflow pipeline is designed to be run on linux and requires nextflow to be installed. 
Dependencies are resolved either via conda or docker images. Support for HPC, docker, singularity, AWS and many other systems are provided via nextflow.

While it is possible to resolve the dependencies using conda for running on macOS, its recommended that this option be used on linux systems for which it has been extensively test.
If running on macOS it recommended that docker be used with the provided image, in which case it is similar to running in a linux environment.

It is also possible to install and run the program on Windows via [wsl](https://docs.microsoft.com/en-us/windows/wsl/install).

### Set up using conda
Instructions for installing nextflow and dependencies via conda
1. Clone the repo
   ```sh
   git clone https://github.com/sid-krish/rhometa_sim.git
   ```
2. Install the conda package manager: [Miniconda download](https://conda.io/en/latest/miniconda.html)
3. Install nextflow
   ```sh
   conda install -c bioconda nextflow
   ```
4. Adjust settings in nextflow.config file, by default it is configured to work with docker with modest resources.
   Disable the use of docker by setting the docker option to false. Disabling the use of container engines will cause conda packages to be used by default:
   ```sh
   docker {
       enabled = false
   }
   ```
5. The pipeline is now ready to run, and all dependencies will be automatically resolved with conda.

### Set up using docker
Instructions for installing nextflow and using the provided docker image for dependencies
1. Clone the repo
   ```sh
    git clone https://github.com/sid-krish/rhometa_sim.git
   ```
2. Install nextflow [Nextflow install](https://www.nextflow.io/index.html#GetStarted)
3. Install docker desktop [Docker install](https://docs.docker.com/desktop/linux/).
4. Adjust settings in nextflow.config file, by default it is configured to work with docker with modest resources.
5. In the sim_gen.nf file comment the lines related to conda, for instance:
   ```
   // conda 'conda-forge::msprime=1.1.1 conda-forge::gsl'
   ```
6. Ensure docker is running.
7. The pipeline is now ready to run, all the required dependencies are present in the docker image, that the pipeline is preconfigured to use.