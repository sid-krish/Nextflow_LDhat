<h1 align="center">LDhat Nextflow Pipeline</h1>
  <p align="center">
    LDhat pipeline configured to work with gene-conversion type recombination

- [About](#about)
- [Getting Started](#getting-started)
  - [Set up using conda](#set-up-using-conda)
  - [Set up using docker](#set-up-using-docker)
- [Quick Start and Output](#quick-start-and-output)
  - [run_ldhat.nf](#run_ldhatnf)
  - [run_theta_est.nf](#run_theta_estnf)
- [Pipeline Options and Advanced Usage](#pipeline-options-and-advanced-usage)

## About
Nextflow pipeline for LDhat, specifically for the LDhat pairwise module, configured to work with gene-conversion type recombination. Based on the journal article https://academic.oup.com/genetics/article/160/3/1231/6052507 and the ldhat package https://github.com/auton1/LDhat. Additionaly, a theta estimate pipeline has been implemented based on equation 1 from the journal article.

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
   git clone https://github.com/sid-krish/Nextflow_LDhat.git
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
    git clone https://github.com/sid-krish/Nextflow_LDhat.git
   ```
2. Install nextflow [Nextflow install](https://www.nextflow.io/index.html#GetStarted)
3. Install docker desktop [Docker install](https://docs.docker.com/desktop/linux/).
4. Adjust settings in nextflow.config file, by default it is configured to work with docker with modest resources.
5. In the nextflow.config file comment the following line:
   ```
   // conda='environment.yaml'
   ```
6. Ensure docker is running.
7. The pipeline is now ready to run, all the required dependencies are present in the docker image, that the pipeline is preconfigured to use.

<!-- QUICK START AND OUTPUT -->
## Quick Start and Output
The following quick start example makes use of the files in [toy_dataset.zip](https://github.com/sid-krish/Nextflow_LDhat/blob/main/toy_dataset.zip).

This example is designed to help ensure that run_ldhat.nf (rho estimation pipeline) and run_theta_est.nf (theta estimation pipeline) are configured properly and to demonstrate a typical workflow.
The toy datasets were simulated using the simulation pipeline [Nextflow_LDhat_Sim](https://github.com/sid-krish/Nextflow_LDhat_Sim).

### run_ldhat.nf
run_ldhat.nf will by default estimate the theta per site and use it for generating the appropriate lookup table for the number of sequences in the fasta file. The lookup table is then used for the rho estimate.

The default settings of the run_ldhat.nf pipeline have been configured to work with the quick start example as is.

For this example the command to run is:
```sh
nextflow run run_ldhat.nf --input_fasta toy_dataset/example.fa
```

The output files will be saved to 'LDhat_Output'. The file ending with pairwise_outfile.txt and processed_results.csv have the rho estimates and file ending with theta_est.csv will have the theta (per site) estimate.

The result in the file ending with processed_results.csv should be:
```
max_rho,max_lk
9.000,-41836138.059
```

The result in the file ending with theta_est.csv should be:
```
0.010264296769556949
```

### run_theta_est.nf
run_theta_est.nf is used to get the theta (per site) estimate only.

For this example the command to run is:
```sh
nextflow run run_theta_est.nf --input_fasta toy_dataset/example.fa 
```

The output files will be saved to 'LDhat_Theta_Output'. The file ending with theta_est.csv will have the theta (per site) estimate. It should have the following results:
```
0.010264296769556949
```

<!-- PIPELINE OPTIONS AND ADVANCED USAGE -->
## Pipeline Options and Advanced Usage


