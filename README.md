# 7Be

This repository contains data and scripts related to the 7Be compound system evaluation for INDEN-LE commitee.

## Table of Contents

* [Overview](#overview)
* [Repository Structure](#repository-structure)

## Overview

The 7Be project focuses on:

* Fetching and processing experimental data from the EXFOR database.
* Preparing input files for nuclear reaction evaluation tools like SAMMY and AZURE2.

## Repository Structure

```
7Be/
├── data/                 # Directory for storing experimental and processed data
├── Compare.ipynb         # Jupyter notebook for comparing Ian's data with EXOFR data
├── Notes.md              # Notes and documentation related to the data transformation with respect to what Ian has used
├── fetch_exfor.py        # Script to fetch data from the EXFOR database (requires FUDGE and X4i installed)
├── nuclear.py            # Utility functions for nuclear data processing
├── prepare_azure.py      # Script to processes data for AZURE2
├── prepare_sammy.py      # Script to processes data for SAMMY
├── ripl2pops.xml         # XML file for RIPL compatibility
├── run.sh                # Shell script to execute the full data processing and simulation pipeline
└── .gitignore            # Specifies files and directories to be ignored by Git
```