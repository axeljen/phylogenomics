# Phylogenomics scripts

This repo contains general parsing and wrapper scripts for some phylogenomics data processing and analyses. 

## Dependencies
Requires Python/3.X (most scripts were latest run with 3.12.3).

Install dependencies with pip:
```bash
pip install -r requirements.txt
```

Some of the scripts are wrappers around other software, which are then required to be installed and in PATH:

windowPhylogenies.py requires either iqtree (tested with 2.2.2.6) or phyml (tested with 3.3.20190321), depending on what options are used.

circularizeAndRotate.py uses blast (2.15.0+).

These scripts were written for specific use-cases that I came across, so they may need modification for other use-cases, but anyone is of course free to use/modify them if they can be useful.