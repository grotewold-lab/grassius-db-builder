# grassius-db-builder
Build the grassius chado database from raw inputs

Output a compressed database snapshot that is compatible with the grassius website


### Requirements

Basic usage requires only Python 3.8 (see usage at the bottom of this page)

Running the full pipeline requires a big specialized system

* Python >= 3.8
* HMMER 3.3.2
* BLAST+ 2.9.0
* BioPerl
* Docker, Docker-compose


### Design Overview 

This repository contains
* input data/metadata
* minimal docker configuration to temporarily host a database
* self-contained scripts to perform individual steps
* pipelines to perform a sequence of steps (build a database from scratch)
* tests to verify the contents of a database

An internet connection is only necessary to download publicly available input data, which can be circumvented (see below).


### Inputs

* a snapshot of a blank chado database
* fasta files
* hmm files
* criteria for assigning families to genes
* lists of names from the old grassius website

This repository contains checksums to ensure the integrity of all input data.

For any inputs that are publicly available, this repository contains just a URL and a checksum (no data). Those inputs will be automatically downloaded into a gitignored folder that may be manually backed up and swapped out.



### steps

* relate gene models between versions of the maize genome (e.g.: GRMZM2G042920 = Zm00001d026317 = Zm00001eb431080)

* identify domains in protein transcripts (e.g.: PF00010 -> GRMZM2G042920_P02)

* assign family names to gene models (e.g.: bHLH -> GRMZM2G042920)

* assign protein names to sets of gene models (e.g.: ZmbHLH7 -> GRMZM2G042920,Zm00001d026317,Zm00001eb431080)


### Usage

clone this repository and install python dependencies

```
git clone https://github.com/grotewold-lab/grassius-db-builder.git
cd grassius-db-builder
pip install -r requirements.txt
```


run the example script which involves preparing all the inputs

```
$ python example.py
 
...
downloading public input "maize_v3_proteins"...
	url: http://ftp.maizegdb.org/MaizeGDB/FTP/B73_RefGen_v3/Zea_mays.AGPv3.22.pep.all.fa.gz
extracting downloaded archive for "maize_v3_proteins"...
checking integrity of input "maize_v3_proteins"...
downloading public input "maize_v4_proteins"...
...
all inputs verified. Recommend backing up this folder:
	/home/tessmero/Documents/github/grassius-db-builder/inputs
```

subsequent runs are fast and do not require an internet connection

```
$ python example.py

checking integrity of input "chado_template"...
checking integrity of input "family_rules"...
checking integrity of input "selfbuild_hmm"...
checking integrity of input "pfam_hmm"...
checking integrity of input "maize_v3_proteins"...
checking integrity of input "maize_v4_proteins"...
checking integrity of input "maize_v5_proteins"...
all inputs verified. Recommend backing up this folder:
	/home/tessmero/Documents/github/grassius-db-builder/inputs
```

### Reset docker environment

The gdb.chado module provides functions to build a new database in a docker container. By default it will attempt to connect to an existing docker container "gdb-chado-container". To test a pipeline from start to finish, the docker environment should be reset:

```
docker container kill gdb-chado-container
docker container prune
```