# grassius-db-builder
Build the grassius chado database from raw inputs

Output a compressed database snapshot that is compatible with the grassius website


### Requirements

Running the full pipeline requires a big specialized system

* Python 3
* HMMER 3.3.2
* BLAST+ 2.9.0
* Docker, Docker-compose


### design overview 

This repository contains
* input data/metadata
* minimal docker configuration to temporarily host a database
* self-contained scripts to perform individual steps
* pipelines to perform a sequence of steps (build a database from scratch)
* tests to verify the contents of a database

An internet connection is only necessary to download publicly available input data, which can be circumvented (see below).


### inputs

* a snapshot of a blank chado database
* fasta files
* hmm files
* criteria for assigning families to genes
* lists of names from the old grassius website

This repository contains checksums to ensure the integrity of all input data.

For any inputs that are publicly available, this repository contains just a URL and a checksum (no data). Those inputs can be automatically downloaded once into a gitignored folder that may be manually backed up and swapped out. 



### steps

* relate gene models between versions of the maize genome (e.g.: GRMZM2G042920 = Zm00001d026317 = Zm00001eb431080)

* identify domains in protein transcripts (e.g.: PF00010 -> GRMZM2G042920_P02)

* assign family names to gene models (e.g.: bHLH -> GRMZM2G042920)

* assign protein names to sets of gene models (e.g.: ZmbHLH7 -> GRMZM2G042920,Zm00001d026317,Zm00001eb431080)
