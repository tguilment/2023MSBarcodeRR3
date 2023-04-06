# Project for this article
Scripts used for the edition of fasta files made with python

## Purpose
Change each label in fasta format and remove unwanted sequences for each fasta files contained in a folder

## Description 
From the website https://www.ncbi.nlm.nih.gov/nuccore/ it is possible to download nucleotide sequences related to a specific gene.

The unwanted sequences are removed if:
1 - The general family name ending by "idae"
2 - The line containing "sp.", "sp", "cf", "cf." or "mitochondrion"
3 - We keep only the 3 longest exemplars of the same type of sequences

### Example
Search "Lumbrineris coi" (direct link at https://www.ncbi.nlm.nih.gov/nuccore/?term=Lumbrineris+coi)
gives 65 results that can be download by going to:
send to: > Complete record > File > FASTA format

The beginning of the file is:
>\>HQ932670.1 Lumbrineris japonica voucher BIOUG<CAN_:BP2010-346 cytochrome oxidase subunit 1 (COI) gene, partial cds; mitochondrial

The script will change every line description to:
>\>HQ932670_Lumbrineris_japonica

Then all the sequences containing in upper or lowercase "sp.", "sp", "cf", "cf." or "mitochondrion" and/or having family name ending by "idae" are removed. Finally, only the 3 longest sequences belonging to the same name are kept.

## Usage

To use this script, you must have Python 3 installed. You can download Python from the [official Python website](https://www.python.org/downloads/).
Once you have Python installed, you can run the script from the command line. Here is an example usage, in the shell or terminal type:

```bash
python Path_to_script/script.py Path_to_folder_to_be_processed
```

Evry FASTA files in the selected folder will be processed and the corresponding files with trimmed and removed sequences will be created in the same folder

## Contributing
If you would like to contribute to this script, please feel free to submit a pull request or write at thomas.guilment@gmail.com.
@author: Thomas GUILMENT
Contact: thomas.guilment@gmail.com

@Contributor: Rannyele Passos Ribeiro

