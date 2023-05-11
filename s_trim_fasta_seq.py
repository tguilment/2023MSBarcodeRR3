#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
                          RPR Institution
                    **  All rights reserved **

PROJECT: Fasta processing toolbox

PURPOSE: Change each label in fasta format and remove unwanted sequences

DESCRIPTION: From the website https://www.ncbi.nlm.nih.gov/nuccore/ it is possible to download nucleotide sequences related to a specific gene.
Example: Search "Lumbrineris coi" (direct link at https://www.ncbi.nlm.nih.gov/nuccore/?term=Lumbrineris+coi)
gives 65 results that can be download by going to:
send to: > Complete record > File > FASTA format

The beginning of the file is:
>HQ932670.1 Lumbrineris japonica voucher BIOUG<CAN_:BP2010-346 cytochrome oxidase subunit 1 (COI) gene, partial cds; mitochondrial

It will be changed to
>HQ932670_Lumbrineris_japonica

A serie of tests are made to remove the unwanted sequences:
1 - The general family name ending by "idae"
2 - The line containing "sp.", "sp", "cf", "cf." or "mitochondrion"
3 - We keep only the 3 longest exemplars of the same type of sequences

HOW TO USE: in the shell or terminal type
python Path_to_script Path_to_folder_to_be_processed

Files with trimmed and removed sequences will be created in the same folder

@author: Thomas GUILMENT
Contact: thomas.guilment@gmail.com
@Contributor: Rannyele Passos Ribeiro
"""

# importing required modules
import os
import click
import numpy as np

# BETTER PYTHONIC WAY TO BE DONE USING FUNCTION
# Function with doc + Tests
# f_clean
# f_sort_alphabetical_order
# f_remove
# f_trim

# No complete mitochondrial genome (mitochondrion) were considered
# Exception for NC_156651

def f_update_file(s_path_filename):
    """
        This function clean the file then start by removing unwanted sequences that contain 
        (in lower or upper case) "sp", "cf" or "mitochondrion" in their name.
        Sequence names ending with "idae" are also removed. 
        Finally, if there are more than 3 sequences associated with the same species name 
        then only the 3 longest sequences are kept.

        Args:
            s_path_filename: Absolute path to the file that will be processed
            

        Returns:
            None

    """

    # Read all the file
    with open(s_path_filename, 'r') as f:
        lines = f.readlines()

    # Clean list (remove potential extra empty line and be sure to have a empty line 
    # before each sequence except the first one)
    l_clean = []

    # First sequence
    b_first_sequence = 1

    for line in lines:

        # Find the lines associated with the symbol '>' as a start
        # If the symbol '>' is found then extract each group of "words"
        if line[0] == '>':

            if b_first_sequence == 1:
                l_clean.append(line)
                b_first_sequence = 0

            else:
                l_clean.append('\n')
                l_clean.append(line)

        else:
            if not (line[0] == '\n') and b_first_sequence == 0:
                l_clean.append(line)

    l_clean.append('\n')

    # List of removed sequences
    l_removed = []

    # Names and sequence size
    l_names = []
    l_seq_size = []

    # Index of names
    l_index_names = []

    # Once the sequences format is cleaned then the selection can be performed
    # Index of lines
    d_line_number = -1

    while d_line_number < len(l_clean)-1:
        d_line_number += 1

        line = l_clean[d_line_number]

        # If the symbol '>' is found then extract each group of "words"
        if line[0] == '>':

            if '_' in line:
                words = line.split('_')
                b_edited = 1
            else:
                words = line.split()
                b_edited = 0

            # Test if one of the 2nd, 3rd or 4th word contain only "sp", "sp.", "cf" or "cf." (passing everything in lowercase for every tests)
            # Boolean to decide if the squence is kept or not
            b_keep = 1

            # Convert the words of interest in lowercase
            list_lowercase = [x.lower() for x in words]  # [1:4]]

            words_to_check = ['sp.', 'sp', 'cf', 'cf.', 'mitochondrion', 'mitochondrion,', 'mitochondrion,\n']
            if any(word in list_lowercase or words[1][-4:] == 'idae' for word in words_to_check):
                b_keep = 0


            # if list_lowercase.count('sp.') > 0 \
            #         or list_lowercase.count('sp') > 0 \
            #         or list_lowercase.count('cf') > 0 \
            #         or list_lowercase.count('cf.') > 0\
            #         or list_lowercase.count('mitochondrion') > 0\
            #         or list_lowercase.count('mitochondrion,') > 0\
            #         or list_lowercase.count('mitochondrion,\n') > 0\
            #         or words[1][-4::] == 'idae':
            #     b_keep = 0

        if b_edited:

            # If the sequence is kept, the name and sequence size are saved
            if b_keep == 1:
                l_names.append(words)

                l_clean[d_line_number] = '_'.join(words)

                d_size_seq = 0

                l_index_names.append(d_line_number)

                while l_clean[d_line_number] != '\n':
                    d_line_number += 1
                    d_size_seq += len(l_clean[d_line_number])

                l_seq_size.append(d_size_seq)

            # Else the sequence is removed and save in the removed list
            else:

                l_removed.append('_'.join(words))

                while l_clean[d_line_number] != '\n':
                    del l_clean[d_line_number]
                    l_removed.append(l_clean[d_line_number])

                del l_clean[d_line_number]

                d_line_number -= 1

        else:
            # If the sequence is kept, the name and sequence size are saved
            if b_keep == 1:
                l_names.append(words)

                if len(words) == 4:
                    l_clean[d_line_number] = '_'.join(
                        [words[0][:-2], words[1], words[2], words[3]]) + '\n'

                    d_size_seq = 0

                    l_index_names.append(d_line_number)

                    while l_clean[d_line_number] != '\n':
                        d_line_number += 1
                        d_size_seq += len(l_clean[d_line_number])

                    l_seq_size.append(d_size_seq)

                else:
                    l_clean[d_line_number] = '_'.join(
                        [words[0], words[1], words[2]]) + '\n'

                    d_size_seq = 0

                    l_index_names.append(d_line_number)

                    while l_clean[d_line_number] != '\n':
                        d_line_number += 1
                        d_size_seq += len(l_clean[d_line_number])

                    l_seq_size.append(d_size_seq)

            # Else the sequence is removed and save in the removed list
            else:

                if len(words) == 4:
                    l_removed.append(
                        '_'.join([words[0][:-2], words[1], words[2], words[3]]) + '\n')

                else:
                    l_removed.append(
                        '_'.join([words[0], words[1], words[2]]) + '\n')

                while l_clean[d_line_number] != '\n':
                    del l_clean[d_line_number]
                    l_removed.append(l_clean[d_line_number])

                del l_clean[d_line_number]

                d_line_number -= 1

    if b_edited:
        # Process and identify the name that are present more than 3 times
        l_aux_names = []
        for name in l_names:
            l_aux_names.append('_'.join(name[-3:-1]))

        # Only keep different names
        l_unique_names = list(set(l_aux_names))

        # For each name count how many examples they are
        l_count_names = []
        for name in l_unique_names:
            l_count_names.append(l_aux_names.count(name))

        # Identify the one that are more than 3 examples:
        l_redondant_seq = []
        l_number_of_seq = []

        # Number of sequence to keep
        d_seq_to_keep = 3

        # Auxiliary sequence (removing sequences by inserting a unique sequences to preserve the previous indexation)
        l_aux_clean = list(l_clean)

        # Better "pythonic way" to code this (should be changed)
        d_count = -1
        for d_size in l_count_names:
            d_count += 1
            if d_size > d_seq_to_keep:
                l_redondant_seq.append(l_unique_names[d_count])
                l_number_of_seq.append(l_count_names[d_count])

        # Identify and only keep the 3 biggest sequences
        # <=> removing the smallest sequences until there are 3 left
        for redondant_name in l_redondant_seq:
            d_count = -1
            l_aux_size_seq = []
            l_aux_index_seq = []
            for name in l_aux_names:
                d_count += 1
                if name == redondant_name:
                    l_aux_size_seq.append(l_seq_size[d_count])
                    l_aux_index_seq.append(d_count)

            # Identify the sequences that need to be removed
            while len(l_aux_size_seq) > d_seq_to_keep:

                d_index_line = l_index_names[l_aux_index_seq[np.argmin(
                    l_aux_size_seq)]]

                for d_aux_index_line in range(d_index_line, l_index_names[np.min([l_aux_index_seq[np.argmin(l_aux_size_seq)]+1, len(l_index_names)-1])]):

                    l_removed.append(l_clean[d_aux_index_line])

                    # Choose a special series of characters to remove later
                    l_aux_clean[d_aux_index_line] = '$!@'

                del l_aux_index_seq[np.argmin(l_aux_size_seq)]
                del l_aux_size_seq[np.argmin(l_aux_size_seq)]

        d_count = -1
        for line in l_aux_clean:
            d_count += 1
            if line == '$!@':
                del l_clean[d_count]
                d_count -= 1

        s_path_filename_updated = s_path_filename[:-6] + '_trimmed' + '.fasta'
        s_path_filename_removed = s_path_filename[:-6] + '_removed' + '.fasta'

        print('Creation of ' + s_path_filename_updated +
              ' and ' + s_path_filename_removed)
        with open(s_path_filename_updated, 'w') as f:
            f.writelines(l_clean)

        with open(s_path_filename_removed, 'w') as f:
            f.writelines(l_removed)

    else:

        # Process and identify the name that are present more than 3 times
        l_aux_names = []
        for name in l_names:
            l_aux_names.append('_'.join(name[1:3]))

        # Only keep different names
        l_unique_names = list(set(l_aux_names))

        # For each name count how many examples they are
        l_count_names = []
        for name in l_unique_names:
            l_count_names.append(l_aux_names.count(name))

        # Identify the one that are more than 3 examples:
        l_redondant_seq = []
        l_number_of_seq = []

        # Number of sequence to keep
        d_seq_to_keep = 3

        # Auxiliary sequence (removing sequences by inserting a unique sequences to preserve the previous indexation)
        l_aux_clean = list(l_clean)

        # Better "pythonic way" to code this (should be changed)
        d_count = -1
        for d_size in l_count_names:
            d_count += 1
            if d_size > d_seq_to_keep:
                l_redondant_seq.append(l_unique_names[d_count])
                l_number_of_seq.append(l_count_names[d_count])

        # Identify and only keep the 3 biggest sequences
        # <=> removing the smallest sequences until there are 3 left
        for redondant_name in l_redondant_seq:
            d_count = -1
            l_aux_size_seq = []
            l_aux_index_seq = []
            for name in l_aux_names:
                d_count += 1
                if name == redondant_name:
                    l_aux_size_seq.append(l_seq_size[d_count])
                    l_aux_index_seq.append(d_count)

            # Identify the sequences that need to be removed
            while len(l_aux_size_seq) > d_seq_to_keep:

                d_index_line = l_index_names[l_aux_index_seq[np.argmin(
                    l_aux_size_seq)]]

                for d_aux_index_line in range(d_index_line, l_index_names[np.min([l_aux_index_seq[np.argmin(l_aux_size_seq)]+1, len(l_index_names)-1])]):

                    l_removed.append(l_clean[d_aux_index_line])

                    # Choose a special series of characters to remove later
                    l_aux_clean[d_aux_index_line] = '$!@'

                del l_aux_index_seq[np.argmin(l_aux_size_seq)]
                del l_aux_size_seq[np.argmin(l_aux_size_seq)]

        d_count = -1
        for line in l_aux_clean:
            d_count += 1
            if line == '$!@':
                del l_clean[d_count]
                d_count -= 1

        s_path_filename_updated = s_path_filename[:-6] + '_trimmed' + '.fasta'
        s_path_filename_removed = s_path_filename[:-6] + '_removed' + '.fasta'

        print('Creation of ' + s_path_filename_updated +
              ' and ' + s_path_filename_removed)
        with open(s_path_filename_updated, 'w') as f:
            f.writelines(l_clean)

        with open(s_path_filename_removed, 'w') as f:
            f.writelines(l_removed)

# %% MAIN


@click.command()
@click.argument('s_path_data', default=os.getcwd(), nargs=1)
@click.option('--f', default='')
def main(s_path_data, f):

    if not (f == ''):
        s_path_data = f

    os.chdir(s_path_data)

    # Read files in folder
    list_of_file = [s_f for s_f in os.listdir(
        s_path_data) if os.path.isfile(os.path.join(s_path_data, s_f))]

    for s_filename in list_of_file:
        if (s_filename[-len(".fasta"):] == ".fasta") and \
                not (s_filename[-len("removed.fasta"):] == "removed.fasta") and \
                not (s_filename[-len("trimmed.fasta"):] == "trimmed.fasta") and \
                not (s_filename[-len("updated.fasta"):] == "updated.fasta"):
            print('Processing ' + s_filename)
            f_update_file(os.path.join(s_path_data, s_filename))
        else:
            print(
                s_filename + " is skipped because it is either already proccessed or not a .fasta file")


if __name__ == '__main__':
    main()
