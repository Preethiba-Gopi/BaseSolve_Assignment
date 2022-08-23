from curses import raw
from email import header
from pydoc import describe
import random
import requests
import pandas as pd


def rev_comp(seq):
    '''
        Reverse_Complement function
    '''
    comp = {'A': 'T', 'T': 'A', 'G': 'C',
            'C': 'G', 'N':'N'}  # complement lookup based on Slide #1
    reversed_seq = seq[::-1]
    rev_comp_seq = ''.join([comp[base] for base in reversed_seq])
    return rev_comp_seq


# Fetching fasta sequence from given url
fasta_seq_res = requests.get(
    'https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/MonkeyPox.fn', verify=False)
raw_sequence = fasta_seq_res.text


# Array of dictionary to store all sequence information
seq_data = []

for each_seq in raw_sequence.split('>')[1:]:
    sequence_data = each_seq.split('\n')
    sequence_identifiers = sequence_data[0]
    seq = ''.join(sequence_data[1:])
    seq_data.append({
        "id": sequence_identifiers.split('|')[3],
        "sequence": seq,
        "desc": sequence_identifiers.split('|')[4],
        "identifier": sequence_identifiers
    })

# dataframe from sequence information
seq_df = pd.DataFrame(data=seq_data)

# dataframe to create BED file with columns(from header in input fasta file) id, start, end and sequence description 
bed_df = pd.DataFrame(data=[i["id"] for i in seq_data], columns=["#id"])
bed_df["start"] = [random.randint(1, 500) for i in seq_df.index] #random start position for extracting subsequence
bed_df["end"] = [random.randint(500, 1000) for i in seq_df.index] #random end position for extracting subsequence
bed_df["desc"] = seq_df["desc"]

# dataframe to BED file
bed_df.to_csv("sample.bed", sep='\t', index=False, header=False)

# Creating output fasta file with processed and modified input data
with open("output.fasta", 'w') as output:
    # Iterating BED dataframe
    for index, bed_info in bed_df.iterrows():
        
        # extracting corresponding original sequence by matching ids between bed and seq dataframe
        extracted_df = seq_df.loc[bed_info["#id"] == seq_df["id"]]
        
        # substring from the random position generated in the BED file 
        extracted_seq = extracted_df["sequence"].iloc[0][int(
            bed_info["start"]):int(bed_info["end"])]

        # Reverse Complement for the extracted sequence
        rev_comp_seq = rev_comp(extracted_seq)
        rev_comp_ls = list(rev_comp_seq)

        # Generating random mutated sequence
        random_mutation_position = random.randint(0, len(rev_comp_seq)-1)
        rev_comp_ls[random_mutation_position] = random.choice(
            ['a', 't', 'g', 'c']) # for better identification, modified bases are in lowercase
        processed_seq = ''.join(rev_comp_ls)

        # Replacing original sequence with modified sequence
        original_Seq = extracted_df["sequence"].iloc[0]
        processed_original_seq = original_Seq[0:int(
            bed_info["start"])-1]+processed_seq+original_Seq[int(bed_info["end"]):]
        _idx = extracted_df["identifier"].iloc[0]

        # writing modified sequence to a multi fasta file
        output.write('>'+_idx+'\n'+processed_original_seq+'\n')
        
