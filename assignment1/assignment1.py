
import pandas as pd


def rev_comp(seq):
    '''
        Reverse_Complement function
    '''
    comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G'} # complement lookup based on Slide #1
    reversed_seq = seq[::-1]
    rev_comp_seq = ''.join([comp[base] for base in reversed_seq])
    return rev_comp_seq


# parse CSV file into DataFrame
input_file = pd.read_csv("problem_sheet.csv")

# Filter for non empty sampleID
filtered = input_file[-pd.isnull(input_file['Unnamed: 2'])]

# extract headers from filtered DataFrame
headers = filtered.head(1).iloc[0]

# all rows other than first row are data
all_rows = filtered.iloc[1:]

# create new DataFrame using filterd data and header
parsed_data = pd.DataFrame(data=all_rows.to_numpy(), columns=headers)

# extracting index2 column
index2 = parsed_data['index2']
index2.name = "Original String"
final_data = pd.DataFrame(data = index2)

final_data['Reverse Complemented String'] = index2.apply(lambda x:rev_comp(x))

print(final_data)
final_data.to_csv('result.csv',index=False)

