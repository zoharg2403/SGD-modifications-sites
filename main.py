import pandas as pd

# 1 - filter 'All SGD modifications' sheen in 'PTM project raw data.xlsx':
#     If for the same ORF, there is a site which appears more than once, and has the same modification,
#     the repetitions is deleted

All_SGD_modifications = pd.read_excel('PTM project raw data.xlsx', sheet_name='All SGD modifications', header=None)
All_SGD_modifications = All_SGD_modifications.drop_duplicates()

# 2 - split column D into 3 columns:
#     This column represents which amino acid has a modification: it contains a letter representing the amino acid
#     and a number representing its position. split this into a column of the amino acid identity (letter) and a
#     column of its position (number).

modifiedAA = All_SGD_modifications[3].apply(lambda v: v[0])
AAposition = All_SGD_modifications[3].apply(lambda v: v[1:])
All_SGD_modifications.insert(4, 'modified AA', modifiedAA)
All_SGD_modifications.insert(5, 'modified AA position', AAposition)

# export to .csv file
All_SGD_modifications.to_csv('All SGD modifications - processed.csv', header=False, index=False)

# 3 - Identify for each AA, where it lies in the protein topology, according to Topologyeast:
#     Columns G-O represent data from the 9 different prediction algorithms.
#     The string of letters represents in which part of the protein each amino acid is: i=inside the cytosol, m=in the membrane, o=outside the cytosol
#     Figure out where the modified site sits:
#     return the letter (i/m/o) that appears in the position corresponding to the modified amino acid

Topologyeast = pd.read_excel('PTM project raw data.xlsx', sheet_name='Topologyeast', usecols='A,G:O')
# create new dataframe 'modified_sites' -
# to fill with the position of the modified amino acid (i/m/o) for each prediction algorithem
modified_sites = All_SGD_modifications[[1, 'modified AA', 'modified AA position']]
modified_sites = modified_sites.rename({1: 'ORF'}, axis=1)
modified_sites[list(Topologyeast)[1:]] = ['None'] * 9
# for each row in 'modified_sites' (each AA modification):
for r in modified_sites.index:
    curORF = modified_sites['ORF'][r]
    cur_mod_pos = int(modified_sites['modified AA position'][r]) - 1  # the position of the modification is indexed strating from 1, python indexes start from 0, hence the '-1'
    curORFpredictions = Topologyeast[Topologyeast['ORF'] == curORF].reset_index(drop=True)
    # fill the letter i/m/o in 'cur_mod_pos' position according to each algorithm
    if not curORFpredictions.empty:
        if cur_mod_pos <= len(curORFpredictions['TMHMM'][0]) - 1:
            modified_sites['TMHMM'][r] = curORFpredictions['TMHMM'][0][cur_mod_pos]
            modified_sites['TOPCONS'][r] = curORFpredictions['TOPCONS'][0][cur_mod_pos]
            modified_sites['OCTOPUS'][r] = curORFpredictions['OCTOPUS'][0][cur_mod_pos]
            modified_sites['Philius'][r] = curORFpredictions['Philius'][0][cur_mod_pos]
            modified_sites['PolyPhobius'][r] = curORFpredictions['PolyPhobius'][0][cur_mod_pos]
            modified_sites['SCAMPI'][r] = curORFpredictions['SCAMPI'][0][cur_mod_pos]
            modified_sites['SPOCTOPUS'][r] = curORFpredictions['SPOCTOPUS'][0][cur_mod_pos]
            modified_sites['HMMtop'][r] = curORFpredictions['HMMtop'][0][cur_mod_pos]
            modified_sites['memsat_svm'][r] = curORFpredictions['memsat_svm'][0][cur_mod_pos]
        elif cur_mod_pos > len(curORFpredictions['TMHMM'][0]) - 1:
            # print the ORFs that the modification position is higher then the length of the prediction
            print(curORF)

# export to .csv file
modified_sites.to_csv('Modified Sites.csv', index=False)

# 4 - If the modified AA is in the transmembrane domain:
# filter the dataset - for each algorithem, modification sites = M or m
modified_sites_M = {'TMHMM': modified_sites[['ORF', 'modified AA position']][modified_sites['TMHMM'].isin(['m', 'M'])],
                    'TOPCONS': modified_sites[['ORF', 'modified AA position']][modified_sites['TOPCONS'].isin(['m', 'M'])],
                    'OCTOPUS': modified_sites[['ORF', 'modified AA position']][modified_sites['OCTOPUS'].isin(['m', 'M'])],
                    'Philius': modified_sites[['ORF', 'modified AA position']][modified_sites['Philius'].isin(['m', 'M'])],
                    'PolyPhobius': modified_sites[['ORF', 'modified AA position']][modified_sites['PolyPhobius'].isin(['m', 'M'])],
                    'SCAMPI': modified_sites[['ORF', 'modified AA position']][modified_sites['SCAMPI'].isin(['m', 'M'])],
                    'SPOCTOPUS': modified_sites[['ORF', 'modified AA position']][modified_sites['SPOCTOPUS'].isin(['m', 'M'])],
                    'HMMtop': modified_sites[['ORF', 'modified AA position']][modified_sites['HMMtop'].isin(['m', 'M'])],
                    'memsat_svm': modified_sites[['ORF', 'modified AA position']][modified_sites['memsat_svm'].isin(['m', 'M'])]}

# add columns to dataframes in 'modified_sites_M':
for k in modified_sites_M:
    v = modified_sites_M[k]
    v[['Transmembrane domain start position', 'Transmembrane domain end position', 'Transmembrane domain length', 'Modified site position in the transmembrane domain']] = ['None'] * 4
    for r in v.index:
        curPrediction = Topologyeast[k][(Topologyeast['ORF'] == v['ORF'][r])].reset_index(drop=True)[0]
        cur_mod_pos = int(v['modified AA position'][r]) - 1  # the position of the modification is indexed strating from 1, python indexes start from 0, hence the '-1'
        # find Transmembrane domain start position:
        start_pos = cur_mod_pos
        while curPrediction[start_pos] == 'M' or curPrediction[start_pos] == 'm':
            start_pos -= 1
        v['Transmembrane domain start position'][r] = start_pos
        # find Transmembrane domain end position:
        end_pos = cur_mod_pos
        while curPrediction[end_pos] == 'M' or curPrediction[end_pos] == 'm':
            end_pos += 1
        v['Transmembrane domain end position'][r] = end_pos
        # calc other parameters:
        v['Transmembrane domain length'][r] = end_pos - start_pos + 1
        v['Modified site position in the transmembrane domain'][r] = cur_mod_pos - start_pos + 1
    # save to .csv file
    v.to_csv('Transmembrane Domain Modified Site - ' + k + ' Algorithm.csv', index=False)
