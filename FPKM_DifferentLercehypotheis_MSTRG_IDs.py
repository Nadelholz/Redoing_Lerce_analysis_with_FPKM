import pandas as pd

#with open("Lerce_Nolerer_genecount", 'r') as Lerce_nolerer_semicol, open("Lerce_Nolerer_genecount_rowsdupl", 'w') as duplicated:
#    for line in Lerce_nolerer_semicol:
#        line=line.strip()
#        if ";" in line:
#            duplicated.write(line + '\n')
#        duplicated.write(line + '\n')

df = pd.read_csv('Lerce_Nolerer_FPKM', sep='\t')

df[["ProteinID_split", "mRNAID_Split"]] = df[["ProteinID", "mRNAID"]].apply(lambda x: x.str.split(';'))

modified_rows = []

for _, row in df.iterrows():
    values_1 = row["ProteinID_split"]
    values_2 = row["mRNAID_Split"]
    if isinstance(values_1, list) and isinstance(values_2, list):
        for i in range(max(len(values_1), len(values_2))):
            new_row = row.copy()
            if i < len(values_1):
                new_row['ProteinID'] = values_1[i]
            if i < len(values_2):
                new_row['mRNAID'] = values_2[i]
            modified_rows.append(new_row)
    else:
        modified_rows.append(row)

final_df = pd.DataFrame(modified_rows)

# Reset the index of the final DataFrame
final_df.reset_index(drop=True, inplace=True)

# Remove the unnecessary columns
final_df.drop(['ProteinID_split', 'mRNAID_Split'], axis=1, inplace=True)

# Save the modified DataFrame to a new file
final_df.to_csv('Lerce_nolerer_FPKM_nolongerhasmeregd_.csv', sep='\t', index=False)







##lerer_lessthanagaric
df = pd.read_csv('Lerce_Lerer_lessthanagaric_FPKM', sep='\t')

df[["ProteinID_split", "mRNAID_Split"]] = df[["ProteinID", "mRNAID"]].apply(lambda x: x.str.split(';'))

modified_rows = []

for _, row in df.iterrows():
    values_1 = row["ProteinID_split"]
    values_2 = row["mRNAID_Split"]
    if isinstance(values_1, list) and isinstance(values_2, list):
        for i in range(max(len(values_1), len(values_2))):
            new_row = row.copy()
            if i < len(values_1):
                new_row['ProteinID'] = values_1[i]
            if i < len(values_2):
                new_row['mRNAID'] = values_2[i]
            modified_rows.append(new_row)
    else:
        modified_rows.append(row)

final_df = pd.DataFrame(modified_rows)

# Reset the index of the final DataFrame
final_df.reset_index(drop=True, inplace=True)

# Remove the unnecessary columns
final_df.drop(['ProteinID_split', 'mRNAID_Split'], axis=1, inplace=True)

# Save the modified DataFrame to a new file
final_df.to_csv('Lerce_Lerer_lessthanagaricFPKM_nolongerhasmeregd_.csv', sep='\t', index=False)



####lerce_Lerernofilt_nolongerhasmerged
df = pd.read_csv('Lerce_lerernofilt_FPKM', sep='\t')

df[["ProteinID_split", "mRNAID_Split"]] = df[["ProteinID", "mRNAID"]].apply(lambda x: x.str.split(';'))

modified_rows = []

for _, row in df.iterrows():
    values_1 = row["ProteinID_split"]
    values_2 = row["mRNAID_Split"]
    if isinstance(values_1, list) and isinstance(values_2, list):
        for i in range(max(len(values_1), len(values_2))):
            new_row = row.copy()
            if i < len(values_1):
                new_row['ProteinID'] = values_1[i]
            if i < len(values_2):
                new_row['mRNAID'] = values_2[i]
            modified_rows.append(new_row)
    else:
        modified_rows.append(row)

final_df = pd.DataFrame(modified_rows)

# Reset the index of the final DataFrame
final_df.reset_index(drop=True, inplace=True)

# Remove the unnecessary columns
final_df.drop(['ProteinID_split', 'mRNAID_Split'], axis=1, inplace=True)

# Save the modified DataFrame to a new file
final_df.to_csv('Lerce_lerernofiltFPKM_nolongerhasmeregd_.csv', sep='\t', index=False)