# Import the necessary libraries
import pandas as pd
import os

# Read in the count matrix data from a tab-separated file
df = pd.read_csv("count_matrix", sep ="\t")

# Print the contents of the dataframe
print(df)

# Remove the first five columns from the dataframe
df = df.drop(df.columns[1:6], axis=1)

# Print the contents of the updated dataframe
print(df)

# Create a list of column names
cols = list(df.columns)

# Remove the first element from the list of column names
cols.pop(0)

# Create an empty list for the new column names
col_names = []

# Iterate through the column names and extract the relevant information
for elem in cols:
    elem = elem.split("/")
    col_names.append(elem[-1])

# Create a list of official names for the columns
official_names = []
for elem in col_names:
    official_names.append(elem.replace(".sam",""))

# Iterate over the column names (starting from the second column) and update them
for i, name in enumerate(df.columns[1:], start=1):
    if i <= len(official_names):
        df = df.rename(columns={name: official_names[i-1]})

# Print the resulting dataframe
print(df)

# Write the resulting dataframe to a tab-separated file
df.to_csv("input_genes_DESeq2.txt", index = False, sep ="\t")
