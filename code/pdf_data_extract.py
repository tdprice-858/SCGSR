'''
This script can read in a pdf and convert it to a pandas dataframe.
I am using this to extract calculation reaction rates
Author: Trevor
Date: Jun 15 2023
'''

from tabula import read_pdf
import numpy as np
from tabulate import tabulate
path_to_pdf = '/Users/tdprice/Desktop/\
cctc201601053_sup_0001_misc_information.pdf'
df = read_pdf(path_to_pdf, pages='all')
extract = False
table = [] #This is where I will store the three df's corresponding to my desired table

# Below the count is for each pandas df that was constructed from the PDF
# For my pdf, the desired data was in a table at the end of the SI.
# So the last three df's were parts of that table

for count, panda_db in enumerate(df):
    if count > 5:
        table.append(panda_db)

# I had to manually include these rows of data. In the pdf, the table, spanned
# three pages and at the start of each page the borderline of that row was missing,
# so tabula did not recognize it and include it within the df
row_12 = [[12, 5.01e-25, 2.50e-22, 4.39e-20, 5.50e-19]]
row_51 = [[51, -7.78e-09, -1.07e-06, -7.47e-05, -5.06e-04]]

# to_numpy converts the panda df to a numpy array
#Had to include columns (labels) of some of the df's because they were actually
# data, not table headers
# The following lines concatenate my rows of data
data = np.concatenate((table[0].to_numpy(),
                       row_12,
                 [table[1].columns.to_numpy()], table[1].to_numpy(),
                       row_51,
                 [table[2].columns.to_numpy()], table[2].to_numpy()),
                 axis=0)

# Convert everything to a float
rates = data.astype(float)
print(rates)

