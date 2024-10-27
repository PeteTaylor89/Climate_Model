# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 13:54:45 2023

@author: Peter Taylor

This script creates a query list where each query is <40000 lines output from cliflo

"""

import pandas as pd

# Read the original CSV file
input_file = 'Z:\\Climate Model Drive Backup_TO SORT\\Climate Model V0.1\\CliFlo Data\\Agent_Numbers_Wind.csv'
df = pd.read_csv(input_file)


# Initialize variables
current_row = []
max_days_per_row = 39999
current_days = 0
output_rows = []

# Iterate through the rows of the DataFrame
for index, row in df.iterrows():
    agent = row['Agent']
    start_date = pd.to_datetime(row['Start_Date'])
    end_date = pd.to_datetime(row['End_Date'])
    days = (end_date - start_date).days
    
    # Check if adding the current agent exceeds the maximum days per row
    if current_days + days <= max_days_per_row:
        current_row.append(agent)
        current_days += days
    else:
        # Start a new row if the maximum days are exceeded
        output_rows.append(current_row)
        current_row = [agent]
        current_days = days

# Add the last row
output_rows.append(current_row)

# Create the output CSV file with no headers
output_file = 'Z:\\Climate Model Drive Backup_TO SORT\\Climate Model V0.1\\CliFlo Data\\Agent_Numbers_Wind_Queried.csv'
with open(output_file, 'w') as file:
    for row in output_rows:
        file.write(','.join(map(str, row)) + '\n')

print(f"New CSV file '{output_file}' with the specified layout has been created.")
