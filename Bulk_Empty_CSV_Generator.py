# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 14:19:34 2023

@author: Peter Taylor

Bulk Empty CSV Generator


"""

import csv

# Number of CSV files to create
num_files = 50

# Create the empty CSV files
for i in range(1, num_files + 1):
    file_name = f'Z:\\Climate Model Drive Backup_TO SORT\\Climate Model V0.1\\CliFlo Data\\RAW_Wind_Max_Gust\\Output_{i}.csv'
    with open(file_name, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        # You can optionally write a header row if needed
        # csv_writer.writerow(['Header1', 'Header2', 'Header3'])
