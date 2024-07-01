# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 21:17:41 2024

@author: PKTAYLOR
"""

import tkinter as tk
import pandas as pd
import os
from tkinter import filedialog, messagebox, ttk

# Function to select CSV files
def select_files():
    file_paths = filedialog.askopenfilenames(filetypes=[("CSV files", "*.csv")])
    if file_paths:
        for file in file_paths:
            process_file(file)

# Function to process the CSV file and display the last date
def process_file(file_path):
    try:
        df = pd.read_csv(file_path)
        if 'Date(NZST)' in df.columns:
            df['Date(NZST)'] = pd.to_datetime(df['Date(NZST)'], dayfirst=True)
            last_date = df['Date(NZST)'].max().strftime('%d/%m/%Y')
            first_date = df['Date(NZST)'].min().strftime('%d/%m/%Y')
            all_dates = pd.date_range(start=df['Date(NZST)'].min(), end=df['Date(NZST)'].max())
            missing_dates = all_dates.difference(df['Date(NZST)'])
            missing_dates_count = len(missing_dates)
            short_file_path = os.path.join(os.path.basename(os.path.dirname(file_path)), os.path.basename(file_path))
            add_summary(file_path, short_file_path, first_date, last_date, missing_dates, missing_dates_count)
        else:
            messagebox.showerror("Error", "No 'Date(NZST)' column found in the CSV file.")
    except Exception as e:
        messagebox.showerror("Error", str(e))

# Function to add summary to the table
def add_summary(full_path, file_path, first_date, last_date, missing_dates, missing_dates_count):
    tree.insert('', 'end', values=(file_path, first_date, last_date, missing_dates_count), tags=(full_path,))

# Function to show missing dates
def show_missing_dates(event):
    item = tree.selection()[0]
    file_path = tree.item(item, 'tags')[0]
    try:
        df = pd.read_csv(file_path)
        df['Date(NZST)'] = pd.to_datetime(df['Date(NZST)'], dayfirst=True)
        all_dates = pd.date_range(start=df['Date(NZST)'].min(), end=df['Date(NZST)'].max())
        missing_dates = all_dates.difference(df['Date(NZST)'])
        missing_dates_str = '\n'.join(missing_dates.strftime('%d/%m/%Y'))
        messagebox.showinfo("Missing Dates", f"Missing dates for {file_path}:\n\n{missing_dates_str}")
    except Exception as e:
        messagebox.showerror("Error", str(e))

# Function to run another script
def run_script():
    global last_dates
    last_dates = [tree.item(child)['values'][2] for child in tree.get_children()]
    # Import and run your other script here, passing the last_dates as needed

# Setting up the main window
root = tk.Tk()
root.title("CSV Date Viewer")

# Setting up the treeview table
columns = ('File Path', 'Earliest Date', 'Latest Date', 'Missing Dates Count')
tree = ttk.Treeview(root, columns=columns, show='headings')
tree.heading('File Path', text='File Path')
tree.heading('Earliest Date', text='Earliest Date')
tree.heading('Latest Date', text='Latest Date')
tree.heading('Missing Dates Count', text='Missing Dates Count')
tree.pack(fill=tk.BOTH, expand=True)

# Bind the treeview to the show_missing_dates function
tree.bind('<Double-1>', show_missing_dates)

# Adding buttons
button_frame = tk.Frame(root)
button_frame.pack(fill=tk.X)

select_button = tk.Button(button_frame, text="Select CSV Files", command=select_files)
select_button.pack(side=tk.LEFT, padx=5, pady=5)

run_button = tk.Button(button_frame, text="Run Script", command=run_script)
run_button.pack(side=tk.RIGHT, padx=5, pady=5)

root.mainloop()
