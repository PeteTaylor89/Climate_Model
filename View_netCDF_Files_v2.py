# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 03:07:43 2024

@author: PKTAYLOR
"""
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 22:50:13 2024

@author: PKTAYLOR
"""
import os
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import ttk, messagebox
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

selected_files = []

index_file_path = 'C:\\Users\\PKTAYLOR\\Documents\\RAW Datasets\\MfE\\ReadMe\\MfE Climate Variables.csv'
index_df = pd.read_csv(index_file_path)

variable_name_to_id = dict(zip(index_df['Climate variable'], index_df['Variable ID']))

def filter_files():
    selected_scenario = scenario_var.get()
    selected_variable_name = variable_var.get()
    selected_variable_id = variable_name_to_id[selected_variable_name]
    selected_period = period_var.get()
    
    filtered_files = []
    for file_name in os.listdir('C:\\Users\\PKTAYLOR\\Documents\\RAW Datasets\\MfE\\netCDF_Files'):
        file_variable_id = file_name.split('_')[0]
        if (selected_scenario in file_name and
            selected_variable_id == file_variable_id and
            selected_period in file_name):
            filtered_files.append(file_name)
    
    # Update the listbox or display area with filtered files
    file_listbox.delete(0, tk.END)  # Clear current list
    for file_name in filtered_files:
        file_listbox.insert(tk.END, file_name)

def preview_data():
    global selected_files
    selected_file_index = file_listbox.curselection()
    if not selected_file_index:
        messagebox.showwarning("No File Selected", "Please select a file from the list to preview.")
        return
    
    selected_files = [file_listbox.get(idx) for idx in selected_file_index]
    
    for file_name in selected_files:
        selected_file_path = os.path.join('C:\\Users\\PKTAYLOR\\Documents\\RAW Datasets\\MfE\\netCDF_Files', file_name)
        
        try:
            dataset = nc.Dataset(selected_file_path, mode='r')
            variable_name = file_name.split('_')[0]  # Extract variable ID from file name
            variable = dataset.variables[variable_name]
            data = variable[0, :, :]  # Extract a 2D slice (first time step)

            # Handling NaN values and coordinate grids
            lats = dataset.variables['latitude'][:]
            lons = dataset.variables['longitude'][:]
            values_array = data[:,:].data

            # Replace extreme outliers with NaN
            max_valid_value = 100000  # Example: adjust based on your data range
            values_array[values_array > max_valid_value] = np.nan

            # Create a new popup window
            popup = tk.Toplevel()
            popup.title(f'Preview: {variable.long_name}')
            
            # Plotting with fixed aspect ratio
            fig, ax = plt.subplots(figsize=(8, 6))
            im = ax.pcolormesh(lons, lats, values_array, cmap='viridis', shading='auto')
            ax.set_title(variable.long_name)
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            cbar = fig.colorbar(im, ax=ax)
            cbar.set_label(variable.units)
    
            # Function to handle mouse click events
            def on_plot_click(event):
                if event.xdata is not None and event.ydata is not None:
                    # Convert clicked coordinates to nearest grid cell index
                    lon_idx = np.abs(lons - event.xdata).argmin()
                    lat_idx = np.abs(lats - event.ydata).argmin()
                    value_at_click = values_array[lat_idx, lon_idx]
            
                    # Display value near the clicked point
                    annotation_text = f'Value: {value_at_click:.2f}\nLon: {lons[lon_idx]:.2f}\nLat: {lats[lat_idx]:.2f}'
                    annotation = ax.annotate(annotation_text, xy=(event.xdata, event.ydata),
                                             xytext=(10, -30), textcoords='offset points', ha='center', va='top',
                                             bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                                             arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
                    ax.figure.canvas.draw_idle()

            # Function to handle mouse scroll events (zoom)
            def on_scroll(event):
                # Check if it's a scroll event and in the plot area
                if event.button == 'up':
                    scale_factor = 1.1
                elif event.button == 'down':
                    scale_factor = 1 / 1.1
                else:
                    scale_factor = 1.0
                
                xlim = ax.get_xlim()
                ylim = ax.get_ylim()
                
                xdata = event.xdata  # x position of the cursor in data coordinates
                ydata = event.ydata  # y position of the cursor in data coordinates
                
                if xdata and ydata:
                    ax.set_xlim([xdata - (xdata - xlim[0]) * scale_factor,
                                 xdata + (xlim[1] - xdata) * scale_factor])
                    ax.set_ylim([ydata - (ydata - ylim[0]) * scale_factor,
                                 ydata + (ylim[1] - ydata) * scale_factor])
                else:
                    ax.set_xlim(xlim[0] * scale_factor, xlim[1] * scale_factor)
                    ax.set_ylim(ylim[0] * scale_factor, ylim[1] * scale_factor)
                
                ax.figure.canvas.draw()

            # Connect the event handlers
            fig.canvas.mpl_connect('button_press_event', on_plot_click)
            fig.canvas.mpl_connect('scroll_event', on_scroll)

            # Adjust aspect ratio and tight layout
            ax.set_aspect('equal')  # Fixed aspect ratio
            fig.tight_layout()

            # Display plot in the popup window using FigureCanvasTkAgg
            canvas = FigureCanvasTkAgg(fig, master=popup)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
            # Close button
            close_button = ttk.Button(popup, text="Close", command=popup.destroy)
            close_button.pack(side=tk.BOTTOM, padx=10, pady=10)
            
        except Exception as e:
            messagebox.showerror("Error", f"Error loading file: {e}")



# Function to compare data from two selected files
def compare_data():
    global selected_files
    selected_files_indices = file_listbox.curselection()
    
    if len(selected_files_indices) != 2:
        messagebox.showwarning("Selection Error", "Please select exactly two files to compare.")
        return
    
    selected_files = [file_listbox.get(selected_files_indices[0]), file_listbox.get(selected_files_indices[1])]
    
    try:
        # Create a new popup window for comparison
        comparison_popup = tk.Toplevel()
        comparison_popup.title("Comparison")

        for idx, file_name in enumerate(selected_files):
            file_path = os.path.join('C:\\Users\\PKTAYLOR\\Documents\\RAW Datasets\\MfE\\netCDF_Files', file_name)
            dataset = nc.Dataset(file_path, mode='r')
            variable_name = file_name.split('_')[0]
            variable = dataset.variables[variable_name]
            data = variable[0, :, :]  # Extract a 2D slice (first time step)

            # Handling NaN values and coordinate grids
            lats = dataset.variables['latitude'][:]
            lons = dataset.variables['longitude'][:]
            values_array = data[:,:].data

            # Replace extreme outliers with NaN
            max_valid_value = 100000  # Example: adjust based on your data range
            values_array[values_array > max_valid_value] = np.nan

            # Plotting with fixed aspect ratio
            fig, ax = plt.subplots(figsize=(8, 6))
            im = ax.pcolormesh(lons, lats, values_array, cmap='viridis', shading='auto')
            ax.set_title(variable.long_name)
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            cbar = fig.colorbar(im, ax=ax)
            cbar.set_label(variable.units)
    
            # Function to handle mouse click events
            def on_plot_click(event):
                if event.xdata is not None and event.ydata is not None:
                    # Convert clicked coordinates to nearest grid cell index
                    lon_idx = np.abs(lons - event.xdata).argmin()
                    lat_idx = np.abs(lats - event.ydata).argmin()
                    value_at_click = values_array[lat_idx, lon_idx]
            
                    # Display value near the clicked point
                    annotation_text = f'Value: {value_at_click:.2f}\nLon: {lons[lon_idx]:.2f}\nLat: {lats[lat_idx]:.2f}'
                    annotation = ax.annotate(annotation_text, xy=(event.xdata, event.ydata),
                                             xytext=(10, -30), textcoords='offset points', ha='center', va='top',
                                             bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                                             arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
                    ax.figure.canvas.draw_idle()
    
            # Connect the event handler
            fig.canvas.mpl_connect('button_press_event', on_plot_click)

            # Adjust aspect ratio and tight layout
            ax.set_aspect('equal')  # Fixed aspect ratio
            fig.tight_layout()

            # Display plot in the comparison popup window using FigureCanvasTkAgg
            canvas = FigureCanvasTkAgg(fig, master=comparison_popup)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
        # Close button for comparison popup
        close_button = ttk.Button(comparison_popup, text="Close", command=comparison_popup.destroy)
        close_button.pack(side=tk.BOTTOM, padx=10, pady=10)

    except Exception as e:
        messagebox.showerror("Error", f"Error loading files: {e}")


def clear_scenario():
    scenario_var.set("")

def clear_variable():
    variable_var.set("")  

def clear_period():
    period_var.set("") 
    
# Create the main window
root = tk.Tk()
root.title("NetCDF File Explorer")
"""root.iconbitmap('path_to_icon_file.ico')"""

# Scenario dropdown
scenario_label = ttk.Label(root, text="Scenario:")
scenario_label.grid(row=0, column=0, padx=10, pady=5, sticky='w')
scenario_var = tk.StringVar(root)
scenario_dropdown = ttk.Combobox(root, textvariable=scenario_var, values=['historical', 'ssp126', 'ssp245', 'ssp370'], width=20)
scenario_dropdown.grid(row=0, column=1, padx=10, pady=5, sticky='w')
scenario_clear_button = ttk.Button(root, text="Clear", command=clear_scenario)
scenario_clear_button.grid(row=0, column=2, padx=10, pady=5, sticky='w')

# Variable dropdown (using climate variable names)
variable_label = ttk.Label(root, text="Climate Variable:")
variable_label.grid(row=1, column=0, padx=10, pady=5, sticky='w')
variable_var = tk.StringVar(root)
variable_dropdown = ttk.Combobox(root, textvariable=variable_var, values=index_df['Climate variable'].tolist(), width=50)
variable_dropdown.grid(row=1, column=1, padx=10, pady=5, sticky='w')
variable_clear_button = ttk.Button(root, text="Clear", command=clear_variable)
variable_clear_button.grid(row=1, column=2, padx=10, pady=5, sticky='w')

# Period dropdown
period_label = ttk.Label(root, text="Period:")
period_label.grid(row=2, column=0, padx=10, pady=5, sticky='w')
period_var = tk.StringVar(root)
period_dropdown = ttk.Combobox(root, textvariable=period_var, values=['bp1986-2005', 'bp1995-2014', 'fp2021-2040', 'fp2041-2060', 'fp2080-2099', 'wl1.5', 'wl2', 'wl3'], width=20)
period_dropdown.grid(row=2, column=1, padx=10, pady=5, sticky='w')
period_clear_button = ttk.Button(root, text="Clear", command=clear_period)
period_clear_button.grid(row=2, column=2, padx=10, pady=5, sticky='w')

# Filter button
filter_button = ttk.Button(root, text="Filter Files", command=filter_files)
filter_button.grid(row=3, column=0, columnspan=3, pady=10)

# Listbox to display filtered files
file_listbox = tk.Listbox(root, selectmode=tk.MULTIPLE, width=100)
file_listbox.grid(row=4, column=0, columnspan=3, padx=10, pady=5)

# Preview button
preview_button = ttk.Button(root, text="Preview", command=preview_data)
preview_button.grid(row=5, column=0, padx=10, pady=10)

# Comparison button
compare_button = ttk.Button(root, text="Compare", command=compare_data)
compare_button.grid(row=5, column=1, padx=10, pady=10)

root.mainloop()
