# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 03:07:43 2024

@author: PKTAYLOR
"""

import os
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import ttk, messagebox
import tkinter.filedialog as filedialog
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
    selected_period = period_var.get()

    filtered_files = []
    for file_name in os.listdir('C:\\Users\\PKTAYLOR\\Documents\\RAW Datasets\\MfE\\netCDF_Files'):
        file_variable_id = file_name.split('_')[0]

        if (selected_scenario in file_name and
            (not selected_variable_name or variable_name_to_id[selected_variable_name] == file_variable_id) and
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
            values_array = data[:, :].data

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
                if event.button == 'down':
                    scale_factor = 1.1
                elif event.button == 'up':
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

            # Initialize previous coordinates for dragging
            drag_data = {'x': None, 'y': None}

            # Function to start dragging
            def on_drag_start(event):
                drag_data['x'] = event.xdata
                drag_data['y'] = event.ydata

            # Function to handle dragging
            def on_drag(event):
                if drag_data['x'] is None or drag_data['y'] is None:
                    return
                if event.inaxes != ax:
                    return
                dx = event.xdata - drag_data['x']
                dy = event.ydata - drag_data['y']
                
                xlim = ax.get_xlim()
                ylim = ax.get_ylim()
                
                ax.set_xlim(xlim[0] - dx, xlim[1] - dx)
                ax.set_ylim(ylim[0] - dy, ylim[1] - dy)
                
                drag_data['x'] = event.xdata
                drag_data['y'] = event.ydata
                
                ax.figure.canvas.draw()

            # Function to end dragging
            def on_drag_end(event):
                drag_data['x'] = None
                drag_data['y'] = None

            # Connect the event handlers
            fig.canvas.mpl_connect('button_press_event', on_plot_click)
            fig.canvas.mpl_connect('scroll_event', on_scroll)
            fig.canvas.mpl_connect('button_press_event', on_drag_start)
            fig.canvas.mpl_connect('motion_notify_event', on_drag)
            fig.canvas.mpl_connect('button_release_event', on_drag_end)

            # Adjust aspect ratio and tight layout
            ax.set_aspect('equal')  # Fixed aspect ratio
            fig.tight_layout()

            # Display plot in the popup window using FigureCanvasTkAgg
            canvas = FigureCanvasTkAgg(fig, master=popup)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

            # Save function
            def save_plot():
                file_path = filedialog.asksaveasfilename(defaultextension=".jpg",
                                                         filetypes=[("JPEG files", "*.jpg"), ("All files", "*.*")])
                if file_path:
                    fig.savefig(file_path, format='jpeg')

            # Save button
            save_button = ttk.Button(popup, text="Save as JPEG", command=save_plot)
            save_button.pack(side=tk.LEFT, padx=10, pady=10)
    
            # Close button
            def close_popup():
                # Reset variables or any other necessary cleanup
                selected_files = []
                popup.destroy()

            close_button = ttk.Button(popup, text="Close", command=close_popup)
            close_button.pack(side=tk.RIGHT, padx=10, pady=10)
            
        except Exception as e:
            messagebox.showerror("Error", f"Error loading file: {e}")

def clear_scenario():
    scenario_var.set("")

def clear_variable():
    variable_var.set("")  

def clear_period():
    period_var.set("") 
    
# Create the main window
root = tk.Tk()
root.title("NZ Climate Explorer")

# Create a notebook for tabs
notebook = ttk.Notebook(root)
notebook.pack(expand=1, fill='both')

# Create the "Preview" tab
preview_frame = ttk.Frame(notebook)
notebook.add(preview_frame, text='Preview')

# Scenario dropdown
scenario_label = ttk.Label(preview_frame, text="Scenario:")
scenario_label.grid(row=0, column=0, padx=10, pady=5, sticky='w')
scenario_var = tk.StringVar(preview_frame)
scenario_dropdown = ttk.Combobox(preview_frame, textvariable=scenario_var, values=['historical', 'ssp126', 'ssp245', 'ssp370'], width=20)
scenario_dropdown.grid(row=0, column=1, padx=10, pady=5, sticky='w')
scenario_clear_button = ttk.Button(preview_frame, text="Clear", command=clear_scenario)
scenario_clear_button.grid(row=0, column=2, padx=10, pady=5, sticky='w')

# Variable dropdown (using climate variable names)
variable_label = ttk.Label(preview_frame, text="Climate Variable:")
variable_label.grid(row=1, column=0, padx=10, pady=5, sticky='w')
variable_var = tk.StringVar(preview_frame)
variable_dropdown = ttk.Combobox(preview_frame, textvariable=variable_var, values=index_df['Climate variable'].tolist(), width=50)
variable_dropdown.grid(row=1, column=1, padx=10, pady=5, sticky='w')
variable_clear_button = ttk.Button(preview_frame, text="Clear", command=clear_variable)
variable_clear_button.grid(row=1, column=2, padx=10, pady=5, sticky='w')

# Period dropdown
period_label = ttk.Label(preview_frame, text="Period:")
period_label.grid(row=2, column=0, padx=10, pady=5, sticky='w')
period_var = tk.StringVar(preview_frame)
period_dropdown = ttk.Combobox(preview_frame, textvariable=period_var, values=['bp1986-2005', 'bp1995-2014', 'fp2021-2040', 'fp2041-2060', 'fp2080-2099', 'wl1.5', 'wl2', 'wl3'], width=20)
period_dropdown.grid(row=2, column=1, padx=10, pady=5, sticky='w')
period_clear_button = ttk.Button(preview_frame, text="Clear", command=clear_period)
period_clear_button.grid(row=2, column=2, padx=10, pady=5, sticky='w')

# Filter button
filter_button = ttk.Button(preview_frame, text="Filter Files", command=filter_files)
filter_button.grid(row=3, column=0, columnspan=3, pady=10)

# Listbox to display filtered files
file_listbox = tk.Listbox(preview_frame, selectmode=tk.MULTIPLE, width=100)
file_listbox.grid(row=4, column=0, columnspan=3, padx=10, pady=5)

# Preview button
preview_button = ttk.Button(preview_frame, text="Preview", command=preview_data)
preview_button.grid(row=5, column=0, padx=10, pady=10)

# Create the "Compare" tab
compare_frame = ttk.Frame(notebook)
notebook.add(compare_frame, text='Compare')

# Variable dropdown (using climate variable names)
compare_variable_label = ttk.Label(compare_frame, text="Climate Variable:")
compare_variable_label.grid(row=0, column=0, padx=10, pady=5, sticky='w')
compare_variable_var = tk.StringVar(compare_frame)
compare_variable_dropdown = ttk.Combobox(compare_frame, textvariable=compare_variable_var, values=index_df['Climate variable'].tolist(), width=50)
compare_variable_dropdown.grid(row=0, column=1, padx=10, pady=5, sticky='w')
compare_variable_clear_button = ttk.Button(compare_frame, text="Clear", command=lambda: compare_variable_var.set(''))
compare_variable_clear_button.grid(row=0, column=2, padx=10, pady=5, sticky='w')

# Period dropdown (focusing on historical datasets)
compare_period_label = ttk.Label(compare_frame, text="Period:")
compare_period_label.grid(row=1, column=0, padx=10, pady=5, sticky='w')
compare_period_var = tk.StringVar(compare_frame)
compare_period_dropdown = ttk.Combobox(compare_frame, textvariable=compare_period_var, values=['bp1986-2005', 'bp1995-2014'], width=20)
compare_period_dropdown.grid(row=1, column=1, padx=10, pady=5, sticky='w')
compare_period_clear_button = ttk.Button(compare_frame, text="Clear", command=lambda: compare_period_var.set(''))
compare_period_clear_button.grid(row=1, column=2, padx=10, pady=5, sticky='w')

# Filter function for Compare tab
def filter_compare_files():
    selected_variable_name = compare_variable_var.get()
    selected_variable_id = variable_name_to_id.get(selected_variable_name, '')
    selected_period = compare_period_var.get()
    filtered_files = []
    for file_name in os.listdir('C:\\Users\\PKTAYLOR\\Documents\\RAW Datasets\\MfE\\netCDF_Files'):
            file_variable_id = file_name.split('_')[0]
            if (selected_variable_id == file_variable_id and 
                selected_period in file_name and 
                'historical' in file_name):
                filtered_files.append(file_name)

    # Update the listbox with filtered files
    compare_file_listbox.delete(0, tk.END)  # Clear current list
    for file_name in filtered_files:
        compare_file_listbox.insert(tk.END, file_name)

# Filter button for Compare tab
compare_filter_button = ttk.Button(compare_frame, text="Filter Files", command=filter_compare_files)
compare_filter_button.grid(row=2, column=0, columnspan=3, pady=10)

# Listbox to display selected files for comparison
compare_file_listbox = tk.Listbox(compare_frame, selectmode=tk.MULTIPLE, width=100)
compare_file_listbox.grid(row=3, column=0, columnspan=3, padx=10, pady=5)

def compare_files():
    global selected_files
    selected_file_index = compare_file_listbox.curselection()
    if len(selected_file_index) != 1:
        messagebox.showwarning("Selection Error", "Please select exactly one file from the list to compare.")
        return
    
    selected_file = compare_file_listbox.get(selected_file_index[0])
    selected_file_path = os.path.join('C:\\Users\\PKTAYLOR\\Documents\\RAW Datasets\\MfE\\netCDF_Files', selected_file)
    
    # Extract period, month, and variable from the base file name
    base_period = next((period for period in ['bp1986-2005', 'bp1995-2014'] if period in selected_file), None)
    base_month = next((month for month in ['ANN', 'DJF', 'MAM', 'JJA', 'SON'] if month in selected_file), None)
    base_variable = next((variable for variable in ['TNn', 'TN', 'TXx', 'TX25','TX30', 'TX', 'T', 'DTR', 'FD', 'GDD5',
                                                    'GDD10', 'CD18', 'HD18', 'PR', 'DD1mm','RR1mm', 'RR25mm', 'R99pVAL',
                                                    'PEDsrad', 'sfcWind', 'Wd10', 'Wd99pVAL', 'hurs', 'rsds'
                                                    ] if variable in selected_file), None)
    
    if not base_period or not base_month or not base_variable:
        messagebox.showerror("File Error", "Selected file does not contain valid period, month, or variable information.")
        return
    
    # Filter files for dropdown
    available_files = []
    for file_name in os.listdir('C:\\Users\\PKTAYLOR\\Documents\\RAW Datasets\\MfE\\netCDF_Files'):
        if ("historical" not in file_name and 
            base_period in file_name and 
            base_month in file_name and
            base_variable == file_name.split('_')[0]):  # Ensure the variable is exact
            available_files.append(file_name)
    
    try:
        dataset = nc.Dataset(selected_file_path, mode='r')
        variable_name = selected_file.split('_')[0]  # Extract variable ID from file name
        variable = dataset.variables[variable_name]
        data = variable[0, :, :]  # Extract a 2D slice (first time step)

        # Handling NaN values and coordinate grids
        lats = dataset.variables['latitude'][:]
        lons = dataset.variables['longitude'][:]
        values_array = data[:, :].data

        # Replace extreme outliers with NaN
        max_valid_value = 100000  # Example: adjust based on your data range
        values_array[values_array > max_valid_value] = np.nan

        # Create a new popup window
        popup = tk.Toplevel()
        popup.title(f'{variable.long_name}')
        
        # Plotting with fixed aspect ratio
        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.pcolormesh(lons, lats, values_array, cmap='viridis', shading='auto')
        ax.set_title(f'{variable.long_name} ({base_month}) - {base_period}')
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
            if event.button == 'down':
                scale_factor = 1.1
            elif event.button == 'up':
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
        
        # Dropdown for available files
        compare_file_label = ttk.Label(popup, text="Compare with:")
        compare_file_label.pack(side=tk.LEFT, padx=10, pady=10)

        compare_file_var = tk.StringVar(popup)
        compare_file_dropdown = ttk.Combobox(popup, textvariable=compare_file_var, values=available_files, width=50)
        compare_file_dropdown.pack(side=tk.LEFT, padx=10, pady=10)

        def compare_and_plot():
            compare_file_name = compare_file_var.get()
            if not compare_file_name:
                messagebox.showerror("File Error", "Please select a file to compare.")
                return

            compare_file_path = os.path.join('C:\\Users\\PKTAYLOR\\Documents\\RAW Datasets\\MfE\\netCDF_Files', compare_file_name)
            
            try:
                compare_dataset = nc.Dataset(compare_file_path, mode='r')
                compare_variable_name = compare_file_name.split('_')[0]  # Extract variable ID from compare file name
                compare_variable = compare_dataset.variables[compare_variable_name]
                compare_data = compare_variable[0, :, :]  # Extract a 2D slice (first time step)
                compare_period = next((cperiod for cperiod in ['fp2021-2040', 'fp2041-2060', 'fp2080-2099', 'wl1.5', 'wl2', 'wl3'] if cperiod in compare_file_name), None)
                compare_scenario = next((scenario for scenario in ['ssp126', 'ssp245', 'ssp370'] if scenario in compare_file_name), None)

                compare_values_array = compare_data[:, :].data

                # Replace extreme outliers with NaN in compare data
                compare_values_array[compare_values_array > max_valid_value] = np.nan

                # Combine data based on variable type
                if base_variable in ['PR', 'R99pVAL', 'sfcWind', 'Wd99pVAL', 'hurs']:
                    combined_data = values_array * (1 + compare_values_array / 100.0)
                else:
                    combined_data = values_array + compare_values_array

                # Determine vmin and vmax for color scale synchronization
                vmin = min(np.nanmin(values_array), np.nanmin(combined_data))
                vmax = max(np.nanmax(values_array), np.nanmax(combined_data))

                # Create a new popup window for combined data plot
                combined_popup = tk.Toplevel()
                combined_popup.title(f'{variable.long_name} + {compare_scenario}, {compare_period}')
                
                # Plotting with fixed aspect ratio and synchronized color scale
                combined_fig, combined_ax = plt.subplots(figsize=(8, 6))
                combined_im = combined_ax.pcolormesh(lons, lats, combined_data, cmap='viridis', shading='auto', vmin=vmin, vmax=vmax)
                combined_ax.set_title(f'{variable.long_name} + {compare_scenario}, {compare_period}')
                combined_ax.set_xlabel('Longitude')
                combined_ax.set_ylabel('Latitude')
                combined_cbar = combined_fig.colorbar(combined_im, ax=combined_ax)
                combined_cbar.set_label(variable.units)

                # Create a new popup window for updated base data plot
                base_updated_popup = tk.Toplevel()
                base_updated_popup.title(f'{variable.long_name}')

                # Plotting updated base data with synchronized color scale
                base_updated_fig, base_updated_ax = plt.subplots(figsize=(8, 6))
                base_updated_im = base_updated_ax.pcolormesh(lons, lats, values_array, cmap='viridis', shading='auto', vmin=vmin, vmax=vmax)
                base_updated_ax.set_title(f'{variable.long_name} ({base_month}) - {base_period}')
                base_updated_ax.set_xlabel('Longitude')
                base_updated_ax.set_ylabel('Latitude')
                base_updated_cbar = base_updated_fig.colorbar(base_updated_im, ax=base_updated_ax)
                base_updated_cbar.set_label(variable.units)

                # Display combined plot in the popup window using FigureCanvasTkAgg
                combined_canvas = FigureCanvasTkAgg(combined_fig, master=combined_popup)
                combined_canvas.draw()
                combined_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

                # Display updated base plot in the popup window using FigureCanvasTkAgg
                base_updated_canvas = FigureCanvasTkAgg(base_updated_fig, master=base_updated_popup)
                base_updated_canvas.draw()
                base_updated_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
                # Adjust aspect ratio and tight layout
                combined_ax.set_aspect('equal')  # Fixed aspect ratio
                base_updated_ax.set_aspect('equal')
                combined_fig.tight_layout()
                base_updated_fig.tight_layout()
                
                # Save combined plot as JPEG
                def save_combined_plot():
                    file_path = filedialog.asksaveasfilename(defaultextension=".jpeg", filetypes=[("JPEG files", "*.jpeg"), ("All files", "*.*")])
                    if file_path:
                        combined_fig.savefig(file_path)

                # Save base updated plot as JPEG
                def save_base_updated_plot():
                    file_path = filedialog.asksaveasfilename(defaultextension=".jpeg", filetypes=[("JPEG files", "*.jpeg"), ("All files", "*.*")])
                    if file_path:
                        base_updated_fig.savefig(file_path)

                # Save Plot button for combined plot
                save_combined_button = ttk.Button(combined_popup, text="Save Plot", command=save_combined_plot)
                save_combined_button.pack(side=tk.RIGHT, padx=10, pady=10)

                # Save Plot button for base updated plot
                save_base_updated_button = ttk.Button(base_updated_popup, text="Save Plot", command=save_base_updated_plot)
                save_base_updated_button.pack(side=tk.RIGHT, padx=10, pady=10)
                
                # Export button for combined plot
                export_combined_button = ttk.Button(combined_popup, text="Export to ArcGIS")
                export_combined_button.pack(side=tk.RIGHT, padx=10, pady=10)

                # Export button for base updated plot
                export_base_updated_button = ttk.Button(base_updated_popup, text="Export to ArcGIS")
                export_base_updated_button.pack(side=tk.RIGHT, padx=10, pady=10)

            except Exception as e:
                messagebox.showerror("Error", f"Error comparing file: {e}")

        # Compare button to initiate comparison
        compare_button = ttk.Button(popup, text="Compare", command=compare_and_plot)
        compare_button.pack(side=tk.RIGHT, padx=10, pady=10)

        # Close button
        def close_popup():
            compare_file_var.set("")
            popup.destroy()
        
        close_button = ttk.Button(popup, text="Close", command=close_popup)
        close_button.pack(side=tk.RIGHT, padx=10, pady=10)
        
    except Exception as e:
        messagebox.showerror("Error", f"Error loading file: {e}")


# Compare button 
compare_button = ttk.Button(compare_frame, text="Select Base Case", command=compare_files)
compare_button.grid(row=4, column=0, columnspan=3, pady=10)

root.mainloop()
