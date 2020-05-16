#Step Heating Experiment Error Propagation

# This program reads data from step heating experiments as 1 or more .csv files,
    # completes diffusivity error propagation calculations, and displays figures.
    # Figures are compiled and may be exported as pdf file given a filepath for the
    # output. The formulas and structure for this program are based on an existing
    # spreadsheet (Ginster, U., & Reiners, P. W., 2018).

# Author: Research assistant for Dr. Marissa Tremblay
    # Isabella Zuffoletti- Purdue Univerity Geology & Geophysics undergraduate, 2022

# Source: Ginster, U., and Peter W. Reiners. "Error Propagation in the Derivation of Noble Gas Diffusion Parameters for Minerals From Step Heating Experiments." Geochemistry, Geophysics, Geosystems 19.10 (2018): 3706-3720.
    # https://doi.org/10.1029/2018GC007531

import tkinter as tk
from matplotlib.backends.backend_pdf import PdfPages
import EVA
import OLS
import WLS
import Arrhenius

def export_csv():
    pass
def Summary_button_fxn(geometry, input_filepath, estimate):
    pass

# All figures saved in figs are exported to pdf
def export_pdf(figs, output_filepath):
    with PdfPages(output_filepath + "/output.pdf") as pdf:
        for fig in figs:
            pdf.savefig(fig)

# Calls Arrhenius.py to create figures
    # Figures added to "figs"
def Diffusivity_button_fxn(geometry, input_filepath):
    if geometry == 1:
        geometry = "infinite sheet"
    elif geometry == 2:
        geometry = "sphere"
    elif geometry == 3:
        geometry = 'infinite cylinder'
    else:
        print('Error: Select geometry')
    arrhenius1, newfigs = Arrhenius.general_arrhenius_calcs(geometry, input_filepath)
    global figs
    for fig in newfigs:
        figs.append(fig)
        fig.show()


# Calls Arrhenius.py & OLS.py to create figures
    # OLS tables and plot added to "figs"
def OLS_button_fxn(geometry, input_filepath):
    if geometry == 1:
        geometry = "infinite sheet"
    elif geometry == 2:
        geometry = "sphere"
    elif geometry == 3:
        geometry = 'infinite cylinder'
    else:
        print('Error: Select geometry')
    arrhenius1, newfigs = Arrhenius.general_arrhenius_calcs(geometry, input_filepath)
    newfigs = OLS.OLS(arrhenius1, geometry, input_filepath)
    global figs
    for fig in newfigs:
        figs.append(fig)
        fig.show()


# Calls Arrhenius.py & WLS.py to create figures
    # WLS tables added to "figs"
def WLS_button_fxn(geometry, input_filepath):
    if geometry == 1:
        geometry = "infinite sheet"
    elif geometry == 2:
        geometry = "sphere"
    elif geometry == 3:
        geometry = 'infinite cylinder'
    else:
        print('Error: Select geometry')
    arrhenius1, newfigs = Arrhenius.general_arrhenius_calcs(geometry, input_filepath)
    newfigs = WLS.WLS(arrhenius1, geometry, input_filepath)
    global figs
    for fig in newfigs:
        figs.append(fig)
        fig.show()


# Calls Arrhenius.py & EVA.py to create figures
    # EVA tables added to "figs"
def EVA_button_fxn(geometry, input_filepath):
    if geometry == 1:
        geometry = "infinite sheet"
    elif geometry == 2:
        geometry = "sphere"
    elif geometry == 3:
        geometry = 'infinite cylinder'
    else:
        print('Error: Select geometry')
    arrhenius1, newfigs = Arrhenius.general_arrhenius_calcs(geometry, input_filepath)
    newfigs = EVA.EVA(arrhenius1, geometry, input_filepath)
    global figs
    for fig in newfigs:
        figs.append(fig)
        fig.show()


# Function creates GUI w/ tkinter
    # Buttons call fxns defined above
def init_GUI():

    # Create GUI window
    root = tk.Tk()
    canvas1 = tk.Canvas(root, width=400, height=375)
    canvas1.pack()
    title_text = tk.Label(root, text='Error Propagation')
    title_text.config(font=("Courier", 20))
    canvas1.create_window(200, 15, window=title_text)


    # -------------------- File path entry box formatting --------------------
    input_filepath_entry = tk.Entry()
    input_filepath_entry.place(x=150, y=30)
    output_filepath_entry = tk.Entry()
    output_filepath_entry.place(x=150, y=60)
    input_text = tk.Label(root, text='Input file path')
    input_text.config(font=("Courier", 14))
    input_text.place(x=14, y=34)
    output_text = tk.Label(root, text='Output file path')
    output_text.config(font=("Courier", 14))
    output_text.place(x=12, y=64)

    # -------------------- Geometry Selection formatting --------------------
    geometry_entry = tk.IntVar()
    geometry_selection_label = tk.Label(root, text = "Select geometry:")
    geometry_selection_label.config(font=("Courier", 14))
    geometry_selection_label.place(x=10, y=90)
    g1 = tk.Radiobutton(root, text = "Infinite sheet", variable=geometry_entry, value = 1)
    g1.config(font=("Courier", 13))
    g1.place(x=20, y=115)
    g2 = tk.Radiobutton(root, text = "Sphere", variable=geometry_entry, value = 2)
    g2.config(font=("Courier", 13))
    g2.place(x=160, y=115)
    g3 = tk.Radiobutton(root, text = "Infinite Cylinder", variable=geometry_entry, value = 3)
    g3.config(font=("Courier", 13))
    g3.place(x=240, y=115)

    # -------------------- Button formatting --------------------

    Diffusivity_button = tk.Button(root, text='Diffusivity',command= lambda:Diffusivity_button_fxn(geometry_entry.get(), input_filepath_entry.get()), fg='purple', height = 3, width = 10)
    Diffusivity_button.config(font=("Courier", 15))
    Diffusivity_button.place(x=275, y=225)
    OLS_button = tk.Button(root, text='OLS',command= lambda:OLS_button_fxn(geometry_entry.get(), input_filepath_entry.get()), fg='blue', height = 3, width = 10)
    OLS_button.config(font=("Courier", 15))
    OLS_button.place(x=25, y=150)
    WLS_button = tk.Button(root, text='WLS',command= lambda:WLS_button_fxn(geometry_entry.get(), input_filepath_entry.get()), fg='dark green', height = 3, width = 10)
    WLS_button.config(font=("Courier", 15))
    WLS_button.place(x=150, y=150)
    EVA_button = tk.Button(root, text='EVA',command= lambda:EVA_button_fxn(geometry_entry.get(), input_filepath_entry.get()), fg='maroon', height = 3, width = 10)
    EVA_button.config(font=("Courier", 15))
    EVA_button.place(x=275, y=150)
    Summary_button = tk.Button(root, text='Summary',command= lambda:Summary_button_fxn(geometry_entry.get(), input_filepath_entry.get(), float(user_estimate_entry.get())), fg='purple', height = 3, width = 10)
    Summary_button.config(font=("Courier", 15))
    Summary_button.place(x=150, y=225)
    CSV_button = tk.Button(root, text='Export CSV',command= lambda:export_csv(), fg='navy', height = 3, width = 10)
    CSV_button.config(font=("Courier", 15))
    CSV_button.place(x=87.5, y=300)
    pdf_button = tk.Button(root, text='Export PDF', command=lambda: export_pdf(figs, output_filepath_entry.get()),fg='navy', height=3, width=10)
    pdf_button.config(font=("Courier", 15))
    pdf_button.place(x=212.5, y=300)

    # -------------------- Estimate entry box and label formatting --------------------
    estimate_text = tk.Label(root, text='Tc Estimate:')
    estimate_text.config(font=("Courier", 14))
    estimate_text.place(x=20, y=225)
    user_estimate_entry = tk.Entry(root, width=10)
    user_estimate_entry.place(x=20, y=250)

def initialize_vars():
    global figs
    figs = []

def main():
    initialize_vars()
    init_GUI()
    tk.mainloop()

if __name__ == '__main__':
    main()

