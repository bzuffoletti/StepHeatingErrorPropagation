import pandas as pd
import numpy as np
import warnings
import matplotlib.pyplot as plt

# Initialize variables
global domain_size
domain_size = 27
global domain_size_uncert
domain_size_uncert = 2

def SUMSQ(x):
    return (sum([i**2 for i in x]))

# Sheet geometry specific equations
def sheet_arrhenius(arrhenius_legend,arrhenius1):
    arrhenius_legend['S'] = '(D/a^2)n  (s¯¹)'
    arrhenius1.loc[arrhenius1.index[0], 'S'] = np.where(arrhenius1.loc[arrhenius1.index[0], 'P'] < 0.6, ((np.pi / 4) * (arrhenius1.loc[arrhenius1.index[0], 'P'] ** 2) / arrhenius1.loc[arrhenius1.index[0], 'H']), (((-4) / (np.pi ** 2)) * np.log((1 - arrhenius1.loc[arrhenius1.index[0], 'P']) * (np.pi ** (2 / 8)))) / arrhenius1.loc[arrhenius1.index[0], 'H'])
    for idx in range(len(arrhenius1)):
        if idx == 0:
            pass
        else:
            arrhenius1.loc[arrhenius1.index[idx], 'S'] = np.where(arrhenius1.loc[arrhenius1.index[idx], 'P'] < 0.6,(np.pi / 4) * (arrhenius1.loc[arrhenius1.index[idx], 'P'] ** 2 - (arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] ** 2)) /arrhenius1.loc[arrhenius1.index[idx], 'H'], (((-4) / (np.pi) ** 2) * np.log((1 -arrhenius1.loc[arrhenius1.index[idx], 'P']) / (1 -arrhenius1.loc[arrhenius1.index[idx - 1], 'P']))) /arrhenius1.loc[arrhenius1.index[idx], 'H'])

    arrhenius_legend['T'] = 'δ(D/a2)(Δtn)  (s¯¹)'
    arrhenius1.loc[arrhenius1.index[0], 'T'] = np.where(arrhenius1.loc[arrhenius1.index[0], 'P'] < 0.6, np.absolute(np.pi * (arrhenius1.loc[arrhenius1.index[0], 'P'] ** 2) * (-1) * arrhenius1.loc[arrhenius1.index[0], 'I'] / 4 / (arrhenius1.loc[arrhenius1.index[0], 'H'] ** 2)), np.absolute((-4) * np.log((1 - arrhenius1.loc[arrhenius1.index[0], 'P']) * (np.pi ** (2 / 8)) * (-1) * arrhenius1.loc[arrhenius1.index[0], 'I'] / (np.pi ** 2) / (arrhenius1.loc[arrhenius1.index[0], 'H'] ** 2))))
    for idx in range(len(arrhenius1)):
        if idx == 0:
            pass
        else:
            arrhenius1.loc[arrhenius1.index[idx], 'T'] = np.where(arrhenius1.loc[arrhenius1.index[idx], 'P'] < 0.6, np.absolute(np.pi * (arrhenius1.loc[arrhenius1.index[idx], 'P'] ** 2 -arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] ** 2) * (-1) * arrhenius1.loc[arrhenius1.index[idx], 'I'] / 4 /arrhenius1.loc[arrhenius1.index[idx], 'H'] ** 2),np.absolute((-4) * np.log((1 - arrhenius1.loc[arrhenius1.index[idx], 'P']) / (1 -arrhenius1.loc[arrhenius1.index[idx - 1], 'P'])) * (-1) * arrhenius1.loc[arrhenius1.index[idx], 'I'] / np.pi ** 2 /arrhenius1.loc[arrhenius1.index[idx], 'H'] ** 2))

    arrhenius_legend['U'] = 'δ(D/a^2)(Xi)  (s¯¹)'
    arrhenius1.loc[arrhenius1.index[0], 'U'] = np.where(arrhenius1.loc[arrhenius1.index[0], 'P'] < 0.6, (np.pi * arrhenius1.loc[arrhenius1.index[0], 'P'] / (2 * arrhenius1.loc[arrhenius1.index[0], 'H'] * TotalXT)) * np.sqrt(((1 - arrhenius1.loc[arrhenius1.index[0], 'P']) * arrhenius1.loc[arrhenius1.index[0], 'K']) ** 2 +arrhenius1.loc[arrhenius1.index[0], 'P'] ** 2 * (SUMSQ(arrhenius1.loc[arrhenius1.index != 0, 'K']) + Final_Degassing_uncert ** 2)), (4 / (np.pi ** 2 * arrhenius1.loc[arrhenius1.index[0], 'H'] * TotalXT)) * np.sqrt(arrhenius1.loc[arrhenius1.index[0], 'K'] ** 2 + (arrhenius1.loc[arrhenius1.index[0], 'P'] / (1 - arrhenius1.loc[arrhenius1.index[0], 'P'])) ** 2 * (SUMSQ(arrhenius1.loc[arrhenius1.index != idx, 'K']) + Final_Degassing_uncert ** 2)))
    for idx in range(len(arrhenius1)):
        if idx == 0:
            pass
        else:
            arrhenius1.loc[arrhenius1.index[idx], 'U'] = np.where(arrhenius1.loc[arrhenius1.index[idx], 'P'] < 0.6, (np.pi / (2 * arrhenius1.loc[arrhenius1.index[idx], 'H'] * TotalXT)) * np.sqrt(((arrhenius1.loc[arrhenius1.index[idx], 'P'] * (1 - arrhenius1.loc[arrhenius1.index[idx], 'P'])) - (arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] * (1 - arrhenius1.loc[arrhenius1.index[idx - 1], 'P']))) ** 2 * SUMSQ(arrhenius1.loc[arrhenius1.index < idx, 'K']) + (arrhenius1.loc[arrhenius1.index[idx], 'P'] * (1 - arrhenius1.loc[arrhenius1.index[idx], 'P']) + (arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] ** 2)) ** 2 * (arrhenius1.loc[arrhenius1.index[idx], 'K'] ** 2) + (arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] ** 2 - arrhenius1.loc[arrhenius1.index[idx], 'P'] ** 2) ** 2 * (SUMSQ(arrhenius1.loc[arrhenius1.index > idx, 'K']) + Final_Degassing_uncert ** 2)), ((4 / np.pi) ** 2 * arrhenius1.loc[arrhenius1.index[idx], 'H'] * TotalXT * np.sqrt((1 + arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] / (1 - arrhenius1.loc[arrhenius1.index[idx - 1], 'P'])) ** 2 * arrhenius1.loc[arrhenius1.index[idx], 'K'] ** 2 + (arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] / (1 - arrhenius1.loc[arrhenius1.index[idx - 1], 'P']) - arrhenius1.loc[arrhenius1.index[idx], 'P'] / (1 - arrhenius1.loc[arrhenius1.index[idx], 'P'])) ** 2 * SUMSQ(arrhenius1.loc[arrhenius1.index > idx, 'K']) + Final_Degassing_uncert ** 2)))
    return arrhenius_legend,arrhenius1

# Sphere geometry specific equations
def sphere_arrhenius(arrhenius_legend,arrhenius1):
    arrhenius_legend['S'] = '(D/a^2)n  (s¯¹)'
    arrhenius1.loc[arrhenius1.index[0], 'S'] = np.where(arrhenius1.loc[arrhenius1.index[0], 'P'] < 0.85, (2 * np.pi - np.pi ** 2 * arrhenius1.loc[arrhenius1.index[0], 'P'] / 3 - 2 * np.pi * np.sqrt(1 - np.pi * arrhenius1.loc[arrhenius1.index[0], 'P'] / 3)) / np.pi ** 2 / arrhenius1.loc[arrhenius1.index[0], 'H'], -1 * np.log((1 - arrhenius1.loc[arrhenius1.index[0], 'P']) * (np.pi ** (2 / 6))) / arrhenius1.loc[arrhenius1.index[0], 'H'])
    for idx in range(len(arrhenius1)):
        if idx == 0:
            pass
        # np.log() is "LN()", while np.log10() is base 10
        else:
            arrhenius1.loc[arrhenius1.index[idx], 'S'] = np.where(arrhenius1.loc[arrhenius1.index[idx], 'P'] < 0.85, ((-1 * np.pi ** 2 * (arrhenius1.loc[arrhenius1.index[idx], 'P'] -arrhenius1.loc[arrhenius1.index[idx - 1], 'P']) / 3 - 2 * np.pi * (np.sqrt(1 - np.pi *arrhenius1.loc[arrhenius1.index[idx], 'P'] / 3) - np.sqrt(1 - np.pi *arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] / 3))) / np.pi ** 2 /arrhenius1.loc[arrhenius1.index[idx], 'H']),((-1) * np.log((1 - arrhenius1.loc[arrhenius1.index[idx], 'P']) / (1 - arrhenius1.loc[arrhenius1.index[idx - 1], 'P'])) / np.pi ** 2 /arrhenius1.loc[arrhenius1.index[idx], 'H']))

    arrhenius_legend['T'] = 'δ(D/a2)(Δtn)  (s¯¹)'
    arrhenius1.loc[arrhenius1.index[0], 'T'] = np.where(arrhenius1.loc[arrhenius1.index[0], 'P'] < 0.85, np.absolute((2 / np.pi -arrhenius1.loc[arrhenius1.index[0], 'P'] / 3 - 2 * np.sqrt(1 - np.pi *arrhenius1.loc[arrhenius1.index[0], 'P'] / 3) / np.pi) * -1 *arrhenius1.loc[arrhenius1.index[0], 'I'] /arrhenius1.loc[arrhenius1.index[0], 'H'] ** 2),np.absolute((-1) * np.log((1 - arrhenius1.loc[arrhenius1.index[0], 'P']) * np.pi ** (2 / 6) * (-1) * arrhenius1.loc[arrhenius1.index[0], 'I'] / np.pi ** 2 / arrhenius1.loc[arrhenius1.index[0], 'H'] ** 2)))
    for idx in range(len(arrhenius1)):
        if idx == 0:
            pass
        else:
            arrhenius1.loc[arrhenius1.index[idx], 'T'] = np.where(arrhenius1.loc[arrhenius1.index[idx], 'P'] < 0.85, np.absolute(((-1)*(arrhenius1.loc[arrhenius1.index[idx], 'P']-arrhenius1.loc[arrhenius1.index[idx-1], 'P'])/3-2*(np.sqrt(1-np.pi*arrhenius1.loc[arrhenius1.index[idx], 'P']/3)-np.sqrt(1-np.pi*arrhenius1.loc[arrhenius1.index[idx-1], 'P']/3))/np.pi)*(-1)*arrhenius1.loc[arrhenius1.index[idx], 'I']/arrhenius1.loc[arrhenius1.index[idx], 'H']**2), np.absolute((-1)*np.log((1-arrhenius1.loc[arrhenius1.index[idx], 'P'])/(1-arrhenius1.loc[arrhenius1.index[idx-1], 'P']))*(-1)*arrhenius1.loc[arrhenius1.index[idx], 'I']/np.pi**2/arrhenius1.loc[arrhenius1.index[idx], 'H']**2))

    arrhenius_legend['U'] = 'δ(D/a^2)(Xi)  (s¯¹)'
    arrhenius1.loc[arrhenius1.index[0], 'U'] = np.where(arrhenius1.loc[arrhenius1.index[0], 'P'] <= 0.85, (1 / (3 * arrhenius1.loc[arrhenius1.index[0], 'H'] * TotalXT)) * ((1 / np.sqrt(1 - np.pi * arrhenius1.loc[arrhenius1.index[0], 'P'] / 3)) - 1) * np.sqrt(((1 - arrhenius1.loc[arrhenius1.index[0], 'P']) * arrhenius1.loc[arrhenius1.index[0], 'K']) ** 2 + (-arrhenius1.loc[arrhenius1.index[0], 'P']) ** 2 * (SUMSQ(arrhenius1.loc[arrhenius1.index != 0, 'K']) + Final_Degassing_uncert ** 2)),(1 / (np.pi ** 2 * arrhenius1.loc[arrhenius1.index[0], 'H'] * TotalXT)) * np.sqrt(arrhenius1.loc[arrhenius1.index[0], 'K'] ** 2 + (arrhenius1.loc[arrhenius1.index[0], 'P'] / (1 - arrhenius1.loc[arrhenius1.index[0], 'P'])) ** 2 * (SUMSQ(arrhenius1.loc[arrhenius1.index != 0, 'K']) + Final_Degassing_uncert ** 2)))
    for idx in range(len(arrhenius1)):
        if idx == 0:
            pass
        else:
            arrhenius1.loc[arrhenius1.index[idx], 'U'] = np.where(arrhenius1.loc[arrhenius1.index[idx], 'P'] <= 0.85, (1 / (3 * arrhenius1.loc[arrhenius1.index[idx], 'H'] * TotalXT)) * np.sqrt(((1 - arrhenius1.loc[arrhenius1.index[idx], 'P']) * (1 / np.sqrt(1 - np.pi * arrhenius1.loc[arrhenius1.index[idx], 'P'] / 3) - 1) + (1 - arrhenius1.loc[arrhenius1.index[idx - 1], 'P']) * (1 - 1 / np.sqrt(1 - np.pi * arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] / 3))) ** 2 * SUMSQ(arrhenius1.loc[arrhenius1.index < idx, 'K']) + ((1 - arrhenius1.loc[arrhenius1.index[idx], 'P']) * (1 / np.sqrt(1 - np.pi * arrhenius1.loc[arrhenius1.index[idx], 'P'] / 3) - 1) - arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] * (1 - 1 / np.sqrt(1 - np.pi * arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] / 3))) ** 2 * arrhenius1.loc[arrhenius1.index[idx], 'K'] ** 2 + (arrhenius1.loc[arrhenius1.index[idx], 'P'] * (1 - 1 / np.sqrt(1 - np.pi *arrhenius1.loc[arrhenius1.index[idx], 'P'] / 3)) -arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] * (1 - 1 / np.sqrt(1 - np.pi *arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] / 3))) ** 2 * (SUMSQ(arrhenius1.loc[arrhenius1.index > idx, 'K']) + Final_Degassing_uncert ** 2)),(1 / (np.pi ** 2 * arrhenius1.loc[arrhenius1.index[idx], 'H'] * TotalXT)) * np.sqrt((1 + arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] / (1 - arrhenius1.loc[arrhenius1.index[idx - 1], 'P'])) ** 2 *arrhenius1.loc[arrhenius1.index[idx], 'K'] ** 2 + (arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] / (1 - arrhenius1.loc[arrhenius1.index[idx - 1], 'P']) -arrhenius1.loc[arrhenius1.index[idx], 'P'] / (1 - arrhenius1.loc[arrhenius1.index[idx], 'P'])) ** 2 * SUMSQ(arrhenius1.loc[arrhenius1.index > idx, 'K']) + Final_Degassing_uncert ** 2))
    return arrhenius_legend, arrhenius1

# Cylinder geometry specific equations
def cylinder_arrhenius(arrhenius_legend,arrhenius1):
    arrhenius_legend['S'] = '(D/a^2)n  (s¯¹)'
    arrhenius1.loc[arrhenius1.index[0], 'S'] = np.where(arrhenius1.loc[arrhenius1.index[0], 'P'] < 0.6, ((np.pi / 4) * (arrhenius1.loc[arrhenius1.index[0], 'P'] ** 2) / arrhenius1.loc[arrhenius1.index[0], 'H']), (((-4) / (np.pi ** 2)) * np.log((1 - arrhenius1.loc[arrhenius1.index[0], 'P']) * (np.pi ** (2 / 8))) / arrhenius1.loc[arrhenius1.index[0], 'H']))
    for idx in range(len(arrhenius1)):
        if idx == 0:
            pass
        else:
            arrhenius1.loc[arrhenius1.index[idx], 'S'] = np.where(arrhenius1.loc[arrhenius1.index[idx], 'P'] > 0.6, ((-1) / (5.78 * arrhenius1.loc[arrhenius1.index[idx], 'H'])) * np.log((1 - arrhenius1.loc[arrhenius1.index[idx], 'P']) / (1 - arrhenius1.loc[arrhenius1.index[idx - 1], 'P'])), (1 / (np.pi ** 2 * arrhenius1.loc[arrhenius1.index[idx], 'H'])) * (-1 * np.pi ** 2 * (arrhenius1.loc[arrhenius1.index[idx], 'P'] - arrhenius1.loc[arrhenius1.index[idx - 1], 'P']) - 4 * np.pi * (np.sqrt(4 - np.pi * arrhenius1.loc[arrhenius1.index[idx], 'P']) - np.sqrt(4 - np.pi * arrhenius1.loc[arrhenius1.index[idx - 1], 'P']))))

    arrhenius_legend['T'] = 'δ(D/a2)(Δtn)  (s¯¹)'
    arrhenius1.loc[arrhenius1.index[0], 'T'] = np.where(arrhenius1.loc[arrhenius1.index[0], 'P'] <= .85, np.absolute((2 / np.pi -arrhenius1.loc[arrhenius1.index[0], 'P'] / 3 - 2 * np.sqrt(1 - np.pi *arrhenius1.loc[arrhenius1.index[0], 'P'] / 3) / np.pi) * (-1) *arrhenius1.loc[arrhenius1.index[0], 'I'] /arrhenius1.loc[arrhenius1.index[0], 'H'] ** 2),np.absolute((-1) * np.log((1 - arrhenius1.loc[arrhenius1.index[0], 'P']) * np.pi ** 2 / 6) * (-1) *arrhenius1.loc[arrhenius1.index[0], 'I'] / np.pi ** 2 /arrhenius1.loc[arrhenius1.index[0], 'H'] ** 2))
    for idx in range(len(arrhenius1)):
        if idx == 0:
            pass
        else:
            arrhenius1.loc[arrhenius1.index[idx], 'T'] = np.where(arrhenius1.loc[arrhenius1.index[idx], 'P'] <= 0.85,np.absolute(((-1) * (arrhenius1.loc[arrhenius1.index[idx], 'P'] -arrhenius1.loc[arrhenius1.index[idx - 1], 'P']) / 3 - 2 * (np.sqrt(1 - np.pi *arrhenius1.loc[arrhenius1.index[idx], 'P'] / 3) - np.sqrt(1 - np.pi * arrhenius1.loc[arrhenius1.index[idx - 1], 'P'] / 3)) / np.pi) * (-1) * arrhenius1.loc[arrhenius1.index[idx], 'I'] /arrhenius1.loc[arrhenius1.index[idx], 'H'] ** 2),np.absolute((-1) * np.log((1 - arrhenius1.loc[arrhenius1.index[idx], 'P']) / (1 -arrhenius1.loc[arrhenius1.index[idx - 1], 'P'])) * (-1) * arrhenius1.loc[arrhenius1.index[idx], 'I'] / np.pi ** 2 /arrhenius1.loc[arrhenius1.index[idx], 'H'] ** 2))

    arrhenius_legend['U'] = 'δ(D/a^2)(Xi)  (s¯¹)'
    arrhenius1.loc[arrhenius1.index[0], 'U'] = np.where(arrhenius1.loc[arrhenius1.index[0], 'P'] > .6, (1 / (5.78 * arrhenius1.loc[arrhenius1.index[0], 'H'] * TotalXT)) * np.sqrt(arrhenius1.loc[arrhenius1.index[0], 'K'] ** 2 + (arrhenius1.loc[arrhenius1.index[0], 'P'] / (1 - arrhenius1.loc[arrhenius1.index[0], 'P'])) ** 2 * (SUMSQ(arrhenius1.loc[arrhenius1.index > idx, 'K']) + Final_Degassing_uncert ** 2)),(1 / (arrhenius1.loc[arrhenius1.index[0], 'H'] * TotalXT)) * (2 / np.sqrt(4 - np.pi * arrhenius1.loc[arrhenius1.index[0], 'P']) - 1) * np.sqrt(((1 -arrhenius1.loc[arrhenius1.index[0], 'P']) *arrhenius1.loc[arrhenius1.index[0], 'K']) ** 2 + (-arrhenius1.loc[arrhenius1.index[0], 'P']) ** 2 * (SUMSQ(arrhenius1.loc[arrhenius1.index > idx, 'K']) + Final_Degassing_uncert ** 2)))
    for idx in range(len(arrhenius1)):
        if idx == 0:
            pass
        else:
            arrhenius1.loc[arrhenius1.index[idx], 'U'] = np.where(arrhenius1.loc[arrhenius1.index[idx], 'P'] > .6,(1/(5.78*arrhenius1.loc[arrhenius1.index[idx], 'H']*TotalXT))*np.sqrt((1+arrhenius1.loc[arrhenius1.index[idx-1], 'P']/(1-arrhenius1.loc[arrhenius1.index[idx-1], 'P']))**2*arrhenius1.loc[arrhenius1.index[idx], 'K']**2 + (arrhenius1.loc[arrhenius1.index[idx-1], 'P']/(1-arrhenius1.loc[arrhenius1.index[idx-1], 'P'])-arrhenius1.loc[arrhenius1.index[idx], 'P']/(1-arrhenius1.loc[arrhenius1.index[idx], 'P']))**2*(SUMSQ(arrhenius1.loc[arrhenius1.index > idx, 'K']) + Final_Degassing_uncert ** 2)), (1/(arrhenius1.loc[arrhenius1.index[idx], 'H']*TotalXT))*np.sqrt(((1-arrhenius1.loc[arrhenius1.index[idx], 'P'])*(2/np.sqrt(4-np.pi*arrhenius1.loc[arrhenius1.index[idx], 'P'])-1) + (1-arrhenius1.loc[arrhenius1.index[idx-1], 'P'])*(1-2/np.sqrt(4-np.pi*arrhenius1.loc[arrhenius1.index[idx-1], 'P'])))**2*SUMSQ(arrhenius1.loc[arrhenius1.index < idx, 'K']) + ((1-arrhenius1.loc[arrhenius1.index[idx], 'P'])*(2/np.sqrt(4-np.pi*arrhenius1.loc[arrhenius1.index[idx], 'P'])-1) + (-arrhenius1.loc[arrhenius1.index[idx-1], 'P'])*(1-2/np.sqrt(4-np.pi*arrhenius1.loc[arrhenius1.index[idx-1], 'P'])))**2*arrhenius1.loc[arrhenius1.index[idx], 'K']**2 + ((-arrhenius1.loc[arrhenius1.index[idx], 'P'])*(2/np.sqrt(4-np.pi*arrhenius1.loc[arrhenius1.index[idx], 'P'])-1) + (-arrhenius1.loc[arrhenius1.index[idx-1], 'P'])*(1-2/np.sqrt(4-np.pi*arrhenius1.loc[arrhenius1.index[idx-1], 'P'])))**2*(SUMSQ(arrhenius1.loc[arrhenius1.index > idx, 'K']) + Final_Degassing_uncert ** 2)))
    return arrhenius_legend, arrhenius1

# Fxn reads csv given file path and initialized dataframe
def parse_data(file_path):

    df1 = pd.read_csv(file_path)
    df1.columns = ['A', 'C', 'D', 'H', 'J', 'K']
    col_headers = [ 'Step', 'T (ºC)', 'T•1σ (ºC)', 'Δt (s)', 'Xn (pcc)','δXn•1σ (pcc)']
    df1_legend = pd.DataFrame(columns = list(df1.columns))
    df1_legend.loc[0] = col_headers
    return df1, df1_legend

# This fxn formats and rounds a dataframe to for output in csv files or pdf tables
def df_round(arrhenius_dataframe):

    # If final degassing step, add last row back into dataframe for pdf output
    # This last row is not calculated until the pdf output is formatted in order
    # to simplify indexing into dataframe during calculations
    global Final_Degassing_step
    global Final_Degassing
    global Final_Degassing_uncert
    if Final_Degassing_step == True:
        df1_finaldegassing = pd.DataFrame(np.repeat("", len(arrhenius_dataframe.columns))[None, :],
                                          columns=arrhenius_dataframe.columns, index=[len(arrhenius_dataframe)])
        df1_finaldegassing["J"] = Final_Degassing
        df1_finaldegassing["K"] = Final_Degassing_uncert
        df1_finaldegassing["L"] = (Final_Degassing_uncert / Final_Degassing) * 100
        df1_finaldegassing["M"] = Final_Degassing / TotalXT
        df1_finaldegassing["N"] = (1 / TotalXT) * (np.sqrt(
            (((1 - (df1_finaldegassing["J"]) / TotalXT) * df1_finaldegassing["K"]) ** 2) + (
                        (df1_finaldegassing["J"] / TotalXT) ** 2) * (SUMSQ(arrhenius_dataframe['K']))))
        df1_finaldegassing["O"] = 100 * (df1_finaldegassing["N"] / df1_finaldegassing["M"])
        df1_finaldegassing["P"] = df1_finaldegassing["M"] + arrhenius_dataframe["P"].iloc[-1]
        df1_finaldegassing.round(2)
        new_df= pd.concat([arrhenius_dataframe, df1_finaldegassing])
    else:
        new_df = arrhenius_dataframe

    # This rounds a copy of the dataframe for creating figures
    # It does not round the dataframe that is used in subsequent calcaulation
    new_df = new_df.applymap(lambda x: round(x, 1) if isinstance(x, (int, float)) else x)

    return new_df

# Makes matplotlib figures to be saved from dataframe
def arrhenius_figs(arrhenius_dataframe, arrhenius_legend, geometry, input_filepath):
    # --------------------------- He Diffusivity Plots ---------------------------
    # Uses unrounded dataframe to plot
    # Initializes 3 rows x 2 columns of graphs as subplots object
    fig0, axs = plt.subplots(2, 3, figsize=(12, 6))
    fig0.subplots_adjust(top=0.88, hspace=.45, wspace=.25, bottom= .11)
    fig0.suptitle('Arrhenius Sheet 4')
    plt.figtext(0.99, 0.01, 'Geometry: ' +str(geometry)+ '\nFilepath: ' + str(input_filepath), horizontalalignment='right')

    # Graph 1 (Top Left): 3He Diffusivity in Infinite Sheet
    axs[0, 0].plot(arrhenius_dataframe['E'], arrhenius_dataframe["S"], c='lightblue', zorder=0)
    axs[0, 0].errorbar(arrhenius_dataframe['E'], arrhenius_dataframe["S"], yerr=arrhenius_dataframe['X'], fmt='none',
                       c='darkblue', capsize=3,
                       zorder=10)
    axs[0, 0].set_title('3He Diffusivity- ' + str(geometry), pad=18)
    axs[0, 0].set_xlabel('10^4/T, K¯¹')
    axs[0, 0].set_ylabel('D/a², s¯¹')

    # Graph 2 (Bottom Left): 3He Diffusivity in Infinite Sheet
    axs[1, 0].errorbar(arrhenius_dataframe['E'], arrhenius_dataframe["S"], yerr=arrhenius_dataframe['X'], fmt='none',
                       c='royalblue', capsize=3,
                       zorder=10)
    axs[1, 0].set_title('3He Diffusivity- ' + str(geometry), pad=18)
    axs[1, 0].set_xlabel('10^4/T, K¯¹')
    axs[1, 0].set_ylabel('D/a², s¯¹')

    # Graph 3 (Top Center): 3He Diffusivity in Infinite Sheet
    axs[0, 1].plot(arrhenius_dataframe['E'], arrhenius_dataframe["Y"], c='lightblue', zorder=0)
    axs[0, 1].errorbar(arrhenius_dataframe['E'], arrhenius_dataframe["Y"], yerr=arrhenius_dataframe['AB'], fmt='none',
                       c='darkblue', capsize=3,
                       zorder=10)
    axs[0, 1].set_title('3He Diffusivity- ' + str(geometry), pad=18)
    axs[0, 1].set_xlabel('10^4/T, K¯¹')
    axs[0, 1].set_ylabel('ln(D/a²), ln(s¯¹)')

    # Graph 4 (Top Right): 3He Diffusivity in Infinite Sheet
    axs[0, 2].plot(arrhenius_dataframe['E'], arrhenius_dataframe["Y"], c='lightblue', zorder=0)
    axs[0, 2].errorbar(arrhenius_dataframe['E'], arrhenius_dataframe["Y"], yerr=arrhenius_dataframe['AB'], fmt='none',
                       c='darkblue', capsize=3,
                       zorder=10)
    axs[0, 2].set_xlim([arrhenius_dataframe['E'].min() - 1, (arrhenius_dataframe['E'].mean() + .5)])
    axs[0, 2].set_title('3He Diffusivity- ' + str(geometry), pad=18)
    axs[0, 2].set_xlabel('10^4/T, K¯¹')
    axs[0, 2].set_ylabel('ln(D/a²), ln(s¯¹)')

    # Graph 5 (Bottom Center): 3He Diffusivity in Infinite Sheet
    axs[1, 1].errorbar(arrhenius_dataframe['E'], arrhenius_dataframe["Y"], yerr=arrhenius_dataframe['AB'], fmt='none',
                       c='royalblue', capsize=3,
                       zorder=10)
    axs[1, 1].set_title('3He Diffusivity- ' + str(geometry), pad=18)
    axs[1, 1].set_xlabel('10^4/T, K¯¹')
    axs[1, 1].set_ylabel('ln(D/a²), ln(s¯¹)')

    # Graph 6 (Bottom Right): Relative SD (1σ) in D/a² vs. T
    axs[1, 2].scatter(arrhenius_dataframe['C'], arrhenius_dataframe["W"], edgecolors='g', marker='d', facecolors='none')
    axs[1, 2].set_title('Relative SD (1σ) in D/a² vs. T', pad=18)
    axs[1, 2].set_xlabel('T, °C')
    axs[1, 2].set_ylabel('Relative 1σ error of (D/a²), %')

    arrhenius1 = df_round(arrhenius_dataframe)

    # --------------------------- Arrhenius: Table 1 ---------------------------
    fig1 = plt.figure(figsize=(8, len(arrhenius1['A']) * .22))
    ax1 = plt.subplot(111)
    ax1.autoscale(enable=True)
    ax1.axis('off')
    c = arrhenius1.loc[:, "A":'I'].shape[1]
    headers = np.array([""] * (c))
    ax1.text(-.0485, .052, 'Step Temperature', fontsize=10)
    ax1.text(.0215, .052, 'Step Duration', fontsize=10)
    table1 = ax1.table(
        cellText=np.vstack([headers, arrhenius_legend.loc[:, "A":'I'], arrhenius1.loc[:, "A":'I'].values]),
        cellColours=[['lightgray'] * c] + [['none'] * c] * (len(arrhenius1["A"]) + 1), bbox=[0, 0, 1, 1])
    table1.auto_set_column_width(col=list(range(c)))
    # Font adjustment
    table1.auto_set_font_size(False)
    table1.set_fontsize(9)
    cells = table1._cells
    for cell in cells:
        if cell[0] == 1:
            table1._cells[cell].set_fontsize(8)
            table1._cells[cell]._loc = "center"
        if cell[0] == 0:
            table1._cells[cell]._loc = "left"
            table1._cells[cell].set_linewidth(0)
    fig1.subplots_adjust(top=0.9, bottom=.08)
    fig1.suptitle('Arrhenius Sheet 1')
    plt.figtext(0.99, 0.01, 'Geometry: ' + str(geometry) + '\nFilepath: ' + str(input_filepath), horizontalalignment='right')

    # --------------------------- Arrhenius: Table 2 ---------------------------
    # Table formatting
    fig2 = plt.figure(figsize=(20, len(arrhenius1['A']) * .22))
    ax2 = plt.subplot(111)
    ax2.autoscale(enable=True)
    ax2.axis('off')
    c = arrhenius1.loc[:, "J":'X'].shape[1] + 1
    col_list = list("A")
    col_list.extend(list(arrhenius1.loc[:, "J":'X'].columns))
    headers = np.array([""] * (c))
    ax2.text(-.052, .049, 'He, pcc \n(cold blk corr.)', fontsize=9.5)
    ax2.text(-.0305, .049, 'Gas fraction', fontsize=9.5)
    ax2.text(-.009, .049, 'Cumultive gas\nfraction', fontsize=9.5)
    ax2.text(.0121, .049, 'Diffusivity divided by \ndiffusion domain size square', fontsize=9.5)
    table2 = ax2.table(cellText=np.vstack([headers, headers, arrhenius_legend[col_list], arrhenius1[col_list].values]),cellColours=[['lightgray'] * c] * 2 + [['none'] * c] * (len(arrhenius1["J"]) + 1),bbox=[0, 0, 1, 1])
    table2.auto_set_column_width(col=list(range(len(arrhenius_legend[col_list]))))
    # Font adjustment
    table2.auto_set_font_size(False)
    table2.set_fontsize(9)
    cells = table2._cells
    for cell in cells:
        if cell[0] == 2:
            table2._cells[cell].set_fontsize(8)
            table2._cells[cell]._loc = "center"
        if cell[0] == 1 or cell[0] == 0:
            table2._cells[cell].set_linewidth(0)
            # If there is a final degassing step, the last row will appear red in the table
        if Final_Degassing_step == True:
            if cell[0] == len(arrhenius1["J"]) + 2:
                table2._cells[cell].get_text().set_color('red')
    fig2.subplots_adjust(top=0.9, bottom=.08)
    fig2.suptitle('Arrhenius Sheet 2')
    plt.figtext(0.99, 0.01, 'Geometry: ' + str(geometry) + '\nFilepath: ' + str(input_filepath),horizontalalignment='right')
    # --------------------------- Arrhenius: Table 3 ---------------------------
    # Table formatting
    fig3 = plt.figure(figsize=(20, len(arrhenius1['A']) * .22))
    ax3 = plt.subplot(111)
    ax3.autoscale(enable=True)
    ax3.axis('off')
    col_list = list("A")
    col_list.extend(list(arrhenius1.loc[:, "Y":'AJ'].columns))
    c = arrhenius1.loc[:, "Y":'AJ'].shape[1] + 1
    ax3.text(-.0525, .052, 'Natural Log of (D/a2)', fontsize=11)
    ax3.text(-.0167, .052, 'Diffusivity', fontsize=11)
    ax3.text(.02, .052, 'Natural Log of (D/a2)', fontsize=11)
    headers = np.array([''] * (c))
    table3 = ax3.table(cellText=np.vstack([headers, arrhenius_legend[col_list], arrhenius1[col_list].values]),
                       cellColours=[['lightgray'] * c] + [['none'] * c] * (len(arrhenius1["Y"]) + 1), bbox=[0, 0, 1, 1])
    table3.auto_set_column_width(col=list(range(len(arrhenius_legend[col_list]))))
    # Font adjustment
    table3.auto_set_font_size(False)
    table3.set_fontsize(9)
    cells = table3._cells
    for cell in table3._cells:
        if cell[0] == 1:
            table3._cells[cell].set_fontsize(8)
            table3._cells[cell]._loc = "center"
        if cell[0] == 0:
            table3._cells[cell]._loc = "left"
            table3._cells[cell].set_linewidth(0)
    fig3.suptitle('Arrhenius Sheet 3')
    fig3.subplots_adjust(top=0.9, bottom=.08)
    plt.figtext(0.99, 0.01, 'Geometry: ' + str(geometry) + '\nFilepath: ' + str(input_filepath), horizontalalignment='right')

    return [fig1, fig2, fig3, fig0]

# Main fxn (calculations used in all geometries)
def general_arrhenius_calcs(geometry, input_filepath):
    # Read csv
    arrhenius1, arrhenius_legend = parse_data(input_filepath)

    # Total gas fraction calculation
    global TotalXT
    TotalXT = arrhenius1['J'].sum()
    TotalXT_uncert = np.sqrt(SUMSQ(arrhenius1['K']))

    # --------------------------- Check if final degassing step ---------------------------
    # removed from dataframe and stored as variables if there is a final
    # degassing step for ease of indexing in following calculations
    global Final_Degassing_uncert
    global Final_Degassing_step
    global Final_Degassing
    if arrhenius1['H'].iloc[-1] < 0:
        Final_Degassing = arrhenius1.loc[arrhenius1.index[-1], 'J']
        Final_Degassing_uncert = arrhenius1.loc[arrhenius1.index[-1], 'K']
        arrhenius1.drop(arrhenius1.tail(1).index, inplace=True)
        Final_Degassing_step = True
    else:
        Final_Degassing_uncert = 0
        Final_Degassing = 0

    # ---------------------------  Step Temperature Calculations ---------------------------
    arrhenius_legend['E'] = '10^4/T (K¯¹)'
    arrhenius1['E'] = 10000 / (arrhenius1['C'] + 273.15)

    arrhenius_legend['F'] = '1σ•10^4/T (K¯¹)'
    arrhenius1['F'] = np.sqrt((((-1.0) * 10.0 ** 4 / (arrhenius1['C'] + 273.15) ** 2) * arrhenius1['D']) ** 2)

    # ---------------------------  Step Duration Calculations ---------------------------
    arrhenius_legend['G'] = 'Δt (h)'
    arrhenius1['G'] = arrhenius1['H'] / 3600

    arrhenius_legend['B'] = 't, cum (h)'
    arrhenius1['B'] = arrhenius1['G'].cumsum()

    arrhenius_legend['I'] = 'δΔtn•1σ (s)'
    arrhenius1['I'] = 30

    # ---------------------------  He, pcc (cold blk corr.) Calculations ---------------------------

    arrhenius_legend['K'] = 'δXn•1σ (pcc)'
    arrhenius1['K'] = np.where(arrhenius1['J'] < .06, arrhenius1['J'] * 0.1,
                               (np.where(arrhenius1['J'] < 0.2, arrhenius1['J'] * 0.06, arrhenius1['J'] * 0.03)))

    arrhenius_legend['L'] = 'Rel 1σ (%)'
    arrhenius1['L'] = (arrhenius1['K'] / arrhenius1['J']) * 100

    # ---------------------------  Gas Fraction Calculations ---------------------------

    arrhenius_legend['M'] = 'Xn/XT'
    arrhenius1['M'] = arrhenius1['J'] / TotalXT

    arrhenius_legend['N'] = 'Xn/XT•1σ'
    for idx in range(len(arrhenius1)):
        arrhenius1.loc[arrhenius1.index[idx], 'N'] = (1 / TotalXT) * (np.sqrt((((1 - (arrhenius1.loc[arrhenius1.index[idx], 'J']) / TotalXT) * arrhenius1.loc[arrhenius1.index[idx], 'K']) ** 2) + ((arrhenius1.loc[arrhenius1.index[idx], 'J'] / TotalXT) ** 2) * (SUMSQ(arrhenius1.loc[arrhenius1.index != idx, 'K']) + (Final_Degassing_uncert ** 2))))

    arrhenius_legend['O'] = 'Rel 1σ (%)'
    arrhenius1['O'] = arrhenius1['N'] / arrhenius1['M'] * 100

    # --------------------------- Cumulative gas fraction Calculations ---------------------------

    arrhenius_legend['P'] = 'fn'
    arrhenius1['P'] = arrhenius1['M'].cumsum()

    arrhenius_legend['Q'] = 'δfn•1σ'
    for idx in range(len(arrhenius1)):
        arrhenius1.loc[arrhenius1.index[idx], 'Q'] = np.sqrt((((1 - arrhenius1.loc[arrhenius1.index[idx], 'P']) / TotalXT) ** 2) * (SUMSQ(arrhenius1.loc[arrhenius1.index <= idx, 'K'])) + (((arrhenius1.loc[arrhenius1.index[idx], 'P'] / TotalXT) ** 2) * (SUMSQ(arrhenius1.loc[arrhenius1.index > idx, 'K']) + Final_Degassing_uncert ** 2)))

    arrhenius_legend['R'] = 'Rel 1σ (%)'
    arrhenius1['R'] = arrhenius1['Q'] / arrhenius1['P'] * 100

    # ----------------- Diffusivity divided by diffusion domain size square Calculations -------------------

    # Call functions for geometry specific calculations
    if geometry == 'infinite sheet':
        sheet_arrhenius(arrhenius_legend,arrhenius1)
    elif geometry == 'sphere':
        sphere_arrhenius(arrhenius_legend,arrhenius1)
    elif geometry == 'infinite cylinder':
        cylinder_arrhenius(arrhenius_legend,arrhenius1)
    else:
        warnings.warn('select geometry',FutureWarning)

    arrhenius_legend['V'] = ' δ(D/a2)(1σ)  (s¯¹)'
    arrhenius1["V"] = np.sqrt(arrhenius1['T'] ** 2 + arrhenius1['U'] ** 2)

    arrhenius_legend['W'] = 'Rel. 1σ (%)'
    arrhenius1["W"] = 100 * (arrhenius1["V"] / arrhenius1["S"])

    arrhenius_legend['X'] = '2σ (s¯¹)'
    arrhenius1["X"] = arrhenius1["V"] * 2

    # --------------------------- Natural Log of (D/a2) Calculations ---------------------------

    arrhenius_legend['Y'] = 'ln(D/a²)n (ln(s¯¹))'
    arrhenius1["Y"] = np.log(arrhenius1["S"])

    arrhenius_legend['Z'] = 'δ(ln(D/a²)•1σ (ln(s¯¹))'
    arrhenius1["Z"] = arrhenius1["V"] / arrhenius1["S"]

    arrhenius_legend['AA'] = 'Rel. 1σ (%)'
    arrhenius1["AA"] = np.absolute(100 * (arrhenius1["Z"] / arrhenius1["Y"]))

    arrhenius_legend['AB'] = '2σ (ln(s¯¹))]'
    arrhenius1["AB"] = 2 * arrhenius1["Z"]

    # --------------------------- Diffusivity Calculations ---------------------------

    arrhenius_legend['AC'] = 'Dn (cm^2/s¯¹)'
    arrhenius1["AC"] = arrhenius1["S"] * (domain_size ** 2) * 10 ** (-8)

    arrhenius_legend['AD'] = 'Dn•1σ (cm^2/s¯¹)'
    arrhenius1["AD"] = arrhenius1["AC"] * np.sqrt(
        (arrhenius1["V"] / arrhenius1["S"]) ** 2 + (2 * (domain_size_uncert / domain_size)) ** 2)

    arrhenius_legend['AE'] = 'Rel 1σ (%)'
    arrhenius1["AE"] = 100 * arrhenius1['AD'] / arrhenius1['AC']

    arrhenius_legend['AF'] = '2σ (ln(s¯¹))'
    arrhenius1["AF"] = 2 * arrhenius1['AD']

    # --------------------------- Natural Log of (D) Calculations ---------------------------

    arrhenius_legend['AG'] = 'Ln(Dn) (Ln(cm2/s¯¹))'
    arrhenius1["AG"] = np.log(arrhenius1["AC"])

    arrhenius_legend['AH'] = '1σ (Ln(cm2/s¯¹))'
    arrhenius1["AH"] = arrhenius1["AD"] / arrhenius1["AC"]

    arrhenius_legend['AI'] = 'Rel 1σ (%)'
    arrhenius1["AI"] = np.absolute(100 * (arrhenius1["AH"] / arrhenius1["AG"]))

    arrhenius_legend['AJ'] = '2σ (Ln(cm2/s¯¹))'
    arrhenius1["AJ"] = 2 * arrhenius1["AH"]

    # Column sorting
        # Sorted by column name and length (mimics A->Z & AA->AZ in excel)
    arrhenius_legend = arrhenius_legend.reindex(sorted(arrhenius_legend.columns), axis=1)
    arrhenius_legend = arrhenius_legend.reindex(sorted(arrhenius_legend.columns, key=len), axis=1)
    arrhenius1 = arrhenius1.reindex(sorted(arrhenius1.columns), axis=1)
    arrhenius1 = arrhenius1.reindex(sorted(arrhenius1.columns, key=len), axis=1)

    figs = arrhenius_figs(arrhenius1, arrhenius_legend, geometry, input_filepath)


    return arrhenius1, figs