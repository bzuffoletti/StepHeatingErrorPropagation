import pandas as pd
import numpy as np
global R_kcal_molK
import matplotlib.pyplot as plt
R_kcal_molK = 0.0019872041
global R_kJ_molK
R_kJ_molK = 0.0083144621
global domain_size
domain_size = 27
global domain_size_uncert
domain_size_uncert = 2

def round_df_tup(df_tup):
    new_df_tup = []
    for df in df_tup:
        df = df.applymap(lambda x: round(x, 1) if isinstance(x, (int, float)) else x)
        new_df_tup.append(df)
    return new_df_tup

def EVA_figs(EVA_init_range1, first_EVA_iter_range1, second_EVA_iter_range1, EVA_init_range2, first_EVA_iter_range2, second_EVA_iter_range2, EVA_init_range3, first_EVA_iter_range3, second_EVA_iter_range3, ranges, input_filepath, geometry):
    fig6a, axs6a = plt.subplots(3, 2,  figsize=(8, 9.5), gridspec_kw={'width_ratios': [1, 1.75]})
    fig6a.suptitle('Effective Variance Approximation (EVA)\nRetrograde steps ' + str(int(ranges[0][0]))+ ':' + str(int(ranges[0][1])))
    axs6a[0, 0].set_title('Initial')
    axs6a[0, 0].table(cellText=EVA_init_range1[1].values, colLabels=EVA_init_range1[1].columns, loc='center')
    tableaa = axs6a[0, 1].table(cellText=EVA_init_range1[0].values, colLabels=EVA_init_range1[0].columns, loc='center')
    axs6a[1, 0].set_title('First Iteration')
    axs6a[1, 0].table(cellText=first_EVA_iter_range1[1].values, colLabels=first_EVA_iter_range1[1].columns,
                      loc='center')
    tableab = axs6a[1, 1].table(cellText=first_EVA_iter_range1[0].values, colLabels=first_EVA_iter_range1[0].columns,
                                loc='center')
    axs6a[2, 0].set_title('Second Iteration')
    axs6a[2, 0].table(cellText=second_EVA_iter_range1[1].values, colLabels=second_EVA_iter_range1[1].columns,
                      loc='center')
    tableac = axs6a[2, 1].table(cellText=second_EVA_iter_range1[0].values, colLabels=second_EVA_iter_range1[0].columns,
                                loc='center')
    tableaa.auto_set_font_size(False)
    tableaa.set_fontsize(6.5)
    tableab.auto_set_font_size(False)
    tableab.set_fontsize(6.5)
    tableac.auto_set_font_size(False)
    tableac.set_fontsize(6.5)
    for array in axs6a:
        for ax in array:
            ax.axis('off')
            ax.autoscale(enable=True)
    fig6a.subplots_adjust(top=0.9, bottom=.08)
    plt.figtext(0.99, 0.01, 'Geometry: ' + str(geometry) + '\nFilepath: ' + str(input_filepath), horizontalalignment='right')


    fig6b, axs6b = plt.subplots(3, 2, figsize=(8, 9.5),
                                gridspec_kw={'width_ratios': [1, 1.75]})
    fig6b.suptitle('Effective Variance Approximation (EVA)\nRetrograde steps ' + str(int(ranges[1][0]))+ ':' + str(int(ranges[1][1])))
    axs6b[0, 0].set_title('Initial')
    axs6b[0, 0].table(cellText=EVA_init_range2[1].values, colLabels=EVA_init_range2[1].columns, loc='center')
    tableba = axs6b[0, 1].table(cellText=EVA_init_range2[0].values, colLabels=EVA_init_range2[0].columns, loc='center')
    axs6b[1, 0].set_title('First Iteration')
    axs6b[1, 0].table(cellText=first_EVA_iter_range2[1].values, colLabels=first_EVA_iter_range2[1].columns,
                      loc='center')
    tablebb = axs6b[1, 1].table(cellText=first_EVA_iter_range2[0].values, colLabels=first_EVA_iter_range2[0].columns,
                                loc='center')
    axs6b[2, 0].set_title('Second Iteration')
    axs6b[2, 0].table(cellText=second_EVA_iter_range2[1].values, colLabels=second_EVA_iter_range2[1].columns,
                      loc='center')
    tablebc = axs6b[2, 1].table(cellText=second_EVA_iter_range2[0].values, colLabels=second_EVA_iter_range2[0].columns,
                                loc='center')
    tableba.auto_set_font_size(False)
    tableba.set_fontsize(6.5)
    tablebb.auto_set_font_size(False)
    tablebb.set_fontsize(6.5)
    tablebc.auto_set_font_size(False)
    tablebc.set_fontsize(6.5)
    for array in axs6b:
        for ax in array:
            ax.axis('off')
            ax.autoscale(enable=True)
    fig6b.subplots_adjust(top=0.9, bottom=.08)
    plt.figtext(0.99, 0.01, 'Geometry: ' + str(geometry) + '\nFilepath: ' + str(input_filepath), horizontalalignment='right')


    fig6c, axs6c = plt.subplots(3, 2, figsize=(8, 9.5),
                                gridspec_kw={'width_ratios': [1, 1.75]})
    fig6c.suptitle('Effective Variance Approximation (EVA)\nRetrograde steps ' + str(int(ranges[0][0])) + ":" + str(int(ranges[0][1])) + " & " + str(int(ranges[1][0])) + ':' + str(int(ranges[1][1])))
    axs6c[0, 0].set_title('Initial')
    axs6c[0, 0].table(cellText=EVA_init_range3[1].values, colLabels=EVA_init_range3[1].columns,
                      loc='center')
    tableca = axs6c[0, 1].table(cellText=EVA_init_range3[0].values, colLabels=EVA_init_range3[0].columns,
                                loc='center')
    axs6c[1, 0].set_title('First Iteration')
    axs6c[1, 0].table(cellText=first_EVA_iter_range3[1].values,
                      colLabels=first_EVA_iter_range3[1].columns, loc='center')
    tablecb = axs6c[1, 1].table(cellText=first_EVA_iter_range3[0].values,
                                colLabels=first_EVA_iter_range3[0].columns, loc='center')
    axs6c[2, 0].set_title('Second Iteration')
    axs6c[2, 0].table(cellText=second_EVA_iter_range3[1].values,
                      colLabels=second_EVA_iter_range3[1].columns, loc='center')
    tablecc = axs6c[2, 1].table(cellText=second_EVA_iter_range3[0].values,
                                colLabels=second_EVA_iter_range3[0].columns, loc='center')
    tableca.auto_set_font_size(False)
    tableca.set_fontsize(6.5)
    tablecb.auto_set_font_size(False)
    tablecb.set_fontsize(6.5)
    tablecc.auto_set_font_size(False)
    tablecc.set_fontsize(6.5)
    for array in axs6c:
        for ax in array:
            ax.axis('off')
            ax.autoscale(enable=True)
    fig6c.subplots_adjust(top=0.9, bottom=.08)
    plt.figtext(0.99, 0.01, 'Geometry: ' + str(geometry) + '\nFilepath: ' + str(input_filepath), horizontalalignment='right')

    return fig6a, fig6b, fig6c

def init_EVA(x,y, step_index, arrhenius1):

    model, cov = np.polyfit(x, y, 1, cov=True)
    p = np.poly1d(model)
    yp = p(x)  # predicted y values based on fit
    m = model[0]

    EVA_iter_df0 = pd.DataFrame()
    EVA_iter_df0[0] = (arrhenius1.loc[step_index, "Z"] ** 2) + (m ** 2) * (arrhenius1.loc[step_index, "F"] ** 2)
    EVA_iter_df0[1] = 1 / EVA_iter_df0[0]
    Xwm = np.dot(arrhenius1.loc[step_index, "E"], EVA_iter_df0[1]) / np.sum(EVA_iter_df0[1])
    Ywm = np.dot(arrhenius1.loc[step_index, "Y"], EVA_iter_df0[1]) / np.sum(EVA_iter_df0[1])
    b = Ywm - m * Xwm
    EVA_iter_df0[2] = arrhenius1.loc[step_index, "E"] - Xwm
    EVA_iter_df0[3] = arrhenius1.loc[step_index, "Y"] - Ywm
    EVA_iter_df0[4] = arrhenius1.loc[step_index, "Y"] - b - m * arrhenius1.loc[step_index, "E"]
    EVA_iter_df0[5] = EVA_iter_df0[1] * EVA_iter_df0[2] * EVA_iter_df0[3]
    EVA_iter_df0[6] = EVA_iter_df0[1] * (EVA_iter_df0[2] ** 2)
    SE2_εw = sum([x * y * z for x, y, z in zip(EVA_iter_df0[1], EVA_iter_df0[4], EVA_iter_df0[4])]) / (
            len(arrhenius1.loc[step_index, "Y"]) - 2)
    mσ = np.sqrt(SE2_εw / sum([x * y * z for x, y, z in zip(EVA_iter_df0[1], EVA_iter_df0[2], EVA_iter_df0[2])]))
    bσ = np.sqrt(SE2_εw / sum(EVA_iter_df0[1]) + (mσ ** 2) * (Xwm ** 2))
    D0 = np.exp(b) * np.sqrt(bσ ** 2 + (2 * domain_size_uncert / domain_size) ** 2)
    EVA_iter_df1 = pd.DataFrame()
    EVA_iter_df1[0] = ["Ea (kcal/mol) =", "Ea (kJ/mol) =", "Do (cm^2/s) =", "Do/a^2 (1/s) =", "", "m=", 'b=', 'Xwm',
                       'Ywm', 'SE^2(εw):']
    EVA_iter_df1[1] = [(-m * 10000 * R_kcal_molK), (-m * 10000 * R_kJ_molK),
                       (np.exp(b) * (.0001 * domain_size) ** 2), (np.exp(b)), '', m, b, Xwm, Ywm, SE2_εw]
    EVA_iter_df1[2] = [(mσ * R_kcal_molK * 10000), (mσ * R_kJ_molK * 10000), D0, (bσ * np.exp(b)), "", mσ, bσ,
                       np.NaN, np.NaN, np.NaN]

    EVA_iter_df1.columns = ['', 'w,uc', '1σ']
    EVA_iter_df0.columns = ['σ2 (yi)-effective', 'wi', 'xi-xwm', 'yi-ywm', 'yi-yi-model', 'Num', 'Denom']

    return EVA_iter_df0, EVA_iter_df1

def iter_EVA(x,y,previous_df, step_index, arrhenius1):
    model, cov = np.polyfit(x, y, 1, cov=True)
    p = np.poly1d(model)
    yp = p(x)  # predicted y values based on fit
    model_m = model[0]
    m = sum(previous_df[0]["Num"]) / sum(previous_df[0]["Denom"])

    EVA_iter_df0 = pd.DataFrame()
    EVA_iter_df0[0] = (arrhenius1.loc[step_index, "Z"] ** 2) + (model_m** 2) * (arrhenius1.loc[step_index, "F"] ** 2)
    EVA_iter_df0[1] = 1 / EVA_iter_df0[0]
    Xwm = np.dot(arrhenius1.loc[step_index, "E"], EVA_iter_df0[1]) / np.sum(EVA_iter_df0[1])
    Ywm = np.dot(arrhenius1.loc[step_index, "Y"], EVA_iter_df0[1]) / np.sum(EVA_iter_df0[1])
    b = Ywm - m * Xwm
    EVA_iter_df0[2] = arrhenius1.loc[step_index, "E"] - Xwm
    EVA_iter_df0[3] = arrhenius1.loc[step_index, "Y"] - Ywm
    EVA_iter_df0[4] = arrhenius1.loc[step_index, "Y"] - b - m * arrhenius1.loc[step_index, "E"]
    EVA_iter_df0[5] = EVA_iter_df0[1] * EVA_iter_df0[2] * EVA_iter_df0[3]
    EVA_iter_df0[6] = EVA_iter_df0[1] * (EVA_iter_df0[2] ** 2)
    SE2_εw = sum([x * y * z for x, y, z in zip(EVA_iter_df0[1], EVA_iter_df0[4], EVA_iter_df0[4])]) / (
                len(arrhenius1.loc[step_index, "Y"]) - 2)
    mσ = np.sqrt(SE2_εw / sum([x * y * z for x, y, z in zip(EVA_iter_df0[1], EVA_iter_df0[2], EVA_iter_df0[2])]))
    bσ = np.sqrt(SE2_εw / sum(EVA_iter_df0[1]) + (mσ ** 2) * (Xwm ** 2))
    D0 = np.exp(b) * ((domain_size * .0001) ** 2)
    D0σ = D0 * np.sqrt((bσ) ** 2 + ((2 * domain_size_uncert / domain_size) ** 2))
    EVA_iter_df1 = pd.DataFrame()
    EVA_iter_df1[0] = ["Ea (kcal/mol) =", "Ea (kJ/mol) =", "Do (cm^2/s) =", "Do/a^2 (1/s) =", "", "m=", 'b=', 'Xwm',
                       'Ywm', 'SE^2(εw):']
    EVA_iter_df1[1] = [(-m * 10000 * R_kcal_molK), (-m * 10000 * R_kJ_molK),
                       (np.exp(b) * (.0001 * domain_size) ** 2), (np.exp(b)), '', m, b, Xwm, Ywm, SE2_εw]
    EVA_iter_df1[2] = [(mσ * R_kcal_molK * 10000), (mσ * R_kJ_molK * 10000), D0σ, (bσ * np.exp(b)), "", mσ, bσ, np.NaN,
                       np.NaN, np.NaN]

    EVA_iter_df1.columns = ['', 'w,uc', '1σ']
    EVA_iter_df0.columns = ['σ2 (yi)-effective', 'wi', 'xi-xwm', 'yi-ywm', 'yi-yi-model', 'Num', 'Denom']

    return EVA_iter_df0, EVA_iter_df1

def make_step_ranges(arrhenius1):
    range1 = [np.floor(len(arrhenius1['A']) / 3), (np.floor(len(arrhenius1['A']) / 3) + 5)]
    range2 = [np.floor(len(arrhenius1['A']) / 3 * 2), (np.floor(len(arrhenius1['A']) / 3 * 2) + 5)]
    range3 = range1 + range2
    return [range1, range2, range3]

def EVA(arrhenius1, input_filepath, geometry):
    ranges = make_step_ranges(arrhenius1)

    x = arrhenius1.loc[ranges[0][0]:ranges[0][1], "Y"]
    y = arrhenius1.loc[ranges[0][0]:ranges[0][1], "E"]
    EVA_init_range1 = init_EVA(x, y, np.r_[ranges[0][0]:ranges[0][1]], arrhenius1)
    first_EVA_iter_range1 = iter_EVA(x,y, EVA_init_range1, np.r_[ranges[0][0]:ranges[0][1]], arrhenius1)
    second_EVA_iter_range1 = iter_EVA(x,y,first_EVA_iter_range1, np.r_[ranges[0][0]:ranges[0][1]], arrhenius1)

    x = arrhenius1.loc[ranges[1][0]:ranges[1][1], "Y"]
    y = arrhenius1.loc[ranges[1][0]:ranges[1][1], "E"]
    EVA_init_range2 = init_EVA(x, y, np.r_[ranges[1][0]:ranges[1][1]], arrhenius1)
    first_EVA_iter_range2 = iter_EVA(x, y, EVA_init_range2,np.r_[ranges[1][0]:ranges[1][1]], arrhenius1)
    second_EVA_iter_range2 = iter_EVA(x, y,first_EVA_iter_range2, np.r_[ranges[1][0]:ranges[1][1]], arrhenius1)

    x = arrhenius1.loc[np.r_[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1]], "Y"]
    y = arrhenius1.loc[np.r_[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1]], "E"]
    EVA_init_range3 = init_EVA(x, y, np.r_[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1]], arrhenius1)
    first_EVA_iter_range3 = iter_EVA(x,y, EVA_init_range3,np.r_[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1]], arrhenius1)
    second_EVA_iter_range3 = iter_EVA(x,y,first_EVA_iter_range3, np.r_[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1]], arrhenius1)

    EVA_init_range1 = round_df_tup(EVA_init_range1)
    first_EVA_iter_range1 = round_df_tup(first_EVA_iter_range1)
    second_EVA_iter_range1 = round_df_tup(second_EVA_iter_range1)
    EVA_init_range2 = round_df_tup(EVA_init_range2)
    first_EVA_iter_range2 = round_df_tup(first_EVA_iter_range2)
    second_EVA_iter_range2 = round_df_tup(second_EVA_iter_range2)
    EVA_init_range3 = round_df_tup(EVA_init_range3)
    first_EVA_iter_range3 = round_df_tup(first_EVA_iter_range3)
    second_EVA_iter_range3 = round_df_tup(second_EVA_iter_range3)

    figs = EVA_figs(EVA_init_range1, first_EVA_iter_range1, second_EVA_iter_range1, EVA_init_range2, first_EVA_iter_range2, second_EVA_iter_range2, EVA_init_range3, first_EVA_iter_range3, second_EVA_iter_range3, ranges, input_filepath, geometry)

    return figs


