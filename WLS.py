import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
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

def make_figures(init_range1, first_iter_range1, second_iter_range1, init_range2, first_iter_range2, second_iter_range2, init_range3, first_iter_range3, second_iter_range3, ranges, geometry, input_filepath):

    fig5a, axs5a = plt.subplots(3, 2, figsize=(8, 9.5), gridspec_kw={'width_ratios': [1, 1.75]}, )
    fig5a.suptitle('Weighted Least Squares (WLS) Regression\nRetrograde steps ' + str(int(ranges[0][0])) +':' + str(int(ranges[0][1])))
    axs5a[0, 0].set_title('Initial')
    axs5a[0, 0].table(cellText=init_range1[1].values, colLabels=init_range1[1].columns, loc='center')
    axs5a[0, 1].table(cellText=init_range1[0].values, colLabels=init_range1[0].columns, loc='center')
    axs5a[1, 0].set_title('First Iteration')
    axs5a[1, 0].table(cellText=first_iter_range1[1].values, colLabels=first_iter_range1[1].columns,loc='center')
    axs5a[1, 1].table(cellText=first_iter_range1[0].values, colLabels=first_iter_range1[0].columns,loc='center')
    axs5a[2, 0].set_title('Second Iteration')
    axs5a[2, 0].table(cellText=second_iter_range1[1].values, colLabels=second_iter_range1[1].columns,loc='center')
    axs5a[2, 1].table(cellText=second_iter_range1[0].values, colLabels=second_iter_range1[0].columns, loc='center')
    for array in axs5a:
        for ax in array:
            ax.axis('off')
            ax.autoscale(enable=True)
    fig5a.subplots_adjust(top=0.9, bottom=.08)
    plt.figtext(0.99, 0.01, 'Geometry: ' + str(geometry) + '\nFilepath: ' + str(input_filepath),horizontalalignment='right')

    fig5b, axs5b = plt.subplots(3, 2, gridspec_kw={'width_ratios': [1, 1.75]},figsize=(8, 9.5))
    fig5b.suptitle('Weighted Least Squares (WLS) Regression\nRetrograde steps '+ str(ranges[1][0]) +':' + str(ranges[1][1]))
    axs5b[0, 0].set_title('Initial')
    axs5b[0, 0].table(cellText=init_range2[1].values, colLabels=init_range2[1].columns, loc='center')
    axs5b[0, 1].table(cellText=init_range2[0].values, colLabels=init_range2[0].columns, loc='center')
    axs5b[1, 0].set_title('First Iteration')
    axs5b[1, 0].table(cellText=first_iter_range2[1].values, colLabels=first_iter_range2[1].columns,loc='center')
    axs5b[1, 1].table(cellText=first_iter_range2[0].values, colLabels=first_iter_range2[0].columns,loc='center')
    axs5b[2, 0].set_title('Second Iteration')
    axs5b[2, 0].table(cellText=second_iter_range2[1].values, colLabels=second_iter_range2[1].columns,loc='center')
    axs5b[2, 1].table(cellText=second_iter_range2[0].values, colLabels=second_iter_range2[0].columns,loc='center')
    for array in axs5b:
        for ax in array:
            ax.axis('off')
            ax.autoscale(enable=True)
    fig5b.subplots_adjust(top=0.9, bottom=.08)
    plt.figtext(0.99, 0.01, 'Geometry: ' + str(geometry) + '\nFilepath: ' + str(input_filepath),horizontalalignment='right')

    fig5c, axs5c = plt.subplots(3, 2, figsize=(8, 9.5), gridspec_kw={'width_ratios': [1, 1.75]})
    fig5c.suptitle('Weighted Least Squares (WLS) Regression\nRetrograde steps '+ str(int(ranges[2][0]))+':'+str(int(ranges[2][1])) +' & ' + str(int(ranges[2][2]))+':'+str(int(ranges[2][3])))
    axs5c[0, 0].set_title('Initial')
    axs5c[0, 0].table(cellText=init_range3[1].values, colLabels=init_range3[1].columns,loc='center')
    axs5c[0, 1].table(cellText=init_range3[0].values, colLabels=init_range3[0].columns,loc='center')
    axs5c[1, 0].set_title('First Iteration')
    axs5c[1, 0].table(cellText=first_iter_range3[1].values,colLabels=first_iter_range3[1].columns, loc='center')
    axs5c[1, 1].table(cellText=first_iter_range3[0].values,colLabels=first_iter_range3[0].columns, loc='center')
    axs5c[2, 0].set_title('Second Iteration')
    axs5c[2, 0].table(cellText=second_iter_range3[1].values,colLabels=second_iter_range3[1].columns, loc='center')
    axs5c[2, 1].table(cellText=second_iter_range3[0].values, colLabels=second_iter_range3[0].columns, loc='center')
    for array in axs5c:
        for ax in array:
            ax.axis('off')
            ax.autoscale(enable=True)
    fig5c.subplots_adjust(top=0.9, bottom=.08)
    plt.figtext(0.99, 0.01, 'Geometry: ' + str(geometry) + '\nFilepath: ' + str(input_filepath),horizontalalignment='right')

    return fig5a, fig5b, fig5c

def init_WLS(x,y, step_index, arrhenius1):
    model, cov = np.polyfit(x, y, 1, cov=True)
    p = np.poly1d(model)
    yp = p(x)  # predicted y values based on fit
    m = model[0]

    WLS_init_df0 = pd.DataFrame()
    WLS_init_df0[0] = 1/(arrhenius1.loc[step_index, "Z"]**2)
    initial_Xwm = np.dot(arrhenius1.loc[step_index, "E"], WLS_init_df0[0]) / np.sum(WLS_init_df0[0])
    initial_Ywm = np.dot(arrhenius1.loc[step_index, "Y"], WLS_init_df0[0]) / np.sum(WLS_init_df0[0])
    b = initial_Ywm - m * initial_Xwm


    WLS_init_df0[1] = arrhenius1.loc[step_index, "E"]- initial_Xwm
    WLS_init_df0[2] = arrhenius1.loc[step_index, "Y"]- initial_Ywm
    WLS_init_df0[3] = arrhenius1.loc[step_index, "Y"] - b - m * arrhenius1.loc[step_index, "E"]
    SE2_εw =  sum([x * y * z for x, y, z in zip(WLS_init_df0[3], WLS_init_df0[3], WLS_init_df0[0])])/ (len(arrhenius1.loc[step_index, "Y"])-2)
    mσ = np.sqrt(SE2_εw / sum([x * y * z for x, y, z in zip(WLS_init_df0[0], WLS_init_df0[3], WLS_init_df0[3])]))
    bσ = np.sqrt(SE2_εw / sum(WLS_init_df0[0]) + (mσ ** 2) * (initial_Xwm ** 2))
    WLS_init_df0[4] = WLS_init_df0[0]* WLS_init_df0[1]* WLS_init_df0[2]
    WLS_init_df0[5] = WLS_init_df0[0]* (WLS_init_df0[1]**2)
    WLS_init_df1 = pd.DataFrame()
    WLS_init_df1[0] = ["Ea (kcal/mol) =", "Ea (kJ/mol) =", "Do (cm^2/s) =", "Do/a^2 (1/s) =", "", "Init m", 'Init b', 'initial Xwm', 'initial Ywm', 'SE^2(εw):']
    WLS_init_df1[1] = [(-10000*R_kcal_molK*m), (-10000*R_kJ_molK*m), (np.exp(b)*(domain_size*.0001)**2), np.exp(b), '', m, b, initial_Xwm,initial_Ywm,SE2_εw]
    WLS_init_df1[2] = [(mσ*R_kcal_molK*10000), (mσ*R_kJ_molK*10000), WLS_init_df1.iloc[2,1]*np.sqrt((bσ)**2+(2*domain_size_uncert/domain_size)**2), bσ*WLS_init_df1.iloc[3,1],'', mσ, bσ, np.NaN, np.NaN, np.NaN]

    WLS_init_df1.columns = ['', 'w,uc', '1σ']
    WLS_init_df0.columns = ['wi', 'xi-xwm', 'yi-ywm', 'yi-yi-model', 'Num', 'Denom']

    return WLS_init_df0, WLS_init_df1

def iter_WLS(previous_df, step_index, arrhenius1):

    m = sum(previous_df[0]['Num'])/sum(previous_df[0]["Denom"])
    WLS_iter_df0 = pd.DataFrame()
    WLS_iter_df0[0] = 1/(arrhenius1.loc[step_index, "Z"]**2)
    Xm = np.dot(arrhenius1.loc[step_index, "E"], WLS_iter_df0[0]) / np.sum(WLS_iter_df0[0])
    Ym = np.dot(arrhenius1.loc[step_index, "Y"], WLS_iter_df0[0]) / np.sum(WLS_iter_df0[0])
    b = Ym - m * Xm

    WLS_iter_df0[1] = arrhenius1.loc[step_index, "E"] - Xm
    WLS_iter_df0[2] = arrhenius1.loc[step_index, "Y"] - Ym
    WLS_iter_df0[3] = arrhenius1.loc[step_index, "Y"] - b - m * arrhenius1.loc[step_index, "E"]
    SE2_εw =  sum([x * y * z for x, y, z in zip(WLS_iter_df0[0], WLS_iter_df0[3], WLS_iter_df0[3])])/ (len(arrhenius1.loc[step_index, "Y"])-2)
    WLS_iter_df0[4] = WLS_iter_df0[0] * WLS_iter_df0[1] * WLS_iter_df0[2]
    WLS_iter_df0[5] = WLS_iter_df0[0] * (WLS_iter_df0[1] ** 2)
    mσ = np.sqrt(SE2_εw/sum([x * y * z for x, y, z in zip(WLS_iter_df0[0], WLS_iter_df0[3], WLS_iter_df0[3])]))
    bσ = np.sqrt(SE2_εw/sum(WLS_iter_df0[0])+(mσ**2)*(Xm**2))

    WLS_iter_df1 = pd.DataFrame()
    WLS_iter_df1[0] = ["Ea (kcal/mol) =", "Ea (kJ/mol) =", "Do (cm^2/s) =", "Do/a^2 (1/s) =", "", 'm=', 'b=','Xwm', 'Ywm', 'SE^2(εw):' ]
    WLS_iter_df1[1] = [(-10000*R_kcal_molK*m), (-10000*R_kJ_molK*m), (np.exp(b)*(domain_size*.0001)**2), np.exp(b), '', m, b, Xm, Ym, SE2_εw,]
    WLS_iter_df1[2] = [(mσ*R_kcal_molK*10000), (mσ*R_kJ_molK*10000), WLS_iter_df1.iloc[2,1]*np.sqrt((bσ)**2+(2*domain_size_uncert/domain_size)**2), bσ*WLS_iter_df1.iloc[3,1],'', mσ, bσ, np.NaN, np.NaN, np.NaN]

    WLS_iter_df1.columns = ['', 'w,uc', '1σ']
    WLS_iter_df0.columns = ['wi', 'xi-xwm', 'yi-ywm', 'yi-yi-model', 'Num', 'Denom']

    return WLS_iter_df0, WLS_iter_df1

def make_step_ranges(arrhenius1):
    range1 = [np.floor(len(arrhenius1['A']) / 3), (np.floor(len(arrhenius1['A']) / 3) + 5)]
    range2 = [np.floor(len(arrhenius1['A']) / 3 * 2), (np.floor(len(arrhenius1['A']) / 3 * 2) + 5)]
    range3 = range1 +range2
    return [range1, range2, range3]

def WLS(arrhenius1, geometry, input_filepath):
    ranges = make_step_ranges(arrhenius1)

    x = arrhenius1.loc[ranges[0][0]:ranges[0][1], "Y"]
    y = arrhenius1.loc[ranges[0][0]:ranges[0][1], "E"]
    WLS_init_range1 = init_WLS(x,y, np.r_[ranges[0][0]:ranges[0][1]], arrhenius1)

    x = arrhenius1.loc[ranges[1][0]:ranges[1][1], "Y"]
    y = arrhenius1.loc[ranges[1][0]:ranges[1][1], "E"]
    WLS_init_range2 = init_WLS(x,y, np.r_[ranges[1][0]:ranges[1][1]], arrhenius1)

    x = arrhenius1.loc[np.r_[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1]], "Y"]
    y = arrhenius1.loc[np.r_[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1]], "E"]
    WLS_init_range3 = init_WLS(x,y, np.r_[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1]], arrhenius1)

    first_WLS_iter_range1 = iter_WLS(WLS_init_range1, np.r_[ranges[0][0]:ranges[0][1]], arrhenius1)
    first_WLS_iter_range2 = iter_WLS(WLS_init_range2, np.r_[ranges[1][0]:ranges[1][1]], arrhenius1)
    first_WLS_iter_range3 = iter_WLS(WLS_init_range3, np.r_[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1]], arrhenius1)

    second_WLS_iter_range1 = iter_WLS(first_WLS_iter_range1, np.r_[ranges[0][0]:ranges[0][1]], arrhenius1)
    second_WLS_iter_range2 = iter_WLS(first_WLS_iter_range2, np.r_[ranges[1][0]:ranges[1][1]], arrhenius1)
    second_WLS_iter_range3 = iter_WLS(first_WLS_iter_range3, np.r_[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1]], arrhenius1)

    init_range1 = round_df_tup(WLS_init_range1)
    first_iter_range1 = round_df_tup(first_WLS_iter_range1)
    second_iter_range1 = round_df_tup(second_WLS_iter_range1)
    init_range2 = round_df_tup(WLS_init_range2)
    first_iter_range2 = round_df_tup(first_WLS_iter_range2)
    second_iter_range2 = round_df_tup(second_WLS_iter_range2)
    init_range3 = round_df_tup(WLS_init_range3)
    first_iter_range3 = round_df_tup(first_WLS_iter_range3)
    second_iter_range3 = round_df_tup(second_WLS_iter_range3)

    figs = make_figures(init_range1, first_iter_range1, second_iter_range1, init_range2, first_iter_range2, second_iter_range2, init_range3, first_iter_range3, second_iter_range3, ranges, geometry, input_filepath)

    return figs