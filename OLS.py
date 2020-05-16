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

def OLS_calculate(x,y,z):
    n = len(y)
    dofreedom = n - 2

    model, cov = np.polyfit(x, y, 1, cov=True)
    p = np.poly1d(model)
    yp = p(x)  # predicted y values based on fit
    slope = model[0]
    intercept = model[1]
    r2 = np.corrcoef(x, y)[0][1] ** 2
    regression_ss = np.sum((yp - np.mean(y)) ** 2)
    residual_ss = np.sum((y - yp) ** 2)
    slope_pm = np.sqrt(residual_ss / (dofreedom * np.sum((x - np.mean(x)) ** 2)))
    intercept_pm = slope_pm * np.sqrt(np.sum(x ** 2) / n)
    s = np.sqrt(residual_ss / dofreedom)
    F = regression_ss / s ** 2

    Ea_kcal_mol = -slope * 10000 * R_kcal_molK
    Ea_kcal_mol_σ = slope_pm * 10000 * R_kcal_molK
    Ea_kcal_mol_rsd = 100 * Ea_kcal_mol_σ / Ea_kcal_mol
    Ea_kj_mol = -slope * 10000 * R_kJ_molK
    Ea_kj_mol_σ = slope_pm * 10000 * R_kJ_molK
    Ea_kj_mol_rsd = 100 * Ea_kj_mol_σ / Ea_kj_mol
    Do_a2 = np.exp(intercept)
    Do_a2_σ = Do_a2 * intercept_pm
    Do_a2_rsd = 100 * Do_a2_σ / Do_a2
    Do = Do_a2 * (domain_size / 10000) ** 2
    Do_σ = Do * np.sqrt((Do_a2_σ / Do_a2) ** 2 + (2 * domain_size_uncert / domain_size) ** 2)
    Do_rsd = 100 * Do_σ / Do

    yi_yi_model2 = (y - intercept - slope * x) ** 2
    σyi2 = 1 / (z ** 2)
    χ2 = round((np.dot(yi_yi_model2, σyi2) / dofreedom), 2)

    df_left = pd.DataFrame()
    df_left[0] = ["Ea (kcal/mol) =", "Ea (kJ/mol) =", "Do (cm^2/s) =", "Do/a^2 (1/s) ="]
    df_left[1] = [Ea_kcal_mol, Ea_kj_mol, Do, Do_a2]
    df_left[2] = [Ea_kcal_mol_σ, Ea_kj_mol_σ, Do_σ, Do_a2_σ]
    df_left[3] = [Ea_kcal_mol_rsd, Ea_kj_mol_rsd, Do_rsd, Do_a2_rsd]
    df_left.columns = ['χ2reduced =', χ2, '1σ', 'RDS %']

    df_right = pd.DataFrame()
    df_right[0] = ['m', 'm 1σ', 'R^2', "F", "regr ss"]
    df_right[1] = [slope, slope_pm, r2, F, regression_ss]
    df_right[2] = ["b", "b 1σ", 'σy', 'degr. freedom', 'residual ss']
    df_right[3] = [intercept, intercept_pm, s, dofreedom, residual_ss]

    return df_left, df_right

def make_figures(df1, df2, df3, df4, ranges, arrhenius1, geometry, input_filepath):
    fig4, axs4 = plt.subplots(4, 2)

    # Initialize of figure with 8 tables
    for array in axs4:
        for ax in array:
            ax.axis('off')
            ax.autoscale(enable=True)
    fig4.suptitle('Ordinary Least Squares (OLS) regression of ln(D/a^2) on (10^4/T)')
    # Create the 8 tables with OLS regression data
    axs4[0, 0].table(cellText=df1[0].values, colLabels=df1[0].columns, loc='center')
    axs4[0, 0].set_title('Retrograde Steps ' + str(int(ranges[0][0]))+':'+str(int(ranges[0][1])), fontsize=6)
    axs4[0, 1].table(cellText=df1[1].values, loc='center')
    axs4[1, 0].table(cellText=df2[0].values, colLabels=df2[0].columns, loc='center')
    axs4[1, 0].set_title('Retrograde Steps ' + str(int(ranges[1][0]))+':'+str(int(ranges[1][1])), fontsize=6)
    axs4[1, 1].table(cellText=df2[1].values, loc='center')
    axs4[2, 0].table(cellText=df3[0].values, colLabels=df3[0].columns, loc='center')
    axs4[2, 0].set_title('Retrograde Steps ' + str(int(ranges[2][0]))+':'+str(int(ranges[2][1])) +' & ' + str(int(ranges[2][2]))+':'+str(int(ranges[2][3])), fontsize=6)
    axs4[2, 1].table(cellText=df3[1].values, loc='center')
    axs4[3, 0].table(cellText=df4[0].values, colLabels=df4[0].columns, loc='center')
    axs4[3, 0].set_title('Retrograde & Prograde  Steps ' + str(int(ranges[3][0]))+':'+ str(int(ranges[3][1])), fontsize=6)
    axs4[3, 1].table(cellText=df4[1].values, loc='center')
    fig4.subplots_adjust(top=0.9, bottom=.1, hspace=.5)
    plt.figtext(0.99, 0.01, 'Geometry: ' + str(geometry) + '\nFilepath: ' + str(input_filepath),horizontalalignment='right')

    plt.style.use('ggplot')
    fig5 = plt.figure(figsize=(8, 8))
    ax5 = plt.subplot(111)
    plt.plot(arrhenius1.loc[ranges[1][0]:ranges[1][1], "E"], arrhenius1.loc[ranges[1][0]:ranges[1][1], "Y"], '-o', color='red')
    plt.plot(arrhenius1.loc[ranges[0][0]:ranges[0][1], "E"], arrhenius1.loc[ranges[0][0]:ranges[0][1], "Y"], '-o', color='green')
    ax5.set_ylim(bottom=ax5.get_ylim()[0] - 6)
    ax5.set_title(str(geometry) + " geometry:")
    ax5.set_xlabel('10^4/T (K^-1)')
    ax5.set_ylabel('lnD/a^2 (ln(1/s))')
    plt.text(13, -16.4,
             "All retrograde steps:\nEa =" + str(df3[0].iloc[1, 1].round(2)) + " kJ/mol \nD0 =" + str(
                 df3[0].iloc[2, 1].round(2)) + ' cm2/s \n')
    plt.text(13, -18.7,
             "Hi-T retrograde steps " +str(int(ranges[1][0]))+':'+str(int(ranges[1][1])) +"\nEa =" + str(df2[0].iloc[1, 1].round(2)) + " kJ/mol \nD0 =" + str(df2[0].iloc[2, 1].round(2)) + ' cm2/s \n', color="red")
    plt.text(13, -21, "Lo-T retrograde steps " + str(int(ranges[1][0]))+':'+str(int(ranges[1][1])) + ":\nEa =" + str(df1[0].iloc[1, 1].round(2)) + " kJ/mol \nD0 =" + str(
                 df1[0].iloc[2, 1].round(2)) + ' cm2/s \n', color='green')
    plt.text(13, -23.3, "Prograde & retrograde steps (" +  str(int(ranges[3][0]))+':'+ str(int(ranges[3][1]))+ '):\nEa =' + str(df4[0].iloc[1, 1].round(2)) + " kJ/mol \nD0 =" + str(df4[0].iloc[2, 1].round(2)) + ' cm2/s \n')
    for x, y, i in zip(arrhenius1.loc[ranges[1][0]:ranges[1][1], "E"], arrhenius1.loc[ranges[1][0]:ranges[1][1], "Y"], range(int(ranges[1][0]+1),int(ranges[1][1]+1))):
        ax5.annotate(i, xy=(x, y), xytext=(x, y + .15))
    for x, y, i in zip(arrhenius1.loc[ranges[0][0]:ranges[0][1], "E"], arrhenius1.loc[ranges[0][0]:ranges[0][1], "Y"], range(int(ranges[1][0]+1),int(ranges[1][1]+1))):
        ax5.annotate(i, xy=(x, y), xytext=(x, y + .15))
    fig5.subplots_adjust(top=0.9, bottom=.08)
    plt.figtext(0.99, 0.01, 'Geometry: ' + str(geometry) + '\nFilepath: ' + str(input_filepath),horizontalalignment='right')


    return [fig4, fig5]

def make_step_ranges(arrhenius1):
    range1 = [np.floor(len(arrhenius1['A']) / 3), (np.floor(len(arrhenius1['A']) / 3) + 5)]
    range2 = [np.floor(len(arrhenius1['A']) / 3 * 2), (np.floor(len(arrhenius1['A']) / 3 * 2) + 5)]
    range3 = range1 +range2
    range4 = [range1[0], range2[1]]
    return [range1, range2, range3, range4]

def OLS(arrhenius1, geometry, input_filepath):
    ranges = make_step_ranges(arrhenius1)

    x = arrhenius1.loc[ranges[0][0]:ranges[0][1], "Y"]
    y = arrhenius1.loc[ranges[0][0]:ranges[0][1], "E"]
    z = arrhenius1.loc[ranges[0][0]:ranges[0][1], "Z"]
    OLS_range1 = OLS_calculate(y, x, z)

    x = arrhenius1.loc[ranges[1][0]:ranges[1][1], "Y"]
    y = arrhenius1.loc[ranges[1][0]:ranges[1][1], "E"]
    z = arrhenius1.loc[ranges[1][0]:ranges[1][1], "Z"]
    OLS_range2 = OLS_calculate(y, x, z)

    x = pd.concat([arrhenius1.loc[ranges[2][0]:ranges[2][1], "Y"], arrhenius1.loc[ranges[2][2]:ranges[2][3], "Y"]])
    y = pd.concat([arrhenius1.loc[ranges[2][0]:ranges[2][1], "E"], arrhenius1.loc[ranges[2][2]:ranges[2][3], "E"]])
    z = pd.concat([arrhenius1.loc[ranges[2][0]:ranges[2][1], "Z"], arrhenius1.loc[ranges[2][2]:ranges[2][3], "Z"]])
    OLS_range3 = OLS_calculate(y, x, z)

    x = arrhenius1.loc[ranges[3][0]:ranges[3][1], "Y"]
    y = arrhenius1.loc[ranges[3][0]:ranges[3][1], "E"]
    z = arrhenius1.loc[ranges[3][0]:ranges[3][1], "Z"]
    OLS_range4 = OLS_calculate(y, x, z)

    df1 = round_df_tup(OLS_range1)
    df2 = round_df_tup(OLS_range2)
    df3 = round_df_tup(OLS_range3)
    df4 = round_df_tup(OLS_range4)

    figs = make_figures(df1, df2, df3, df4, ranges, arrhenius1, geometry, input_filepath)

    return figs

