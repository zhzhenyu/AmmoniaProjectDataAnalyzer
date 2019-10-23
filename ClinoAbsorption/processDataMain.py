"""
this code processes absorption/desorption csv file
by first reading the data
second building a dataframe using time as index column and relavant dependent columns
third calculating important metrics including breakthrough capacity and absorbed/desorbed quantities
fourth plot
"""

import csv
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

def builddf(filename, path):
    """
    @summary: read csv file and return a dataframe
    @param filename: csv filename
    @param path: csv file path in sub folder Data
    :return: df
    """
    rows = []

    with open(os.path.join(path, filename)) as csvfile:
        csvreader = csv.reader(csvfile)

        for row in csvreader:
            rows.append(row)

    # build dataframe
    df_dict = {'time': 0, 'nh3PPM': 0, 'temperature': 0, 'n2flow': 0, 'nh3flow': 0}
    df_list = []
    for row in rows[5:]:
        rowsplit = row[0].split()
        nh3PPM = float(rowsplit[1])
        temperature = float(rowsplit[2])
        n2flow = float(rowsplit[3]) * 1000  # [slpm to sccm]
        nh3flow = float(rowsplit[4])  # [sccm]

        # calculate time
        timesplit = rowsplit[0].split(':')
        timesplit = [float(i) for i in timesplit]
        time = timesplit[0] * 24 * 3600 + timesplit[1] * 3600 + timesplit[2] * 60 + timesplit[
            3]  # [day:hour:min:sec] to secs

        df_dict.update(
            {'time': time, 'nh3PPM': nh3PPM, 'temperature': temperature, 'n2flow': n2flow, 'nh3flow': nh3flow})
        df_list.append(df_dict.copy())

    df = pd.DataFrame(df_list)
    df['time'] = df['time'] - df.iloc[0]['time']  # reorigin the dataframe
    df.set_index('time', inplace=True)  # use time as index col

    return df

def normdf(df):
    """summary: normalize nh3PPM and temperature"""
    df['temperature'] = df['temperature'] / df['temperature'].max()
    df['nh3PPM'] = df['nh3PPM'] / df['nh3PPM'].max()
    df['nh3flow'] = df['nh3flow'] / df['nh3flow'].max()
    return df

def movingAveragedf(df, window):
    """summary: take the moving average of df and backfill the data"""
    df['nh3PPM'] = df['nh3PPM'].rolling(window).mean()
    df.bfill(inplace=True)
    return df

def plot(df, tstart, colList, title, xlabel, ylabel, legendList, fontname, fontsize):
    ax = df.loc[tstart:, colList].plot(fontsize=fontsize)
    ax.set_xlabel(xlabel, fontname=fontname, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontname=fontname, fontsize=fontsize)
    ax.set_title(title, fontname=fontname, fontsize=fontsize)
    for tick in ax.get_xticklabels(): tick.set_fontname(fontname)
    for tick in ax.get_yticklabels(): tick.set_fontname(fontname)
    ax.legend(legendList)
    plt.show()

def plotdataframe(tstart, plotdf, plotdfNormAbsorb, plotdfNormDesorb, plotdfMA):
    if plotdf:
        colList = ['nh3PPM']
        title, xlabel, ylabel, fontname, fontsize = filename.split('_')[0], 'Time (s)', 'NH$_3$ PPM', 'Arial', 13
        legendList = ['NH$_3$ PPM']
        plot(df, tstart, colList, title, xlabel, ylabel, legendList, fontname, fontsize)

    if plotdfNormAbsorb:
        colList = ['nh3PPM', 'temperature', 'nh3flow']
        title, xlabel, ylabel, fontname, fontsize = filename.split('_')[0], 'Time (s)', 'Normalized NH$_3$ PPM/Temperature/Flowrate', 'Arial', 14
        legendList = ['NH$_3$ PPM', 'Temperature', 'NH$_3$ Flowrate']
        plot(df_norm, tstart, colList, title, xlabel, ylabel, legendList, fontname, fontsize)

    if plotdfNormDesorb:
        colList = ['nh3PPM', 'temperature']
        title, xlabel, ylabel, fontname, fontsize = filename.split('_')[0], 'Time (s)', 'Normalized NH$_3$ PPM/Temperature/Flowrate', 'Arial', 14
        legendList = ['NH$_3$ PPM', 'Temperature']
        plot(df_norm, tstart, colList, title, xlabel, ylabel, legendList, fontname, fontsize)

    if plotdfMA:
        colList = ['nh3PPM']
        title, xlabel, ylabel, fontname, fontsize = filename.split('_')[0], 'Time (s)', 'Normalized NH$_3$ PPM', 'Arial', 14
        legendList = ['NH$_3$ PPM']
        plot(df_ma, tstart, colList, title, xlabel, ylabel, legendList, fontname, fontsize)

def breakthroughCapacity(df):
    t_start = df[df['nh3flow'] > 4].index.tolist()[0] # when the NH3 flowrate is above 4 sccm
    t_end = df.loc[ (df['nh3PPM'] > 20) & (df['nh3flow'] > 4) ].index.tolist()[0] # when the NH3 PPM is above 20
    t_blank = 6
    # print(df.loc[t_start:t_end, :])
    nh3_molarrate = 500*0.000000745*10800*0.000001 # [mol/s]
    nh3_g = nh3_molarrate*(t_end - t_start - t_blank)*17 # [g]
    return  nh3_g

def totalAbsorb(df):
    t_start = df[df['nh3flow'] > 4].index.tolist()[0] # when the NH3 flowrate is above 4 sccm
    t_end = df[df['nh3PPM'] > 20].index.tolist()[0] # when the NH3 PPM is above 20
    t_blank = 6
    nh3_PPM = df.loc[t_start:t_end, 'nh3PPM'].tolist()
    nh3_molarrate = sum([500*0.000000745*(10800-i)*0.000001 for i in nh3_PPM]) # [mol/s]
    nh3_g = nh3_molarrate*(t_end - t_start - t_blank)*17 # [g]
    return  nh3_g


if __name__ == "__main__":
    # Read in csv file
    filename = "ClinoTestDesorp5_VCR_10-21-2019.csv"
    path = "Z:\AmmoniaSynthesis\ReactionAbsorption\PythonScript\Data"
    window = 20
    m_absorbent = 0.2035

    """build dataframe"""
    df = builddf(filename, path)
    df_norm = normdf(df.copy())
    df_ma = movingAveragedf(df.copy(), window)

    """Plot dataframe"""
    t_start = df.index.tolist()[0]
    plotdf, plotdfNormAbsorb, plotdfNormDesorb, plotdfMA = 1, 0, 1, 0
    if 'Absorp' in filename.split('_')[0]:
        t_start = df[df['nh3flow'] > 4].index.tolist()[0] - 12.0 # when the NH3 flowrate is above 4 sccm
        plotdf, plotdfNormAbsorb, plotdfNormDesorb, plotdfMA = 0, 1, 0, 0
    plotdataframe(t_start, plotdf, plotdfNormAbsorb, plotdfNormDesorb, plotdfMA)

    """calculate breakthrough capacity and total absorbed/desorbed
       max NH3 PPM and NH3 flow delta_t = 6 s
       max Temperature and NH3 flow delta_t = 4 s
       derived from ClinoTestAbsorb2
    """
    if 'Absorp' in filename.split('_')[0]:
        m_nh3 = breakthroughCapacity(df.copy())
        weightper = m_nh3/m_absorbent * 100.00 # %
        print('Breakthrough capacity is {}%'.format(round(weightper, 2)))

    elif 'Desorp' in filename.split('_')[0]:
        print('Desorb')

