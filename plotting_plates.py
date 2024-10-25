# import statements
import pandas as pd
from pandas import DataFrame
from numpy import nan
import csv
from matplotlib import pyplot as plt
import os

# read the assay files and load the data into a data frame
def load_df(filename):
    f = open(filename)
    data_raw = f.read()
    f.close()

    data1 = data_raw.split('\n')
    data2 = data1[1:]
    header = data1[0]
    indices = header.split('\t')
    data3 = [i.split('\t') for i in data2]
    data = DataFrame(data3, columns = indices)
    data = data.replace('', nan)
    data = data.replace('None', nan)
    data = data.dropna()
    df = data.iloc[:150] # trim to consider about 24 hrs
    for col in df: # type conversions for all columns
        if col == 'Time' or col == 'T 500':
            continue
        df.loc[:,col] = df[col].astype(float)
    return df 

# works for 384 (16x24) welled plate
def graph_df(df, red = []):
    fig, ax = plt.subplots(nrows = 16, ncols = 24, sharex = True, sharey = True, figsize = (20, 20))
    cols = list(df.columns[2:])
    for col in df.columns[2:]: # exclude the non-well columns
        # print(col, end=' ')
        index = cols.index(col)
        ax[index//24][index%24].plot(df.index, df[col])
        ax[index//24][index%24].get_xaxis().set_visible(False)
        ax[index//24][index%24].get_yaxis().set_visible(False)
    
    for col in red:
        index = cols.index(col)
        ax[index//24][index%24].set_facecolor('xkcd:light pastel green')
        ax[index//24][index%24].set_title(col)

# inspect 4 wells
def inspect_well(df, w1, w2, w3, w4):
    fig, ax = plt.subplots(nrows = 2, ncols = 2, sharex = True, sharey = True, figsize = (4,4))
    ax[0,0].plot(df.index, df[w1])
    ax[0,0].get_xaxis().set_visible(False)
    ax[0,0].get_yaxis().set_visible(False)
    ax[0,1].plot(df.index, df[w2])
    ax[0,1].get_xaxis().set_visible(False)
    ax[0,1].get_yaxis().set_visible(False)
    ax[1,0].plot(df.index, df[w3])
    ax[1,0].get_xaxis().set_visible(False)
    ax[1,0].get_yaxis().set_visible(False)
    ax[1,1].plot(df.index, df[w4])
    ax[1,1].get_xaxis().set_visible(False)
    ax[1,1].get_yaxis().set_visible(False)