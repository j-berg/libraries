#!/usr/bin/python

"""
Statistical plotting package of common routines
Author: Jordan Berg
Version: 0.0.1

***Genes must be rows in dataframes***

Add to path for importability in python as:
import my_functions as mf

To do:
Documentation with required formats for datatypes
Test on sample datasets -- create a 10*10 df for testing
finish checking fig_titles, etc if not None
"""

# Import required libraries
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress


# Data format type validation
check_format(data, type):

    #Check input data format types are correct
    if isinstance(data, type) == False:
        print(type + 'type required')
        sys.exit(1)


# RPM normalization
def rpm_convert(df):

    #Check inputs
    check_format(df, pd.DataFrame)

    #Normalize dataframe by sample RPM metric
    df_rpm = df / (df.sum() / 1e6)
    df = df_rpm

    return df


# Translation efficiency normalization
def te_convert(df, log2=True, samples):

    #Remind user of correct formatting scheme for proper FUNCTIONS
    print('')

    #Check inputs
    check_format(df, pd.DataFrame)
    check_format(log10, bool)
    check_format(samples, list)

    #Prepare blank dataframe based on user input dataframe
    df_te = df
    df_te = df_te[df_te.columns[[]]]

    #Initialize division scheme for fp/rna
    y = 0
    z = 1

    #Perform translation efficiency calculations
    for x in samples:
        df_te[x] = pd.Series(df_te.index)
        df_te[x] = df[df.columns[y]]/df[df.columns[z]]
        y = y + 2
        z = z + 2

    #Perform log2 scaling of data
    if log2 == True:
        df_te = np.log2(df_te)

    df = df_te

    return df


# First column zero normalization
def zero_convert(df):

    #Check inputs
    check_format(df, pd.DataFrame)

    #Scale to first column zero
    df_zero = df.sub(df[df.columns[0]], axis=0)

    df = df_zero

    return df


# Count cut-off trimmming
def cutoff(df, min_counts):

    #Check inputs
    check_format(df, pd.DataFrame)
    check_format(min_counts, int)

    #Trim dataframe of rows with
    df = df.T[df.T.columns[df.T.min() > min_counts]]
    df = df.T

    return df


# Heatmap-based scaling (mean=0, stdev=1)
def heat_scaling(df):

    #Check inputs
    check_format(df, pd.DataFrame)

    # Drop proteins with NA values
    df = df.dropna(axis=0)

    # Scale by gene rows
    dfT = df.T #code requires to scale by column, so have to do a swap-a-roo
    df_norm_col=(dfT -dfT.mean())/dfT.std()
    df_norm_col = df_norm_col.T

    df = df_norm_col

    return df


# Clustered heatmap
def basic_heat(df, cluster_title=None, scaled=True, col_clust=False, x_ax=4, y_ax=10):

    #Check inputs
    check_format(df, pd.DataFrame)
    check_format(cluster_title, str)
    check_format(scaled, bool)
    check_format(col_clust, bool)
    check_format(x_ax, int)
    check_format(y_ax, int)

    if scaled == True:
        heat_scaling(df)

    #Generate clustermap
    ax = sns.clustermap(df, cmap="RdBu_r",center=0,
                   xticklabels=True,linewidths=.1, linecolor='#DCDCDC',
                   col_cluster=col_clust,figsize=(x_ax,y_ax))

    #Add title
    ax = plt.title(cluster_title)

    #Save heatmap
    plt.savefig(cluster_title + '.png')


# Basic heatmap
def basic_heat(df, heat_title=None, scaled=True):

    #Check inputs
    check_format(df, pd.DataFrame)
    check_format(heat_title, str)
    check_format(scaled, bool)

    if scaled == True:
        heat_scaling(df)

    #Generate heatmap
    ax = sns.heatmap(df, cmap="RdBu_r", center=0,
                xticklabels=True, linewidths=.1,
                linecolor='#DCDCDC')

    #Add title
    ax = plt.title(heat_title)

    #Save heatmap
    plt.savefig(heat_title + '.png')



# Volcano Plot




# Violin Plot




# Boxplot -- option for list of genes to pull from dataframe
def boxplot(df, samples, genes=None, box_title=None):

    #Check inputs
    check_format(df, pd.DataFrame)
    check_format(samples, list)
    if genes is not None:
        check_format(genes, list)
    if box_title is not None:
        check_format(box_title, str)

    #Unstack dataframe for plotting
    df_unstacked = df.unstack().reset_index()
    df_unstacked = df_unstacked.rename(columns={df_unstacked.columns[0]: 'condition', df_unstacked.columns[2]: 'log2TE'})

    #Prepare subset list of genes to plot if provided
    if samples is not None:
        df_pex = df_unstacked[df_unstacked.columns[1].isin(genes)]
        df_unstacked = df_pex

    #Plot data
    ax = sns.set_style("whitegrid")
    ax = sns.boxplot(x="condition", y="log2TE", data=df_unstacked, order=samples, width=.35)

    #Add title
    ax = plt.title(box_title)

    #Save heatmap
    plt.savefig(box_title + '.png')


# Linear Regression -- all genes against one (interest_gene)
def big_linreg(df, interest_gene, dataset_name):

    #Check inputs
    check_format(df, pd.DataFrame)
    check_format(gene, str)
    check_format(dataset_name, str)

    #Convert gene of interest metrics to list
    interest = df.loc[interest_gene].values.tolist()

    lm_interest=[]

    for row in df.iterrows():
        index, data = row
        gene = df.loc[index].values.tolist()
        if len(gene) is not len(interest):
            continue
        else:
            slope, intercept, r_value, p_value, std_err = linregress(gene, interest)
            lm_interest.append([index, slope, intercept, r_value, (r_value**2), p_value, std_err])

    df_lm_interest = pd.DataFrame(lm_interest, columns=['gene', 'slope', 'intercept', 'r_value', 'r_squared', 'p_value', 'std_err'])

    # Save linear modeling metrics to .csv
    df_lm_interest.to_csv('lm_' + interest_gene + '_' + dateset_name + '.csv',sep=',')


# Linear Regression -- targets
def target_linreg(df, sample1, sample2, log10=True, fig_size=6):

    #Check inputs
    check_format(df, pd.DataFrame)
    check_format(sample1, str)
    check_format(sample2, str)
    check_format(log10, bool)
    check_format(fig_size, int)

    #Pass samples to lists
    sampl1 = df[sample1].values.tolist()
    sampl2 = df[sample2].values.tolist()

    #Scale samples on log10 if user-required
    if log10 == True:
        sampl1 = np.log10(sampl1)
        sampl2 = np.log10(sampl2)

    #Run linear regression on user input samples and collect stats
    slope, intercept, r_value, p_value, std_err = linregress(sampl1, sampl2)

    #Create dataframe of user input sample lists for graphing
    df_graph = pd.DataFrame({sample1:sampl1p, sample2:sampl2p})

    #Plot scatter of two genes of interest
    ax = df_graph.plot.scatter(x=sample1,y=sample2, c='Black', figsize=(fig_size,fig_size), s=0.5)
    ax.grid(False) #Hide gridlines

    #Plot trend line
    z = np.polyfit(sampl1, sampl2, 1)
    p = np.poly1d(z)
    ax = plt.plot(sampl1,p(sampl1),"purple")

    #Calculate r_sq with conserved sign for anticorrelation
    if r_value < 0:
        r_sq = -1 * (r_value**2)
    else:
        r_sq = r_value**2

    #Add title with r_squared value
    ax = plt.title("R$^2$ = " + str(round(r_sq, 2)))

    #Save graph
    plt.savefig(sample1 + '_' + sample2 + '_linregress.png')
