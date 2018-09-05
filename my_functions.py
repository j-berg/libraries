#!/usr/bin/python

"""
Statistical plotting package of common routines
Author: Jordan Berg
Version: 0.0.1

To do:
Test on sample datasets -- create a 10*10 df for testing
finish checking fig_titles, etc if not None
Add user-given subsetting of heatmap/clustermap
PCA
MLM
"""

# Import required libraries
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
from scipy.stats import linregress


#Print help info
def help():

    #Print library info
    print('Statistical plotting package of common routines\nAuthor: Jordan Berg\nVersion: 0.0.1\n')

    #Print usage info
    print('Usage:\nAdd to $PATH for libary importability in python as:\nimport my_functions as mf\n')

    #Print formatting instructions
    print('Required file formatting:\n***Genes must be rows in dataframe***')
    print("***Samples must be in order in columns in dataframe***\nCorrect file formatting is the user's responsibility")
    print('Pay attention to function format specification\n')

    #Print function options
    print('Functions:')
    print('Read in csv file to pandas DataFrame: read_df(csv_file_name, output_df_name)')
    print('RPM normalize ribosome profiling dataset: rpm_convert(df)')
    print('Normalize ribosome profiling dataset for translation efficiency (sample order must be sample1_fp, sample1_rna, ...): rpm_convert(df)')
    print('Normalize all sample columns to first sample: zero_convert(df)')
    print('Implement minimum counts cutoff to sequencing count table: cutoff(df, min_counts)')
    print('Normalize gene rows to mean=0, stdev=1 for heatmaps: mean_scaling(df)')
    print('Basic clustermap with gene clustering: basic_cluster(df, cluster_title=None, scaled=True, col_clust=False, x_ax={int}, y_ax={int})')
    print('Basic heatmap: basic_heat(df, heat_title=None, scaled=True)')
    print('Volcano plot (base samples (base_name) must be first columns in df, followed by comparison samples (comp_names)): volcano(df, comp_name={str}, base_name={str}, highlights{list}, figure_title={str})')
    print('violinswarm(df, samples, genes=None, title)')
    print('boxplot(df, samples, genes=None, title)')
    print('Simple Linear Regression, Iterates over all genes against your gene of interest: total_linreg(df, interest_gene, dataset_name)')
    print('Simple Linear Regression, Runs and plots: target_linreg(df, interest_gene, dataset_name)')


# Data format type validation
def check_format(data, type):

    #Check input data format types are correct
    if isinstance(data, type) == False:
        print(type + 'type required')
        sys.exit(1)


# Auto read DataFrame
def read_df(file_name):

    #Check inputs
    check_format(file_name, str)

    df_name = pd.read_csv(file_name, sep=",", index_col=0)

    return df_name


# RPM normalization
def rpm_convert(df):

    #Check inputs
    check_format(df, pd.DataFrame)

    #Normalize dataframe by sample RPM metric
    df_rpm = df / (df.sum() / 1e6)

    return df_rpm


# Translation efficiency normalization
def te_convert(df, samples, log2=True):

    #Check inputs
    check_format(df, pd.DataFrame)
    check_format(samples, list)
    check_format(log2, bool)

    #Scale all DataFrame values for division where x = x + 1e-7
    df += 1e-7

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
        df_log = np.log2(df_te)

    return df_log


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
def mean_scaling(df):

    #Check inputs
    check_format(df, pd.DataFrame)

    # Drop proteins with NA values
    df = df.dropna(axis=0)

    # Scale by gene rows
    dfT = df.T #code requires to scale by column, so have to do a swap-a-roo
    df_norm_col=(dfT -dfT.mean())/dfT.std()
    df_norm_col = df_norm_col.T

    return df_norm_col


# Clustered heatmap
def basic_cluster(df, cluster_title=None, scaled=True, col_clust=False, x_ax=4, y_ax=10):

    #Check inputs
    check_format(df, pd.DataFrame)
    check_format(cluster_title, str)
    check_format(scaled, bool)
    check_format(col_clust, bool)
    check_format(x_ax, int)
    check_format(y_ax, int)

    if scaled == True:
        df = mean_scaling(df)
    else:
        # Drop proteins with NA values
        df = df.dropna(axis=0)

    #Generate clustermap
    ax = sns.clustermap(df, cmap="RdBu_r",center=0,
                   xticklabels=True,linewidths=.1, linecolor='#DCDCDC',
                   col_cluster=col_clust,figsize=(x_ax,y_ax))

    #Add title
    ax = plt.title(cluster_title)

    #Save heatmap
    plt.savefig(cluster_title + '.png')


# Basic heatmap
def basic_heat(df, heat_title, scaled=True):

    #Check inputs
    check_format(df, pd.DataFrame)
    check_format(heat_title, str)
    check_format(scaled, bool)

    if scaled == True:
        df = mean_scaling(df)
    else:
        # Drop proteins with NA values
        df = df.dropna(axis=0)

    #Generate heatmap
    ax = sns.heatmap(df, cmap="RdBu_r", center=0,
                xticklabels=True, linewidths=.1,
                linecolor='#DCDCDC')

    #Add title
    ax = plt.title(heat_title)

    #Save heatmap
    plt.savefig(heat_title + '.png')


# Volcano Plot
def volcano(df, comp_name, base_name, highlights, figure_title):

    #Check datatypes
    check_format(df, pd.DataFrame)
    check_format(comp_name, str)
    check_format(base_name, str)
    check_format(highlights, list)
    check_format(figure_title, str)

    #Drop rows with NaN, na, etc.
    df = df.dropna()
    #df_volcano = df_volcano.drop('nan')

    #Make df of sample averages
    #<----------need to do this row by row
    df_avg = {comp_name: df.filter(regex=comp_name).mean(), base_name: df.filter(regex=base_name).mean()}

    #Calculate fold change using average
    df_volcano = np.log2(df_avg[comp_name] / df_avg[base_name])
    df_volcano = pd.DataFrame(df_volcano, columns=['log2 Fold Change'])
    df_volcano['-log10 P-Value'] = ''

    #Calculate number of samples for each name
    comp_num = int(df.filter(regex=comp_name, axis=1).count(axis=1).mean())
    base_num = int(df.filter(regex=base_name, axis=1).count(axis=1).mean())

    # Calculate p-value using 1-way ANOVA with replicates and append to df_...
    for row in df.iterrows():
        index, data = row
        comp_row = data[base_num:(base_num+comp_num)].values.tolist() #<----------how to slice this?
        base_row = data[0:base_num].values.tolist() #<----------how to slice this?

        # Append p_value to df_volcano
        statistic, p_value = stats.f_oneway(comp_row, base_row)
        df_volcano.loc[index,'-log10 P-Value'] = float(-1 * (np.log10(p_value)))

    # Scatter plot log2FoldChange vs. log10P-Value
    ax = sns.scatterplot(x='log2 Fold Change', y='-log10 P-Value', data=df_volcano, color='LightGray')
    ax = sns.scatterplot(x='log2 Fold Change', y='-log10 P-Value', data=df_volcano.loc[highlights], color='DarkRed')

    #Assign title to plot
    ax = ax.set_title(figure_title)

    #Save volcano plot
    plt.savefig(figure_title + '_volcano.png')


# Seaborns plot prep
def sns_prep(df, samples):

    #Unstack dataframe for plotting
    df_unstacked = df.unstack().reset_index()
    df_unstacked = df_unstacked.rename(columns={df_unstacked.columns[0]: 'condition', df_unstacked.columns[2]: 'expr'})

    return df_unstacked


# Violin Plot -- subset vs general
def violinswarm(df, samples, title, genes=None):

    #Check inputs
    check_format(df, pd.DataFrame)
    check_format(samples, list)
    check_format(title, str)
    if genes is not None:
        check_format(genes, list)

    #Run prep on DataFrame
    df_unstacked = sns_prep(df, samples)

    #Create DataFrame with subset genes for swarmplotting
    swarm_target = df_unstacked.loc[(df_unstacked.ix[:,1].isin(genes))]

    #Plot data, overlayings subset swarm plot with total violin plot
    sns.set_style('whitegrid')
    swarm_plot = sns.violinplot(x='condition', y='expr', data=df_unstacked, inner=None, order=samples)
    swarm_plot = sns.swarmplot(x='condition',y='expr',data=swarm_target,color='black',size=7, order=samples)

    #Save violinswarm plot
    fig = swarm_plot.get_figure()
    fig.savefig(title + '_violinswarm.png')


# Boxplot -- option for list of genes to pull from dataframe
def boxplot(df, samples, title, genes=None):

    #Check inputs
    check_format(df, pd.DataFrame)
    check_format(samples, list)
    check_format(title, str)
    if genes is not None:
        check_format(genes, list)

    #Run prep on DataFrame
    df_unstacked = sns_prep(df, samples)

    #Prepare subset list of genes to plot if provided
    if genes is not None:
        df_genes = df_unstacked[df_unstacked.ix[:,1].isin(genes)]
        df_unstacked = df_genes

    #Plot data
    ax = sns.set_style('whitegrid')
    ax = sns.boxplot(x='condition', y='expr', data=df_unstacked, order=samples, width=.35)

    #Add title
    ax = plt.title(title)

    #Save heatmap
    plt.savefig(title + '_boxplot.png')


# Linear Regression -- all genes against one (interest_gene)
def total_linreg(df, interest_gene, output_name):

    #Check inputs
    check_format(df, pd.DataFrame)
    check_format(interest_gene, str)
    check_format(output_name, str)

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
    df_lm_interest.to_csv('lm_' + interest_gene + '_' + output_name + '.csv',sep=',')


# Linear Regression -- targets
def target_linreg(df, gene_x, gene_y, log10=True, fig_size=6):

    #Check inputs
    check_format(df, pd.DataFrame)
    check_format(gene_x, str)
    check_format(gene_y, str)
    check_format(log10, bool)
    check_format(fig_size, int)

    #Pass samples to lists
    sampl1 = df.loc[gene_x].values.tolist()
    sampl2 = df.loc[gene_y].values.tolist()

    #Scale samples on log10 if user-required <----- Maybe do this after LM?
    if log10 == True:
        sampl1 = np.log10(sampl1)
        sampl2 = np.log10(sampl2)

    #Run linear regression on user input samples and collect stats
    slope, intercept, r_value, p_value, std_err = linregress(sampl1, sampl2)

    #Create dataframe of user input sample lists for graphing
    df_graph = pd.DataFrame({gene_x:sampl1, gene_y:sampl2})

    #Plot scatter of two genes of interest
    ax = df_graph.plot.scatter(x=gene_x,y=gene_y, c='Black', figsize=(fig_size,fig_size), s=0.5)
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
    plt.savefig(gene_x + '_' + gene_y + '_linregress.png')
