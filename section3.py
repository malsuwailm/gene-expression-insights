import pandas as pd
from scipy.stats import ttest_ind
from scipy import stats
import numpy as np

def sample_rows(df_, sample_size, dataset_n, sampling_times):
    df_result = pd.DataFrame()
    for i in range(0, sampling_times):
        sample = df_.sample(n=sample_size)
        abundance = sample.value_counts()
        abundance = abundance.sort_index()
        while(len(abundance) <= 14):
            temp = pd.Series([0])
            abundance = pd.concat([abundance, temp])

        df_result['(' + str(dataset_n) + ')'+' Abundance ' + str(i + 1)] = abundance.values

    return df_result

def paired_ttest(df_1, df_2):
    result = stats.ttest_ind(df_1, df_2,axis=1)
    return result.statistic, result.pvalue

def foldchange(df_1, df_2):
    df = pd.DataFrame()
    df['mean1'] = df_1.mean(axis=1)
    df['mean2'] = df_2.mean(axis=1)

    foldchange = []
    for index, row in df.iterrows():
        foldchange.append(row['mean2']/row['mean1'])

    return df['mean1'], df['mean2'], foldchange

def significant_genes(df):
    result_list = []
    sig_genes = {}

    for index, row in df.iterrows():
        if float(row['P-Value']) <= 0.05:
            result_list.append(True)
            sig_genes.update({row['Gene'] : row['Fold Change']})
        else:
            result_list.append(False)

    df['Significant Genes'] = result_list
    return sig_genes

df_1 = pd.read_csv('/home/themokona/development/genomics/assignment_5/data/data1.txt')
df_2 = pd.read_csv('/home/themokona/development/genomics/assignment_5/data/data2.txt')

#5000 samples, 3 times
sample_1 = sample_rows(df_1, sample_size=5000, dataset_n=1, sampling_times=3)
sample_2 = sample_rows(df_2, sample_size=5000, dataset_n=2, sampling_times=3)

#get data
st, p_value = paired_ttest(sample_1, sample_2)
mean1, mean2, fold = foldchange(sample_1, sample_2)

#combine into one table
df = pd.DataFrame()
df = pd.concat([sample_1, sample_2], axis=1)

gene_number = []
for i in range(1,16):
    gene_number.append(i)

df['Gene'] = gene_number
df['(1) Mean'] = mean1
df['(2) Mean'] = mean2
df['Fold Change'] = fold
df['P-Value'] = p_value
df['Statistic'] = st

sig_genes = significant_genes(df)

print('Significant Genes:\n')
for gene in sig_genes:
    print(int(gene), '->', str(round(sig_genes[gene], 2)))


import plotly.graph_objects as go
import plotly.express as px

def log_values(df):
    df['Log2FC'] = np.log2(df['Fold Change'])
    df['-LogPVal'] = np.log10(df['P-Value']) * (-1)

def volcano_plot(df, sample_size, sampling_times):
    fig = px.scatter(df, x='Log2FC', y='-LogPVal', color='Significant Genes', text='Gene')
    fig.update_traces(textposition='top center', marker_size=20)
    fig.update_layout(title_text='Volcano plot for Gene Expression using ' + str(sample_size) + ' rows sampled ' + str(sampling_times) + ' times')
    fig.show()


log_values(df)
volcano_plot(df, sample_size=5000, sampling_times=3)

#5000 samples, 10 times
sample_1 = sample_rows(df_1, sample_size=5000, dataset_n=1, sampling_times=10)
sample_2 = sample_rows(df_2, sample_size=5000, dataset_n=2, sampling_times=10)

#get data
st, p_value = paired_ttest(sample_1, sample_2)
mean1, mean2, fold = foldchange(sample_1, sample_2)

#combine into one table
df = pd.DataFrame()
df = pd.concat([sample_1, sample_2], axis=1)

gene_number = []
for i in range(1,16):
    gene_number.append(i)

df['Gene'] = gene_number
df['(1) Mean'] = mean1
df['(2) Mean'] = mean2
df['Fold Change'] = fold
df['P-Value'] = p_value
df['Statistic'] = st

sig_genes = significant_genes(df)

print('Significant Genes:\n')
for gene in sig_genes:
    print(int(gene), '->', str(round(sig_genes[gene], 2)))


log_values(df)
volcano_plot(df, sample_size=5000, sampling_times=10)

#50000 samples, 3 times
sample_1 = sample_rows(df_1, sample_size=50000, dataset_n=1, sampling_times=3)
sample_2 = sample_rows(df_2, sample_size=50000, dataset_n=2, sampling_times=3)

#get data
st, p_value = paired_ttest(sample_1, sample_2)
mean1, mean2, fold = foldchange(sample_1, sample_2)

#combine into one table
df = pd.DataFrame()
df = pd.concat([sample_1, sample_2], axis=1)

gene_number = []
for i in range(1,16):
    gene_number.append(i)

df['Gene'] = gene_number
df['(1) Mean'] = mean1
df['(2) Mean'] = mean2
df['Fold Change'] = fold
df['P-Value'] = p_value
df['Statistic'] = st

sig_genes = significant_genes(df)

print('Significant Genes:\n')
for gene in sig_genes:
    print(int(gene), '->', str(round(sig_genes[gene], 2)))

log_values(df)
volcano_plot(df, sample_size=50000, sampling_times=3)


#50000 samples, 10 times
sample_1 = sample_rows(df_1, sample_size=50000, dataset_n=1, sampling_times=10)
sample_2 = sample_rows(df_2, sample_size=50000, dataset_n=2, sampling_times=10)

#get data
st, p_value = paired_ttest(sample_1, sample_2)
mean1, mean2, fold = foldchange(sample_1, sample_2)

#combine into one table
df = pd.DataFrame()
df = pd.concat([sample_1, sample_2], axis=1)

gene_number = []
for i in range(1,16):
    gene_number.append(i)

df['Gene'] = gene_number
df['(1) Mean'] = mean1
df['(2) Mean'] = mean2
df['Fold Change'] = fold
df['P-Value'] = p_value
df['Statistic'] = st

sig_genes = significant_genes(df)

print('Significant Genes:\n')
for gene in sig_genes:
    print(int(gene), '->', str(round(sig_genes[gene], 2)))

log_values(df)
volcano_plot(df, sample_size=50000, sampling_times=10)

#Entire data set, 3 times
sample_1 = sample_rows(df_1, sample_size=len(df_1.index), dataset_n=1, sampling_times=3)
sample_2 = sample_rows(df_2, sample_size=len(df_2.index), dataset_n=2, sampling_times=3)

#get data
st, p_value = paired_ttest(sample_1, sample_2)
mean1, mean2, fold = foldchange(sample_1, sample_2)

#combine into one table
df = pd.DataFrame()
df = pd.concat([sample_1, sample_2], axis=1)

gene_number = []
for i in range(1,16):
    gene_number.append(i)

df['Gene'] = gene_number
df['(1) Mean'] = mean1
df['(2) Mean'] = mean2
df['Fold Change'] = fold
df['P-Value'] = p_value
df['Statistic'] = st

sig_genes = significant_genes(df)

print('Significant Genes:\n')
for gene in sig_genes:
    print(int(gene), '->', str(round(sig_genes[gene], 2)))


#Entire data set, 10 times
sample_1 = sample_rows(df_1, sample_size=len(df_1.index), dataset_n=1, sampling_times=10)
sample_2 = sample_rows(df_2, sample_size=len(df_2.index), dataset_n=2, sampling_times=10)

#get data
st, p_value = paired_ttest(sample_1, sample_2)
mean1, mean2, fold = foldchange(sample_1, sample_2)

#combine into one table
df = pd.DataFrame()
df = pd.concat([sample_1, sample_2], axis=1)

gene_number = []
for i in range(1,16):
    gene_number.append(i)

df['Gene'] = gene_number
df['(1) Mean'] = mean1
df['(2) Mean'] = mean2
df['Fold Change'] = fold
df['P-Value'] = p_value
df['Statistic'] = st

sig_genes = significant_genes(df)

print('Significant Genes:\n')
for gene in sig_genes:
    print(int(gene), '->', str(round(sig_genes[gene], 2)))


