import pandas as pd
import itertools

def sample_rows(df_, sample_size_):
    df_result = pd.DataFrame()
    for i in range(0,10):
        sample = df_.sample(n=sample_size_)
        abundance = sample.value_counts()
        abundance = abundance.sort_index()
        while(len(abundance) <= 14):
            temp = pd.Series([0])
            abundance = pd.concat([abundance, temp])

        df_result['Abundance ' + str(i + 1)] = abundance.values
    
    return df_result

def plot_gene(df_):
    df_['mean'] = df_.mean(axis=1)
    df_['variance'] = df_.var(axis=1)

    colors = itertools.cycle(["r", "b", "g"])
    ax = df_.plot.scatter(x='mean',
                                y='variance',
                                color=next(colors))




#load data
df = pd.read_csv('/home/moneera/development/genomics/assignment_5/data/data1.txt')
#sample 
result = sample_rows(df, 500)
#plot mean vs variance
plot_gene(result)


#sample 
result = sample_rows(df, 5000)
#plot mean vs variance
plot_gene(result)

#sample 
result = sample_rows(df, 50000)
#plot mean vs variance
plot_gene(result)





