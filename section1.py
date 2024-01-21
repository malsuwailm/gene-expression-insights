#Approach 1 - General Cluster Visualization

from sklearn.preprocessing import StandardScaler
import seaborn as sns
from sklearn.cluster import KMeans
sns.set(style='whitegrid')

df = pd.read_csv('data/expression.txt', delimiter='\t', index_col=0)
scaler = StandardScaler()
segmentation_std = scaler.fit_transform(df)

kmeans = KMeans(n_clusters=3, init='k-means++', random_state = 42)
kmeans.fit(segmentation_std)

df_segm_kmeans = pd.DataFrame(segmentation_std)
df_segm_kmeans['Segment K-means'] = kmeans.labels_

df_segm_kmeans['Segment'] = df_segm_kmeans['Segment K-means'].map({0:'1st cluster',
                                                                            1: '2nd cluster',
                                                                            2: '3rd cluster'})

sns.scatterplot(data=pd.DataFrame(segmentation_std))


#Approach 2 - Determine genes with most expression

import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
from tslearn.clustering import TimeSeriesKMeans

df = pd.read_csv('data/expression.txt', delimiter='\t')
df = df.dropna().reset_index().drop('index',axis=1)

#encode gategory to number gene
def encode(x):
    x = x.split('_')
    x = x[1]
    return int(x)

df['name'] = df['name'].apply(encode)
df= df.dropna().reset_index().drop('index', axis=1)

#vizualize
genes_list = df.columns.tolist()[1:]
plt.plot(df['name'],df[genes_list[4]])

#clustering
undersample_data = df.loc[np.linspace(df.index.min(),df.index.max(),100).astype(int)]
data_array = np.array(df.T.drop('name').values)

model = TimeSeriesKMeans(n_clusters=3, metric="dtw", max_iter=10)
model.fit(data_array)
expression_list = undersample_data.T.drop('name').index.tolist()

#apply fitted model to dataset
y = model.predict(data_array)

x = undersample_data.name

plt.figure(figsize=(20,20))
k_dict = {'1':0,'2':0,'3':0,'4':1,'5':1,'6':1,'7':2,'8':2,'9':2}
colors = ['navy']*3+['darkorange']*3+['k']*3
Names = ['Class 0']*3+['Class 1']*3+['Class 2']*3
for j in range(1,10):
    plt.subplot(3,3,j)
    k = np.random.choice(np.where(y==k_dict[str(j)])[0])
    plt.plot(x,y,'.',color=colors[j-1])
    plt.ylabel('Expression',fontsize=20)
    plt.xlabel('Gene',fontsize=20)
    plt.title('Gene=%s, Class = %s'%(expression_list[k],Names[j-1]),fontsize=20)
    plt.ylim(data_array.min(),data_array.max())

plt.figure(figsize=(20,20))
k_dict = {'1':0,'2':0,'3':0,'4':1,'5':1,'6':1,'7':2,'8':2,'9':2}
colors = ['navy']*3+['darkorange']*3+['k']*3
Names = ['Class 0']*3+['Class 1']*3+['Class 2']*3
for j in range(1,10):
    plt.subplot(3,3,j)
    k = np.random.choice(np.where(y==k_dict[str(j)])[0])
    plt.hist(data_array[k],color=colors[j-1])
    plt.ylabel('Expression',fontsize=20)
    plt.xlabel('Gene',fontsize=20)
    plt.title('Gene=%s, Class = %s'%(expression_list[k],Names[j-1]),fontsize=20)
    plt.xlim(data_array.min(),data_array.max())


from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import seaborn as sns
sns.set(style='whitegrid')

df = pd.read_csv('data/expression.txt', delimiter='\t', index_col=0)
scaler = StandardScaler()
segmentation_std = scaler.fit_transform(df)

pcamodel = PCA()
pca = pcamodel.fit(segmentation_std)
scores_pca = pca.transform(segmentation_std)

#kmeans clustering with PCA
#fitting kmeas using transformed data from PCA

#implement
kmeans_pca = KMeans(n_clusters=3, init='k-means++', random_state = 42)

#fit data with k-means pca model
kmeans_pca.fit(scores_pca)

#analyze result

#create new dataframe with original features and now with PCA scores and clusters
df_segm_pca_kmeans = pd.concat([df.reset_index(drop=True), pd.DataFrame(scores_pca)], axis=1)
df_segm_pca_kmeans.columns.values[-3:] = ['Component 1', 'Component 2', 'Component 3']

#last column with pca k-means clustering labels
df_segm_pca_kmeans['Segment K-means PCA'] = kmeans_pca.labels_

#map 3 clusters to segments
df_segm_pca_kmeans['Segment'] = df_segm_pca_kmeans['Segment K-means PCA'].map({0:'1st cluster',
                                                                            1: '2nd cluster',
                                                                            2: '3rd cluster'})

#visualize first 2 components
x_axis = df_segm_pca_kmeans['Component 2']
y_axis = df_segm_pca_kmeans['Component 1']

sns.scatterplot(data=df_segm_pca_kmeans, x='Component 2', y='Component 1', hue = 'Segment')
plt.title('Clusters based on PCA Components')
plt.show()


df = pd.read_csv('data/expression.txt', delimiter='\t')
df['name'] = df['name'].apply(encode)

df['Cluster'] = kmeans.labels_
df = df.sort_values(by='Cluster')

dataplot = sns.heatmap(df.corr())

plt.show()

from sklearn.manifold import TSNE

df = pd.read_csv('data/expression.txt', delimiter='\t', index_col=0)

tsne = TSNE(n_components=3, verbose=1, random_state=42)
z = tsne.fit_transform(df) 

sns.scatterplot(data=z).set(title="Gene Expresssion T-SNE projection") 

from umap import UMAP
import plotly.express as px

umap_2d = UMAP(n_components=2, init='random', random_state=0)

proj_2d = umap_2d.fit_transform(df)

fig_2d = px.scatter(
    proj_2d, x=0, y=1)

fig_2d.show()
