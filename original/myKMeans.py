from sklearn.cluster import KMeans
import numpy as np

def executaKMeans(cluster_atual, qtd_clusters, centroid):
    x = []
    # print (cluster_atual.cluster)
    for k in cluster_atual.cluster:
        # print (k.cluster)
        aux = np.array(k.histo.T[1:2][0])
        for i in range(len(aux)):
            x.append([i,aux[i]])

    kmeans = KMeans(n_clusters=2, random_state=0)
    kmeans.fit(x)
    aux = np.array(centroid.histo.T[1:2][0])
    pre = []
    for i in range(len(aux)):
        pre.append([i,aux[i]])    
    p = kmeans.predict(pre)
    return p