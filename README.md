# clustering Kmeans (PyKmClust)
* Nome: Lucas Barbosa Rocha
* Disciplina: Inteligência Artificial
* Trabalho: Implementar um clustering para sequências de DNA utilizando a ideia do K-means.
* Contato: lucas.lb.rocha@gmail.com
## Objetivo: 
Implementar a ideia do trabalho: Application of k-means clustering algorithm in grouping the DNA sequences of hepatitis B virus (HBV). O foco foi implementar (fazer funcionar) um clusteriazador com a ideia de K-means e utilizar para comparar com a quantidade de clusters gerado pelo outro clusterizador utilizando MeanShift. 
Eu implementei a ideia do artigo. Não foi possível (no momento) implementar utilizando a biblioteca do sklearn para k-means. O artigo entregou um fluxograma do funcionamento do k-means dele e funcionou da forma que ele explicou e utilizando as fórmulas propostas no artigo. Mas no final, é a ideia do k-means sendo utilizada. Eu implementei em python. A ideia é definir k clusters, definir os centrois e colocar as sequências mais próximas do centroid no cluser correto, assim como atualizar o centroid após inserirmor uma sequência dentro do cluster.

## Problema: 
Cuidado com o tamanho do arquivo de entrada. A utilização da memória não está otimizada, ou seja, estou criando várias variáveis e o consumo de memória está alta.

