# ###################################################################################
# Nome: Lucas Barbosa Rocha
# Disciplina: Inteligência Artificial
# Trabalho: Implementar um clustering para sequências de DNA 
#           utilizando KMeans.
# Contato: lucas.lb.rocha@gmail.com
# Git: Lucasbarbosarocha
#
# Classe main
# Objetivo: implementar o trabalho: MeShClust: an intelligent tool for clustering DNA
#           sequences. O foco foi implementar (fazer funcionar). Os autores utilizaram
#           c++ e fizeram implementações na mão de muitos pontos. Eu trouxe para o python
#           e tentei utilizar as bibliotecas do python. Alguns pontos não estão 
#           apresentando os mesmos resultados porque ainda estou em processo de
#           interpretação de como os autores estão fazendo no trabalho deles.
# Problema: Cuidado com o tamanho do arquivo de entrada. A utilização da memória não 
#           está otimizada, ou seja, estou criando várias variáveis e o consumo de 
#           memória está alta.
# ###################################################################################

# ###################################################################################
# importações
# ###################################################################################
import os
import requests
import numpy as np
from histograma import *
from cluster import *
from myKMeans import *
import random
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn import linear_model
from sklearn.linear_model import LogisticRegression
from myKMeans import *
from sklearn.metrics import accuracy_score, jaccard_score
from sklearn.linear_model import SGDClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

# Calculando a distância de uma sequência para o centroid
def equacao1(centroid, sequencia):
    a = centroid.histo.T[1:2][0]
    b = sequencia.histo.T[1:2][0]

    soma = 0
    for i in range(len(a)):
        soma = soma + ((a[i] - b[i])**2)
    return np.sqrt(soma)

# Calculando a distância de n sequências  paara um centroid e determinado o cluster ideal
def equacao2(clusters, sequencia, meio):
    m = 10000
    p = 0
    for i in clusters:
        aux = equacao1(i.centroid, sequencia)
        if (aux < m):
            m = np.minimum(m, aux)
            p = p + 1
    return m, (p - 1)

# Atualizando os valores do centroid
def equacao3(clusters, c, sequencia_old):
    aux = 1
    qtd_membros = 0
    for j in clusters[c].cluster:
        qtd_membros = qtd_membros + 1
        if (aux == 1):
            soma = j.histo.T[1:2][0]
            aux = 2
        else:
            soma = soma + j.histo.T[1:2][0]
    soma = soma / len(clusters[c].cluster)
    clusters[c].centroid.histo.T[1:2] = soma

# ###################################################################################
# main()
# ###################################################################################
# Abrindo arquivo
nomeArquivo = "hbv.fasta"
arquivo = open(nomeArquivo)
k = 2
kmer = 3
# Convertendo 'todas' as sequências para histograma
sequencias = []
print ("### Convertendo sequências para k-mer histograma.")
qtd_sequencias = 0
comprimentos = []
while True:
    nome = arquivo.readline()
    if len(nome) == 0:
        break
    sequencia = arquivo.readline()
    histo = sequenceToHistograma(nome, sequencia)
    comprimentos.append(histo.length)
    sequencias.append(histo)
    qtd_sequencias = qtd_sequencias + 1
# Fechando arquivo
arquivo.close()
print ("==> ",qtd_sequencias, " Sequência convertidas!")

if (k > qtd_sequencias):
    print("k (k-means) <= qtd_sequencias")

comprimentos.sort()
meio = (comprimentos[0] + comprimentos[len(comprimentos)-1]) / 2
comprimentos_menores = []
comprimentos_maiores = []
for i in comprimentos:
    if (i < meio):
        comprimentos_menores.append(i)
    else:
        comprimentos_maiores.append(i)

# Melhorando o algoritmo, no original as sequências são random
# forçar o centroid ser menor e maior em outros clusteres é melhor para o algoritmo
clusters = []
for i in (range(k)):
    if (i % 2 != 0): #max
        tam = int(len(comprimentos_maiores)/2)
        if (len(comprimentos_maiores) > 0):
            comp = comprimentos_maiores[tam]
            comprimentos_maiores.remove(comp)
    else: # min
        tam = int(len(comprimentos_menores)/2)
        if (len(comprimentos_menores) > 0):
            comp = comprimentos_menores[tam]
            comprimentos_menores.remove(comp)
    pos = 0
    for j in range(len(sequencias)):
        if comp == sequencias[j].length:
            pos = j 
    # pos = (random.randint(0,qtd_sequencias))
    # pos = 0
    cluster_aux = Cluster(sequencias[pos], sequencias[pos])
    sequencias.remove(sequencias[pos])
    clusters.append(cluster_aux)

while True:
    if (len(sequencias) == 0):
        break
    menor, p = equacao2(clusters,sequencias[0], meio)
    clusters[p].insereCluster(sequencias[0])
    equacao3(clusters, p, sequencias[0])

    if (len(sequencias)> 0):
        sequencias.remove(sequencias[0])


# imprimindo resultados
print("### Escrevendo no arquivo de saída.")
auxNome = "outputk"+str(k)+".clstr"
arquivoSaida = open(auxNome, 'w')
qtd_cluster = 0
for i in clusters:
    # print ('>Cluster ', qtd_cluster)
    e = ">Cluster " + str(qtd_cluster)
    arquivoSaida.write(e+"\n")
    aux = 0
    for j in i.cluster:
        if (j.nome == i.centroid.nome):
            # print("%d %dnt, %s *\n" %(aux, j.length, j.nome[0:len(j.nome)-1]), end='') 
            e = str(aux) + " " + str(j.length)+"nt, " + j.nome[0:len(j.nome)-1] + " *"
            arquivoSaida.write(str(e)+"\n")
        else:
            # print("%d %dnt, %s" %(aux, j.length, j.nome), end='') 
            e = str(aux) + " " + str(j.length) +"nt, "+ j.nome[0:len(j.nome)-1]
            arquivoSaida.write(str(e)+"\n")
        aux = aux + 1
    qtd_cluster = qtd_cluster + 1
print("==> Done!")
print ("==>", auxNome, " criado!")
arquivoSaida.close()
