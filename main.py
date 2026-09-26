# ###################################################################################
# Nome: Lucas Barbosa Rocha
# Disciplina: Inteligência Artificial
# Trabalho: Implementar um clustering para sequências de DNA utilizando KMeans.
# Contato: lucas.lb.rocha@gmail.com
# Git: Lucasbarbosarocha
#
# Objetivo: implementar a ideia do trabalho "Application of k-means clustering
#           algorithm in grouping the DNA sequences of hepatitis B virus (HBV)".
#           Cada sequência vira um vetor de frequência de k-mers e o k-means agrupa
#           esses vetores. A saída segue o formato .clstr do CD-HIT.
#
# Uso: python main.py hbv.fasta -k 2
# ###################################################################################

import argparse
import time

from kmer import carregar, escrever_clstr
from kmeans import kmeans, representantes


def main():
    parser = argparse.ArgumentParser(description="Clustering de sequências de DNA com k-means.")
    parser.add_argument("entrada", nargs="?", default="hbv.fasta", help="arquivo FASTA")
    parser.add_argument("-k", "--clusters", type=int, default=2, help="quantidade de clusters")
    parser.add_argument("--kmer", type=int, default=3, help="tamanho do k-mer (memória: 4^k floats por sequência)")
    parser.add_argument("--n-init", type=int, default=10, help="execuções com sementes diferentes")
    parser.add_argument("--max-iter", type=int, default=300)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("-o", "--saida", help="arquivo .clstr (padrão: outputk<k>.clstr)")
    args = parser.parse_args()
    saida = args.saida or f"outputk{args.clusters}.clstr"

    inicio = time.perf_counter()
    print("### Convertendo sequências para vetores de k-mers.")
    nomes, comprimentos, X = carregar(args.entrada, args.kmer)
    print(f"==> {len(nomes)} sequências convertidas ({X.nbytes / 2**20:.2f} MB).")

    print(f"### Rodando k-means com k={args.clusters}.")
    labels, centros, inercia, iteracoes = kmeans(
        X, args.clusters, n_init=args.n_init, max_iter=args.max_iter, seed=args.seed)
    print(f"==> Convergiu em {iteracoes} iterações, inércia {inercia:.6f}.")

    print("### Escrevendo no arquivo de saída.")
    escrever_clstr(saida, nomes, comprimentos, labels, representantes(X, labels, centros))
    print(f"==> {saida} criado! ({time.perf_counter() - inicio:.2f}s)")


if __name__ == "__main__":
    main()
