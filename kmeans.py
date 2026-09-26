# ###################################################################################
# Módulo kmeans
# Objetivo: k-means (Lloyd) sobre os vetores de k-mers.
#   1. Inicialização k-means++ (centróides iniciais espalhados, não dependem do
#      comprimento das sequências nem da ordem do arquivo).
#   2. Repete atribuição (centróide mais próximo) e atualização (média dos membros)
#      até os centróides pararem de mudar.
#   3. Roda n_init vezes com sementes diferentes e fica com a menor inércia.
# Os centróides são vetores próprios: nenhuma sequência é alterada.
# ###################################################################################

import numpy as np


def distancias2(X, normas2, centros):
    """Distância euclidiana ao quadrado de cada sequência para cada centróide, (N, k).

    Usa |x - c|^2 = |x|^2 - 2 x.c + |c|^2 para não criar uma matriz (N, k, 4^k).
    O produto é feito no dtype de X para não criar uma cópia float64 de X, e o resto
    é feito in-place para existir só uma matriz (N, k) em memória.
    """
    d2 = (X @ centros.T.astype(X.dtype)).astype(np.float64)
    d2 *= -2
    d2 += normas2[:, None]
    d2 += (centros * centros).sum(axis=1)[None, :]
    return np.maximum(d2, 0, out=d2)


def kmeans_pp(X, normas2, k, rng):
    n = len(X)
    centros = np.empty((k, X.shape[1]), dtype=np.float64)
    centros[0] = X[rng.integers(n)]
    d2 = distancias2(X, normas2, centros[:1])[:, 0]
    for c in range(1, k):
        total = d2.sum()
        i = rng.choice(n, p=d2 / total) if total > 0 else rng.integers(n)
        centros[c] = X[i]
        d2 = np.minimum(d2, distancias2(X, normas2, centros[c:c + 1])[:, 0])
    return centros


def _lloyd(X, normas2, centros, max_iter, tol):
    k = len(centros)
    labels = None
    for iteracao in range(1, max_iter + 1):
        d2 = distancias2(X, normas2, centros)
        anteriores, labels = labels, d2.argmin(axis=1)
        if anteriores is not None and np.array_equal(anteriores, labels):
            break  # ninguém mudou de cluster: convergiu

        # Soma dos membros de cada cluster, uma coluna (k-mer) por vez
        novos = np.column_stack([np.bincount(labels, weights=X[:, j], minlength=k)
                                 for j in range(X.shape[1])])
        tamanhos = np.bincount(labels, minlength=k)
        vazios = tamanhos == 0
        novos[~vazios] /= tamanhos[~vazios, None]
        if vazios.any():
            # Cluster vazio recebe a sequência mais distante do seu centróide
            mais_longe = np.argsort(d2[np.arange(len(X)), labels])[::-1]
            novos[vazios] = X[mais_longe[:vazios.sum()]]

        deslocamento = ((novos - centros) ** 2).sum()
        centros = novos
        if deslocamento <= tol:
            break

    d2 = distancias2(X, normas2, centros)
    labels = d2.argmin(axis=1)
    inercia = d2[np.arange(len(X)), labels].sum()
    return labels, centros, inercia, iteracao


def kmeans(X, k, n_init=10, max_iter=300, tol=1e-10, seed=None):
    """Retorna (labels, centros, inercia, iteracoes) da melhor das n_init execuções."""
    if not 1 <= k <= len(X):
        raise ValueError(f"k deve estar entre 1 e {len(X)} (quantidade de sequências)")
    rng = np.random.default_rng(seed)
    normas2 = np.einsum("ij,ij->i", X, X, dtype=np.float64)
    melhor = None
    for _ in range(n_init):
        centros = kmeans_pp(X, normas2, k, rng)
        resultado = _lloyd(X, normas2, centros, max_iter, tol)
        if melhor is None or resultado[2] < melhor[2]:
            melhor = resultado
    return melhor


def representantes(X, labels, centros):
    """Para cada cluster, a sequência real mais próxima do centróide (marcada com '*')."""
    d2 = distancias2(X, np.einsum("ij,ij->i", X, X, dtype=np.float64), centros)
    rep = np.full(len(centros), -1, dtype=np.int64)
    for c in range(len(centros)):
        membros = np.flatnonzero(labels == c)
        if membros.size:
            rep[c] = membros[d2[membros, c].argmin()]
    return rep
