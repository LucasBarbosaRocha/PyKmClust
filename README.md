# Clustering K-means (PyKmClust)
* Nome: Lucas Barbosa Rocha
* Disciplina: Inteligência Artificial
* Trabalho: Implementar um clustering para sequências de DNA utilizando a ideia do K-means.
* Contato: lucas.lb.rocha@gmail.com

Projeto irmão: **PyMeShClust** (mesmo problema, resolvido com MeanShift). Os dois foram feitos para serem comparados.

## A ideia
Baseado no artigo *Application of k-means clustering algorithm in grouping the DNA sequences of hepatitis B virus (HBV)*.

Agrupar sequências de DNA alinhando todas contra todas é caro: são O(n²) alinhamentos. A alternativa é comparar sequências **sem alinhar** (*alignment-free*), pela composição de k-mers:

1. **Sequência → vetor.** Cada sequência vira um vetor com a frequência de cada k-mer (palavra de tamanho k). Para k=3 são 64 posições: AAA, AAC, ..., TTT. O vetor é normalizado (soma 1), para que sequências de tamanhos diferentes sejam comparáveis.
   ```
   ACGTAC  ->  ACG, CGT, GTA, TAC  ->  [0, ..., 0.25 (ACG), ..., 0.25 (CGT), ...]
   ```
2. **K-means.** Você escolhe quantos grupos quer (k). O algoritmo:
   * escolhe k centróides iniciais espalhados (k-means++);
   * coloca cada sequência no cluster do centróide mais próximo (distância euclidiana);
   * recalcula cada centróide como a média dos vetores dos seus membros;
   * repete os dois passos anteriores até nenhuma sequência mudar de cluster;
   * roda 10 vezes com sementes diferentes e fica com o resultado de menor inércia (soma das distâncias ao centróide).
3. **Saída.** Um arquivo `.clstr` no formato do CD-HIT; a sequência mais próxima do centróide é marcada com `*`.

A diferença central para o PyMeShClust: aqui **você define o número de clusters**; lá, o número de clusters surge de um limiar de similaridade.

## Como usar

### Instalação
A única dependência é o `numpy`. Em Ubuntu/Debian o `pip` não instala no Python do sistema (erro `externally-managed-environment`), então use uma das opções:

```bash
# Opção 1: numpy do sistema
sudo apt install python3-numpy

# Opção 2: ambiente virtual (precisa do pacote python3-venv)
sudo apt install python3.12-venv
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
source .venv/bin/activate
```

### Rodar o clustering
```bash
python3 main.py hbv.fasta -k 2
python3 main.py sequencias.fasta -k 8 --kmer 4 -o saida.clstr
```

| Parâmetro | Padrão | Descrição |
|---|---|---|
| `entrada` | `hbv.fasta` | arquivo FASTA (sequência em uma ou várias linhas) |
| `-k` | 2 | quantidade de clusters |
| `--kmer` | 3 | tamanho do k-mer |
| `--n-init` | 10 | execuções com sementes diferentes |
| `--max-iter` | 300 | máximo de iterações por execução |
| `--seed` | 0 | semente aleatória |
| `-o` | `outputk<k>.clstr` | arquivo de saída |

### Avaliar um resultado
`avaliar.py` compara um `.clstr` com um rótulo verdadeiro extraído do cabeçalho das sequências:

```bash
python3 avaliar.py outputk2.clstr --rotulo regiao
```

| Rótulo | Exemplo de cabeçalho | Classe |
|---|---|---|
| `genotipo` | `HBV genotype C DNA, complete genome` | C |
| `regiao` | `HBV genotype F S gene ..., partial cds` | gene S parcial |
| `especie` | `Klebsiella pneumoniae strain SIKP041` | Klebsiella pneumoniae |
| `genero` | `Klebsiella pneumoniae strain SIKP041` | Klebsiella |

Sequências sem rótulo (os reads `ERR599000...` do `sequencias.fasta` e 4 sequências do HBV sem genótipo) são agrupadas normalmente, mas ficam fora da avaliação.

Métricas (1 = perfeito):
* **ARI** (Adjusted Rand Index): concordância entre os pares de sequências, corrigida pelo acaso. ~0 = agrupamento aleatório. É a métrica principal.
* **NMI**: informação mútua normalizada entre clusters e classes.
* **Pureza**: fração das sequências que estão no cluster dominado pela sua classe. Sobe sozinha quando há muitos clusters pequenos, então olhe sempre junto com o ARI.

ARI e NMI foram conferidos contra o scikit-learn.

### Reproduzir os experimentos
```bash
python3 experimentos.py   # gera RESULTADOS.md e resultados.csv
```

## Resultados
Tabelas completas em [RESULTADOS.md](RESULTADOS.md) e [resultados.csv](resultados.csv). Valores em ARI.

**Dados:**
* `hbv.fasta`: 16 sequências de HBV, com 7 genótipos e 4 regiões do genoma (genoma completo, gene S, gene P, P/preS/S).
* `sequencias.fasta`: 500 sequências, sendo 237 contigs de 13 espécies/8 gêneros de bactérias e 263 reads curtos (~100nt) sem rótulo.

### Efeito do tamanho do k-mer (k = nº de classes)
| k-mer | HBV região (k=4) | HBV genótipo (k=7) | gênero, todas (k=8) | gênero, só rotuladas (k=8) | matriz p/ 1 milhão de seqs |
|---|---|---|---|---|---|
| 3 | **0.801** | 0.127 | 0.205 | 0.213 | 256 MB |
| 4 | **0.801** | **0.300** | **0.225** | 0.322 | 1 GB |
| 5 | **0.801** | 0.124 | 0.131 | **0.329** | 4 GB |
| 6 | **0.801** | 0.124 | 0.000 | 0.245 | 16 GB |

* **k-mer 4 é o melhor compromisso.** Com k-mer 5 e 6 a memória multiplica por 4 a cada passo e a qualidade não melhora.
* **Os reads curtos atrapalham o k-means.** Com k-mer 6, um read de 100nt tem ~95 k-mers espalhados por 4096 posições; o vetor fica esparso e distante de tudo. Os reads viram outliers e "roubam" centróides. Com k=8, todas as sequências rotuladas caem num único cluster (ARI 0). Sem os reads, o ARI volta a 0.245.

### Comparação PyKmClust × PyMeShClust × versão original
Melhor configuração de cada método (entre parênteses: k-mer e parâmetro).

| Tarefa | PyKmClust (k-means) | PyMeShClust | Original 2020 |
|---|---|---|---|
| HBV região | 0.801 (3, k=4) | **0.966** (3, s=0.95, 6 clusters) | 0.867 (k-means k=2) |
| HBV genótipo | 0.300 (4, k=7) | **0.552** (6, s=0.8, 13 clusters) | 0.025 (k-means k=2) |
| Bactérias, espécie (todas) | **0.223** (4, k=8) | 0.198 (3, s=0.8, 247 clusters) | 0.033 (MeShClust s=0.95) |
| Bactérias, gênero (só rotuladas) | **0.329** (5, k=8) | 0.194 (3, s=0.8, 25 clusters) | - |

### O que os números mostram
* **A composição de k-mers separa bem a região do genoma, e mal o genótipo.** Genótipos do HBV diferem em ~8% dos nucleotídeos, o que quase não muda a frequência de k-mers curtos. Já um trecho do gene S e um genoma completo têm composições bem diferentes.
* **O original de 2020 acertava a região com k=2 por acaso.** A região no HBV quase coincide com o comprimento (genoma completo × parcial), e o histograma com bug agrupava justamente por comprimento. No genótipo o original fica em ~0 (0.025).
* **Nas bactérias, nenhum método passa de ~0.33.** Os contigs têm ~800nt e vêm de regiões diferentes do genoma. A "assinatura genômica" por k-mers costuma precisar de trechos de vários kb para separar espécies.
* **O k-means vai melhor nas bactérias, o PyMeShClust no HBV.** O k-means força exatamente k grupos; o PyMeShClust fragmenta as bactérias em muitos clusters pequenos (pureza alta, ARI baixo).
* **Cuidado com o HBV:** são só 16 sequências; uma sequência trocada de cluster muda bastante o ARI.

## Arquivos
* `main.py`: linha de comando do clustering.
* `kmer.py`: leitura de FASTA, vetores de k-mers e escrita do `.clstr`.
* `kmeans.py`: k-means (k-means++ e iterações de Lloyd).
* `avaliar.py`: rótulos, ARI, NMI e pureza.
* `experimentos.py`: gera `RESULTADOS.md` e `resultados.csv`.
* `original/`: versão original de 2020, mantida como referência.

## Memória
Todas as sequências ficam em uma matriz numpy `(N, 4^k)` float32, pré-alocada. Cada k-mer a mais multiplica a memória por 4:

| Sequências | k=3 | k=4 | k=5 | k=6 |
|---|---|---|---|---|
| 100 mil | 25 MB | 100 MB | 400 MB | 1.6 GB |
| 1 milhão | 256 MB | 1 GB | 4 GB | 16 GB |

Com 200 mil sequências de 800nt (arquivo de 167MB), k-mer 3 e k=20: ~265MB de pico de RSS e ~18s.

## Mudanças em relação à versão original (2020)
O problema de memória descrito na época vinha de:
* cada sequência guardava uma tabela `(74, 4)` quando só uma coluna era usada;
* o `khmer` criava tabelas, threads e um arquivo temporário para cada sequência;
* o centróide era recalculado somando o cluster inteiro a cada inserção (O(n²)).

Além disso:
* O histograma era a **distribuição de abundância** do khmer (quantos k-mers aparecem 1x, 2x, ...), não a contagem de cada k-mer. As linhas com zero eram puladas, desalinhando as posições entre sequências, e na prática os clusters separavam as sequências pelo comprimento. Agora é o vetor de frequência de k-mers normalizado.
* A escolha do cluster mais próximo devolvia o índice errado para k ≥ 3. Agora usa `argmin`.
* O centróide era a própria sequência semente e tinha o histograma sobrescrito. Agora é um vetor separado.
* Havia uma única passada. Agora repete até convergir.
* Os centróides iniciais eram escolhidos pelo comprimento (só funcionava para k=2). Agora usa k-means++.
* O comprimento contava o `\n` (1nt a mais) e o leitor só aceitava FASTA com a sequência em uma linha.
* A inércia final foi conferida contra o `KMeans` do scikit-learn.
