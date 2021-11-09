## Experiments from the manuscript

### Prerequisites

* python3

### Experiments

(more to come)

### Sequences

In order to track sequences generated under different models, sequence pairs
from the same model were, in most cases, assigned names with a specific fruit.
They have no connection in reality to the given fruit (e.g. they are *not*
reads sequenced from that fruit). 

Sequence pairs generated from a randomly chosen segment of E.coli, under a
mutation model. There is one pair for each mutation rate (.1%, .5%, 1%, 5%,
and 10%).

```bash 
ecoli.K12.dupfree.K=16.L=10K.K=16.mutation_model.first_pair.fa.gz
```

Unrelated sequence pairs generated as random duplicate-free walks. There is
one pair for each target Jaccard (10%, 20%, ... , 90%).

```bash 
plantain.K=7.first_pair.fa.gz
plantain.K=8.first_pair.fa.gz
```

A randomly chosen 1Kbp segment from human (hg38) chromosome 20.

```bash 
hg38.chr20.L=1K.fa.gz
```

(more to come)
