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

#### ecoli.K12.dupfree.K=16.L=10K.K=16.mutation_model.first_pair.fa.gz

Related sequence pairs generated from a randomly chosen segment of E.coli,
under a mutation model. There is one pair for each mutation rate (.1%, .5%, 1%,
5%, and 10%). The first sequence in each pair is the same.

#### tangerine.ecoli.K12.L=10K.R=\*.mutation_model.fa.gz

Sequences (*not* pairs) drawn from randomly chosen segments of E.coli and
subjected to a Bernoulli mutation model. For each mutation rate (1% 5% 10%)
there are one hundred sequences. The same one hundred E. coli segments were
selected for each mutation rate.

#### lemon.K=16.R=10%.fa.gz 

Related sequence pairs generated from a mutation model applied to a random
duplicate-free sequence. For each length (100, 1000, 2000, ... , 9000, 10000)
there are fifty pairs (each with no relation to the other pairs) for mutation
rate 10%.

#### lemonB.K=16.R=10%.fa.gz 

Related sequence pairs generated from a mutation model applied to a random
duplicate-free sequence. There are fifty pairs (each with no relation to the
other pairs) for mutation rate 10%.

#### plantain.K=7.first_pair.fa.gz and plantain.K=8.first_pair.fa.gz

Unrelated sequence pairs generated as duplicate-free random walks. There is
one pair for each target Jaccard (10%, 20%, ... , 90%).

#### hg38.chr20.L=1K.fa.gz

A randomly chosen 1Kbp segment from human (hg38) chromosome 20.

(more to come)
