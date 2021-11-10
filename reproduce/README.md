## Experiments from the manuscript

### Prerequisites

* python3
* gzip
* mashmap (from https://github.com/marbl/MashMap)
* modified mashmap (from https://github.com/medvedevgroup/MashMap), named as
mashmapW1

### Experiments

Most of the experiment involve running jaccard_correction_test.py on a
collection of sequence pairs. The output of jaccard_correction_test is a
tab-delimited table with a row for each pair or each group of pairs, as shown
below. Labeled columns provide details and results beyond what is reflected in
the manuscript's figures and tables. The sequence files are described later in
this document. 

```bash 
#nameA                     nameA                      replicates w   k  length.nt |a|   I(A,B) U(A,B) J(A,B)   I(A,B;w) U(A,B;w) J(A,B;w) D(A,B;w) Jd(A,B;w) C(A,B;w)  Bias(A,B;w) J(A,B;w)-J(A,B) I(A,B;w)-C(A,B;w)
LEMON_L10015_K16_R10%_1_0  LEMON_L10015_K16_R10%_1_1  50         100 16 10015     10000 1714   18286  0.093733 23.160   373.420  0.062117 1687     0.092120  22.612257 -0.031567   -0.031616       0.547743
LEMON_L10015_K16_R10%_2_0  LEMON_L10015_K16_R10%_2_1  50         100 16 10015     10000 1805   18195  0.099203 23.160   371.780  0.062371 1802     0.099022  24.096816 -0.034236   -0.036832       -0.936816
LEMON_L10015_K16_R10%_3_0  LEMON_L10015_K16_R10%_3_1  50         100 16 10015     10000 1945   18055  0.107726 25.620   368.860  0.069591 1919     0.106134  26.194298 -0.035308   -0.038135       -0.574298
LEMON_L10015_K16_R10%_4_0  LEMON_L10015_K16_R10%_4_1  50         100 16 10015     10000 1921   18079  0.106256 23.960   369.240  0.065008 1901     0.105033  25.798101 -0.035354   -0.041248       -1.838101
LEMON_L10015_K16_R10%_5_0  LEMON_L10015_K16_R10%_5_1  50         100 16 10015     10000 1936   18064  0.107174 26.780   370.340  0.072444 1903     0.105156  25.606760 -0.036029   -0.034731       1.173240
```

#### Figure 3

For the unrelated pairs:

```bash 
gzip -dc plantain.K=7.first_pair.fa.gz \
  | jaccard_correction_test.py K=7 W=20 --replicates=50 --prng=20210908A \
  > plantain.K=7.W=20.first_pair.dat

gzip -dc plantain.K=8.first_pair.fa.gz \
  | jaccard_correction_test.py K=8 W=20 --replicates=50 --prng=20210908A \
  > plantain.K=8.W=20.first_pair.dat
```

For the related pairs:

```bash 
gzip -dc ecoli.K12.dupfree.K=16.L=10K.K=16.mutation_model.first_pair.fa.gz \
  | jaccard_correction_test.py K=16 W=20 --replicates=50 --prng=20210908A \
  > ecoli.K12.dupfree.K=16.L=10K.K=16.mutation_model.first_pair.W=20.dat

gzip -dc ecoli.K12.dupfree.K=16.L=10K.K=16.mutation_model.first_pair.fa.gz \
  | jaccard_correction_test.py K=16 W=200 --replicates=50 --prng=20210908A \
  > ecoli.K12.dupfree.K=16.L=10K.K=16.mutation_model.first_pair.W=200.dat
```

#### Figure 4

```bash 
cat hg38.chr20.fa \
  | sliding_jaccard.py --minimizers:local data/hg38.chr20.L=1K.fa.gz K=16 W=200 \
  > hg38.chr20.L=1K.K=16.W=200.local.dat.gz
```

#### Figure 5, Figure 6B

```bash 
declare -a windowSizes=("20" "100" "200" "300" "400" "500" "600" "700" "800" "900" "1000")
for W in "${windowSizes[@]}" ; do
  gzip -dc lemonB.K=16.R=10%.fa.gz \
    | jaccard_correction_test.py K=16 W=${W} \
      --replicates=50 --prng=20210908A \
    > lemonB.K=16.R=10%.W=${W}.dat
  done
```

#### Figure 6A

```bash 
gzip -dc lemon.K=16.R=10%.fa.gz \
  | jaccard_correction_test.py K=16 W=20 --replicates=50 --prng=20210908A \
  > lemon.K=16.R=10%.W=20.dat
```

#### Table 1

```bash 
declare -a errorRates=("1" "5" "10")
for R in "${errorRates[@]}" ; do
  # modified mashmap, with W=1
  mashmapW1 \
    -r ${ecoliDir}/ecoli.K12.fa.gz \
    -q tangerine.ecoli.K12.L=${L}.R=${R}%.mutation_model.fa.gz \
    -o temp/tangerine.ecoli.K12.L=${L}.W=1.R=${R}%.mutation_model.mashmap \
    2> tangerine.ecoli.K12.L=${L}.W=1.R=${R}%.mutation_model.mashmap.log \
    > tangerine.ecoli.K12.L=${L}.W=1.R=${R}%.mutation_model.mashmap.out
  #
  cat temp/tangerine.ecoli.K12.L=${L}.W=1.R=${R}%.mutation_model.mashmap \
    | awk 'BEGIN { print "#query q.len q.start q.end strand reference r.len r.start r.end identity.mash identity.truth" }
                 { print $0,identity; }' identity=$((100-R)) \
    > tangerine.ecoli.K12.L=${L}.W=1.R=${R}%.mutation_model.mashmap.dat
  #
  # regular mashmap, with default W=200
  mashmap \
    -r ${ecoliDir}/ecoli.K12.fa.gz \
    -q tangerine.ecoli.K12.L=${L}.R=${R}%.mutation_model.fa.gz \
    -o temp/tangerine.ecoli.K12.L=${L}.W=200.R=${R}%.mutation_model.mashmap \
    2> tangerine.ecoli.K12.L=${L}.W=200.R=${R}%.mutation_model.mashmap.log \
    > tangerine.ecoli.K12.L=${L}.W=200.R=${R}%.mutation_model.mashmap.out
  #
  cat temp/tangerine.ecoli.K12.L=${L}.W=200.R=${R}%.mutation_model.mashmap \
    | awk 'BEGIN { print "#query q.len q.start q.end strand reference r.len r.start r.end identity.mash identity.truth" }
                 { print $0,identity; }' identity=$((100-R)) \
    > tangerine.ecoli.K12.L=${L}.W=200.R=${R}%.mutation_model.mashmap.dat
  done
```

#### Table S2

```bash 
:> hash_specs
echo "minimap2_noncanon    minimap2    --noncanonical --noinhibit:correction" >> hash_specs
echo "murmurhash3_noncanon murmurhash3 --noncanonical --inhibit:correction"   >> hash_specs
echo "splitmix64_noncanon  splitmix64  --noncanonical --inhibit:correction"   >> hash_specs

cat hash_specs \
  | while read hashName hashFunc hashCanon inhibition ; do
      gzip -dc data/plantain.K=7.first_pair.fa.gz \
        | jaccard_correction_test.py K=7 W=20 \
          --hash=${hashFunc}.0 ${hashCanon} --replicates=50 --prng=20210908A \
          ${inhibition} \
        | awk '/^#/ { $3="hash "$3; print $0; }
              !/^#/ { $3=hashName" "$3; print $0; }' hashName=${hashName} \
        > results/plantain.K=7.W=20.first_pair.hash=${hashName}.dat
      #
      gzip -dc data/plantain.K=8.first_pair.fa.gz \
        | jaccard_correction_test.py K=8 W=20 \
          --hash=${hashFunc}.0 ${hashCanon} --replicates=50 --prng=20210908A \
          ${inhibition} \
        | awk '/^#/ { $3="hash "$3; print $0; }
              !/^#/ { $3=hashName" "$3; print $0; }' hashName=${hashName} \
        > results/plantain.K=8.W=20.first_pair.hash=${hashName}.dat
      done
```

#### Table S3

```bash 
:> hash_specs
echo "minimap2_noncanon    minimap2    --noncanonical --noinhibit:correction" >> hash_specs
echo "murmurhash3_noncanon murmurhash3 --noncanonical --inhibit:correction"   >> hash_specs
echo "splitmix64_noncanon  splitmix64  --noncanonical --inhibit:correction"   >> hash_specs

cat hash_specs \
  | while read hashName hashFunc hashCanon inhibition ; do
      echo "=== ${hashName} ==="
      gzip -dc data/ecoli.K12.dupfree.K=16.L=10K.K=16.mutation_model.first_pair.fa.gz \
        | jaccard_correction_test.py K=16 W=20 \
          --hash=${hashFunc}.0 ${hashCanon} --replicates=50 --prng=20210908A \
          ${inhibition} \
          --progress=1 \
        | awk '/^#/ { $3="hash "$3; print $0; }
              !/^#/ { $3=hashName" "$3; print $0; }' hashName=${hashName} \
        | line_up_columns \
        > results/ecoli.K12.dupfree.K=16.L=10K.K=16.mutation_model.first_pair.W=20.hash=${hashName}.dat
      done
```

### Sequences

In order to track random sequences generated under different models, sequence
pairs from the same model were, in most cases, assigned names with a specific
fruit. They have no connection in reality to the given fruit (e.g. they are
*not* reads sequenced from that fruit).

Where sequences were drawn from E.coli, these were taken from Escherichia coli
strain K-12 substrain MG1655, https://www.ncbi.nlm.nih.gov/nuccore/U00096.

The Genome Reference Consortium's Human Build 38, also known as GRCh38 or hg38,
was also used. Specifically chromsome 20, which is available at
https://www.ncbi.nlm.nih.gov/nuccore/CM000682.2. In the experiments above, the
filename is hg38.chr20.fa.

Sequence pairs are stored in tandem in the following files. That is, the first
two sequences in the file are one pair, the next two sequences are the second
pair, and so on.

#### ecoli.K12.dupfree.K=16.L=10K.K=16.mutation_model.first_pair.fa.gz

Related sequence pairs generated from a randomly chosen segment of E.coli,
under a Poisson mutation model. There is one pair for each mutation rate (.1%,
.5%, 1%, 5%, and 10%). The first sequence in each pair is the same.

#### tangerine.ecoli.K12.L=10K.R=\*.mutation_model.fa.gz

Sequences (*not* pairs) drawn from randomly chosen segments of E.coli and
subjected to a Bernoulli mutation model. For each mutation rate (1%, 5%, and
10%) there are one hundred sequences. The same one hundred E. coli segments
were selected for each mutation rate.

#### lemon.K=16.R=10%.fa.gz

Related sequence pairs generated from a Poisson mutation model applied to a
random duplicate-free sequence. For each length (100, 1000, 2000, ... , 9000,
10000) there are fifty pairs (each with no relation to the other pairs) for
mutation rate 10%.

#### lemonB.K=16.R=10%.fa.gz

Related sequence pairs generated from a Poisson mutation model applied to a
random duplicate-free sequence. There are fifty pairs (each with no relation to
the other pairs) for mutation rate 10%.

#### plantain.K=\*.first_pair.fa.gz

Unrelated sequence pairs generated as independent duplicate-free random walks.
There is one pair for each target Jaccard (10%, 20%, ... , 90%).

#### hg38.chr20.L=1K.fa.gz

A randomly chosen 1Kbp segment from human (hg38) chromosome 20.

