# Background

The original WASTER in ASTER software package requires more than 48 GB to filter out SNP sites, and the default k-size is 19 (k21 requires ~512GB memory)

# Aim

1. reduce memory usage to lower than 16GB
2. expand k-size to 21
3. enable mismatch at kmer flanks (to some extent)
4. parallelization without lock

# Implementation

## 1. initial `FilterTable`

Randomly select shared kmers for SNP calling while reducing memory cost to 8GB [`char [16][1<<29]`].

1. 16 samples. Transfer a k21 to 2-bit encoding, split left flank (n), right flank (m), central base (c);
2. sort n and m by lexicographical order to decide swap them or not;
3. hash n and m to get t, b, and r value;
4. store hash values: `FilterTable[t][b] = r` AND `t == fileorder`.
5. parallel: 16 threads, each assigned to an FilterTable[fileorder] array. For each file input, allocate 64MB memory for buffer.

\* Why `t==fileorder`: Assume each sample have a% shared kmers and (1-a%) individual kmers (which are noises for us). For each file, `t == fileorder` compresses a% and (1-a%) to 1/16 in the `FilterTable`. However, the a% can be fully recovered (since each file shares them), while (1-a%) is truly reduced to 1/16.

## 2. `EsTable`

Roughly estimate a kmer's frequency in a group (4 samples) to judge whether it is shared or not.

1. `blockEsTable`: project the first index of `EsTable` to fileorder;
2. read sample sequences, hash the kmers and look up t, b, and r in `FilterTable` to state. 
3. parallelization: also 16 threads

## 3. Compress: `kMerInfoTable`

1. A bucket every 64 b indices;
2. store kmers of the most frequency.

Remaining: 128M kmers -> max 128 Mbase SNPs 

## 4. Mismatch processing: `largeKMers`

1. Call out full-length k31 from sequences;
2. Judge whether a k21 is valid by: central base (c) heterogeneity; k31 flank mismatch. Keep a k21 if its corresponding k31s have less than 2 mismatch. Discard k21s if not, or if central base is polymorphic;
3. consensus k31 flanks.
4. parallelization: 16 threads.