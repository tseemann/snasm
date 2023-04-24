[![License: GPLv3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-3.0.en.html)
![Don't judge me](https://img.shields.io/badge/Language-Perl_5-steelblue.svg)

# snasm

SNP-based alignments from assemblies

## Introduction

When you have FAATQ, you use SNippy.
If you also have aseemblies, you hack
then into Snippy with the `--ctgs` 
option. If you only have assemblies,
you can now use `snasm`.

## Installation

Because `snasm` is undergoing development
it is not in Bioconda yet, so installation
is still somewhat manual:
```
git clone https://github.com/tseemann/snasm.git

conda create -n snasm \
  any2fasta minimap2 bedtools snp-sites\
  'bcftools>=1.17' 'samtools>=1.17' \
  seqkit csvtk goalign gotree iqtree \
  perl-bioperl-core perl-file-slurp \
  perl-file-which perl-text-csv

conda activate snasm

export PATH=$PATH:$PWD/snasm/bin

snasm -h
```

## Quick Start

```
# Align all FASTA files, use first one as reference
% snasm -d outdir asms/*.fasta

# Provide a list of filenames in a file
% snasm -d outdir -i genomes.fofn

# Choose a specific reference
% snasm -d outdir -r ref.gbk asms/*.fa

# Don't use all the CPU cores
% snasm -j 8 ...

# We handles lots of file types
% snasm in.gbk in.fa in.fa.gz in.gb.bz2 in.embl.xz in.gff3

```
 
## Input

We accept any file that `any2fasta` 
can extract sequence data from.
This includes Genbank, FASTA, GFA,
EMBL or GFF3.

To obtain mutation consequences, the
reference genome needs to be in 
Gebbank format.

Any file can be optionally compressed
with GZIP, BZIP2, LZMA, or XZ.

## Output

Each genome will have a folder named 
after it, contain the sequence (fna),
variants (vcf, vcf.gz), the alignments
(paf, bed) and the unaligned regions (bed).

The overall output folder will have
summary output files, including a
multisampel VCF, a reference-based
multiple sequence alignment (afa, aln),
a coverage report (cov), and a list
of sample names (ids).


## Etymology

The name `snasm` is a contraction of the
words "SNP" and "Assembly".

## Citation

Seemann T, SNP based alignments from assemblies  **Github** https://github.com/tseemann/snasm

## Feedback

Please submit via the [Github Issues page](https://github.com/tseemann/snasm/issues)

## Licence

[GPL v3](https://raw.githubusercontent.com/tseemann/snasm/master/LICENSE)

## Author

* Torsten Seemann
* Web: https://tseemann.github.io/
* Twitter: [@torstenseemann](https://twitter.com/torstenseemann)
* Blog: [The Genome Factory](https://thegenomefactory.blogspot.com/)

