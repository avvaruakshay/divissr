# DiviSSR

[![Build](https://img.shields.io/badge/Build-passing-brightgreen.svg)]()
[![License](https://img.shields.io/badge/Licence-MIT-blue.svg)]()

DiviSSR is a DNA tandem repeat identification tool. Tandem repeats are important 
genomic sequences which have functional and evolutionary significance.

DiviSSR is scripted in C++. 

### Working principle

DiviSSR captures a mathematical property exhibited by binary number representations of DNA tandem repeats. In the binary number or the 2-bit format of DNA sequences, nucleotides are represented by 4 different combinations of 2 binary digits (bits). The binary number representations of DNA tandem repeats conform to a unique division rule. DiviSSR scans the genome in windows, converts each window sequence to a binary number and checks if the number qualifies the division rule. The detailed explanation of the algorithm is below and an overview of the method is depicted in the figure below.

<img src="./divissr/lib/4097_figure1_v2.jpg" alt=""/>
A) 2-bit representation for each nucleotide. B) Example repeat sequence of ACT motif with 2 complete units and a partial suffix of 2 bases. C) Binary representation of the repeat sequence generates a number denoted by S. A divisor (D) is built based on the length of the repeat, the motif and the partial suffix. The division of S with D yields a quotient which is equal to binary representation of the motif and the remainder equals binary representation of partial suffix.
<br><br>


### Installation

```bash
$ pip install divissr
```

### Usage

The help message and available options can be accessed using
```bash
$ divissr -h 
```
which gives the following output
```
usage: divissr [-h] -i <FILE> [-o <FILE>] [-m <INT>] [-M <INT>] [-l <INT>] [-a] [-g <FILE>]
                 [--anno-format <STR>] [--gene-key <STR>] [--up-promoter <INT>]
                 [--down-promoter <INT>]

Required arguments:
  -i, --input	<FILE>	Input sequence file. Fasta format.

Optional arguments:
  -o, --output			<FILE>	Output file name. Default: Input file name + _divissr.tsv
  -m, --min-motif-size	<INT>	Minimum size of a repeat motif in bp. Default: 1
  -M, --max-motif-size	<INT>	Maximum size of a repeat motif in bp. Default: 6
  -l, --min-length		<INT>	Cutoff repeat length. Default: 2*M.
 								Should at least be twice of maximum motif size.
  -a, --analyse					Generate a summary HTML report.

Compound repeat arguments:
  --compound            Report compound repeats. The output of compound repeats is a separate file with the suffix ".compound".
  -d <INT>, --comp-dist <INT>
                        Maximum distance between individual repeats of compound repeat. Use negative to denote overlap. Default: 0

Annotation arguments:
  -g, --annotate		<FILE>	Genomic feature file to annotate repeats w.r.t genes.
  								Both GFF and GTF can be processed.
  --anno-format			<STR>   Format of genomic feature file. 
  								Valid inputs: GFF, GTF. Default: GFF
  --gene-key			<STR>   Attribute used as unique identifier for gene name.
  								The default identifier is "gene". 
  --up-promoter			<INT>   Upstream distance(bp) from TSS to be considered as
  								promoter region. Default: 1000
  --down-promoter		<INT>   Downstream distance(bp) from TSS to be considered as
  								promoter region. Default: 1000
```

### `-i or --input`
**Expects:** *STRING (to be used as filename)*<br>
**Default:** *None*<br>
This is the only required argument for the program. The input file must be a 
valid FASTA file. 

### `-o or --output`
**Expects:** *STRING (to be used as filename)*<br>
**Default:** *Input Filename + _divissr.tsv (see below)*<br>
If this option is not provided, the default output filename will be the same as the input filename, with its extension replaced with '_divissr.tsv'. For example, if the input filename is `my_seq.fa`, the default output filename will be `my_seq.fa_divissr.tsv`. If the input filename does not have any extension, `_divissr.tsv` will be appended to the filename. Please note that even in the case of no identified SSRs, the output file is still created (therefore overwriting any previous file of the same name) but with no content in the file.
#### Output for fasta
The output is a tab-delimited file, with one SSR record per line. 
The output columns follow the [BED](https://genome.ucsc.edu/FAQ/FAQformat.html) format. The details of the columns are given below:

| S.No | Column | Description |
|:----:| ------ | ----------- |
| 1 | Chromosome | Chromosome or Sequence Name as specified by the first word in the FASTA header |
| 2 | Repeat Start | 0-based start position of SSR in the Chromosome |
| 3 | Repeat Stop | End position of SSR in the Chromosome |
| 4 | Repeat Class | Class of repeat as grouped by their cyclical variations |
| 5 | Repeat Length | Total length of identified repeat in nt |
| 6 | Repeat Strand | Strand of SSR based on their cyclical variation |
| 7 | Motif Number | Number of times the base motif is repeated |
| 8 | Actual Repeat | Starting sequence of the SSR irrespective of Repeat class and strand|

An example output showing some of the largest repeats from *Drosophila melanogaster* is given below
```
X       22012826  22014795  ACTGGG  1969    -       328     TCCCAG
2RHet   591337    591966    AATACT  629     -       104     ATTAGT
4       1042143   1042690   AAATAT  547     +       91      AAATAT
2RHet   598244    598789    AATACT  545     -       90      AGTATT
XHet    122       663       AGAT    541     +       135     GATA
X       22422335  22422827  AGAT    492     +       123     GATA
3R      975265    975710    AAAT    445     -       111     TTAT
X       15442288  15442724  ACAGAT  436     +       72      ACAGAT
2L      22086818  22087152  AATACT  334     -       55      TATTAG
YHet    137144    137466    AAGAC   322     -       64      CTTGT
```

### `-m or --min-motif-size`
**Expects:** *INTEGER*<br>
**Default:** *1*<br>
Minimum length of motifs to be considered. By default, divissr ignores redundant 
motifs. For example, a stretch of 12 A's is considered a monomer repeat of 12 
A's rather than a dimer repeat of 6 AA's. 

### `-M or --max-motif-size`
**Expects:** *INTEGER*<br>
**Default:** *6*<br>
Maximum length of motifs to be considered. Setting a large value of `-M` has a 
non-trivial effect on both the runtime and memory usage of divissr.

### `-l or --min-length`
**Expects:** *INTEGER*<br>
**Default:** *2* * *M*<br>
Minimum length cut-off to be considered when finding an SSR. The same cut-off 
will apply for SSRs of all motif lengths, even if the motif length is not a 
divisor of this value. In such cases, SSRs that end with a partial motif are 
also picked if they pass the length cut-off. This value should be at least twice
of the maximum motif size.

### `-a or --analyze`
**Expects:** *None*<br>
**Default:** *False*<br>
In addition to the default tab-separated output, DiviSSR can also generate a fully
interactive HTML report for easy downstream analysis of the repeat data. The 
filename will be the same prefix as that of the main output. For example, if the
input filename was my_seq.fa, the analysis report will be my_seq_divissr.html. An 
example HTML report, generated from the repeat data of Homo sapiens (build hg19),
can be accessed here (Right click -> Save As).

### `--compound`
**Expects:** *None*<br>
**Default:** *False*<br>
This is flag which when set to true reports all compound repeats. Compound repeats
are repeats which are either overlapping or separated by a small gap. Compound repeats
are reported in a separate file with the extension ".compound". The maximum distance
between two individual repeats of a compound repeat can be set using the option `-d`.
The compoud repeat has four columns with first four denoting the sequence, start, 
end and length of the repeat. The fifth column denotes the repeat classes of individual
repeats of the compound repeat and the number of different cyclical variations that
have occured as compound repeat. The sixth column is the strand of the individual
repeats. The seventh column is denotes the actual motifs of individual repeats, 
repeat length and the distance between the individual repeats. Below is an example
output reporting compound repeats.

```
chr1	94808	94847	(AAAAG)1(A)1	39	-|-	(CTTTT)20|D0|(T)23
chr1	113871	113896	(AGAGGG)1(AAAGAG)1	25	-|-	(CCCTCT)13|D0|(TCTCTT)12
chr1	130955	130977	(AAAAAC)1(AAAAC)1	22	+|+	(CAAAAA)12|D0|(AAAAC)13
chr1	147433	147475	(AAC)1(AAAAAC)1(AAC)1	42	+|+|+	(CAA)21|D0|(AACAAA)17|D-2|(AAC)14
chr1	156944	156967	(AAAAAT)2	23	+|+	(TAAAAA)12|D0|(AAAAAT)14
chr1	174902	174957	(A)1(AG)1(AAAG)1	55	+|+|+	(A)31|D0|(AG)13|D-2|(AGAA)15
chr1	180444	180472	(ACCCTG)1(AACCCT)1	28	+|+	(CCCTGA)16|D0|(TAACCC)12
chr1	180504	180537	(ACCCTC)1(AACCCT)1	33	+|+	(CCTCAC)20|D0|(CCCTAA)15
chr1	180552	180590	(AACCCT)2(AACCC)1	38	+|+|+	(AACCCT)14|D0|(AACCCT)14|D-2|(AACCC)12
chr1	180748	180912	(AACCCT)4	164	+|+|+|+	(CCCTAA)55|D0|(ACCCTA)44|D-2|(AACCCT)21|D-7|(AACCCT)44
```

### `-d or --comp-dist`
**Expects:** *INTEGER*<br>
**Default:** *0*<br>
Maximum distance between two individual repeats of a compound repeat. To report repeats overlapping by X bp use negative numbers i.e., -X. The default is 0bp. 

### `-g or --annotate`
**Expects:** *FILE*<br>
**Default:** *None*<br>
Input a genomic feature file to annotate the repeats in the genomic context. 
DiviSSR accepts both GFF and GTF format genomic feature files. Each repeat is 
annotated w.r.t the closest gene and classified either as Genic, Exonic, 
Intronic and Intergenic according to the position of the repeat. Besides this, 
the repeat is also checked if it falls in the promoter region of the gene. 
Annotation adds 7 columns to the default divissr output which already consist 8 
columns.

| S.No | Column | Description |
|:----:| ------ | ----------- |
| 9 | Gene name | Name of the closest gene |
| 10 | Gene Start | Start position of gene in the Chromosome |
| 11 | Gene Stop | End position of gene in the Chromosome |
| 12 | Strand | The strand orientation of the gene |
| 13 | Genomic annotation | Annotation of the repeat w.r.t to the gene. Possible annotations are {Genic, Exonic, Intronic, Intergenic} |
| 14 | Promoter annotation | If repeat falls in the promoter region of the closest gene. The default promoter region is 1Kb upstream and downstream of TSS. |
| 15 | Distance from TSS | Distance of the repeat from the TSS of the gene. |

### `--anno-format`
**Expects:** *STRING*<br>
**Default:** *GFF*<br>
Option to specify the format of the input genomic feature file. Accepted inputs 
are GFF or GTF. More details about the GFF and GTF formats can be found 
[here](https://asia.ensembl.org/info/website/upload/gff.html).

### `--gene-key`
**Expects:** *STRING*<br>
**Default:** *gene*<br>
The attribute key used for the name of the gene in the GFF/GTF file. In the 
below example GFF file, we have the location of a gene and it's mRNA and exon 
locations. The last column of the file specifies attributes associated with each
feature, like ID, Parent, gene etc. DiviSSR uses on of the attribute to identify 
the gene and also it's exons. In th below example the key "gene" can be used to 
identify gene and the exons of the gene as they have the same gene name. Please 
check your GFF/GTF file for a robust attribute key which can identify all genes 
and their corresponding exons. We are actively working on better annotation 
where we can identify genes and their exons based on the ID and Parent.

```
# Sample GFF
NC_004354.4	RefSeq	gene	124370	126714	.	-	.	ID=gene1;Name=CG17636;gbkey=Gene;gene=CG17636;gene_biotype=protein_coding;gene_synonym=DmelCG17636,EG:23E12.1;
NC_004354.4	RefSeq	mRNA	124370	126714	.	-	.	ID=rna1;Parent=gene1;Name=NM_001103384.3;gbkey=mRNA;gene=CG17636;transcript_id=NM_001103384.3
NC_004354.4	RefSeq	exon	126626	126714	.	-	.	ID=id13;Parent=rna1;gbkey=mRNA;gene=CG17636;transcript_id=NM_001103384.3
NC_004354.4	RefSeq	exon	125495	126259	.	-	.	ID=id14;Parent=rna1;gbkey=mRNA;gene=CG17636;transcript_id=NM_001103384.3
```

### `--up-promoter`
**Expects:** *INT*<br>
**Default:** *1000*<br>
Upstream distance(bp) from the TSS of the gene to be considered as promoter 
region. Default 1000.

### `--down-promoter`
**Expects:** *INT*<br>
**Default:** *1000*<br>
Downstream distance(bp) from the TSS of the gene to be considered as promoter 
region. Default 1000.

## Change log
### v 0.1.1

### v 0.1.2
- Solved compound repeat error for gzipped fasta input

1. Handling error of non ACGTN nucleotides.

## Contact
For queries or suggestions, please contact:

Akshay Kumar Avvaru - <avvaru@ccmb.res.in><br>
Divya Tej Sowpati - <tej@ccmb.res.in>