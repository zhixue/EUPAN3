# EUPAN3

## 1. Introduction
EUPAN3 is a EUkaryote Pan-genome ANalysis toolkit for 3rd generation sequencing.


## 2. Install
### Dependency
* [Quast](https://github.com/ablab/quast) (>=v5.0.2)
* [Minimap2](https://github.com/lh3/minimap2) (>=v2.17)
* [Cd-hit](https://github.com/weizhongli/cdhit) (>=v4.8.1)
* [Gclust](https://github.com/niu-lab/gclust) (>=v1.0)
* [Blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST) (>=2.11.0)
* [Bedtools](https://github.com/arq5x/bedtools2) (>=2.29.2)
### 

## 3. Usage
```
usage: eupan3.py [-h] [-v]
                 {assemsta,unalnbseq,mergebseq,rmredundant,rmctm,fastasta,ptpg,genecov,elecov}
                 ...

EUPAN3: EUkaryote Pan-genome ANalysis toolkit for 3rd generation sequencing
Version 0.1.0

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

sub commands:
  {assemsta,unalnbseq,mergebseq,rmredundant,rmctm,fastasta,ptpg,genecov,elecov}
    assemsta            Use Quast to get assemblies statistics
    unalnbseq           Get partial/full unaligned block sequences
    mergebseq           Merge the block sequences from many files
    rmredundant         Cluster the sequences and remove redundant sequences
    rmctm               Detect and discard the contaminated sequences
    fastasta            Calculate statistics of fasta file (DNA)
    ptpg                Extract the longest transcript elements from genes
    genecov             Compute gene coverage and depth
    elecov              Compute genomic element coverage and depth
```

### 3.1 assemsta
```
usage: eupan3.py assemsta [-h] -a <assembly.fa> -r <reference.fa> -o
                            <output_directory> [-q <quast.py>] [-t <int>]
                            [-l <int>] [-g <reference.gff>] [-i <int>]
                            [--large_genome]

optional arguments:
  -h, --help            show this help message and exit
  -a <assembly.fa>, --assembly_path <assembly.fa>
                        Path of contigs or scaffolds
  -r <reference.fa>, --reference_path <reference.fa>
                        Path of reference genome
  -o <output_directory>, --output_dir <output_directory>
                        Path of output of Quast
  -q <quast.py>, --quast_path <quast.py>
                        Path of quast.py (default: quast.py in $PATH)

quast parameters:
  -t <int>, --thread <int>
                        number of threads using minimap2 in Quast (default: 1)
  -l <int>, --length_threshold <int>
                        Min length of contigs/scaffolds to consider (default:
                        500, if --large_genome, default: 3000)
  -g <reference.gff>, --gff_path <reference.gff>
                        Min length of sequences to consider (default: None)
  -i <int>, --identity <int>
                        Min alignment identity of sequences (choices:
                        80/90/95. default: 90)
  --large_genome        Use optimal parameters for evaluation of large genomes
                        (typically > 100 Mbp)
```
### 3.2 unalnbseq
```
usage: eupan3.py unalnbseq [-h] -a <assembly.fa> -u <unaligned.info> -o
                             <output.fa> [-l <int>] [-k <str>] [-s <str>]
                             [-rr <reference.fa>] [-ri <int>] [-rc <int>]
                             [-rm <minimap2_path>] [-rt <int>] [-rd <str>]

optional arguments:
  -h, --help            show this help message and exit
  -a <assembly.fa>, --assembly_path <assembly.fa>
                        Path of contigs or scaffolds
  -u <unaligned.info>, --unaln_path <unaligned.info>
                        Path of unaligned table (at quast_output/contigs_repor
                        ts/contigs_report_xxx.unaligned.info)
  -o <output.fa>, --output_path <output.fa>
                        Path of output unaln sequences
  -l <int>, --length_filter <int>
                        Min length of sub sequences to consider (default: 500)
  -k <str>, --kind_unaln <str>
                        Use full/partial/all unaligned sequences (choices:
                        full/partial/all. default: all)
  -s <str>, --sample_tag <str>
                        Add sample tag before each contig/scaffold (default:
                        None, e.g. -s Sample1 , "_" will be used, Chr1 ->
                        Sample1_Chr1)

realign parameters:
  -rr <reference.fa>, --realign_reference <reference.fa>
                        Using minimap2 to realign to reference
                        genome/mitochondrion/plastid and drop high similar
                        unaligned sub sequences
  -ri <int>, --realign_identity <int>
                        [Only use when -rr on] Min alignment identity of
                        sequences in realign step (choices: 80/90/95. default:
                        90)
  -rc <int>, --realign_coverage <int>
                        [Only use when -rr on] Min alignment coverage of
                        sequences (default: 80)
  -rm <minimap2_path>, --realign_minimap2 <minimap2_path>
                        [Only use when -rr on] Path of minimap2 (default:
                        minimap2 in $PATH)
  -rt <int>, --realign_thread <int>
                        [Only use when -rr on] Number of threads using
                        minimap2 (default: 1)
  -rd <str>, --realign_dir <str>
                        [Only use when -rr on] Temp directory to realign
                        (default: temp_dir)
```

### 3.3 mergebseq
```
usage: eupan3.py mergebseq [-h] [-o <output.fa>] [--allow_samectg]
                           input_fa [input_fa ...]

positional arguments:
  input_fa              Path of input fasta

optional arguments:
  -h, --help            show this help message and exit
  -o <output.fa>, --output <output.fa>
                        Path of output fasta
  --allow_samectg       Allow the sample contig/scaffold name
```

### 3.4 rmredundant
```
usage: eupan3.py rmredundant [-h] -i <input.fa> [-o <output_dir>] [-c <int>]
                             [-m <cluter_method_path>] [-t <int>]

optional arguments:
  -h, --help            show this help message and exit
  -i <input.fa>, --input_fa <input.fa>
                        Path of input fasta
  -o <output_dir>, --output_dir <output_dir>
                        Path of output directory
  -c <int>, --sequence_identity <int>
                        Sequence identity threshold (default: 90)
  -m <cluter_method_path>, --method_path <cluter_method_path>
                        Path of cluster method, support gclust/cd-hit-
                        est/blastn (default: gclust in $PATH)
  -t <int>, --thread <int>
                        Number of threads when clustering (default: 1)
```

### 3.5 rmctm
```
usage: eupan3.py rmctm [-h] -i <input.fa> [-o <output_dir>] [--rewrite] -nt
                       <nt_path> -at <accession2taxid_path> -rl
                       <rankedlineage_path> [-wl <str>] [-ai <int>]
                       [-al <int>] [--force_remove] [-m <blastn_path>]
                       [-e <float>] [-t <int>]

optional arguments:
  -h, --help            show this help message and exit
  -i <input.fa>, --input_fa <input.fa>
                        Path of input fasta
  -o <output_dir>, --output_dir <output_dir>
                        Path of output directory
  --rewrite             Rewrite the blastout in output directory

database parameters:
  -nt <nt_path>, --nt <nt_path>
                        Path of nt (download and decompress from
                        https://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz)
  -at <accession2taxid_path>, --at <accession2taxid_path>
                        Path of accession2taxid (download and decompress from
                        https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/
                        nucl_gb.accession2taxid.gz)
  -rl <rankedlineage_path>, --rl <rankedlineage_path>
                        Path of rankedlineage (download and decompress from ht
                        tps://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_ta
                        xdump.tar.gz)

filter parameters:
  -wl <str>, --white_list <str>
                        legal tax, others will be considered as contamination
                        (default: Viridiplantae)
  -ai <int>, --alignment_identity <int>
                        Alignment identity threshold (default: 90)
  -al <int>, --alignment_length <int>
                        Alignment length threshold (default: 100)
  --force_remove        Force remove all contamination even if same regions of
                        sequences map to white list

blastn parameters:
  -m <blastn_path>, --blastn_path <blastn_path>
                        Path of blastn (default: blastn in $PATH)
  -e <float>, --evalue <float>
                        E-value threshold of blastn (default: 1e-5)
  -t <int>, --thread <int>
                        Number of threads when blastn, recommend less than 4
                        (default: 1)
```

### 3.6 fastasta
```
usage: eupan3.py fastasta [-h] -i <input.fa> -o <output.fa> [-l <int>]

optional arguments:
  -h, --help            show this help message and exit
  -i <input.fa>, --input_path <input.fa>
                        Path of input fasta file
  -o <output.fa>, --output_path <output.fa>
                        Path of output fasta file
  -l <int>, --length_filter <int>
                        Min length of sequences to consider (default: 500)
```

### 3.7 ptpg
```
usage: eupan3.py ptpg [-h] -i <input.gff/gtf> [-r <str>] -o <output.gff/gtf>

optional arguments:
  -h, --help            show this help message and exit
  -i <input.gff/gtf>, --input <input.gff/gtf>
                        Path of input gff/gtf
  -r <str>, --region <str>
                        CDS or exon (default: CDS)
  -o <output.gff/gtf>, --output <output.gff/gtf>
                        Path of output gff/gtf
```

### 3.8 genecov
```
usage: eupan3.py genecov [-h] -a <input.gff/gtf> -b <input.bed> [-r <str>] -o
                         <output.cov> -n <str> [-m <int>]

optional arguments:
  -h, --help            show this help message and exit
  -a <input.gff/gtf>, --annotation <input.gff/gtf>
                        Path of input gff/gtf
  -b <input.bed>, --bed <input.bed>
                        bed of of coverage from bedtools
  -r <str>, --region <str>
                        CDS or exon (default: CDS)
  -o <output.cov>, --output <output.cov>
                        Path of output cov
  -n <str>, --sample_name <str>
                        Name of sample
  -m <int>, --min_depth <int>
                        Min depth (default: 1)
```

### 3.9 elecov
```
usage: eupan3.py elecov [-h] -a <anotation.bed> -b <input.bed> -o <output.cov> -n
                 <str> [-m <int>]

optional arguments:
  -h, --help            show this help message and exit
  -a <anotation.bed>, --annotation <anotation.bed>
                        Path of anotation.bed
  -b <input.bed>, --bed <input.bed>
                        bed of of coverage from bedtools
  -o <output.cov>, --output <output.cov>
                        Path of output cov
  -n <str>, --sample_name <str>
                        Name of sample
  -m <int>, --min_depth <int>
                        Min depth (default: 1)
```