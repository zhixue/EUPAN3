# EUPAN3

## 1. Introduction
EUPAN3 is a EUkaryote Pan-genome ANalysis toolkit for 3rd generation sequencing.


## 2. Install
### Dependency
* [Quast](https://github.com/ablab/quast) (>=v5.0.2)
* [Minimap2](https://github.com/lh3/minimap2) (>=v2.17)
* [Cd-hit](https://github.com/weizhongli/cdhit) (>=v4.8.1)
* [Gclust](https://github.com/niu-lab/gclust) (>=v1.0)
* [Blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (>=2.11.0)
### 

## 3. Usage
```
usage: eupan3.py [-h] [-v] {assemsta,unalnsseq} ...

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

sub commands:
  {assemsta,unalnsseq}
    assemsta            Use Quast to get assembled contigs/scaffolds to
                        reference to check the statistics (The script will
                        call QUAST program, so the directory where quast.py
                        locates is needed.)
    unalnsseq           Get partial/full/all unaligned sub sequences of
                        contigs/scaffolds from Quast
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
### 3.2 unalnsseq
```
usage: eupan3.py unalnsseq [-h] -a <assembly.fa> -u <unaligned.info> -o
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