# collapse_rescue
Method for rescuing collapsed contigs.

## Installation
- install `collapse_rescue`
    ```bash
    git clone https://github.com/wangyibin/collapse_rescue.git
    cd collapse_rescue
    chmod +x bin/* 
    pip install -r requirements.txt
    ```
- configure  
`vim ~/.bash_profile`
    ```bash
    export PATH=$HOME/software/collapse_rescue/bin:$PATH
    export PATH=$HOME/software/ALLHiC_adjuster/bin:$PATH
    export PATH=$HOME/software/ALLHiC/bin:$PATH
    export PATH=$HOME/software/popCNV/bin:$PATH
    ```
## Dependencies
1. python packages
    - rich
    - joblib
    - numpy
    - pandas
    - pyfaidx
    - pytools

2. softwares
    - [ALLHiC](https://github.com/tangerzhang/ALLHiC)
    - [ALLHiC_adjuster](https://github.com/sc-zhang/ALLHiC_adjuster)
    - [popCNV](https://github.com/sc-zhang/popCNV)

## Methods
1. calculate the copy number of each contig using `popCNV`
    ```bash
    mkdir read_depth && cd read_depth
    mosdepth -t 10 -b 1000 LAp LAp.contig.bam
    cd ..
    popCNV -g LAp.contig.fasta -s 1000 -r read_depth/ -b bamfile/ -l contig.bed -w wrk_dir --group group.list  --sample sample.list --wild 0 -t 10
    ```
    Alignment must use contig-level assembly as reference.
2. rescue
    ```bash
    agp2cluster.py LAp.agp > LAp.clusters.txt
    dup_collapsed_contigs.py Allele.ctg.table 06.genes.round.cn LAp.clusters.txt 8 -o output
    ```
    `Allele.ctg.table` generated from `ALLHiC`  
    `06.genes.round.cn` generated from `popCNV`  
3. optimize
    ```bash
    convert_agp_to_tour.py LAp.agp tour
    dup_collapsed_optimize.py collapsed.rescued.txt LAp.pairs.txt LAp.clm tour -o new_tour > collapsed.optimize.txt
    ```
    `LAp.pairs.txt` and `LAp.clm` were generated from `allhic extract`
4. build
    ```bash
    dup_collapsed_fasta.py collapsed.optimize.txt LAp.contig.fasta > contig.dup.fasta
    cd new_tour
    ALLHiC_build contig.dup.fasta
    ```
    Results are the `groups.agp` and `groups.asm.fasta`.   
    Sequences ID with `_d` are the collapsed contigs.
## Citing 