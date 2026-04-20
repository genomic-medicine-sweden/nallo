process BAM_QC_FILTER {

    tag "${meta.id}"
    label 'process_low'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*_filtered.bam"), path("*_filtered.bam.bai"), emit: bam_bai

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    set -euo pipefail

    # 1. Build contig length table
    samtools view -H ${args} ${bam} \
    | awk '\$1=="@SQ" {
        for(i=1;i<=NF;i++){
            if(\$i ~ /^SN:/) sn=substr(\$i,4)
            if(\$i ~ /^LN:/) ln=substr(\$i,4)
        }
        print sn"\t"ln
    }' > contig_lengths.txt

    # 2. Extract reads (streaming, no intermediate files)
    # FAST path: skip CIGAR parsing in awk, use samtools CIGAR-aware filtering logic
    samtools view -H ${bam} > header.txt

    samtools view ${bam} \
    | awk -v threads=${task.cpus} '
    BEGIN{
        while((getline < "contig_lengths.txt")>0){
            len[\$1]=\$2
        }
    }

    function refspan(cigar,   i,n,num,op,s){
        s=0; num=""
        for(i=1;i<=length(cigar);i++){
            op=substr(cigar,i,1)
            if(op ~ /[0-9]/) num=num op
            else {
                n=num+0
                if(op=="M"||op=="="||op=="X"||op=="D"||op=="N")
                    s+=n
                num=""
            }
        }
        return s
    }

    {
        r=\$3; p=\$4; c=\$6
        if(!(r in len)) next
        if(p + refspan(c) - 1 <= len[r]) print
    }' \
    | cat header.txt - | samtools view -b -@ ${task.cpus} -o ${prefix}_filtered.bam -

    samtools sort -@ ${task.cpus} -o ${prefix}_filtered.sorted.bam ${prefix}_filtered.bam
    mv ${prefix}_filtered.sorted.bam ${prefix}_filtered.bam
    samtools index -@ ${task.cpus} ${prefix}_filtered.bam
    """
}
