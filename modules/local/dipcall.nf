process DIPCALL {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::dipcall=0.3"
    container "quay.io/biocontainers/dipcall:0.3--hdfd78af_0"

    // This is bad but I don't know where it's going wrong, a test dataset is really needed for this process
    // Could change to a subworkflow as well...
    shell '/bin/bash', '-uo pipefail'

    // "To make proper calls on sex chromosomes, users should hard mask PARs on the reference chrY."

    input:
    tuple val(meta), path(haplotype_1), path(haplotype_2), val(sex)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(mmi)
    path(par)

    output:
    tuple val(meta), path("*.dip.vcf.gz") , emit: variant_calls
    tuple val(meta), path("*.dip.bed")    , emit: confident_regions
    tuple val(meta), path("*.paf.gz")     , emit: paf
    tuple val(meta), path("*.sam.gz")     , emit: sam
    tuple val(meta), path("*.var.gz")     , emit: var
    tuple val(meta), path("*.bam")        , emit: bam
    tuple val(meta), path("*.hap1.bed")   , emit: hap1_bed
    tuple val(meta), path("*.hap2.bed")   , emit: hap2_bed
    tuple val(meta), path("*.pair.vcf.gz"), emit: pair
    tuple val(meta), path("*.tmp")        , emit: tmp, optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if(sex == 2) {
    """
    minimap2 ${args} -c --paf-no-hit -xasm5 --cs -t ${task.cpus} ${mmi} ${haplotype_1} 2> ${prefix}.hap1.paf.gz.log | gzip > ${prefix}.hap1.paf.gz
    minimap2 ${args} -c --paf-no-hit -xasm5 --cs -t ${task.cpus} ${mmi} ${haplotype_2} 2> ${prefix}.hap2.paf.gz.log | gzip > ${prefix}.hap2.paf.gz

    minimap2 ${args} -a -xasm5 --cs -t ${task.cpus} ${mmi} ${haplotype_1} 2> ${prefix}.hap1.sam.gz.log | gzip > ${prefix}.hap1.sam.gz
    minimap2 ${args} -a -xasm5 --cs -t ${task.cpus} ${mmi} ${haplotype_2} 2> ${prefix}.hap2.sam.gz.log | gzip > ${prefix}.hap2.sam.gz

    gzip -dc ${prefix}.hap1.paf.gz | sort -k6,6 -k8,8n | k8 /usr/local/bin/paftools.js call - 2> ${prefix}.hap1.var.gz.vst | gzip > ${prefix}.hap1.var.gz
    gzip -dc ${prefix}.hap2.paf.gz | sort -k6,6 -k8,8n | k8 /usr/local/bin/paftools.js call - 2> ${prefix}.hap2.var.gz.vst | gzip > ${prefix}.hap2.var.gz

    k8 /usr/local/bin/dipcall-aux.js samflt ${prefix}.hap1.sam.gz | samtools sort --threads ${task.cpus} -o ${prefix}.hap1.bam -
    k8 /usr/local/bin/dipcall-aux.js samflt ${prefix}.hap2.sam.gz | samtools sort --threads ${task.cpus} -o ${prefix}.hap2.bam -

    gzip -dc ${prefix}.hap1.var.gz | grep ^R | cut -f2- > ${prefix}.hap1.bed
    gzip -dc ${prefix}.hap2.var.gz | grep ^R | cut -f2- > ${prefix}.hap2.bed

    htsbox pileup -q5 -evcf ${reference} ${prefix}.hap1.bam ${prefix}.hap2.bam | htsbox bgzip > ${prefix}.pair.vcf.gz

    bedtk isec -m ${prefix}.hap1.bed ${prefix}.hap2.bed > ${prefix}.dip.bed

    k8 /usr/local/bin/dipcall-aux.js vcfpair -s ${prefix} ${prefix}.pair.vcf.gz | htsbox bgzip > ${prefix}.dip.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dipcall: \$(dipcall-aux.js version)
    END_VERSIONS
    """
    } else if (sex == 1) {
    """
    minimap2 ${args} -c --paf-no-hit -xasm5  -t ${task.cpus} ${mmi} ${haplotype_1} 2> ${prefix}.hap1.paf.gz.log | gzip > ${prefix}.hap1.paf.gz
    minimap2 ${args} -c --paf-no-hit -xasm5  -t ${task.cpus} ${mmi} ${haplotype_2} 2> ${prefix}.hap2.paf.gz.log | gzip > ${prefix}.hap2.paf.gz

    minimap2 ${args} -a -xasm5  -t ${task.cpus} ${mmi} ${haplotype_1} 2> ${prefix}.hap1.sam.gz.log | gzip > ${prefix}.hap1.sam.gz
    minimap2 ${args} -a -xasm5  -t ${task.cpus} ${mmi} ${haplotype_2} 2> ${prefix}.hap2.sam.gz.log | gzip > ${prefix}.hap2.sam.gz

    gzip -dc ${prefix}.hap1.paf.gz | sort -k6,6 -k8,8n | k8 /usr/local/bin/paftools.js call - 2> ${prefix}.hap1.var.gz.vst | gzip > ${prefix}.hap1.var.gz
    gzip -dc ${prefix}.hap2.paf.gz | sort -k6,6 -k8,8n | k8 /usr/local/bin/paftools.js call - 2> ${prefix}.hap2.var.gz.vst | gzip > ${prefix}.hap2.var.gz

    k8 /usr/local/bin/dipcall-aux.js samflt ${prefix}.hap1.sam.gz | samtools sort --threads ${task.cpus} -o ${prefix}.hap1.bam -
    k8 /usr/local/bin/dipcall-aux.js samflt ${prefix}.hap2.sam.gz | samtools sort --threads ${task.cpus} -o ${prefix}.hap2.bam -

    gzip -dc ${prefix}.hap1.var.gz | grep ^R | cut -f2- > ${prefix}.hap1.bed
    gzip -dc ${prefix}.hap2.var.gz | grep ^R | cut -f2- > ${prefix}.hap2.bed

    htsbox pileup -q5 -evcf ${reference} ${prefix}.hap1.bam ${prefix}.hap2.bam | htsbox bgzip > ${prefix}.pair.vcf.gz

    bedtk isec ${prefix}.hap1.bed ${prefix}.hap2.bed | egrep -v '^(chr)?[XY]' > ${prefix}.tmp.bed
    bedtk isec ${prefix}.hap1.bed ${prefix}.hap2.bed | bedtk isec ${par} >> ${prefix}.tmp.bed
    bedtk sub ${prefix}.hap2.bed ${prefix}.hap1.bed | egrep '^(chr)?X' | bedtk sub - ${par} >> ${prefix}.tmp.bed
    bedtk sub ${prefix}.hap1.bed ${prefix}.hap2.bed | egrep '^(chr)?Y' >> ${prefix}.tmp.bed
    bedtk sort ${prefix}.tmp.bed > ${prefix}.dip.bed

    k8 /usr/local/bin/dipcall-aux.js vcfpair -s ${meta.id} -p ${par} ${prefix}.pair.vcf.gz | htsbox bgzip > ${prefix}.dip.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dipcall: \$(dipcall-aux.js version)
    END_VERSIONS
    """
    } else {
        error "Sex needs to be either '1' or '2'"
    }
}
