process ADD_FOUND_IN_TAG {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf, stageAs: "?/*"), path(index, stageAs: "?/*")
    val(variant_caller)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*.tbi")                    , emit: tbi, optional: true
    tuple val(meta), path("*.csi")                    , emit: csi, optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args2.contains("--output-type b") || args2.contains("-Ob") ? "bcf.gz" :
                    args2.contains("--output-type u") || args2.contains("-Ou") ? "bcf" :
                    args2.contains("--output-type z") || args2.contains("-Oz") ? "vcf.gz" :
                    args2.contains("--output-type v") || args2.contains("-Ov") ? "vcf" :
                    "vcf"

    """
    bcftools view \\
        $args \\
        --threads $task.cpus \\
        $vcf |\\
    awk '
        # Print all existing header lines
        /^##/ {
            print;
            next;
        }

        # Add a INFO header line before the #CHROM line, then print the #CHROM line
        /^#CHROM/ {
            printf("##INFO=<ID=FOUND_IN,Number=1,Type=String,Description=\\"Program that called the variant\\">\\n");
            print;
            next;
        }

        # Then for all other lines, append FOUND_IN=variant_caller to the INFO field (column 8)
        {
            OFS="\\t";
            \$8 = sprintf("%s;FOUND_IN=${variant_caller}", \$8 );
            print;
        }
    ' |\\
    bcftools view \\
        $args2 \\
        --threads $task.cpus \\
        --output ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        busybox-awk: \$(awk 2>&1 | sed -n 's/.*v\\([^ ]*\\) (.*/\\1/p')
    END_VERSIONS
    """

    stub:
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args2.contains("--output-type b") || args2.contains("-Ob") ? "bcf.gz" :
                    args2.contains("--output-type u") || args2.contains("-Ou") ? "bcf" :
                    args2.contains("--output-type z") || args2.contains("-Oz") ? "vcf.gz" :
                    args2.contains("--output-type v") || args2.contains("-Ov") ? "vcf" :
                    "vcf"
    def index = args2.contains("--write-index=tbi") || args2.contains("-W=tbi") ? "tbi" :
                args2.contains("--write-index=csi") || args2.contains("-W=csi") ? "csi" :
                args2.contains("--write-index") || args2.contains("-W") ? "csi" :
                ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && index.matches("csi|tbi") ? "touch ${prefix}.${extension}.${index}" : ""

    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        busybox-awk: \$(awk 2>&1 | sed -n 's/.*v\\([^ ]*\\) (.*/\\1/p')
    END_VERSIONS
    """
}
