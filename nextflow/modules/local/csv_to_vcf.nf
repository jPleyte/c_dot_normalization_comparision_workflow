/*
 * Convert csv variant list to vcf
 */
process CSV_TO_VCF {
    publishDir "${params.outdir}/intermediate", mode: 'symlink'

    input:
    tuple val(meta), path(csv) 
    path src_dir
    
    output:
    tuple val(meta), path("${meta.id}_variants_sorted.vcf"), emit: vcf
    
    script:
    """
    export PYTHONPATH="${src_dir}:\${PYTHONPATH:-}"
    python ${src_dir}/nomenclature_comparison/transform/csv_to_vcf.py --in $csv --out ${meta.id}_variants_unsorted.vcf
    vcf-sort ${meta.id}_variants_unsorted.vcf > ${meta.id}_variants_sorted.vcf
    """
}