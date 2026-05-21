/*
Convert the VEP output (tsv)  to a csv with common headings to be used for comparison
*/
process NORMALIZE_VEP_TSV {
    publishDir "${params.outdir}/normalize", mode: 'symlink'

    input:
    tuple val(meta), path(vep_csv) 
    path src_dir

    output:
    tuple val(meta), path("${meta.id}_vep_${meta.ref_mode}_normalized.csv"), emit: normalized_vep

    script:
    """
    export PYTHONPATH="${src_dir}:\${PYTHONPATH:-}"
    python ${src_dir}/nomenclature_comparison/transform/normalize_vep_tsv.py --in $vep_csv --out ${meta.id}_vep_${meta.ref_mode}_normalized.csv
    """
}    