/*
Convert the Transcript Effects output (json) to a csv with common headings to be used for comparison
*/
process NORMALIZE_TFX {
    publishDir "${params.outdir}/normalize", mode: 'symlink'

    input:
    tuple val(meta), path(tfx_json) 
    path src_dir

    output:
    tuple val(meta), path("${meta.id}_${meta.tfx_id}_tfx_normalized.csv"), emit: normalized_tfx

    script:
    """
    export PYTHONPATH="${src_dir}:\${PYTHONPATH:-}"
    python ${src_dir}/nomenclature_comparison/transform/normalize_tfx.py --in $tfx_json --out ${meta.id}_${meta.tfx_id}_tfx_normalized.csv
    """
}