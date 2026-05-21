/*
Convert a csv file to tab delimitted avinput file
*/
process CSV_TO_AVINPUT {
    input:
    tuple val(meta), path(csv) 
    path src_dir

    output:
    tuple val(meta), path("${meta.id}.avinput"), emit: avinput
    
    script:
    """
    export PYTHONPATH="${src_dir}:\${PYTHONPATH:-}"
    python ${src_dir}/nomenclature_comparison/transform/csv_to_avinput.py --in $csv --out ${meta.id}.avinput
    """
}