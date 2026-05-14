/*
Convert a csv file to tab delimitted avinput file
*/
process CSV_TO_JSON {
    input:
    tuple val(meta), path(csv) 
    path tfx_src_dir

    output:
    tuple val(meta), path("${meta.id}_out.json"), emit: json_file
    
    script:
    """
    export PYTHONPATH="${tfx_src_dir}:\${PYTHONPATH:-}"
    python ${tfx_src_dir}/edu/ohsu/compbio/txeff/util/csv_to_json.py --csv $csv --json ${meta.id}_out.json --variants
    """
}