process TX_EFF_CONTROL {
    input:
    tuple val(meta), path(variants_json), path(annovar_variant_function), path(annovar_exonic_variant_function)            
    path refseq_ccds_map
    path reference_fasta
    path seqrepo_dir
    val uta_db_url
    val number_of_threads    
    path tfx_src_dir

    output:
    tuple val(meta), path("${meta.id}_tfx_out.json"), emit: tfx
    
    script:
    """
    export PYTHONPATH="${tfx_src_dir}:\${PYTHONPATH:-}"
    export HGVS_SEQREPO_DIR=${seqrepo_dir}
    export UTA_DB_URL=${uta_db_url}

    python ${tfx_src_dir}/edu/ohsu/compbio/txeff/tx_eff_control.py               \\
       --refseq_ccds_map ${refseq_ccds_map}                                      \\
       --variants ${variants_json}                                               \\
       --reference_fasta ${reference_fasta}                                      \\
       --annovar_variant_function ${annovar_variant_function}                    \\
       --annovar_exonic_variant_function ${annovar_exonic_variant_function}      \\
       --threads ${number_of_threads}                                            \\
       --out ${meta.id}_tfx_out.json
    """
}