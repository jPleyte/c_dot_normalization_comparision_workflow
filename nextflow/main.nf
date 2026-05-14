/* 
The C Dot Normalization Comparison Workflow starts with a list of variants 
and sends those variants to any number of tools that generate the predicted 
effects cDNA and protein transcripts, and produces a spreadsheet showing 
the transcript effects from each tool. 
*/
include { RUN_TRANSCRIPT_EFFECTS } from './subworkflows/local/run_transcript_effects'

workflow {
    main:
        log.info "Variant source: ${params.variant_source_file_csv}"
        log.info "tfx source: ${params.compbio_cgd_tx_eff_src_dir}"

        ch_variants = channel.fromPath(params.variant_source_file_csv)
            .map { 
                    file -> def meta = [ id: file.baseName, source: params.variant_source_file_csv ]
                    return [ meta, file ] 
            }

        // ANNOVAR installation directory
        annovar_dir = file(params.annovar_dir, checkIfExists: true)

        // Python src directories
        local_src_dir = channel.value(file(params.local_python_src_dir, checkIfExists: true))
        tfx_source_dir = channel.value(file(params.compbio_cgd_tx_eff_src_dir, checkIfExists: true))
        

        // Transcript Effects configuration 
        refseq_ccds_map = channel.value(file(params.refseq_ccds_map, checkIfExists: true))
        reference_fasta = channel.value(file(params.reference_fasta_file, checkIfExists: true))
        seqrepo_dir = channel.value(file(params.seq_repo_dir, checkIfExists: true))
        uta_db_url = channel.value(params.uta_db_url ?: error("The uta_db_url parameter is required"))
        tfx_threads = channel.value(params.number_of_tfx_threads ?: 3)

        RUN_TRANSCRIPT_EFFECTS(ch_variants, 
                               refseq_ccds_map, 
                               reference_fasta,
                               seqrepo_dir,
                               uta_db_url,
                               tfx_threads,
                               annovar_dir,
                               local_src_dir,
                               tfx_source_dir)
}    