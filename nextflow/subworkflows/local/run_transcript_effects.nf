include { CSV_TO_AVINPUT } from '../../modules/local/csv_to_avinput'
include { ANNOVAR_ANNOTATE_VARIATION } from '../../modules/local/annovar_annotate_variation'
include { CSV_TO_JSON } from '../../modules/local/csv_to_json'
include { TX_EFF_CONTROL  } from '../../modules/local/tx_eff_control'

workflow RUN_TRANSCRIPT_EFFECTS {

    take:
        ch_tfx_inputs   // Receives: [ meta, path(variants), path(tfx_source_dir) ]
        refseq_ccds_map
        reference_fasta
        seqrepo_dir
        uta_db_url
        number_of_tfx_threads
        annovar_dir     
        local_src_dir   
        
    main:
        // Convert the variant csv Annovar avinput format
        ch_for_avinput = ch_tfx_inputs.map { meta, variants, tfx_dir -> [ meta, variants ] }
        CSV_TO_AVINPUT(ch_for_avinput, local_src_dir)

        // Run Annovar             
        ANNOVAR_ANNOTATE_VARIATION(CSV_TO_AVINPUT.out.avinput, annovar_dir)

        // Convert csv to json
        CSV_TO_JSON(ch_tfx_inputs)

        // Extract the tfx src dir
        ch_tfx_dir_lookup = ch_tfx_inputs.map { meta, variant_file, tfx_dir -> [ meta, tfx_dir ] }

        // JOIN the two streams by their 'meta' key
        // This creates a channel that emits: [ meta, variants_json_file, [variant_function, exonic_function] ]
        //ch_combined_inputs = CSV_TO_JSON.out.json_file.join(ANNOVAR_ANNOTATE_VARIATION.out.annovar_results)
        ch_files_joined = CSV_TO_JSON.out.json_file.join(ANNOVAR_ANNOTATE_VARIATION.out.annovar_results)
        ch_combined_inputs = ch_files_joined.join(ch_tfx_dir_lookup)

        ch_final_control_inputs = ch_combined_inputs.map { meta, json, var_func, exon_func, tfx_dir ->
            // Adjust the order here to match TX_EFF_CONTROL's exact input declaration!
            return [ meta, json, var_func, exon_func, tfx_dir ]
        }
        // Run Transcript Effects     
        TX_EFF_CONTROL(
            ch_final_control_inputs,
            refseq_ccds_map,
            reference_fasta,
            seqrepo_dir,
            uta_db_url,
            number_of_tfx_threads,
        )
    
    emit:
    tfx_results = TX_EFF_CONTROL.out.tfx
}