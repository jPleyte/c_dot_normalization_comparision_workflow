include { CSV_TO_AVINPUT } from '../../modules/local/csv_to_avinput'
include { ANNOVAR_ANNOTATE_VARIATION } from '../../modules/local/annovar_annotate_variation'
include { CSV_TO_JSON } from '../../modules/local/csv_to_json'
include { TX_EFF_CONTROL  } from '../../modules/local/tx_eff_control'

workflow RUN_TRANSCRIPT_EFFECTS {

    take:
        variants_ch // This is a tuple: [ val(meta), path(varaints_csv_file) ]
        refseq_ccds_map
        reference_fasta
        seqrepo_dir
        uta_db_url
        number_of_tfx_threads
        annovar_dir     
        local_src_dir   
        tfx_source_dir    

    main:
        // Convert the variant csv Annovar avinput format
        CSV_TO_AVINPUT(variants_ch, local_src_dir)

        // Run Annovar             
        ANNOVAR_ANNOTATE_VARIATION(CSV_TO_AVINPUT.out.avinput, annovar_dir)

        // Convert csv to json
        CSV_TO_JSON(variants_ch, tfx_source_dir)

        // JOIN the two streams by their 'meta' key
        // This creates a channel that emits: [ meta, variants_json_file, [variant_function, exonic_function] ]
        ch_combined_inputs = CSV_TO_JSON.out.json_file.join(ANNOVAR_ANNOTATE_VARIATION.out.annovar_results)

        // Run Transcript Effects 
        TX_EFF_CONTROL(
            ch_combined_inputs,
            refseq_ccds_map,
            reference_fasta,
            seqrepo_dir,
            uta_db_url,
            number_of_tfx_threads,
            tfx_source_dir
        )
    
    emit:
    tfx_results = TX_EFF_CONTROL.out.tfx
}