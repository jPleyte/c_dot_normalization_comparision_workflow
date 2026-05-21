/* 
The C Dot Normalization Comparison Workflow starts with a list of variants 
and sends those variants to any number of tools that generate the predicted 
effects cDNA and protein transcripts, and produces a spreadsheet showing 
the transcript effects from each tool. 
*/
include { RUN_TRANSCRIPT_EFFECTS } from './subworkflows/local/run_transcript_effects'
include { CSV_TO_VCF } from './modules/local/csv_to_vcf'
include { VEP } from './modules/local/vep'
include { NORMALIZE_TFX } from './modules/local/normalize_tfx'
include { NORMALIZE_VEP_VCF } from './modules/local/normalize_vep_vcf'
include { COMPARE_RESULTS } from './modules/local/compare_results'

workflow {
    main:
        log.info "Variant source: ${params.variant_source_file_csv}"
        log.info "${params.compbio_cgd_tx_eff_src_label_one}: ${params.compbio_cgd_tx_eff_src_dir_one}"
        log.info "${params.compbio_cgd_tx_eff_src_label_two}: ${params.compbio_cgd_tx_eff_src_dir_two}"

        ch_variants = channel.fromPath(params.variant_source_file_csv)
            .map { 
                    file -> def meta = [ id: file.baseName, source: params.variant_source_file_csv ]
                    return [ meta, file ] 
            }

        // ANNOVAR installation directory
        annovar_dir = file(params.annovar_dir, checkIfExists: true)

        // Python src directories
        local_src_dir = channel.value(file(params.local_python_src_dir, checkIfExists: true))        
        
        // Transcript Effects configuration 
        refseq_ccds_map = channel.value(file(params.refseq_ccds_map, checkIfExists: true))
        reference_fasta = channel.value(file(params.reference_fasta_file, checkIfExists: true))
        seqrepo_dir = channel.value(file(params.seq_repo_dir, checkIfExists: true))
        uta_db_url = channel.value(params.uta_db_url ?: error("The uta_db_url parameter is required"))
        tfx_threads = channel.value(params.number_of_tfx_threads ?: 3)

        ch_tfx_sources = channel.of(
            [ [id: params.compbio_cgd_tx_eff_src_label_one], file(params.compbio_cgd_tx_eff_src_dir_one, checkIfExists: true) ],
            [ [id: params.compbio_cgd_tx_eff_src_label_two], file(params.compbio_cgd_tx_eff_src_dir_two, checkIfExists: true) ])

        // Combine the two inputs: variant list(s) and tfx source directories
        // This creates a channel that looks like 
        // 1. [ variants, [id: 'run_one'], dir_one ]
        // 2. [ variants, [id: 'run_two'], dir_two ]
        // ch_tfx_inputs = ch_variants.combine(ch_tfx_sources)
        //     .map { variants, meta, tfx_source_dir -> [ meta, variants, tfx_source_dir ] }
        ch_tfx_inputs = ch_variants
            .combine(ch_tfx_sources)
            .map { variant_meta, variant_file, tfx_meta, tfx_dir -> 
                // Merge the TFX label into your variant metadata map 
                // so you don't lose track of which directory is being run!
                def unified_meta = variant_meta + [ tfx_id: tfx_meta.id ]
        
                    // 2. Return the clean 3-item tuple your subworkflow expects
                    return [ unified_meta, variant_file, tfx_dir ]
                }
        
        RUN_TRANSCRIPT_EFFECTS(ch_tfx_inputs, 
                                refseq_ccds_map, 
                                reference_fasta,
                                seqrepo_dir,
                                uta_db_url,
                                tfx_threads,
                                annovar_dir,
                                local_src_dir)

        NORMALIZE_TFX(RUN_TRANSCRIPT_EFFECTS.out.tfx_results, local_src_dir)
        
        ch_tfx_labeled = NORMALIZE_TFX.out.normalized_tfx        
        .map { meta, file -> 
            // I'm tring "${meta.tfx_id}" but maybe ${meta.id}_${meta.tfx_id} would be better 
            return [ meta, "tfx_${meta.tfx_id}", file ] 
        }

        // VEP Configuration 
        vep_bin = channel.value(file(params.vep_bin, checkIfExists: true))
        vep_data_dir = channel.value(file(params.vep_data_dir, checkIfExists: true))
        vep_fasta = channel.value(file(params.vep_fasta_file, checkIfExists: true))
        
        // We want to run VEP twice, using each of these flags so a channel is created that combines the input file and the parameter.
        def vep_mode_labels = [
            "--use_transcript_ref": "refseq",
            "--use_given_ref": "fasta"
        ]        
        ch_vep_reference_modes = channel.fromList(vep_mode_labels.keySet().toList())

        CSV_TO_VCF(ch_variants, local_src_dir)

        ch_vep_inputs = CSV_TO_VCF.out.vcf
            .combine(ch_vep_reference_modes)
            .map { meta, vcf, flag ->
                def new_meta = meta.clone()
                new_meta.ref_mode = vep_mode_labels[flag] // Adds 'ref_mode' field with value 'refseq' or 'fasta'
                new_meta.ref_flag = flag
                return [ new_meta, vcf ]
            }
        
        VEP(ch_vep_inputs,
            vep_bin,
            vep_data_dir,
            vep_fasta
        )        

        // There are two VEP files, this debug statements let's us confirm it.
        VEP.out.vep_results.view { meta, file -> "DEBUG: Sending VEP ref_mode=${meta.ref_mode} for variants=${meta.id}" }
        NORMALIZE_VEP_VCF(VEP.out.vep_results, local_src_dir)

        ch_vep_labeled = NORMALIZE_VEP_VCF.out.normalized_vep
            .map { meta, file -> [ meta, "vep_${meta.ref_mode}", file ] }

        // Mix all 4 files into one channel 
        ch_all_labeled = ch_tfx_labeled.mix(ch_vep_labeled)

        // Your groupTuple block remains identical and catches all 4!
        ch_comparison_inputs = ch_all_labeled
            .map { meta, label, file -> [ meta.id, meta, label, file ] }
            .groupTuple(by: 0)
            .map { id, metas, labels, files -> [ metas[0], labels, files ] }
        
        // Log the variant list(s) and annotation tools. 
        ch_comparison_inputs.view { meta, labels, files -> "DEBUG: ${meta.id} was processed by ${labels}" }

        COMPARE_RESULTS(ch_comparison_inputs, local_src_dir)
}