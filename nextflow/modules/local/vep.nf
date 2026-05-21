/*
 * Run VEP on a variant list in VCF format 
 * The flag parameter can be "--use_transcript_ref" or "--use_given_ref":
 *   --use_transcript_ref (The Default): VEP will use the sequence from the transcript record. 
 *                                         This is usually what you want for HGVS nomenclature because it ensures the "c." and "p." coordinates make sense.
 *   --use_given_ref: This would force VEP to use your FASTA sequence instead.
 * 
 * Notice that the second and third parameters "vep_bin" and "vep_data_dir" have "stageAs" aliases.  This is because the vep bin and vep data directory are both named "vep" so when nextflow creates symlinks to them in the run directory they conflict with eachother. 
*/
process VEP {
    publishDir "${params.outdir}/vep", mode: 'symlink'

    input:    
    tuple val(meta), path(variants_vcf)
    path vep_bin, stageAs: 'vep_executable'
    path vep_data_dir, stageAs: 'vep_data_dir'
    path vep_fasta
    
    output:
    // tuple val(meta), path("${meta.id}_vep_${ref_label}_output.tsv"), emit: vep_results
    tuple val(meta), path("${meta.id}_vep_${meta.ref_mode}_out.vcf"), emit: vep_results

    script:
    assert meta.ref_flag == "--use_transcript_ref" || meta.ref_flag == "--use_given_ref" : "Invalid VEP flag: ${meta.ref_flag}. Must be --use_transcript_ref or --use_given_ref"    
    """
    ./vep_executable             \\
    --dir vep_data_dir           \\
    --dir_cache vep_data_dir     \\
    --fasta ${vep_fasta}         \\
    ${meta.ref_flag}             \\
    --offline                    \\
    --merged                     \\
    --species homo_sapiens       \\
    --assembly GRCh37            \\
    --hgvs                       \\
    --hgvsg                      \\
    --ccds                       \\
    --symbol                     \\
    --numbers                    \\
    --biotype                    \\
    --vcf                        \\
    --protein                    \\
    --shift_3prime 1             \\
    -i ${variants_vcf}           \\
    -o ${meta.id}_vep_${meta.ref_mode}_out.vcf
    """
}