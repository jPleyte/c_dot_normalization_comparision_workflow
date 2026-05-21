/*
Run annotate_variation.pl to generate annovar.variant_function and annovar.exonic_variant_function
*/
process ANNOVAR_ANNOTATE_VARIATION {
    publishDir "${params.outdir}/annovar", mode: 'symlink'

    input:
    tuple val(meta), path(avinput) 
    path annovar_dir

    output:
    tuple val(meta), path("${meta.id}_annovar.variant_function"), path("${meta.id}_annovar.exonic_variant_function"), emit: annovar_results
    
    script:
    """
    ${annovar_dir}/annotate_variation.pl -geneanno              \\
                                    -out ${meta.id}_annovar     \\
                                    -build hg19                 \\
                                    ${avinput}                  \\
                                    ${annovar_dir}/humandb/     \\
                                    --transcript_function       \\
                                    --hgvs                      \\
                                    --separate                  \\
                                    --dbtype                    \\
                                    refGeneWithVer              \\
                                    --otherinfo                 \\
                                    --splicing_threshold 5      \\
                                    --exonicsplicing
    """
}