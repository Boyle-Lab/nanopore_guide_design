//params.bedfile = "/path/to/bedfile"
//params.flank = 500
//params.distance = 4000
//params.output = "/path/to/output"
//params.threads = 8

process parallel_chopchop {
    conda '/home/crone/anaconda3/envs/chopchop'

    input:
    path bedfile
    val flank
    val min_distance
    val max_distance
    path output
    val threads
    path chopchop_path
    val run_chopchop_path

    output:
    val("process_complete"), emit: control_1

    script:
    """
    parallel -j "$threads" -d "\n" $run_chopchop_path {} "$min_distance" "$max_distance" "$flank" "$output" "$chopchop_path" :::: "$bedfile"
    """
}

process filter_chopchop{
    input:
    path output
    val distance
    val min_gc
    val max_gc
    val max_self_complementarity
    val min_efficiency_score
    val max_mm0
    val max_mm1
    val max_mm2
    val max_mm3
    val filter_chopchop_path
    val control_1

    output:
    val("process_complete"), emit: control_2

    script:
    """
    python $filter_chopchop_path --path $output --distance $distance --min_gc $min_gc --max_gc $max_gc --max_self_complementarity $max_self_complementarity \
    --min_efficiency_score $min_efficiency_score --max_mm0 $max_mm0 --max_mm1 $max_mm1 --max_mm2 $max_mm2 --max_mm3 $max_mm3
    """
}

process parallel_crispron{
    conda '/home/crone/anaconda3/envs/crispron'

    input:
    path output
    val threads
    val distance
    val run_crispron_path
    val bedtools_path
    val crispron_path
    val control_2

    output:
    val("process_complete"), emit: control_3

    script:
    """
    ls "$output"/*.upstream.bed| parallel -j "$threads" $run_crispron_path {} "$output" "$distance" "$bedtools_path" "$crispron_path"
    ls "$output"/*.downstream.bed | parallel -j "$threads" $run_crispron_path {} "$output" "$distance" "$bedtools_path" "$crispron_path"

    cat "$output"/*high_scoring_guides.*bed > "$output"/output_guides.bed

    rm "$output"/*.fasta
    rm "$output"/*bp.csv
    """

}

process filter_crispron{
    input:
    path output
    val distance
    val score_threshold
    val filter_crispron_path
    val control_3

    output:
    val("process_complete"), emit: control_4

    script:
    """
    python $filter_crispron_path --path $output --distance $distance --score_threshold $score_threshold
    """

}

process select_candidates{
    input:
    val targets_bed
    val distance
    val score_threshold
    val select_candidate_guides_path
    path output
    val control_4

    script:
    """
    python $select_candidate_guides_path --targets-bed $targets_bed --distance $distance --score-threshold $score_threshold --output $output
    """

}

workflow guide_design{
    parallel_chopchop(params.bedfile, params.flank, params.min_distance, params.max_distance,
                params.output, params.threads, params.chopchop_path, params.run_chopchop_path)
    filter_chopchop(params.output, params.min_distance, params.min_gc, params.max_gc, params.max_self_complementarity, params.min_efficiency_score, 
                    params.max_mm0, params.max_mm1, params.max_mm2, params.max_mm3, params.filter_chopchop_path, parallel_chopchop.out.control_1.collect())
    parallel_crispron(params.output, params.threads, params.min_distance, params.run_crispron_path, params.bedtools_path, params.crispron_path, filter_chopchop.out.control_2.collect())
    filter_crispron(params.output, params.min_distance, params.score_threshold, params.filter_crispron_path, parallel_crispron.out.control_3.collect())
    //select_candidates(params.targets_bed, params.distance, params.score_threshold, params.select_candidate_guides_path, params.output, filter_crispron.out.control_4.collect())
}

workflow{
    guide_design()
}