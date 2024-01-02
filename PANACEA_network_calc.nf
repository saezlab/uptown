params.iterations = 10
params.cutoff = 6
params.network
params.source_file
params.tf_file
params.scripts_dir = projectDir
params.dataset

// set random seed
java.util.Random rand = new java.util.Random(42)

process get_bio_contexts {
    input:
    path tf_file
    
    output:
    path 'bio_contexts.txt'
    
    script:
    """
    #!/bin/bash
    # remove header and get the index names of the tf_file
    tail -n +2 $tf_file | cut -f1 > bio_contexts.txt
    """
}

process random_idgen {
    input:
    val iterations
    
    output:
    path 'random_ids.txt'
    
    script:
    """
    #!/bin/bash
    # create random identifiers and emit them individually
    echo real
    for i in {1..${iterations}}
    do
        if [ \$i -eq 1 ]
        then
            echo real
        fi
        echo \$RANDOM | md5sum | cut -c1-6
    done > random_ids.txt
    """
}

process panacea_network_calc{
    publishDir "$params.scripts_dir/results_$params.dataset", mode: 'move'
    input:
    path network
    path source_file
    path tf_file
    each random_label
    each bio_context

    // optional output
    output:
    path '*', optional: true

    script:
    """
    #!/bin/bash
    python3 $params.scripts_dir/PANACEA_analysis_cluster.py -n $network -s $source_file -t $tf_file -r $random_label -tf 20 -a 6 -b $bio_context -d $params.dataset
    """
}

workflow {
    Channel.fromPath(params.network)
        .set{network}
    Channel.fromPath(params.source_file)
        .set{source_file}
    Channel.fromPath(params.tf_file)
        .set{tf_file}
    random_idgen(params.iterations)
        .splitText()
        // remove the end of line character
        .map{it.trim()}
        .set{random_ids}
    
    get_bio_contexts(tf_file)
        .splitText()
        // remove the end of line character
        .map{it.trim()}
        // .view()
        .set{bio_contexts}
    
    panacea_network_calc(network, source_file, tf_file, random_ids, bio_contexts)
        // .view()
        .set{panacea_results}
}

