params.iterations = 10
params.network
params.dirpath
params.tf_file
params.scripts_dir = projectDir
params.offtargets


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
    publishDir "$params.scripts_dir/results", mode: 'copy'
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

process panacea_network_eval{
    input:
    path network
    path dirpath
    path offtargets
    path tf_file
    each random_label
    each bio_context

    output:
    path "PANACEA_${random_label}_${bio_context}_graphdata_df.csv", optional: true

    script:
    """
    #!/bin/bash
    python3 $params.scripts_dir/PANACEA_evaluation_cluster.py -n $network -d $dirpath -r $random_label -o $offtargets -b $bio_context -t $tf_file
    echo "PANACEA_${random_label}_${bio_context}_graphdata_df.csv"
    """
}

process concatenate_results {
    publishDir "$params.scripts_dir/eval", mode: 'move'
    input:
    list panacea_results

    output:
    path 'all_results.csv'

    script:
    """
    #!/bin/bash
    # Use the header from the first file
    head -n 1 \$(ls -1 ${panacea_results}[0]) > all_results.csv
    
    # Concatenate all the files, skipping the header
    for file in ${panacea_results}/*
    do
        tail -n +2 \$file >> all_results.csv
    done
    """
}


workflow {
    Channel.fromPath(params.network)
        .set{network}
    Channel.fromPath(params.dirpath)
        .set{dirpath}
    Channel.fromPath(params.tf_file)
        .set{tf_file}
    Channel.fromPath(params.offtargets)
        .set{offtargets}
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
    
    panacea_network_eval(network, dirpath, offtargets, tf_file, random_ids, bio_contexts)
        // .view()
        .collectFile(keepHeader: true, skip:1, name: 'panacea_graphdata_results.csv', storeDir: './')
        .set{panacea_results}
    // python3 /home/victo/tfm/vp-mthesis/PANACEA_evaluation_cluster.py -n network_collectri.sif -d results -r real -o panacea_offtargets.tsv -b HSTS_LAPATINIB -t tf_activity_results.tsv
}
