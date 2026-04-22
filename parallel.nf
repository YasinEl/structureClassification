#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/code"
DATA_FOLDER = "$baseDir/data"
CLASSYFIRE_TOOL = "./Classyfire_2024.jar"

params.structure_csv = "$DATA_FOLDER/classyfire_test.csv"
params.outdir        = "$baseDir/nf_output"
params.chunks        = 44

// ── 1. Split process ──────────────────────────────────────────────────────────
process prepareAndSplit {
    conda "$baseDir/envs/py_env.yml"

    input:
    path csv_file
    val  num_chunks

    output:
    path "chunk_*.csv", emit: chunks  // one file per chunk → one channel item each

    script:
    """
    # Clean/prepare the input
    python ${TOOL_FOLDER}/prepare_for_classyfire_local.py \
        ${csv_file} \
        classyfire_clean_input.csv

    # Split into N chunks, preserving the header in each
    total_lines=\$(tail -n +2 classyfire_clean_input.csv | wc -l)
    lines_per_chunk=\$(( (\$total_lines + ${num_chunks} - 1) / ${num_chunks} ))

    head -1 classyfire_clean_input.csv > header.csv
    tail -n +2 classyfire_clean_input.csv | split -l \$lines_per_chunk -d -a 2 - chunk_body_

    for file in chunk_body_*; do
        cat header.csv \$file > \${file}.csv
        rm \$file
    done
    rm header.csv
    """
}

// ── 2. Per-chunk process (runs N times, fully in parallel) ────────────────────
process runClassyfireLocal {
    conda    "$baseDir/envs/py_env.yml"
    tag      "$chunk"
    maxRetries 2

    // Increase memory on retry — classic nf-core pattern
    memory   { 4.GB * task.attempt }
    time     { 2.hour * task.attempt }

    errorStrategy { task.exitStatus in [130, 137, 140] ? 'retry' : 'terminate' }
    // 130=OOM killed, 137=SIGKILL, 140=timeout

    input:
    path chunk

    output:
    path "${chunk.baseName}.tsv", emit: chunk_results

    script:
    """
    java -jar ${CLASSYFIRE_TOOL} ${chunk} ${chunk.baseName}.tsv
    """
}

// ── 3. Merge all chunk outputs back into one directory ────────────────────────
process mergeResults {
    conda "$baseDir/envs/py_env.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    path tsv_files

    output:
    path "Classyfire_merged.tsv"

    script:
    """
    # Get header from any file
    head -1 \$(ls *.tsv | sort | head -1) > Classyfire_merged.tsv
    # Append all data rows, sorted by filename for reproducibility  
    for f in \$(ls chunk_*.tsv | sort); do
        tail -n +2 \$f >> Classyfire_merged.tsv
    done
    """
}

// ── Workflow ──────────────────────────────────────────────────────────────────
workflow {
    csv_ch = Channel.fromPath(params.structure_csv)

    // Split once → emits a *list* of chunk files
    split_ch = prepareAndSplit(csv_ch, params.chunks)

    // .flatten() is the key: turns [chunk_aa.csv, chunk_ab.csv, ...]
    // into individual channel items → Nextflow spawns one process per item
    chunk_ch = split_ch.chunks.flatten()

    // Each chunk runs as an independent, parallel job
    results_ch = runClassyfireLocal(chunk_ch)

    // Collect all results, then merge once everything is done
    mergeResults(results_ch.chunk_results.collect())
}