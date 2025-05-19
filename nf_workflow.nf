#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/code"
DATA_FOLDER = "$baseDir/data"
CLASSYFIRE_TOOL = "$baseDir/classyfire/Classyfire_2024.jar"



params.structure_csv = "$baseDir/data/gnps_cleaned.csv"

params.output_directory_Classyfire = "$baseDir/nf_output/Classyfire"
params.output_directory_Classyfire_local = "$baseDir/nf_output/Classyfire_local.tsv"
params.output_directory_Npclassifier = "$baseDir/nf_output/Npclassifier"
params.output_directory_ChemInfoService = "$baseDir/nf_output/ChemInfoService"

params.log_Classyfire = "$baseDir/nf_output/Classyfire.log"
params.log_Classyfire_local = "$baseDir/nf_output/Classyfire_local.log"
params.log_Npclassifier = "$baseDir/nf_output/Npclassifier.log"
params.log_ChemInfoService = "$baseDir/nf_output/ChemInfoService.log"

process getClassyfire {
    conda "$baseDir/envs/py_env.yml"

    input:
    val x

    script:
    """
    python $TOOL_FOLDER/batchprocess_classyfire.py \
    $params.structure_csv \
    --outdir $params.output_directory_Classyfire \
    --logfile $params.log_Classyfire \
    --rate-limit 0.3
    """
}

process getClassyfire_local {
    conda "$baseDir/envs/py_env.yml"

    input:
    val x

    script:
    """
    python $TOOL_FOLDER/prepare_for_classyfire_local.py \
    $params.structure_csv \
    classyfire_input.csv \
    --exclude-tsv $params.log_Classyfire_local

    java -jar $CLASSYFIRE_TOOL \
    classyfire_input.csv \
    $params.output_directory_Classyfire_local

    python $TOOL_FOLDER/classyfire_local_logging.py \
    $params.output_directory_Classyfire_local \
    $params.log_Classyfire_local
    """
}

process getNpclassifier {
    conda "$baseDir/envs/py_env.yml"

    input:
    val x

    script:
    """
    python $TOOL_FOLDER/batchprocess_npclassifier.py \
    $params.structure_csv \
    --outdir $params.output_directory_Npclassifier \
    --logfile $params.log_Npclassifier \
    --rate-limit 100
    """
}

process getChemInfoService {
    conda "$baseDir/envs/py_env.yml"

    input:
    val x

    script:
    """
    python $TOOL_FOLDER/batchprocess_chemInfoService.py \
    $params.structure_csv \
    --outdir $params.output_directory_ChemInfoService \
    --logfile $params.log_ChemInfoService \
    --rate-limit 10
    """
}



workflow {

    getClassyfire(1)
    getNpclassifier(1)
    getChemInfoService(1)
    getClassyfire_local(1)
    
}
