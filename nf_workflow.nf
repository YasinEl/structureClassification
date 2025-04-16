#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/code"
DATA_FOLDER = "$baseDir/data"



params.structure_csv = "$baseDir/data/gnps_cleaned.csv"

params.output_directory_Classyfire = "$baseDir/nf_output/Classyfire"
params.output_directory_Npclassifier = "$baseDir/nf_output/Npclassifier"
params.output_directory_ChemInfoService = "$baseDir/nf_output/ChemInfoService"

params.log_Classyfire = "$baseDir/nf_output/Classyfire.log"
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
    --rate-limit 10
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
    
}
