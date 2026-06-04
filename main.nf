// inputs to the function
params.fasta = params.fasta
params.template_pdb = params.template_pdb
params.outdir = params.outdir

// Split FASTA
process SPLIT_FASTA {
    publishDir "results/split_sequences", mode: 'copy'
    
    input:
    path multi_fasta

    output:
    path "*.split.fa"

    script:
    """
    awk '/^>/{f=substr(\$0,2); gsub(/[^a-zA-Z0-9_]/,"_",f); f=f ".split.fa"} {print > f}' ${multi_fasta}
    """
}

// Homology Modeling
process BUILD_HOMOLOGY_MODEL {
    publishDir "results/models", mode: 'copy'
    container 'vulcan:latest'

    input:
    path gisaid_fasta
    path template_pdb

    output:
    path "${gisaid_fasta.baseName}_model.pdb"

    script:
    """
    python ${projectDir}/scripts/run_modeller.py --sequence ${gisaid_fasta} --template ${template_pdb}
    """
}

// Data Cleaning
process CLEAN_PDB {
    publishDir "results/cleaned", mode: 'copy'
    container 'vulcan:latest'

    input:
    path pdb_file

    output:
    path "${pdb_file.baseName}_clean.pdb"

    script:
    """
    export PYTHONPATH=.:\${PYTHONPATH:-}
    python ${projectDir}/scripts/clean_pdb.py ${pdb_file} ${pdb_file.baseName}_clean.pdb
    """
}

// FoldX Repair
process REPAIR_PDB {
    publishDir "results/repair", mode: 'copy'
    container 'vulcan:latest'
    
    input:
    path model_pdb

    output:
    path "*_Repair.pdb"

    script:
    """
    foldx --command=RepairPDB --pdb=${model_pdb}
    """
}

// Stability Analysis
process CALC_STABILITY {
    publishDir "results/stability", mode: 'copy'
    container 'vulcan:latest'

    input:
    path repaired_pdb

    output:
    path "*.fxout"

    script:
    """
    foldx --command=Stability --pdb=${repaired_pdb}
    """
}

// Hamming Distance
process CALC_HAMMING {
    publishDir "results/hamming", mode: 'copy'
    container 'vulcan:latest'

    input:
    path gisaid_fasta
    path template_pdb

    output:
    path "${gisaid_fasta.baseName}_dist.txt"

    script:
    """
    python ${projectDir}/scripts/calc_hamming.py ${template_pdb} ${gisaid_fasta} > ${gisaid_fasta.baseName}_dist.txt
    """
}

process RUN_ESM2 {
    publishDir "results/esm2", mode: 'copy'
    container 'vulcan:latest'

    input:
    path split_fasta

    output:
    path "*.perplexity", emit: scores
    path "*.embeddings.pt", emit: tensors

    script:
    """
    python ${projectDir}/scripts/run_esm2.py ${split_fasta}
    """
}

process GENERATE_GRAPH {
    publishDir "results/graphs", mode: 'copy'
    container 'vulcan:latest'

    input:
    tuple val(sample_id), path(pdb_file), path(esm_embeddings)

    output:
    path "${sample_id}.graph.pt"

    script:
    """
    python ${projectDir}/scripts/pdb_to_graph.py ${pdb_file} ${esm_embeddings} ${sample_id}.graph.pt
    """
}

workflow {
    if (!params.fasta) { // Check for required input
        error "Missing required FASTA input. Run with: --fasta data/[your-input].fasta"
    }
    if (!params.template_pdb) { // Check for required input
        error "Missing required template PDB input. Run with: --template_pdb data/[your-template].pdb"
    }
    if (!params.outdir) { // Check for required input
        error "Missing required output directory input. Run with: --outdir results/[your-output-dir]"
    }   

    bundle_ch = Channel.fromPath(params.fasta)
    template_pdb = file(params.template_pdb)

    individual_fasta_ch = SPLIT_FASTA(bundle_ch).flatten()

    // Parallel Branch
    CALC_HAMMING(individual_fasta_ch, template_pdb)

    // Main Modeling Branch
    models_ch = BUILD_HOMOLOGY_MODEL(individual_fasta_ch, template_pdb)
    cleaned_ch = CLEAN_PDB(models_ch)
    repaired_ch = REPAIR_PDB(cleaned_ch)
    
    // Final Scoring
    CALC_STABILITY(repaired_ch)

    // Machine Learning Branch
    esm2_out = RUN_ESM2(individual_fasta_ch)
    
    cleaned_keyed = cleaned_ch.map { file -> [file.baseName.replace('_model_clean', ''), file] }
    esm_keyed = esm2_out.tensors.map { file -> [file.baseName.replace('.embeddings', ''), file] }
    paired_ch = cleaned_keyed.join(esm_keyed)

    GENERATE_GRAPH(paired_ch)
}
