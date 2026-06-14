// === Nextflow pipeline for analysis of H3N2 HA sequences ===

// Split FASTA for protein sequences
process SPLIT_PROTEIN_FASTA {
    publishDir "${params.outdir}/split_protein_sequences", mode: 'copy'

    input:
    path multi_fasta

    output:
    path "*.split.fa"

    script:
    """
    awk '/^>/{f=substr(\$0,2); gsub(/[^a-zA-Z0-9_]/,"_",f); f=f ".split.fa"} {print > f}' ${multi_fasta}
    """
}

// Split FASTA for nucleotide sequences
process SPLIT_NUCLEOTIDE_FASTA {
    publishDir "${params.outdir}/split_nucleotide_sequences", mode: 'copy'

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
    publishDir "${params.outdir}/models", mode: 'copy'
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
    publishDir "${params.outdir}/cleaned", mode: 'copy'
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
    publishDir "${params.outdir}/repair", mode: 'copy'
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
    publishDir "${params.outdir}/stability", mode: 'copy'
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
    publishDir "${params.outdir}/hamming", mode: 'copy'
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

process CALC_CPG {
    publishDir "${params.outdir}/cpg", mode: 'copy'
    container 'vulcan:latest'

    input:
    path nucleotide_fasta

    output:
    path "${nucleotide_fasta.baseName}_cpg.txt"

    script:
    """
    python ${projectDir}/scripts/calc_cpg.py ${nucleotide_fasta} > ${nucleotide_fasta.baseName}_cpg.txt
    """
}

process RUN_ESM2 {
    publishDir "${params.outdir}/esm2", mode: 'copy'
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
    publishDir "${params.outdir}/graphs", mode: 'copy'
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
    // Check for required inputs
    if (!params.protein) { 
        error "Missing required protein input. Run with: --protein data/[your-protein].fasta"
    }
    if (!params.nucleotide) {
        error "Missing required nucleotide input. Run with: --nucleotide data/[your-nucleotide].fasta"
    }
    if (!params.template_pdb) {
        error "Missing required template PDB input. Run with: --template_pdb data/[your-template].pdb"
    }
    if (!params.outdir) {
        error "Missing required output directory input. Run with: --outdir [your-output-dir]"
    }   

    protein_bundle_ch = Channel.fromPath(params.protein)
    nucleotide_bundle_ch = Channel.fromPath(params.nucleotide)
    template_pdb = file(params.template_pdb)

    individual_fasta_ch = SPLIT_PROTEIN_FASTA(protein_bundle_ch).flatten()
    individual_nucleotide_ch = SPLIT_NUCLEOTIDE_FASTA(nucleotide_bundle_ch).flatten()

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

    CALC_CPG(individual_nucleotide_ch)
}
