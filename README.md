# NanoPhage

Nextflow pipeline for Nanopore phage genome assembly and characterization.

## Pipeline overview

```mermaid
graph LR
    %% Nodos Principales
    Start((fa:fa-play)) --> A
    
    A[fa:fa-cloud-download Google Drive] -- gdown --> B[fa:fa-file-code Raw FASTQ]
    
    B --> C{QC & Filter}
    C -- NanoPlot --> QC[fa:fa-chart-bar Read QC]
    C -- Chopper --> D[fa:fa-filter Filtered Reads]
    
    D -- Flye --> E[fa:fa-dna Assembly]
    
    E --> F[fa:fa-microscope Assembly QC]
    F --- F1[QUAST]
    F --- F2[minimap2]
    F --- F3[CheckV]
    
    E --> G[fa:fa-project-diagram Characterization]
    
    %% Herramientas de Caracterización en abanico
    G --> G1[geNomad]
    G --> G2[BLASTn]
    G --> G3[BACPHLIP]
    G --> G4[tRNAscan-SE]
    G --> G5[ABRicate]
    G --> G6[Pharokka]
    
    %% Unión al Reporte
    F1 & F2 & F3 & G1 & G2 & G3 & G4 & G5 & G6 --> M
    
    M[[fa:fa-file-alt HTML REPORT]]

    %% Estilos Mejorados
    style Start fill:#333,stroke:#fff,color:#fff
    style A fill:#4285F4,stroke:#fff,color:#fff
    style B fill:#f8f9fa,stroke:#4285F4
    style E fill:#27ae60,stroke:#fff,color:#fff,stroke-width:2px
    style G fill:#f39c12,stroke:#fff,color:#fff
    style M fill:#e67e22,stroke:#fff,color:#fff,stroke-width:3px
    
    %% Colores suaves para herramientas
    style G1,G2,G3,G4,G5,G6 fill:#fffef0,stroke:#f39c12
    style F1,F2,F3 fill:#f0fff4,stroke:#27ae60
```

## Quick start

```bash
git clone https://github.com/wgiovanni13/nanophage.git
cd nanophage
nextflow run main.nf
```

All tools are installed automatically via Conda. Databases (CheckV, geNomad, Pharokka) are downloaded during execution.

## Requirements

- Nextflow ≥ 24.10
- Conda or Mamba
- Internet access (for Google Drive download, remote BLAST, and database downloads)

## Parameters

| Parameter    | Default | Description                     |
|--------------|---------|---------------------------------|
| `--file_id`  | (set)   | Google Drive file ID            |
| `--filename` | (set)   | Output filename for raw reads   |
| `--outDir`   | results | Output directory                |
| `--min_len`  | 500     | Minimum read length             |
| `--min_qual` | 10      | Minimum quality score (Phred)   |

## Output

```
results/
├── 00_raw/                    Downloaded FASTQ
├── 01_read_qc/                NanoPlot stats and plots
├── 02_filtered/               Filtered reads and report
├── 03_assembly/               Flye assembly and info
├── 04_assembly_qc/            QUAST, coverage, CheckV
├── 05_characterization/       geNomad, BLAST, BACPHLIP, tRNAs, AMR, Pharokka
└── 06_report/                 nanophage_report.html
```

## Cluster configuration

Edit `nextflow.config` to match your HPC:

```nextflow
process {
    executor = 'slurm'    // or 'local', 'pbs', 'sge' depending on your HPC queue's names
    queue    = 'm'        // or 'batch' depending on your HPC queue's names
    qos      = 'medium'   // or 'm' depending on your HPC queue's names
}
```

## Author

Wagner Giovanni Guzman Mendez
