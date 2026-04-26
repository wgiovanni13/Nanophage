# NanoPhage

Nextflow pipeline for Nanopore phage genome assembly and characterization.

## Pipeline overview

```mermaid
graph TD
    %% Define Styles
    classDef main fill:#fff,stroke:#333,stroke-width:1px,rx:5,ry:5;
    classDef process fill:#e8f4fd,stroke:#2980b9,stroke-width:1.5px,color:#333;
    classDef data fill:#e8f8e8,stroke:#27ae60,stroke-width:1.5px,color:#333,rx:20,ry:20;
    classDef final fill:#fef3e2,stroke:#e67e22,stroke-width:2px,color:#333,font-weight:bold;
    classDef tool fill:#f9f9f9,stroke:#999,stroke-width:1px,stroke-dasharray: 5 5,color:#555;

    %% Main Subgraphs
    subgraph Input_Stage ["1. Input Data"]
        direction LR
        A[fa:fa-cloud-download Google Drive ID] -->|gdown| B(fa:fa-file-code Raw FASTQ)
    end

    subgraph QC_Stage ["2. Read QC & Filtering"]
        direction LR
        B -->|NanoPlot| C{fa:fa-bar-chart Read QC}
        B -->|Chopper| D(fa:fa-filter Filtered Reads)
    end

    subgraph Assembly_Stage ["3. De Novo Assembly"]
        direction LR
        D -->|Flye| E(fa:fa-dna Assembly Contigs)
    end

    subgraph Assembly_QC_Stage ["4. Assembly QC"]
        direction LR
        E -->|QUAST + minimap2 + CheckV| F(fa:fa-check-circle Assembly Quality)
    end

    subgraph Analysis_Stage ["5. Phage Characterization"]
        direction LR
        E --> Analysis
        
        subgraph Analysis_Tools
            direction TB
            genomad[geNomad: Viral ID & Taxonomy]
            blast[BLASTn: Closest Relative]
            bacphlip[BACPHLIP: Lifestyle]
            trna[tRNAscan-SE: tRNAs]
            abricate[ABRicate: AMR Screening]
            pharokka[Pharokka: Annotation]
        end
    end

    subgraph Output_Stage ["6. Integrated Report"]
        Analysis_Tools --> M((fa:fa-file-signature Final HTML Report))
    end

    %% Apply Classes
    class A,M main;
    class C,F,genomad,blast,bacphlip,trna,abricate,pharokka tool;
    class B,D,E data;
    class M final;
    class Input_Stage,QC_Stage,Assembly_Stage,Assembly_QC_Stage,Analysis_Stage,Output_Stage process;

    %% FontAwesome icons require a compatible font/rendering environment. 
    %% If icons don't render, remove 'fa:fa-...' part from labels.
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
