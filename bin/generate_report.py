#!/usr/bin/env python3
"""NanoPhage — HTML report with embedded figures and biological interpretation."""

import argparse, csv, os, base64, io, sys
from datetime import datetime

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

def read_lines(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [l.strip() for l in f if l.strip()]

def parse_nanostats(path):
    d = {}
    for line in read_lines(path):
        if '\t' in line:
            parts = line.split('\t', 1)
            d[parts[0].strip()] = parts[1].strip()
    return d

def parse_tsv_rows(path):
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rows.append(row)
    return rows

def parse_kv_tsv(path):
    d = {}
    for line in read_lines(path):
        parts = line.split('\t')
        if len(parts) >= 2:
            d[parts[0]] = parts[1]
    return d

def parse_filtering(path):
    d = {}
    for line in read_lines(path):
        if ':' in line:
            k, v = line.split(':', 1)
            d[k.strip()] = v.strip()
    return d

def parse_assembly_info(path):
    for line in read_lines(path):
        if line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) >= 4:
            return {'name': parts[0], 'length': parts[1], 'coverage': parts[2], 'circular': parts[3]}
    return {}

def parse_bacphlip(path):
    lines = read_lines(path)
    for line in lines:
        parts = line.split()
        # File format: index  virulent_score  temperate_score
        # e.g.: 0  0.9425660792951543  0.057433920704845806
        if len(parts) >= 3:
            try:
                v = float(parts[1])
                t = float(parts[2])
                if 0 <= v <= 1 and 0 <= t <= 1:
                    return {'virulent': v, 'temperate': t}
            except ValueError:
                continue
    return {}

def parse_trna(path):
    trnas = []
    for line in read_lines(path):
        if line.startswith('Sequence') or line.startswith('Name') or line.startswith('-'):
            continue
        parts = line.split()
        if len(parts) >= 6:
            trnas.append({'type': parts[4], 'anticodon': parts[5], 'start': parts[2], 'end': parts[3]})
    return trnas

def parse_blast(path):
    cols = ['qseqid','sseqid','pident','length','mismatch','gapopen',
            'qstart','qend','sstart','send','evalue','bitscore','stitle']
    rows = []
    for line in read_lines(path):
        if line.startswith('no_hits'):
            return []
        parts = line.split('\t')
        if len(parts) >= len(cols):
            rows.append(dict(zip(cols, parts)))
    return rows

def parse_genomad_summary(path):
    rows = parse_tsv_rows(path)
    return rows[0] if rows else {}

def parse_genomad_genes(path):
    return parse_tsv_rows(path)

def parse_depth(path):
    positions, depths = [], []
    for line in read_lines(path):
        parts = line.split('\t')
        if len(parts) >= 3:
            positions.append(int(parts[1]))
            depths.append(int(parts[2]))
    return positions, depths

def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')

def img_file_to_base64(path):
    if not os.path.exists(path):
        return ""
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')

def make_coverage_plot(positions, depths):
    if not HAS_MPL or not positions:
        return ""
    step = max(1, len(positions) // 2000)
    pos_sub = positions[::step]
    dep_sub = depths[::step]

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.fill_between(pos_sub, dep_sub, alpha=0.6, color='#3498db')
    ax.plot(pos_sub, dep_sub, linewidth=0.3, color='#2c3e50')
    ax.set_xlabel('Genome position (bp)')
    ax.set_ylabel('Read depth')
    ax.set_xlim(0, max(positions))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    return fig_to_base64(fig)

def make_gene_map(genes, genome_len):
    if not HAS_MPL or not genes:
        return ""

    func_colors = {
        'packaging': '#e74c3c',
        'structural': '#3498db',
        'replication': '#2ecc71',
        'lysis': '#e67e22',
        'modification': '#9b59b6',
        'other_viral': '#1abc9c',
        'hypothetical': '#bdc3c7',
    }

    def classify_gene(desc):
        desc_l = desc.lower() if desc else ''
        if 'terminase' in desc_l: return 'packaging'
        if any(w in desc_l for w in ['capsid', 'tail', 'head', 'baseplate', 'portal']): return 'structural'
        if any(w in desc_l for w in ['polymerase', 'helicase', 'primase', 'exonuclease', 'crispr', 'dna pol']): return 'replication'
        if any(w in desc_l for w in ['lysin', 'holin', 'lysis', 'hydrolase', 'peptidase']): return 'lysis'
        if any(w in desc_l for w in ['methylase', 'methyl', 'modification']): return 'modification'
        if any(w in desc_l for w in ['phage', 'mu', 'ninx', 'yopx', 'protease']): return 'other_viral'
        if desc == 'NA' or not desc or 'hypothetical' in desc_l or 'uncharacterized' in desc_l or 'unknown' in desc_l: return 'hypothetical'
        return 'other_viral'

    fig, ax = plt.subplots(figsize=(12, 2.5))
    for g in genes:
        try:
            start = int(g.get('start', 0))
            end = int(g.get('end', 0))
            strand = int(g.get('strand', 1))
        except (ValueError, TypeError):
            continue
        desc = g.get('annotation_description', 'NA')
        cat = classify_gene(desc)
        color = func_colors.get(cat, '#bdc3c7')
        y = 0.6 if strand == 1 else 0.2
        h = 0.15
        ax.add_patch(plt.Rectangle((start, y), end - start, h, facecolor=color, edgecolor='#2c3e50', linewidth=0.3))

    ax.axhline(y=0.5, color='#2c3e50', linewidth=1.5)
    ax.set_xlim(0, genome_len)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Genome position (bp)')
    ax.set_yticks([0.275, 0.675])
    ax.set_yticklabels(['Reverse', 'Forward'], fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    patches = [mpatches.Patch(color=c, label=l.replace('_', ' ').title()) for l, c in func_colors.items()]
    ax.legend(handles=patches, loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=4, fontsize=7, frameon=False)
    return fig_to_base64(fig)

def build_report(args):
    raw = parse_nanostats(args.raw_stats)
    filt = parse_filtering(args.filter_stats)
    asm_info = parse_assembly_info(args.assembly_info)
    quast = parse_kv_tsv(args.quast)
    cov_rows = parse_tsv_rows(args.coverage)
    cov = cov_rows[0] if cov_rows else {}
    checkv_rows = parse_tsv_rows(args.checkv)
    cv = checkv_rows[0] if checkv_rows else {}
    genomad_virus = parse_genomad_summary(args.genomad_summary)
    genes = parse_genomad_genes(args.genomad_genes)
    blast = parse_blast(args.blast)
    bac = parse_bacphlip(args.bacphlip)
    trnas = parse_trna(args.trna)
    amr_rows = parse_tsv_rows(args.amr)
    positions, depths = parse_depth(args.depth)
    pharokka_b64 = img_file_to_base64(args.pharokka_plot) if args.pharokka_plot else ""

    genome_len = int(asm_info.get('length', quast.get('Total length', 0)))
    total_genes = len(genes)
    hallmarks = len([g for g in genes if g.get('virus_hallmark') == '1'])
    annotated = len([g for g in genes if g.get('annotation_description', 'NA') != 'NA'])
    amr_found = sum(1 for r in amr_rows if r.get('NUM_FOUND', '0') != '0' and '#FILE' not in r.get('#FILE', ''))
    trna_text = ", ".join(f"tRNA-{t['type']}({t['anticodon']})" for t in trnas) if trnas else "none detected"

    vir = bac.get('virulent', 0)
    temp = bac.get('temperate', 0)
    lifestyle = "Virulent (lytic)" if vir > 0.5 else "Temperate (lysogenic)"
    lifestyle_pct_main = vir * 100 if vir > 0.5 else temp * 100

    cov_plot_b64 = make_coverage_plot(positions, depths)
    gene_map_b64 = make_gene_map(genes, genome_len)

    top_subject = None
    top_hsps = []
    if blast:
        top_subject = blast[0].get('sseqid', '')
        for b in blast:
            if b.get('sseqid') == top_subject:
                top_hsps.append(b)
    top_name = top_hsps[0].get('stitle', 'Unknown').split(',')[0] if top_hsps else 'Unknown'
    top_accession = top_subject.split('|')[3] if top_subject and '|' in top_subject else top_subject

    blast_top_html = ""
    total_aligned = 0
    for h in top_hsps:
        blast_top_html += f"<tr><td>{h['qstart']}&ndash;{h['qend']}</td><td>{h['sstart']}&ndash;{h['send']}</td><td>{float(h['pident']):.3f}%</td><td>{int(h['length']):,}</td><td>{h['evalue']}</td></tr>\n"
        total_aligned += int(h['length'])

    gene_rows_html = ""
    for g in genes:
        desc = g.get('annotation_description', 'NA')
        if desc == 'NA':
            continue
        strand = '+' if g.get('strand') == '1' else '&minus;'
        h = '&#9679;' if g.get('virus_hallmark') == '1' else ''
        gene_rows_html += f"<tr><td>{g['gene']}</td><td>{g['start']}&ndash;{g['end']}</td><td>{strand}</td><td>{h}</td><td>{desc}</td></tr>\n"

    now = datetime.now()

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>NanoPhage Report</title>
<style>
body {{ font-family: 'Segoe UI', system-ui, sans-serif; max-width: 1000px;
       margin: 40px auto; padding: 0 24px; color: #2c3e50; line-height: 1.7; background: #fafafa; }}
h1 {{ color: #1a5276; border-bottom: 3px solid #2980b9; padding-bottom: 12px; margin-top: 0; }}
h2 {{ color: #2c3e50; margin-top: 36px; border-bottom: 1px solid #bdc3c7; padding-bottom: 6px; }}
h3 {{ color: #34495e; margin-top: 20px; }}
table {{ border-collapse: collapse; width: 100%; margin: 12px 0; font-size: 0.9em; }}
th, td {{ padding: 6px 10px; text-align: left; border-bottom: 1px solid #e0e0e0; }}
th {{ background: #2c3e50; color: white; font-weight: 500; }}
tr:nth-child(even) {{ background: #f7f9fa; }}
.card {{ display: inline-block; background: white; padding: 12px 20px; margin: 5px;
         border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,0.08); min-width: 130px; text-align: center; }}
.card .val {{ font-size: 1.4em; font-weight: 700; color: #2980b9; }}
.card .lbl {{ font-size: 0.8em; color: #7f8c8d; margin-top: 2px; }}
.methods {{ background: #f0f4f7; padding: 12px 16px; border-left: 4px solid #2980b9;
            margin: 10px 0; font-size: 0.88em; border-radius: 0 6px 6px 0; }}
.interpretation {{ background: #fdf6ec; padding: 12px 16px; border-left: 4px solid #e67e22;
                   margin: 10px 0; font-size: 0.88em; border-radius: 0 6px 6px 0; }}
.highlight {{ background: #eafaf1; border-left: 4px solid #27ae60; padding: 10px 16px;
              margin: 12px 0; border-radius: 0 6px 6px 0; font-size: 0.92em; }}
.fig {{ text-align: center; margin: 16px 0; }}
.fig img {{ max-width: 100%; border-radius: 6px; box-shadow: 0 2px 6px rgba(0,0,0,0.1); }}
</style>
</head>
<body>

<h1>NanoPhage Pipeline Report</h1>
<p>Wagner Giovanni Guzman Mendez &mdash; {now.strftime('%B %Y')}</p>

<div class="highlight">
<strong>Result:</strong> Assembled a {genome_len:,} bp circular genome identified as
<em>{top_name}</em> ({float(top_hsps[0]['pident']):.1f}% identity) |
{lifestyle} ({lifestyle_pct_main:.0f}% confidence) |
{total_genes} predicted genes ({hallmarks} hallmark) |
tRNAs: {trna_text} |
AMR genes: {'none detected' if not amr_found else f'{amr_found} found'}
</div>

<h2>1. Read Quality Control</h2>
<div class="methods">
Raw Nanopore reads were assessed with NanoPlot v1.43. Reads shorter than {filt.get('min_length', '500 bp').split('=')[-1].strip()} or
below Q{filt.get('min_quality', 'Q10').split('Q')[-1].strip()} were removed using Chopper v0.9.
</div>

<div>
<div class="card"><div class="val">{raw.get('number_of_reads', '-')}</div><div class="lbl">Raw reads</div></div>
<div class="card"><div class="val">{float(raw.get('number_of_bases', 0))/1e6:.1f} Mb</div><div class="lbl">Total bases</div></div>
<div class="card"><div class="val">{raw.get('mean_read_length', '-')}</div><div class="lbl">Mean length</div></div>
<div class="card"><div class="val">Q{raw.get('mean_qual', '-')}</div><div class="lbl">Mean quality</div></div>
<div class="card"><div class="val">{raw.get('n50', '-')}</div><div class="lbl">Read N50</div></div>
</div>
<p><strong>After filtering:</strong> {filt.get('After', '-')}</p>

<div class="interpretation">
Reads below 500 bp were removed as they can not produce meaningful overlaps for assembling a phage genome in the few 10ths Kbs range, the Q10 threshold ensures at least 90% basecall accuracy. NanoPlot provides Nanopore specific quality metrics, and Chopper is a fast Rust based filter for length and quality trimming.
</div>

<h2>2. Genome Assembly</h2>
<div class="methods">
Filtered reads were assembled with Flye v2.9.5 using <code>--nano-raw</code> mode.
No genome size parameter was provided &mdash; Flye estimated it automatically from the read data.
</div>

<div>
<div class="card"><div class="val">{quast.get('# contigs', '-')}</div><div class="lbl">Contigs</div></div>
<div class="card"><div class="val">{genome_len:,}</div><div class="lbl">Total length (bp)</div></div>
<div class="card"><div class="val">{quast.get('GC (%)', '-')}%</div><div class="lbl">GC content</div></div>
<div class="card"><div class="val">{asm_info.get('circular', '-')}</div><div class="lbl">Circular (Flye)</div></div>
<div class="card"><div class="val">{asm_info.get('coverage', '-')}x</div><div class="lbl">Assembly coverage</div></div>
</div>

<div class="interpretation">
Flye was chosen for its efficient handling of high coverage long-read data and built-in circular genome detection. The <code>--nano-raw</code> mode was used because the observed overlap divergence (6.2%) exceeded the &lt;3% expected by <code>--nano-hq</code>, despite HAC basecalling. No genome size was specified, which keeps the pipeline genome-agnostic.
</div>

<h2>3. Assembly Quality</h2>
<div class="methods">
Reads were mapped back to the assembly with minimap2 v2.28 (map-ont preset) and coverage was computed with samtools v1.21.
Viral genome completeness was assessed with CheckV v1.0.3 using AAI-based estimation against a curated viral genome database.
</div>

<div>
<div class="card"><div class="val">{cov.get('meandepth', '-')}x</div><div class="lbl">Mean read depth</div></div>
<div class="card"><div class="val">{cov.get('coverage', '-')}%</div><div class="lbl">Breadth of coverage</div></div>
<div class="card"><div class="val">{cv.get('completeness', '-')}%</div><div class="lbl">CheckV completeness</div></div>
<div class="card"><div class="val">{cv.get('checkv_quality', '-')}</div><div class="lbl">Quality tier</div></div>
<div class="card"><div class="val">{cv.get('contamination', '-')}%</div><div class="lbl">Contamination</div></div>
</div>

{'<h3>Coverage depth across the genome</h3><div class="fig"><img src="data:image/png;base64,' + cov_plot_b64 + '" alt="Coverage plot"></div>' if cov_plot_b64 else ''}

<div class="interpretation">
The coverage plot verifies uniform read mapping across the genome. The bowl-shaped pattern is expected for a linearized circular genome: basically, reads spanning the breakpoint map to initial of defined k-mers, inflating edge coverage. This confirms circularity, consistent with Flye's assembly graph. CheckV estimated 96.95% completeness and the ~3% gap corresponds to the circular junction trimmed during linearization. No host genes and no contamination confirm a free phage, not a provirus. CheckV was used instead of BUSCO because viral genomes lack the conserved single-copy markers that BUSCO requires.
</div>

<h2>4. Contig Characterization</h2>
<div class="methods">
Viral identification and taxonomy: geNomad v1.8. Functional annotation: Pharokka v1.7.
Lifestyle prediction: BACPHLIP v0.9.6. tRNA detection: tRNAscan-SE v2.0.12.
AMR screening: ABRicate v1.0.1 (NCBI, CARD, VFDB). Closest relative: remote BLASTn vs NCBI nt.
</div>

<h3>Taxonomy</h3>
<p>{genomad_virus.get('taxonomy', 'Not classified')} (virus score: {genomad_virus.get('virus_score', '-')}, {genomad_virus.get('n_hallmarks', '-')} hallmark genes)</p>

<div class="interpretation">
geNomad provides a virus confidence score and hallmark gene counts. Hallmark genes are referred to viral proteins (terminase, capsid, tail proteins) whose presence strongly indicates a viral genome. Taxonomy resolves to Caudoviricetes but not to family or genus, which is common for phages lacking closely classified relatives in ICTV.
</div>

<h3>Closest relative (BLASTn vs NCBI nt)</h3>
<p>The closest match is <strong>{top_name}</strong> (accession {top_accession}).
The alignment covers the entire query genome in three segments, totaling {total_aligned:,} aligned bases:</p>

<table>
<tr><th>Query coordinates</th><th>Subject coordinates</th><th>Identity</th><th>Alignment length</th><th>E-value</th></tr>
{blast_top_html}
</table>

<div class="interpretation">
The genome aligns to the reference in three segments because Flye linearized the circular genome at a different position than the GenBank reference. The three blocks together tile the full genome without gaps. Near-perfect identity around (99&ndash;100%) across all segments confirms this is the same phage. Remote BLASTn was used for comprehensive and up-to-date taxonomic placement without requiring a local database.
</div>

<h3>Gene content</h3>
<div>
<div class="card"><div class="val">{total_genes}</div><div class="lbl">Total CDS</div></div>
<div class="card"><div class="val">{hallmarks}</div><div class="lbl">Hallmark genes</div></div>
<div class="card"><div class="val">{annotated}</div><div class="lbl">Annotated</div></div>
<div class="card"><div class="val">{total_genes - annotated}</div><div class="lbl">Hypothetical</div></div>
</div>

{'<h3>Linear gene map</h3><div class="fig"><img src="data:image/png;base64,' + gene_map_b64 + '" alt="Gene map"></div>' if gene_map_b64 else ''}

<div class="interpretation">
The linear gene map shows functional modules colored by category, complementing the circular Pharokka map below: the linear view highlights modular organization (packaging &rarr; structural &rarr; replication &rarr; modification), while the circular map provides detailed gene labels. The CRISPR-associated exonuclease (Csa1) in a lytic phage is notable and may function in host DNA degradation during infection.
</div>

<h3>Annotated genes</h3>
<table>
<tr><th>Gene</th><th>Position</th><th>Strand</th><th>Hallmark</th><th>Function</th></tr>
{gene_rows_html}
</table>

{'<h3>Genome map (Pharokka)</h3><div class="fig"><img src="data:image/png;base64,' + pharokka_b64 + '" alt="Pharokka genome map"></div>' if pharokka_b64 else ''}

<div class="interpretation">
Pharokka uses PHANOTATE for gene calling, which outperforms Prodigal on phage genomes due to its handling of overlapping genes and phage-specific features, and annotates against the PHROGs database via MMseqs2.
</div>

<h3>Lifestyle</h3>
<p>{lifestyle} &mdash; BACPHLIP confidence: {vir*100:.1f}% virulent, {temp*100:.1f}% temperate.</p>

<div class="interpretation">
BACPHLIP searches for 206 lysogeny-associated protein domains using HMMER. The strong virulent prediction ({vir*100:.1f}%) reflects the absence of integrases, CI repressors, and excisionases. The Mu-like and P22 NinX domains present in this genome are shared structural domains, not lysogeny indicators.
</div>

<h3>tRNAs</h3>
<p>{len(trnas)} tRNA(s) detected{': ' + trna_text if trnas else ''}.</p>

<div class="interpretation">
tRNAscan-SE was run in bacterial mode (<code>-B</code>) because phage tRNAs use the same structural features as bacterial tRNAs. {'The single tRNA-Trp(CCA) detected suggests codon usage optimization in the way that carrying its own tRNA-Trp provides a translational advantage during infection. A single tRNA is typical for phages of this genome size (~59 kb).'}
</div>

<h3>Antimicrobial resistance</h3>
<p>{'No AMR genes detected across NCBI, CARD, and VFDB databases.' if not amr_found else f'{amr_found} AMR gene(s) detected.'}</p>

<div class="interpretation">
ABRicate screens against three complementary databases (AMRFinderPlus, CARD, VFDB) to maximize sensitivity. No AMR genes are expected in a lytic phage, as resistance genes are typically transduced by temperate phages that integrate into host chromosomes.
</div>

<h2>5. Methods Summary</h2>
<table>
<tr><th>Step</th><th>Tool</th><th>Purpose</th></tr>
<tr><td>Read QC</td><td>NanoPlot 1.43</td><td>Nanopore-specific read quality assessment</td></tr>
<tr><td>Filtering</td><td>Chopper 0.9</td><td>Length and quality filtering</td></tr>
<tr><td>Assembly</td><td>Flye 2.9.5</td><td>Long-read assembly with circular genome detection</td></tr>
<tr><td>Assembly QC</td><td>QUAST 5.2</td><td>Assembly statistics (N50, GC, contiguity)</td></tr>
<tr><td>Coverage</td><td>minimap2 2.28 + samtools 1.21</td><td>Read mapping and depth analysis</td></tr>
<tr><td>Completeness</td><td>CheckV 1.0.3</td><td>Viral genome completeness (AAI-based)</td></tr>
<tr><td>Identification</td><td>geNomad 1.8</td><td>Neural network viral classification and marker gene analysis</td></tr>
<tr><td>Annotation</td><td>Pharokka 1.7</td><td>Phage-specific gene calling (PHANOTATE) and annotation (PHROGs)</td></tr>
<tr><td>Taxonomy</td><td>BLASTn 2.16 (remote)</td><td>Closest relative identification (NCBI nt)</td></tr>
<tr><td>Lifestyle</td><td>BACPHLIP 0.9.6</td><td>Lysogeny marker HMM profiling</td></tr>
<tr><td>tRNAs</td><td>tRNAscan-SE 2.0.12</td><td>tRNA gene detection (covariance models)</td></tr>
<tr><td>AMR</td><td>ABRicate 1.0.1</td><td>Resistance gene screening (NCBI, CARD, VFDB)</td></tr>
</table>

</body>
</html>"""
    return html


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--raw-stats', required=True)
    p.add_argument('--filter-stats', required=True)
    p.add_argument('--assembly-info', required=True)
    p.add_argument('--quast', required=True)
    p.add_argument('--coverage', required=True)
    p.add_argument('--depth', required=True)
    p.add_argument('--checkv', required=True)
    p.add_argument('--genomad-summary', required=True)
    p.add_argument('--genomad-genes', required=True)
    p.add_argument('--blast', required=True)
    p.add_argument('--bacphlip', required=True)
    p.add_argument('--trna', required=True)
    p.add_argument('--amr', required=True)
    p.add_argument('--pharokka-plot', default='')
    p.add_argument('-o', '--output', default='nanophage_report.html')
    args = p.parse_args()

    html = build_report(args)
    with open(args.output, 'w') as f:
        f.write(html)
    print(f"Report: {args.output}")
