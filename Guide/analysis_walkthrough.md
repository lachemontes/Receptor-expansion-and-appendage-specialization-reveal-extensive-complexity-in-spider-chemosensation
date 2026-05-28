# Step-by-Step Analysis Guide — Spider Chemosensation (*Argiope bruennichi*)

### By Zaide Montes-Ortiz

## Table of Contents

1. [Data](#1-data)
   - 1.1 [Genome](#11-genome)
   - 1.2 [RNAseq](#12-rnaseq)
2. [Transcriptome Assembly](#2-transcriptome-assembly)
   - 2.1 [Quality Control — FastQC &amp; MultiQC](#21-quality-control--fastqc--multiqc)
   - 2.2 [Adapter Trimming — TrimGalore &amp; Cutadapt](#22-adapter-trimming--trimgalore--cutadapt)
   - 2.3 [Taxonomic classification — Kraken2](#23-taxonomic-classification--kraken2)
   - 2.4 [Genome-guided Assembly — HISAT2 → SAM/BAM → StringTie2](#24-genome-guided-assembly--hisat2--sambam--stringtie2)
   - 2.5 [De novo Assembly — Trinity → TransDecoder → CD-HIT → BUSCO](#25-de-novo-assembly--trinity--transdecoder--cd-hit--busco)
3. [CR Gene Annotation](#3-cr-gene-annotation)
   - 3.1 [InterProScan](#31-interproscan)
   - 3.2 [BLAST searches](#32-blast-searches)
   - 3.3 [CLANS clustering](#33-clans-clustering)
4. [Chromosome Mapping](#4-chromosome-mapping)
   - 4.1 [Miniprot](#41-miniprot)
   - 4.2 [MG2C visualization](#42-mg2c-visualization)
5. [Phylogenetic Analysis](#5-phylogenetic-analysis)
   - 5.1 [Sequence extraction with Seqtk](#51-sequence-extraction-with-seqtk)
   - 5.2 [Alignment — MAFFT &amp; TrimAl](#52-alignment--mafft--trimal)
   - 5.3 [Tree inference — IQ-TREE2](#53-tree-inference--iq-tree2)
   - 5.4 [Tree visualization — iTOL](#54-tree-visualization--itol)
6. [Expression Analysis](#6-expression-analysis)
   - 6.1 [Kallisto — Index &amp; Quantification](#61-kallisto--index--quantification)
   - 6.2 [Abundance matrix](#62-abundance-matrix)
   - 6.3 [Differential expression — DESeq2](#63-differential-expression--deseq2)

---

## 1. Data

### 1.1 Genome

We used the chromosome-level genome assembly of *A. bruennichi* — the first of its kind in the order Araneae. The assembly was generated using 21.8× PacBio sequencing, polished with 19.8× Illumina paired-end data, and scaffolded with Hi-C proximity ligation, resulting in an N50 scaffold size of 124 Mb and 98.4% complete arthropod BUSCOs.

| GenBank         | Name         | Level      | Date     |
| --------------- | ------------ | ---------- | -------- |
| GCA_015342795.1 | ASM1534279v1 | Chromosome | Nov 2020 |
| GCA_947563725.1 | qqArgBrue1.1 | Chromosome | 2021     |

> 📌 More info: https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_015342795.1/
> 📌 Taxonomy: https://www.ncbi.nlm.nih.gov/data-hub/taxonomy/tree/?taxon=94029

### 1.2 RNAseq

Six male and six female adult *A. bruennichi* were collected from Kargower Weg, Waren, Germany (July 2022). Males were exposed to females prior to dissection to boost pheromone receptor gene expression. Four tissue types were dissected independently for each individual: the **distal** and **proximal segments of the first leg**, **pedipalps**, and **mouthparts**. Each sample was sequenced with four technical replicates and two paired-end libraries, yielding **48 libraries per tissue per sex**. In total, **384 libraries** were generated across all tissues and conditions. Libraries were prepared with NEBNext Ultra II Directional RNA Library Prep Kit and sequenced on NextSeq High-Throughput (Illumina).

| **Tissue**           | **Condition** | **Libraries** | **Sequence data yield**         |
| -------------------------- | ------------------- | ------------------- | ------------------------------------- |
| **1st leg distal**   | **Male**      | **48**        | **—**                          |
| **1st leg proximal** | **Male**      | **48**        | **—**                          |
| **1st leg distal**   | **Female**    | **48**        | **—**                          |
| **1st leg proximal** | **Female**    | **48**        | **192.53 GB total (legs)**      |
| **Pedipalps**        | **Male**      | **48**        | **—**                          |
| **Pedipalps**        | **Female**    | **48**        | **38.08 GB total (pedipalps)**  |
| **Mouthparts**       | **Male**      | **48**        | **—**                          |
| **Mouthparts**       | **Female**    | **48**        | **62.19 GB total (mouthparts)** |

**The distal and proximal parts of the first legs were dissected separately, as they are thought to be involved in contact chemosensation and olfaction respectively **[Talukder et al., 2025](https://www.pnas.org/doi/10.1073/pnas.2415468121)**, and may therefore express CR genes differentially.**

## 2. Transcriptome Assembly

Two strategies were used in parallel to maximize transcript recovery: a **genome-guided assembly** with StringTie2 and a **de novo assembly** with Trinity. Both outputs were subsequently used as complementary resources for CR gene identification.

---

### 2.1 Quality Control — FastQC & MultiQC

Quality control was performed on raw reads using FastQC v0.11.9 and MultiQC v1.12. Jobs were submitted as SLURM arrays to process all libraries in parallel.

**Legs (Uppmax):**

```bash
#!/bin/bash
#SBATCH -A snic2022-5-454
#SBATCH -p core -n 8
#SBATCH -t 12:00:00
#SBATCH -J Fastqc
#SBATCH --array=1-96
#SBATCH --output=fastqc_%A_%a.out
#SBATCH --error=fastqc_%A_%a.err

module load bioinfo-tools
module load FastQC/0.11.9
module load MultiQC/1.12

input_file=$(sed -n "$SLURM_ARRAY_TASK_ID p" List_fastq_files.txt)
fastqc $input_file
multiqc .
```

**Mouthparts & Pedipalps (Dardel):**

```bash
#!/bin/bash
#SBATCH -A naiss2024-5-647
#SBATCH -p shared
#SBATCH --cpus-per-task=6
#SBATCH -t 12:00:00
#SBATCH -J Fastqc
#SBATCH --array=1-104
#SBATCH --output=/cfs/klemming/projects/supr/naiss2023-23-109/Analysis/quality/logs/fastqc_%A_%a.out
#SBATCH --error=/cfs/klemming/projects/supr/naiss2023-23-109/Analysis/quality/logs/fastqc_%A_%a.err

input_dir="/cfs/klemming/projects/supr/naiss2023-23-109/Data/Fastq_mouthparts"
output_dir="/cfs/klemming/projects/supr/naiss2023-23-109/Analysis/quality"
list_file="/cfs/klemming/projects/supr/naiss2023-23-109/Analysis/List_files/fastqc_mouth.txt"

mkdir -p $output_dir/logs
input_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $list_file)

if [[ -f "$input_dir/$input_file" ]]; then
    fastqc -o $output_dir "$input_dir/$input_file"
fi

if ls $output_dir/*.zip 1> /dev/null 2>&1; then
    multiqc $output_dir
fi
```

---

### 2.2 Adapter Trimming — TrimGalore & Cutadapt

Low-quality regions and Illumina adapters were removed with TrimGalore v0.6.1. For mouthparts and pedipalps, an additional Cutadapt step was run to remove poly-A and poly-G tails — a second FastQC/MultiQC round was done afterwards to confirm their removal.

**Legs (Uppmax):**

```bash
#!/bin/bash
#SBATCH -A snic2022-5-454
#SBATCH -p node
#SBATCH -t 48:00:00
#SBATCH -J Trimgalore
#SBATCH --array=1-96

module load bioinfo-tools
module load TrimGalore/0.6.1
module load FastQC/0.11.9
module load cutadapt/2.1

input_dir="/proj/naiss2023-23-109/Spiders_Project/Data/RNAseq/Maleleg/.../Fastq"
output_dir="/proj/naiss2023-23-109/Spiders_Project/Data/RNAseq/Maleleg/.../Fastq/Trimgalore3"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID} p" List.txt)
r1="${input_dir}/${sample}_R1_001.fastq.gz"
r2="${input_dir}/${sample}_R2_001.fastq.gz"

trim_galore -j 4 --illumina --paired "${r1}" "${r2}" --fastqc -o "${output_dir}"
```

**Mouthparts (Dardel) — TrimGalore + Cutadapt for poly-A/G removal:**

```bash
#!/bin/bash
#SBATCH -A naiss2024-5-647
#SBATCH -p shared
#SBATCH -c 6
#SBATCH -t 72:00:00
#SBATCH -J TrimGalore_mouthp
#SBATCH --array=1-104

input_dir="/cfs/klemming/projects/supr/naiss2023-23-109/Data/Fastq_mouthparts"
output_dir="/cfs/klemming/projects/supr/naiss2023-23-109/Analysis/trimgalore/mouthparts"
list_file="/cfs/klemming/projects/supr/naiss2023-23-109/Analysis/List_files/fastqc_mouth.txt"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $list_file)
r1="${input_dir}/${sample}_R1_001.fastq.gz"
r2="${input_dir}/${sample}_R2_001.fastq.gz"

trim_galore -j 6 --illumina --paired "$r1" "$r2" --fastqc -o "$output_dir"
```

```bash
# Cutadapt — poly-G and poly-A removal
cutadapt -j 6 -a "G{100}" -a "A{100}" -A "G{100}" -A "A{100}" \
  -o "${r1_output}" -p "${r2_output}" "${r1_input}" "${r2_input}"
```

---

### 2.3 Taxonomic classification — Kraken2

Trimmed reads were classified with Kraken2 to detect contamination.

```bash
#!/bin/bash
#SBATCH -A naiss2024-5-647
#SBATCH -p shared
#SBATCH -c 16
#SBATCH --mem=60GB
#SBATCH -t 7-00:00:00
#SBATCH -J kraken2_mouth
#SBATCH --array=1-96

source ~/miniconda3/etc/profile.d/conda.sh
conda activate Quality

input_dir="/cfs/klemming/projects/supr/naiss2023-23-109/Analysis/trimgalore/mouthparts"
output_dir="/cfs/klemming/projects/supr/naiss2023-23-109/Analysis/kraken2/mouthparts"
db_kraken="/cfs/klemming/projects/databases/kraken2_standard"
list_file="/cfs/klemming/projects/supr/naiss2023-23-109/Analysis/List_files/kraken2/mouthparts.txt"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$list_file")
r1="${input_dir}/${sample}_R1_001_val_1.fq.gz"
r2="${input_dir}/${sample}_R2_001_val_2.fq.gz"

kraken2 --db "$db_kraken" --threads 16 --paired \
    --report "${output_dir}/${sample}_kraken_report.txt" \
    --output "${output_dir}/${sample}_kraken_output.txt" "$r1" "$r2"
```

---

### 2.4 Genome-guided Assembly — HISAT2 → SAM/BAM → StringTie2

Trimmed reads were aligned to the reference genome (GCA_947563725.1) with HISAT2 v2.2.1. SAM files were converted and sorted with Samtools v1.9, then assembled into transcripts with StringTie2 v2.1.4.

**HISAT2 mapping:**

```bash
#!/bin/bash
#SBATCH -A snic2022-5-454
#SBATCH -p node
#SBATCH -t 12:00:00
#SBATCH -J Fem_HISAT_Mapping
#SBATCH --array=1-96

module load bioinfo-tools HISAT2/2.2.1

input_dir="/proj/naiss2023-23-109/Spiders_Project/Data/RNAseq/Femaleleg/.../Trimgalore3"
output_dir="/proj/naiss2023-23-109/Spiders_Project/Analysis/HISAT/Out2_Female_Leg"
genome_index="/proj/snic2022-23-541/Spider_project/Analysis/HISAT/Genome_Index"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID} p" List.txt)
r1="${input_dir}/${sample}_R1_001_val_1.fq.gz"
r2="${input_dir}/${sample}_R2_001_val_2.fq.gz"

hisat2 --dta -p 8 -x "$genome_index" -1 "${r1}" -2 "${r2}" -S "${output_dir}/${sample}.sam"

# Merge .err files to extract alignment rates
tail -n +1 *.err > concatenated_filename_output.txt
```

**Convert HISAT2 .err files to Excel (Python):**

```python
import pandas as pd

with open('concatenated_Ticks_output.txt') as f:
    data = f.read().split('== ')

data = data[1:]
results = []
for section in data:
    lines = section.strip().split('\n')
    filename = lines[0].strip()
    reads = lines[1].split('reads')[0].strip()
    overall_alignment_rate = lines[-1].split()[0].strip()
    results.append((filename, reads, overall_alignment_rate))

df = pd.DataFrame(results, columns=['Filename', 'Reads', 'Overall Alignment Rate'])
df.to_excel('output_Ticks.xlsx', index=False)
```

**SAM → sorted BAM:**

```bash
module load bioinfo-tools
module load samtools/1.9

for sam_file in *.sam; do
    bam_file=${sam_file/.sam/.bam}
    samtools view -Sb $sam_file > $bam_file
    sorted_bam=${bam_file/.bam/_sorted.bam}
    samtools sort $bam_file > $sorted_bam
    samtools index $sorted_bam
done

samtools merge merged.bam *_sorted.bam
samtools index merged.bam
```

**StringTie2 assembly:**

```bash
stringtie merged.bam -p 4 -o transcripts.gtf -A gene_abundances.txt
```

**GTF → FASTA:**

```bash
# Convert GTF to GFF3
gffread transcripts.gtf -o transcripts.gff3

# Extract transcript sequences
gffread transcripts.gff3 -g genome.fasta -w transcripts.fasta
```

---

### 2.5 De novo Assembly — Trinity → TransDecoder → CD-HIT → BUSCO

Trinity v2.9.1 was run independently for each appendage. Assembled transcriptomes were processed to predict ORFs (TransDecoder), remove redundancy (CD-HIT), and assess completeness (BUSCO v5.2.2 with `arachnida_odb10`). Trinity assemblies showed higher completeness (95.9–98.7% BUSCOs) than genome-guided assemblies (74.9–78.7%), so both were retained as complementary resources.

**Concatenate Trinity outputs and tag headers:**

```bash
# Concatenate all Trinity fasta files
cat Trinity_*/Trinity.fasta* > Trinity_filename.fasta

# Add library name to fasta header
awk '/^>/ {print $1 "_S413_S2 " $2, $3;next}1' Trinity_fastafile.fasta > Trinity_file_name_fix.fasta
```

**TransDecoder — predict ORFs:**

```bash
#!/bin/bash
#SBATCH -A snic2022-5-454
#SBATCH -p core -n 12
#SBATCH -t 00:00:00
#SBATCH -J Transdecoder

module load bioinfo-tools
module load TransDecoder/5.7.0

TransDecoder.LongOrfs -t transcripts.fasta
TransDecoder.Predict -t transcripts.fasta
```

**Clean headers and remove stop codons (`*`):**

```bash
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' \
  < Trinity_merge.fasta.transdecoder.pep | sed "1d" > transdecoder_No_redundancy.fasta

awk '/^>/ {print $1 " " ;next}1' transdecoder_No_redundancy.fasta | \
  sed 's/*//g' > transdecoder_No_redundancy_fix.fasta
```

**Optional — filter sequences by length before CD-HIT:**

```python
def filter_sequences(input_file, output_file, max_length):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        seq_id = ''
        sequence = ''
        for line in f_in:
            if line.startswith('>'):
                if sequence and len(sequence) <= max_length:
                    f_out.write(seq_id + sequence + '\n')
                seq_id = line
                sequence = ''
            else:
                sequence += line.strip()
        if sequence and len(sequence) <= max_length:
            f_out.write(seq_id + sequence + '\n')

filter_sequences('transdecoder_singleLine_fix.fasta', 'filtered_sequences.fasta', max_length=2000)
```

**CD-HIT — remove redundancy at 95% identity:**

```bash
#!/bin/sh
#SBATCH -A snic2022-5-454
#SBATCH -p core -n 16
#SBATCH -t 72:00:00
#SBATCH -J CD-HIT

cd-hit -i filtered_sequences.fasta \
       -o CD_HIT_Spider_Transcriptome_deNovo_95.fasta \
       -c 0.95 -n 5 -M 0 -T 0
```

**BUSCO — assembly completeness:**

```bash
#!/bin/sh
#SBATCH -A snic2022-5-454
#SBATCH -p core -n 16
#SBATCH -t 3-00:00:00
#SBATCH -J BUSCO_arachnida_odb10

module load augustus/3.4.0
source $AUGUSTUS_CONFIG_COPY
export AUGUSTUS_CONFIG_PATH=/proj/naiss2023-23-109/Spiders_Project/Analysis/BUSCO/augustus_config

busco -i CD_HIT_Spider_Transcriptome_deNovo_95.fasta \
      -l arachnida_odb10 -o busco_output -m proteins
```

---

## 3. CR Gene Annotation

### 3.1 InterProScan

InterProScan v5 was used to assign functional domains to all candidate proteins. For IRs/iGluRs, we required at least one of: `PF00060`, `PF10613`, or `PF01094`. For GRs, we looked for `PF08395` or descriptors such as "gustatory receptor" and "taste receptor". Databases used: Gene3D, ProSitePatterns, PANTHER, CDD, Pfam, Phobius, SUPERFAMILY, TMHMM.

```bash
#!/bin/bash
#SBATCH -A snic2022-5-454
#SBATCH -p core -n 12
#SBATCH -t 2-00:00:00
#SBATCH -J InterPro_DB

module load bioinfo-tools InterProScan/5.52-86.0

interproscan.sh -i CD_HIT_Spider_Transcriptome.fasta \
  -goterms -t p -f TSV,XML \
  -appl Gene3D,ProSitePatterns,PANTHER,CDD,Pfam,Phobius,SUPERFAMILY,TMHMM \
  -b output_InterPro
```

**Extract hits matching chemosensory terms (Python, run in Jupyter):**

```python
# Step 1 — extract lines matching InterPro hit list
with open('Hit_List_uniq.txt', 'r') as f1:
    items = set(f1.read().splitlines())

with open('gene_abundances.txt', 'r') as f2, open('hits.txt', 'w') as f3:
    header = f2.readline()
    f3.write(header)
    for line in f2:
        if line.startswith(tuple(items)):
            f3.write(line)

# Step 2 — get unique chromosome names from hits
with open('hits.txt', 'r') as f4:
    refs = set()
    next(f4)
    for line in f4:
        cols = line.split()
        refs.add(cols[2])
    uniq_refs = sorted(list(refs))

print(uniq_refs)
```

**Extract chromosomes from genome with Seqtk:**

```bash
module load bioinfo-tools
module load seqtk/1.2-r101

seqtk subseq GCA_015342795.1_genomic.fna Chromosme_List.txt > Chromosome_subset.fasta

# Verify number of extracted chromosomes
grep ">" Chromosome_subset.fasta | wc
```

---

### 3.2 BLAST searches

A custom BLAST database was built from the curated chemosensory receptor gene set and used in reciprocal searches. Multiple strategies were applied: `blastp` against the transcriptome, `blastn` against the genome, `tblastn` for protein-to-nucleotide searches, and `blastp` against the *D. melanogaster* proteome and NCBI nr database for GRs with no initial hits.

**Build BLAST database:**

```bash
makeblastdb -in CD_HIT_Spider_Transcriptome.fasta \
  -title Spiders -dbtype prot -out Spiders_DB -parse_seqids
```

**BLASTp — IRs/iGluRs:**

```bash
#!/bin/sh
#SBATCH -A snic2022-5-454
#SBATCH -p core -n 12
#SBATCH -t 46:00:00
#SBATCH -J Blastp_IR

module load bioinfo-tools blast/2.9.0+

blastp -query 329_IRs_iGluRs_for_spider_blast.fasta \
  -db Spiders_DB -num_threads 4 -evalue 1e-5 -outfmt 6 -max_target_seqs 5
```

**BLASTn — confirm identity against genome:**

```bash
makeblastdb -in GCA_015342795.1_genomic.fna -title Spider_DB \
  -dbtype nucl -out Spider_DB -parse_seqids

blastn -query IR_Predicted_Spiders_CDS.fasta \
  -db Spider_DB -num_threads 12 -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident nident length qstart qend sstart send evalue mismatch gapopen" \
  -max_target_seqs 3
```

**tBLASTn:**

```bash
tblastn -query IR_StringTie_Trinity_Candidates.fasta \
  -db Spider_DB.nhr -num_threads 12 -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident nident length qstart qend sstart send evalue mismatch gapopen" \
  -max_target_seqs 5

# Filter results by identity > 80%
awk '$3>80.0{print $0}' blast_output.txt > tBLASTn_IR_80_Identity.txt
```

**BLASTp GRs with no hits — against Dmel proteome:**

```bash
makeblastdb -in Dmel_proteome.fasta -title Dmel -dbtype prot -out Dmel_DB -parse_seqids

blastp -query GR_Predicted_Spiders_No_Blast_Hits_Fix.fasta \
  -db Dmel_DB -num_threads 12 -evalue 1e-5 -outfmt 6 -max_target_seqs 5
```

**BLASTp GRs — against NCBI nr (Insecta):**

```bash
blastp -query GR_Predicted_Spiders_No_Blast_Hits.fasta \
  -db nr -out blastp_results.txt -outfmt 6 -taxids 50557
```

**Diagnostic BLAST — check if GRs are in StringTie but unannotated:**

```bash
# If hits >90% identity, >80% coverage: GRs are in StringTie under STRG.XXXX names
# → update annotation, no need to add sequences
# If no hits: StringTie genuinely did not assemble them
# → concatenate Trinity GRs + StringTie CDS → CD-HIT → re-quantify with Kallisto

blastn -query trinity_GR_cds.fasta \
       -subject stringtie2_cds.fasta \
       -outfmt "6 qseqid sseqid pident length qlen slen" \
       -evalue 1e-10 > GR_diagnostic.txt
```

---

### 3.3 CLANS clustering

An all-to-all BLASTp network was built with CLANS (E-value cutoff 0.01) to assess whether putative *A. bruennichi* GRs cluster with *Drosophila* ORs, GRs, and PHTF proteins. This confirmed that spider GR candidates group close to the *Drosophila* GR cluster, supporting their annotation.

---

## 4. Chromosome Mapping

### 4.1 Miniprot

Candidate protein sequences were mapped onto the reference genome (GCA_947563725.1) with Miniprot to obtain chromosomal coordinates. This was done independently for IRs/iGluRs and GRs.

**Get chromosome lengths:**

```bash
samtools faidx GCA_947563725.1_genomic.fna
cut -f1,2 GCA_947563725.1_genomic.fna.fai > chromosome_lengths.txt
```

**Extract candidate sequences with Seqtk:**

```bash
conda activate Seqkt

# IRs
seqtk subseq CD_HIT_Spider_Transcriptome.fasta IR_Chromosome_map.txt > \
  IR_Predicted_Spiders_Chromosome_Map.fasta

# GRs — extract the 84 expressed GRs
seqtk subseq GR_Predicted_Spiders_AA.fasta GR_List_81.txt > GR_Predicted_Spiders_AA_82.fasta

# Clean headers
awk '{print $1;next}1' GR_Predicted_Spiders_AA_82.fasta | sed 's/*//g' > \
  GR_Predicted_Spiders_AA_82_Fix.fasta
```

**Run Miniprot:**

```bash
# IRs
./miniprot -t16 GCA_947563725.1_genomic.fna IR_Candidates_transcripts_fix.fasta \
  --gff-only > output_miniprot_IR

# GRs
./miniprot -t16 GCA_947563725.1_genomic.fna GR_Predicted_Spiders_AA_82_Fix.fasta \
  --gff > output_miniprot_GR
```

---

### 4.2 MG2C visualization

Miniprot GFF output and chromosome lengths were used as input to [MG2C](http://mg2c.iask.in/mg2c_v2.0/) for visualization. IRs/iGluRs are shown in purple and GRs in red. Sex chromosomes (X1, X2) are indicated in red. Color maps were prepared in Excel.

> Reference: Chao et al. (2021) *Molecular Horticulture*, https://doi.org/10.1186/s43897-021-00020-x

---

## 5. Phylogenetic Analysis

### 5.1 Sequence extraction with Seqtk

Candidate sequences from StringTie2 and Trinity were extracted and concatenated for phylogenetic dataset construction.

```bash
# Extract IRs from de novo assembly
seqtk subseq CD_HIT_Spider_Transcriptome_deNovo_95.fasta IR_List_2024.txt > \
  IR_Candidates_2024_Trinity.fasta

# Concatenate StringTie + Trinity candidates
cat IR_Predicted_Spiders_111.fasta IR_Trinity_consensus_10.fasta > ConcatenateIR_data.fasta

# Extract final sequence set
seqtk subseq ConcatenateIR_data.fasta List_to_Extract.txt > IR_Trinity_Stringtie_49.fasta

# Remove duplicates and short sequences (<200 aa)
cd-hit -i Spiders_Phylogeny_2024.fasta -o CD_1.0_Spiders_Phylogeny_2024.fasta \
  -c 1.0 -n 5 -l 200 -M 0 -T 0

# Final concatenation
cat CD_1.0_Spiders_Phylogeny_2024.fasta IR_Trinity_Stringtie_49.fasta > Spiders_phylogeny.fasta
```

---

### 5.2 Alignment for iGluR and IRs — MAFFT & TrimAl

Multiple sequence alignment was performed with MAFFT v7. Gapped regions at the alignment ends were edited with Jalview. TrimAl was used to remove poorly aligned columns.

```bash
# MAFFT alignment
mafft --auto Spiders_phylogeny.fasta > Alignment_Spiders_phylogeny.fasta

# TrimAl — remove gapped columns
trimal -in Alignment_Spiders_phylogeny \
       -out Alignment_Spiders_phylogeny_Trimal \
       -st 0 -gt 0.9 -cons 60
```

> ⚠️ After TrimAl, verify that key glutamate-binding residues are preserved:
>
> - S1 region: **R** at position ~551
> - S2 first half: **S or T** at position ~742
> - S2 second half: **D or E** at position ~787

---

### 5.3 Tree inference — IQ-TREE2

Maximum-likelihood trees were inferred with IQ-TREE2, using ModelFinder to determine the best substitution model. Bootstrap support was assessed with 1000 ultrafast bootstrap replicates.

```bash
#!/bin/sh
#SBATCH -A naiss2023-5-461
#SBATCH -p core -n 16
#SBATCH -t 7-00:00:00
#SBATCH -J Spider_Phylo

# conda activate IQtree

# Determine best-fit model and run tree
iqtree -s Alignment_Spiders_phylogeny_Trimal -m MFP -bb 1000
```

**Best-fit models obtained:**

| Analysis                 | BIC model     |
| ------------------------ | ------------- |
| IR/iGluR alignment 1 & 2 | Q.pfam+I+R7   |
| IR/iGluR alignment 3     | Q.pfam+F+I+R8 |

---

### 5.4 Tree visualization — iTOL

Trees were visualized and annotated with [iTOL v6](https://itol.embl.de/). Colored outer bands indicate receptor subfamilies. Black dots indicate antennal IRs; red dots indicate differentially expressed receptors.

> Reference: Letunic & Bork (2021) *Nucleic Acids Research*, https://doi.org/10.1093/nar/gkab301

---

## 6. Expression Analysis

### 6.1 Kallisto — Index & Quantification

Transcript abundance was quantified with Kallisto v0.48.0. The index was built from the deduplicated CDS file generated by TransDecoder and CD-HIT. All 192 libraries were quantified in a SLURM array job.

**Build index:**

```bash
kallisto index -i transcripts_Single_CDS.idx transcripts_Single_CDS.fasta
```

**Quantification — Legs (Uppmax):**

```bash
#!/bin/bash
#SBATCH -A snic2022-5-454
#SBATCH -p node
#SBATCH -t 5-00:00:00
#SBATCH -J Kallisto_spiders
#SBATCH --array=1-192

module load bioinfo-tools kallisto/0.48.0

input_dir="/proj/naiss2023-23-109/Spiders_Project/Data/RNAseq/Trimmed_Reads_All"
output_dir="/proj/naiss2023-23-109/Spiders_Project/Analysis/Kallisto/Expression_dir_Single_CDS_Libs"
index="/proj/naiss2023-23-109/Spiders_Project/Analysis/Kallisto/transcripts_Single_CDS.idx"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID} p" List.txt)
R1="${input_dir}/${sample}_R1_001_val_1.fq.gz"
R2="${input_dir}/${sample}_R2_001_val_2.fq.gz"

kallisto quant -i "$index" -o "${output_dir}/${sample}" \
  --threads 12 --fr-stranded "$R1" "$R2"
```

**Quantification — Mouthparts (Dardel):**

```bash
#!/bin/bash
#SBATCH -A lu2025-7-42
#SBATCH --cpus-per-task=12
#SBATCH --mem=48GB
#SBATCH -t 1-00:00:00
#SBATCH -J mKallistoQuant
#SBATCH --array=1-96

input_dir="/home/zaidemo/lu2025-17-17/Spider/data/trimmomatic/mouthparts"
output_dir="/home/zaidemo/lu2025-17-17/Spider/analysis/kallisto/mouthparts"
index="/home/zaidemo/lu2025-17-17/Spider/analysis/kallisto/mouthparts/index"
list="/home/zaidemo/lu2025-17-17/Spider/analysis/kallisto/mouthparts/list.txt"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$list")
R1="${input_dir}/${sample}_R1_paired.fastq.gz"
R2="${input_dir}/${sample}_R2_paired.fastq.gz"

mkdir -p "${output_dir}/${sample}"
kallisto quant -i "$index" -o "${output_dir}/${sample}" \
  --threads 12 --rf-stranded "$R1" "$R2"
```

**Rename abundance files with directory name:**

```bash
for dir in */; do
    dir_name=$(basename "$dir")
    if [ -f "$dir/abundance.tsv" ]; then
        mv "$dir/abundance.tsv" "$dir/abundance_${dir_name}.tsv"
    fi
done
```

---

### 6.2 Abundance matrix

TPM values from all Kallisto outputs were aggregated into a single matrix for downstream analysis.

```python
import os

input_dir = os.getcwd()
input_files = [f for f in os.listdir(input_dir) if f.endswith(".tsv")]
input_file_paths = [os.path.join(input_dir, f) for f in input_files]

tpm_dict = {}
for input_file_path in input_file_paths:
    with open(input_file_path, "r") as input_file:
        next(input_file)  # skip header
        for line in input_file:
            values = line.strip().split("\t")
            if len(values) < 5:
                continue
            target_id = values[0]
            tpm = values[4]
            if target_id not in tpm_dict:
                tpm_dict[target_id] = []
            tpm_dict[target_id].append(tpm)

with open("tpm_abundance_all_samples.tsv", "w") as output_file:
    output_file.write("target_id\t" + "\t".join(input_files) + "\n")
    for target_id, tpm_list in tpm_dict.items():
        output_file.write(target_id + "\t" + "\t".join(tpm_list) + "\n")
```

> A custom script `KallistoToDeseq2.py` was also used to aggregate estimated counts (`est_counts`) into a single file per transcript across all samples, for input into DESeq2.

---

### 6.3 Differential expression — DESeq2

Differential expression analysis was conducted in R using DESeq2 v1.42.1. The **proximal segment of the male leg** was set as the reference level, based on the hypothesis that the abundance of wall-pore sensilla in that region drives expression of olfactory chemosensory genes.

Four key contrasts were tested:

1. Male distal vs. male proximal
2. Female distal vs. female proximal
3. Male distal vs. female distal
4. Male proximal vs. female proximal

**Expression level categories (averaged TPM):**

| Category              | TPM range |
| --------------------- | --------- |
| Very low              | < 1       |
| Low                   | 1 – 2    |
| Moderate              | 2 – 5    |
| Abundant              | 5 – 20   |
| Highly abundant       | 20 – 50  |
| Very highly expressed | > 50      |

Results include log2 fold changes, p-values, and adjusted p-values. Significant genes were visualized as z-score heatmaps and volcano plots.

---

*Last updated: 2026 — Montes-Ortiz et al.*
