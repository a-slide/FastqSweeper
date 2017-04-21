# FastqSweeper  0.1

**Simple python 2.7 pipeline to clean fastq files with cutadapt and bwa substractive mapping**

---

**Creation : 2015/03/08**

**Last update : 2017/04/21**

---

## Principle

FastqSweeper use Cutadapt to trim reads based on quality, presence of N, and adapters and then perform a step of subtractive mapping with bwa Mem. BWA stream is parsed directly to separate the properly mapped reads that are written in a bam file, from the non-properly mapped reads that are written in a fastq file ready to be realign for further analysis.

A report containing the options used and reads counts depending of their category is also generated.

## Dependencies

* Python 2.7
* Cutadapt
* BWA-Mem
* pysam (python package including samtools)

## Usage

Usage: FastqSweeper.py -i INDEX -1 FASTQ_R1 [-2 FASTQ_R2][...]

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INDEX, --index=INDEX

    Path to the bwa index basename
  -1 FASTQ_R1, --fastq_R1=FASTQ_R1
     Path to the fastq file
  -2 FASTQ_R2, --fastq_R2=FASTQ_R2
    [facultative] Path to the pair fastq file if paired end
  -r, --run
    [facultative] Run command lines else the sotware run in demo mode (default: False)
  -t THREAD, --thread=THREAD
    [facultative] Number of thread to use (default: 1)
  -b BWA_OPT, --bwa_opt=BWA_OPT
    [facultative] bwa options for the mapping step (facultative and quoted) (default: -M)
  -c CUTADAPT_OPT, --cutadapt_opt=CUTADAPT_OPT
    [facultative] cutadapt options for the qc step (facultative and quoted) (default: -m 25 -q 30,30 --trim-n)
  -a ADAPTER, --adapter=ADAPTER
    [facultative] Path to a fasta file containing adapters to be 3' trimmed
  -m MIN_MAPQ, --min_mapq=MIN_MAPQ
    [facultative] Minimal mapq quality to be considered mapped (default: 0)
  -s MIN_MATCH_SIZE, --min_match_size=MIN_MATCH_SIZE
    [facultative] Minimal size of match to be considered mapped (default: 0)
  --skip_cutadapt
    [facultative] Skip the cutadapt step (default: False)

## Authors and Contact

Adrien Leger - 2016

Enright's lab

EMBL EBI

* <aleg@ebi.ac.uk>
* [Github](https://github.com/a-slide)
* [Website](https://a-slide.github.io/)