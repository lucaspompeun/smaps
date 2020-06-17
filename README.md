<p align="center">
    <img src="https://img.shields.io/github/issues/lucaspompeun/smaps" />
    <img src="https://img.shields.io/github/stars/lucaspompeun/smaps" />
    <img src="https://img.shields.io/github/forks/lucaspompeun/smaps" />
    <img src="https://img.shields.io/github/license/lucaspompeun/smaps" /
</p>

# Smaps
A pipeline based on unmapped reads, which combines different tools for extends contigs and closing gaps present in the genomes.

## Requirements
- Linux (64-bit and 32-bit with slightly limited functionality).

- [Python3](https://www.python.org/downloads/) or higher.

- [Java](https://www.java.com/download/) version 1.8.0_201 or higher.

- [Quast](http://bioinf.spbau.ru/quast) in `PATH` variable.

- [Perl](https://www.perl.org/get.html) language installed.

- [Prokka](https://github.com/tseemann/prokka) in `PATH` variable

## Installation
You need to install git

```sh
$ sudo apt install git
```

Clone smaps repository to your machine

```
$ cd && git clone https://github.com/lucaspompeun/smaps.git
```

Creat a symbolic link to Smaps
```
$ sudo chmod 777 smaps/ && ln -s ~/smaps/smaps.py /usr/local/bin/smaps
```

With all dependecies satisfied you can simple run on terminal

```
$ smaps -h
usage: smaps.py [-h] -read1 READ1 [-read2 READ2] -output OUTPUT [-gff GFF] [-reference REFERENCE] [-sspace SSPACE] [-minreads MINREADS] [--version]

Smaps A tool to extends contigs to reduce gaps with unmapped reads

optional arguments:
  -h, --help            show this help message and exit
  -read1 READ1          first fastq file (required)
  -read2 READ2          second fastq file (optional)
  -output OUTPUT        output for results files (please, give the full path for the output)
  -gff GFF              gff file (optional)
  -reference REFERENCE  reference file (optional)
  -sspace SSPACE        number of runs of sspace software (optional, default=5)
  -minreads MINREADS    minimum number of unmapped reads to extend a base (optional, default=5)
  --version             show program's version number and exit

```

> Please, always submit the input files with pathway to the file.

