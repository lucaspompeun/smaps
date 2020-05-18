<p align="center">
    <img src="https://img.shields.io/github/issues/lucaspompeun/smaps" />
    <img src="https://img.shields.io/github/stars/lucaspompeun/smaps" />
    <img src="https://img.shields.io/github/forks/lucaspompeun/smaps" />
    <img src="https://img.shields.io/github/license/lucaspompeun/smaps" /
</p>

# Smaps
A pipeline based on unmapped reads, which combines different tools for extends contigs and closing gaps present in the genomes.

## Requirements
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) in `PATH` variable

- [Python3](https://www.python.org/downloads/) or higher

- [Java](https://www.java.com/download/) version 1.8.0_201 or higher

- [Quast](http://bioinf.spbau.ru/quast) in `PATH` variable

- [Perl](https://www.perl.org/get.html) language installed

- [Prokka](https://github.com/tseemann/prokka) in `PATH` variable

## Installation
You need to install git

```sh
sudo apt install git
```

Clone smaps repository to your machine

```
cd && git clone https://github.com/lucaspompeun/smaps.git
```

Creat a symbolic link to Smaps
```
sudo chmod 777 smaps/ && ln -s ~/smaps/smaps.py /usr/local/bin/smaps
```

With all dependecies satisfied you can simple run on terminal

```
smaps
```

> Please, always submit the input files with pathway to the file.

