#!/usr/bin/perl 
#       gaa.pl		Version 1.0
#
#       Copyright 2010, 2011 Guohui Yao <gyao at wustl dot edu>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

use warnings;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);

BEGIN{
    my $prog_path = abs_path($0);
       $prog_path = readlink($0) if readlink($0);
    $lib_dir = dirname $prog_path;
    push @INC, $lib_dir;
}
use GAA::UniqueGenerator;
use GAA::Simplifier;	
use GAA::Filter;
use GAA::Validator;
use GAA::Planner;
use GAA::Builder;
use GAA::Gap;
use GAA::Patch;
use GAA::CE;

#VARIABLE SETUP
use vars qw/%opt/;

my $usage="
	-.pl	--target|-t 		<target.fa>		 : target consensus file, the contig names should be in the format of >Contig1.1
	    	--query|-q 		<query.fa> 		 : query  consensus file, the contig names should be in the format of >Contig1.1
	    	
	The following options are optional.
	    	--target-qual|-T 	<target.qual.fa> 	 : target quality file
	    	--query-qual|-Q 	<query.qual.fa> 	 : query quality file
		--target_gap|-g         <target gap file>	 : e.g. gap.txt for Newbler
		--query_gap|-G          <query gap file>	 : e.g. gap.txt for Newbler
	    	--match|-m 		<alignment.psl>          : blat target.fa query.fa -fastMap query_vs_target.psl
	    	
		--pair_status_target|r	<454PairStatus.txt>	 : generated from GSMapper
		--pair_status_query|R	<454PairStatus.txt>	 : generated from GSMapper
	    	--ce                    <ce.psl>          	 : ce status psl file

	    	--output-dir|-o 	<output_dir=pwd> 	 : default is merged_assembly under current directory
            	--tip-size|-p 		<tip_size=90>		 : tip tolerance 
            	--score-cutoff-lower|-l	<score_cutoff_lower=30>	 : confident unique contig, if score < 30
            	--score-cutoff-upper|-u	<score-cutoff_upper=100> : confident match,         if score > 100
            	--contig-size-cutoff|-s	<contig_size_cutoff=100> : keep unique contigs of length > 100
	Using -m option would shortcut the alignment step if alignment.psl alread exists.
";

#GATHER INPUT
&GetOptions(
	"target|t=s"         		=> \$opt{target},
	"query|q=s"          		=> \$opt{query},
	"target-qual|T=s"           	=> \$opt{target_qual},
	"query-qual|Q=s"            	=> \$opt{query_qual},
	"target-gap|g=s"            	=> \$opt{target_gap},
	"query-gap|G=s"            	=> \$opt{query_gap},
	"match|m=s"			=> \$opt{match},
	
	"pair_status_target|r=s"	=> \$opt{pair_status_target},
	"pair_status_query|R=s"		=> \$opt{pair_status_query},
	"ce=s"				=> \$opt{ce},
	
	"output-dir|o=s"            	=> \$opt{output_dir},
	"score-cutoff-lower|l=i"	=> \$opt{score_cutoff_lower},
	"score-cutoff-upper|u=i"	=> \$opt{score_cutoff_upper},
	"contig-size-cutoff|s=i"	=> \$opt{contig_size_cutoff},
	"tip-size|p=i"			=> \$opt{tip_size},
		
	"no-unique|U"	 		=> \$opt{no_unique},
	"no-simplifier|S"		=> \$opt{no_simplifier}
);

my $cmdline = join " ", $0, @ARGV;
$cmdline .= " -t $opt{target}" if defined $opt{target}; 
$cmdline .= " -q $opt{query}"  if defined $opt{query}; 
$cmdline .= " -T $opt{target_qual}" if defined $opt{target_qual};
$cmdline .= " -Q $opt{query_qual}"  if defined $opt{query_qual};
$cmdline .= " -g $opt{target_gap}"  if defined $opt{target_gap};
$cmdline .= " -G $opt{query_gap}"   if defined $opt{query_gap};
$cmdline .= " -m $opt{match}"       if defined $opt{match};
$cmdline .= " -r $opt{pair_status_target}" if defined $opt{pair_status_target};
$cmdline .= " -R $opt{pair_status_query}"  if defined $opt{pair_status_query};
$cmdline .= " --ce $opt{ce}"               if defined $opt{ce};
$cmdline .= " -o $opt{output_dir}"  if defined $opt{output_dir};
$cmdline .= " -l $opt{score_cutoff_lower}" if defined $opt{score_cutoff_lower};
$cmdline .= " -u $opt{score_cutoff_upper}" if defined $opt{score_cutoff_upper};
$cmdline .= " -s $opt{contig_size_cutoff}" if defined $opt{contig_size_cutoff};
$cmdline .= " -p $opt{tip_size}"           if defined $opt{tip_size};
$cmdline .= " -U $opt{no_unique}"          if defined $opt{no_unique};
$cmdline .= " -S $opt{no_simplifier}"      if defined $opt{no_simplifier};


# Verify input
die "$usage\n" unless defined $opt{target} and defined $opt{query};
die "target file doesn't exist!\n$usage" unless -e $opt{target};
die "query  file doesn't exist!\n$usage" unless -e $opt{query};

# Default parameters
$opt{tip_size}           = 90  unless defined $opt{tip_size};
$opt{score_cutoff_lower} = 30  unless defined $opt{score_cutoff_lower};
$opt{score_cutoff_upper} = 200 unless defined $opt{score_cutoff_upper};
$opt{contig_size_cutoff} = 100 unless defined $opt{contig_size_cutoff};
$opt{merged_prefix} = basename($opt{target}) . "_n_" . basename($opt{query}); 	# -.fa -.qual

# Absolute path
$opt{target}      = abs_path $opt{target};
$opt{query}       = abs_path $opt{query};
$opt{match}       = abs_path $opt{match}       if defined $opt{match};
$opt{target_qual} = abs_path $opt{target_qual} if defined $opt{target_qual};
$opt{query_qual}  = abs_path $opt{query_qual}  if defined $opt{query_qual};
$opt{target_gap}  = abs_path $opt{target_gap}  if defined $opt{target_gap};
$opt{query_gap}   = abs_path $opt{query_gap}   if defined $opt{query_gap};
$opt{ce}          = abs_path $opt{ce}          if defined $opt{ce};
$opt{pair_status_target} = abs_path $opt{pair_status_target} if defined $opt{pair_status_target};
$opt{pair_status_query}  = abs_path $opt{pair_status_query}  if defined $opt{pair_status_query};

# Output directory
$opt{output_dir} = $ENV{'PWD'} . "/merged_assembly"  unless defined $opt{output_dir};
`mkdir -p $opt{output_dir}` unless -d $opt{output_dir};
chdir $opt{output_dir};
open LOG, ">log" or die "Cannot write to log file\n";
print LOG $cmdline, "\n";
close LOG;

print STDERR "\n", `date`, ">>> GAA <<<\n\n";

# Get started with Blat Aligner
if(defined $opt{match}){
	die "$opt{match} dosn't exist!$usage" unless -e $opt{match};
	print STDERR "Escape mapping, proceeding to next step!\n"; }
else{
	$opt{match} = "match.psl";
	my $blat_cmd = "blat -fastMap -threads=24 $opt{target} $opt{query} $opt{match}";
	print STDERR `date`, "$blat_cmd\n";
	system $blat_cmd;
	die "Blat Failed!\n$usage" unless -e $opt{match};
}

# Simplify
my $simplifier = GAA::Simplifier->new(
	match_file   => $opt{match},
	score_cutoff => $opt{score_cutoff_upper}				# get reliable match
);
$simplifier->simplify unless $opt{no_simplifier};
undef $simplifier;


#################################################
# Filter
my $filter = GAA::Filter->new(
	match_file  => "1match.unique",
	target_file => $opt{target},
	tip_size    => $opt{tip_size}
);
my $merged_track = $filter->filt;

# Validator
my $validator = GAA::Validator->new(
#	track_file  => "2track",
	track       => $merged_track,
	target_file => $opt{target},
);
my $track = $validator->validate;

# Unique contig generator
my $unique = GAA::UniqueGenerator->new(
	match_file    => $opt{match},
	query_file    => $opt{query},	
	score_cutoff  => $opt{score_cutoff_lower},
	length_cutoff => $opt{contig_size_cutoff},
#	first_scaf    => 0,
#	unique_file   => $opt{merged_prefix}.".unique.fa",
);
my ($unique_scaf, $unique_length) = $unique->generate_unique_scaf_n_length unless $opt{no_unique};

# Planner
my $planner = GAA::Planner->new(
#	track_file  => "v.track",
	track       => $track,
	target_file => $opt{target},
	query_file  => $opt{query},
	unique_scaf => $unique_scaf
);
my $plan = $planner->generate_plan;
undef $unique;
undef $filter;
undef $planner;
#################################################

# Plan Reviser
if(defined $opt{ce}){
	my $patch = GAA::Patch->new(
		plan_file => "3floor_plan",
		ce_file   => $opt{ce},
	);
	$patch->patch;
	undef $patch;
}elsif(defined $opt{pair_status_target} and -e $opt{pair_status_target} and
   defined $opt{pair_status_query}  and -e $opt{pair_status_query})
{
	print STDERR "Run break module: ...\n";
	my $ce = GAA::CE->new(
		match_file   => $opt{match},
		pair_status_target => $opt{pair_status_target},
		pair_status_query => $opt{pair_status_query}
	);
	$ce->ce;	# output ce.psl

	my $patch = GAA::Patch->new(
		plan_file => "3floor_plan",
		ce_file   => "ce.psl",
	);
	$patch->patch;
	undef $patch;
}

# Builder
my $builder = GAA::Builder->new(
	plan_file 	 => "3floor_plan",
	target_file 	 => $opt{target},
	query_file       => $opt{query},
	target_qual_file => $opt{target_qual},
	query_qual_file  => $opt{query_qual},
	
	merged_prefix	 => $opt{merged_prefix}
);
$builder->build;

# Gap
if(defined $opt{target_gap}){
	my $gap = GAA::Gap->new(
		plan_file   => "3floor_plan",
		gap_t_file  => $opt{target_gap},
		target_file => $opt{target},
		merged_file => $opt{merged_prefix}.".fa",
		qual_file => $opt{merged_prefix}.".qual",
	);
	$gap->generate_gap;
}

