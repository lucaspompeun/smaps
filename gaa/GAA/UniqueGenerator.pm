#	class GAA::UniqueGenerator
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


package GAA::UniqueGenerator;
use File::Basename;

=h
class UniqueGenerator{
	#is => 'Genome::Model::Tools::GAA',
	has => [
    		query_file => { is => "file name", default => "", doc => "target assembly" },
    		match_file => { is => "file name", default => "", doc => "all match.psl"},
    		first_scaf => { is => "number", default => 0, doc => "the first scaf of unique contigs"}
		score_cutoff  = 30;
		length_cutoff = 100;
		unique_file   = query.uni;
    	]		
}
=cut

sub new{
	my $class = shift;
	my $self  = {@_};
	bless($self, $class);
	
	print STDERR ">>>\n>>> P12. Unique Generator\n>>>\n";
	
	die "UniqueGenerator->new(query_file =>, match_file =>, first_scaf => 0, score_cutoff => 30, contig-length-cutoff => 100, unique_file => query.unique\n" 
	unless -e $self->{query_file} and -e $self->{match_file};

	$self->score_cutoff(30)         unless defined $self->score_cutoff;
	$self->length_cutoff(100)       unless defined $self->length_cutoff;
	$self->unique_file("unique.fa") unless defined $self->unique_file;
	$self->first_scaf(0)            unless defined $self->first_scaf;
	
	return $self;
}


#accessor methods
sub query_file{  $_[0]->{query_file}  = $_[1] if defined $_[1]; $_[0]->{query_file}; }
sub match_file{  $_[0]->{match_file}  = $_[1] if defined $_[1]; $_[0]->{match_file}; }
sub unique_file{ $_[0]->{unique_file} = $_[1] if defined $_[1]; $_[0]->{unique_file}; }

sub query{       $_[0]->{query}       = $_[1] if defined $_[1]; $_[0]->{query}; }
sub match_score{ $_[0]->{match_score} = $_[1] if defined $_[1]; $_[0]->{match_score}; }
sub unique{      $_[0]->{unique}      = $_[1] if defined $_[1]; $_[0]->{unique}; }

sub first_scaf{    $_[0]->{first_scaf}    = $_[1] if defined $_[1]; $_[0]->{first_scaf}; }
sub score_cutoff{  $_[0]->{score_cutoff}  = $_[1] if defined $_[1]; $_[0]->{score_cutoff}; }
sub length_cutoff{ $_[0]->{length_cutoff} = $_[1] if defined $_[1]; $_[0]->{length_cutoff}; }

sub scaf{              $_[0]->{scaf}              = $_[1] if defined $_[1]; $_[0]->{scaf}; }
sub scaf_length{       $_[0]->{scaf_length}       = $_[1] if defined $_[1]; $_[0]->{scaf_length}; }
sub unique_floor_plan{ $_[0]->{unique_floor_plan} = $_[0] if defined $_[1]; $_[0]->{unique_floor_plan}; }


sub generate_unique_scaf_n_length{
	my $self = shift;
	
	$self->read_query;			# get all contigs
	$self->get_match_score;			# get non-unique contigs
	$self->get_unique_contig;		# remove non-unique & short
	$self->set_scaf_n_length;		# get unique scaf
	$self->dump_query_unique_contigs;

	print STDERR `date`;
	($self->scaf, $self->scaf_length);
}


sub generate_unique_contigs{
	my $self = shift;
	
	$self->generate_unique_scaf_n_length;
	
	$self->sort_scaf;
	$self->dump_contigs;

	print STDERR `date`;
	($self->scaf, $self->scaf_length);
}


sub read_query{
	my $self = shift;

	$self->prompt("Read query:\t", $self->{query_file});
	
	open IN, $self->{query_file} or die "$!:query\n";
	my $name;
	while(<IN>){
		chomp;
		if(/^>/){ 
			$_ =~ s/^>//; 
			($name) = split /\s/, $_;
		}else{ 
			$self->{query}{$name} .= $_;
		} 	
	}
	print STDERR "# query contigs:     \t", scalar keys %{$self->{query}}, "\n";
}


sub get_match_score{
	my $self = shift;
	
	$self->prompt("Get match score:\t", $self->{match_file});
	
	open IN, $self->{match_file} or die "$!:match_file\n";;
#	while(<IN>){ last if /^-/; } 		# scap psl header	
	while (<IN>){
		next unless /^\d/;
		chomp;
		my @parts = split;
		my ($matches, $q_name) = @parts[0, 9];
		$self->{match_score}{$q_name} = $matches if $matches >= $self->{score_cutoff} and $matches > $self->{match_score}{$q_name};
	}
	close IN;
	
	print STDERR "# non-unique contigs:\t", scalar keys %{$self->{match_score}}, "\n";
}


sub get_unique_contig{
	my $self = shift;
	
	for(keys %{$self->{match_score}}){
		(defined $self->{query}{$_}) ? (delete $self->{query}{$_}) : (die "$_ is not in query\n");
	}
	print STDERR "# unique contigs:    \t", scalar keys %{$self->{query}}, "\n";
	
	for(keys %{$self->{query}}){
		delete $self->{query}{$_} if length($self->{query}{$_}) < $self->{length_cutoff};
	}
	print STDERR "# unique-short contigs:\t", scalar keys %{$self->{query}}, "\n";
}




#################
### Scaffolding
#################

# Q$s -> [+-]QContig$s.$c
sub set_scaf_n_length{
	my $self = shift;	
	
	for(sort by_pcap_name keys %{$self->{query}}){
		my ($s, $c) = $_ =~ /^Contig(\d+)\.(\d+)/;
		push @{$self->{scaf}{"Q$s"}}, "+QContig$s.$c";
		$self->{scaf_length}{"Q$s"} += length $self->{query}{$_};
	}
	print STDERR "# unique scaffolds:  \t", scalar keys %{$self->{scaf}}, "\n";
}


sub sort_scaf{
	my $self = shift;
		
	my $i = $self->first_scaf;
	my $l = $self->scaf_length;
	for my $s (reverse sort { $l->{$a} <=> $l->{$b} } keys %{$self->{scaf}}){
		my $j = 0;
		for my $c (@{$self->{scaf}{$s}}){
			$j++;
			my $ori = "Contig" . substr($c, 1);
			my $ord = "Contig$i.$j";
			$self->{unique}{$ord} = $self->{query}{$ori};
			$self->{unique_floor_plan}{$ord} = "q".$ori;
			$self->{contig_length}{$ord} = length $self->{query}{$ori};
		}
		$i++;
	}
	$self->{last_scaf} = $i - 1;
	print STDERR "Last unique scaffold:\t", $self->{last_scaf}, "\n";
	
	open TRK, ">P12floor_plan.unique" or die "$!:>>track\n";
	print TRK "$_\t$self->{unique_floor_plan}{$_}\t$self->{contig_length}{$_}\n" for sort by_pcap_name keys %{$self->{unique_floor_plan}};
}

# [>+-]Contig0.1
sub by_pcap_name{
	my ($ai, $af) = $a =~ /^>?[+-]?Contig(\d+)\.(\d+)/;
	my ($bi, $bf) = $b =~ /^>?[+-]?Contig(\d+)\.(\d+)/;
	
	if($ai == $bi){
		$af <=> $bf;
	}else{
		$ai <=> $bi;
	}
}

sub dump_contigs{
	my $self = shift;
	
	open OUT, ">$self->{unique_file}" or die "$!:>unique_file\n";
	for(sort by_pcap_name keys %{$self->{unique}}){
		my $len = length $self->{unique}{$_}; 
		
		print OUT ">$_\t$len\n";
		for(my $i = 0; $i < $len; $i += 60){
			print OUT substr($self->{unique}{$_}, $i, 60), "\n";
		}
	}
	print STDERR "unique contigs file: $self->{unique_file}\n";
}

sub dump_query_unique_contigs{
	my $self = shift;
	
	open OUT, ">$self->{unique_file}" or die "$!:>unique_file\n";
	for(sort by_pcap_name keys %{$self->{query}}){
		my $len = length $self->{query}{$_}; 
		
		print OUT ">$_\t$len\n";
		for(my $i = 0; $i < $len; $i += 60){
			print OUT substr($self->{query}{$_}, $i, 60), "\n";
		}
	}
	print STDERR "unique contigs file:\t", $self->{unique_file}, "\n";
}

sub prompt{
	my $self = shift;
	my $date = `date`;
	my ($time) = $date =~ /(\d\d:\d\d:\d\d)/;
	print STDERR join("", @_), "\n$time\n";
}

# PSL format:
# matches	mismatches	rep_matches	n_count		q_num_insert	q_base_insert	t_num_insert	t_base_insert	strand		q_name
# q_length	q_start		$q_end		13t_name	t_length	t_start		t_end		block_count	block_sizes	q_starts
# t_starts	q_seqs		t_seqs

1;

