#	class GAA::Planner
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


package GAA::Planner;
use List::Util qw(sum);

sub new{
	my $class = shift;
	my $self  = {@_};
	bless($self, $class);
	
	print STDERR ">>>\n>>> 3. Planner\n>>>\n";
	die "Planner->new(track_file => 2track, target_file =>, query_file =>)\n" 
	unless (-e $self->track_file or defined $self->track) and -e $self->target_file and -e $self->query_file;

	return $self;
}

#accessor methods
sub target_file{ $_[0]->{target_file} = $_[1] if defined $_[1]; $_[0]->{target_file}; }
sub query_file{  $_[0]->{query_file}  = $_[1] if defined $_[1]; $_[0]->{query_file}; }
sub track_file{  $_[0]->{track_file}  = $_[1] if defined $_[1]; $_[0]->{track_file}; }

sub target{ $_[0]->{target} = $_[1] if defined $_[1]; $_[0]->{target}; }
sub query{  $_[0]->{query}  = $_[1] if defined $_[1]; $_[0]->{query}; }
sub track{  $_[0]->{track}  = $_[1] if defined $_[1]; $_[0]->{track}; }
sub track_domina{$_[0]->{track_domina}= $_[1] if defined $_[1]; $_[0]->{track_domina}; }

sub scaf{        $_[0]->{scaf}        = $_[1] if defined $_[1]; $_[0]->{scaf}; }
sub scaf_degree{ $_[0]->{scaf_degree} = $_[1] if defined $_[1]; $_[0]->{scaf_degree}; }
sub scaf_length{ $_[0]->{scaf_length} = $_[1] if defined $_[1]; $_[0]->{scaf_length}; }
sub unique_scaf{ $_[0]->{unique_scaf} = $_[1] if defined $_[1]; $_[0]->{unique_scaf}; }
sub sorted_scaf{ $_[0]->{sorted_scaf} = $_[1] if defined $_[1]; $_[0]->{sorted_scaf}; }

sub scaf_adj{ $_[0]->{scaf_adj} = $_[1] if defined $_[1]; $_[0]->{scaf_adj}; }
sub scaf_link{$_[0]->{scaf_link}= $_[1] if defined $_[1]; $_[0]->{scaf_link};}

sub plan{ $_[0]->{plan} = $_[1] if defined $_[1]; $_[0]->{plan}; }	# sorted name => track


sub generate_plan{
	my $self = shift;
	
	$self->set_target;
	$self->set_scaf;		# format: [+-]\d+ => [+-][TQM]Contig
	$self->set_track unless defined $self->track;
	$self->set_query;

	$self->set_scaf_adj_n_domina;
	$self->set_scaf_link;
	$self->merge_scaf;
	$self->add_unique_scaf if defined $self->unique_scaf;
	$self->sort_scaf;
	$self->prompt(">>> Done Planner <<<");
	
	$self->plan;
}


############
### Read
############

sub set_target{
	my $self = shift;
	
	$self->prompt("Read Target:\t", $self->target_file);
	open IN, $self->target_file or die "$! target_file\n";
	my $name;
	while(<IN>){
		chomp;
		if(/^>/){ 
			$_ =~ s/>//; 
			($_) = split /\s/, $_;
			$name = $_;
		}else{ 
			$self->{target}{$name} .= $_;
		} 	
	}
}
# scaffold range: [0 -]
# contig range [1 -]
sub set_scaf{
	my $self = shift;

	$self->prompt("Set Scaffold");

	for(sort {$self->by_pcap_name} keys %{$self->target}){
		my ($s, $c) = $_ =~ /Contig(\d+)\.(\d+)/; 
		$self->{scaf}{$s}[$c-1] = "+T$_";			#distinguish old contig from merged contig by if have "Contig"
	}

	for(keys %{$self->scaf}){
		$self->{scaf_degree}[$_] = @{$self->scaf->{$_}};
	}
	
	my $n_scaf = scalar keys %{$self->scaf};
	my $n_cntg = sum(0, @{$self->scaf_degree});
	print STDERR "# scaffolds\t$n_scaf\n# contigs\t$n_cntg\n# average\t", $n_cntg/$n_scaf, "\n";

	open OUT, ">3.scaffold.degree" or die $!, " 3.scaffold.degree\n";
	print OUT $_, "\n" for @{$self->scaf_degree};
	open OUT, ">3.scaffold.target" or die "$! 3.scaffold.target\n";
	for( sort {abs $a <=> abs $b} keys %{$self->scaf}){
		print OUT "$_\n", join( "\n", @{$self->scaf->{$_}} ), "\n";
	}
}

sub set_track{
	my $self = shift;

	$self->prompt("Read Track\t", $self->track_file);

	open IN, $self->track_file or die "$! track_file\n";
	while(<IN>){
		my ($t) = split /\s+/, $_;
		$self->{track}{$t} = $self->get_track_length($t);
	}

}

sub set_query{
	my $self = shift;
	
	$self->prompt("Read Query\t", $self->query_file);
	open IN, $self->query_file or die "$! query_file\n";
	my $name;
	while(<IN>){
		chomp;
		if(/^>/){ 
			$_ =~ s/>//; 
			($_) = split /\s/, $_;
			$name = $_;
		}else{ 
			$self->{query}{$name} .= $_;
		} 	
	}
}


#####################
# Scaffold Adjacency
#####################

sub set_scaf_adj_n_domina{
	my $self = shift;

	$self->prompt("Remove Merged Contig =>Scaffold Adj, Track Domina");

	for(keys %{$self->track}){
		$self->set_scaf_adj_n_domina_of($_);
	}
	
	open OUT, ">3.track_domina" or die "$! 3.track_domina\n";
	for(sort {$self->by_pcap_name} keys %{$self->{track_domina}}){
		print OUT "$_\t", $self->track_domina->{$_}, "\n";
	}
	open OUT, ">3.scaffold.adj" or die "$! 3.scaffold.adj\n";
	for(sort {(abs $a)+(abs $self->scaf_adj->{$a}) <=> (abs $b)+(abs $self->scaf_adj->{$b})} keys %{$self->{scaf_adj}}){
		print OUT "$_\t", $self->scaf_adj->{$_}, "\n";
	}

	$self->clean_up_scaf;
	my $n_scaf = scalar keys %{$self->scaf};
	my $n_cntg = 0;
	$n_cntg += @{$self->scaf->{$_}} for keys %{$self->scaf};
	print STDERR "# scaffolds\t$n_scaf\n# contigs\t$n_cntg\n# average\t", $n_cntg/$n_scaf, "\n";

	open OUT, ">3.scaffold.merged_contig" or die "$! 3.scaffold.merged_contig\n";
	for( sort {abs $a <=> abs $b} keys %{$self->scaf}){
		print OUT "$_\n", join( "\n", @{$self->scaf->{$_}} ), "\n";
	}
}
sub set_scaf_adj_n_domina_of{
	my ($self, $track) = @_;

	my @osc = $self->get_track_osc($track);

	my ($o, $s, $c) = @osc[0..2];
	my $domina = "$o"."MContig$s.$c";
	$self->{track_domina}{$domina} = $track;
	if($o eq "+"){
		$self->scaf->{$s}[$c-1] = $domina;
	}else{
		$self->scaf->{$s}[$c-1] = $self->negative($domina);
	}

	if(@osc == 6){
		my $os1 = $osc[0].$osc[1];
		my $os2 = $osc[3].$osc[4];
		$self->{scaf_adj}{$os1} = $os2;
		$self->{scaf_adj}{$self->negative($os2)} = $self->negative($os1);
	}
	
}
sub sife_closed_scaf{
	my ($self, $track) = @_;

    #remove duplicate
	my @t = $track =~ /([+-]TContig\d+\.\d+)/g;
	my $tt = join("", @t);
	$tt =~ s/([+-]TContig\d+\.\d+)\1+/$1/g;
	my @osc = $tt =~ /([+-])TContig(\d+)\.(\d+)/g;

	my %open;
	for(my $i = 0; $i < @osc; $i += 3){
		my ($o, $s, $c) = @osc[$i..$i+2];
		$self->{scaf}{$s}[$c-1] = 0;		# remove merged contigs

		my @occurs = $tt =~ /Contig$s\.\d+/g;
		if(@occurs < $self->{scaf_degree}[$s]){
			$open{$s}{ort} = $o;
			$open{$s}{min} = $c if not defined $open{$s}{min} or $c < $open{$s}{min};
			$open{$s}{max} = $c if not defined $open{$s}{max} or $c > $open{$s}{max};
		}
	}
	unless(scalar keys %open){
		$open{$osc[1]}{ort} = $osc[0];
		$open{$osc[1]}{min} = $osc[2];
	}
	%open;
}
sub get_track_osc{
	my ($self, $track) = @_;
	my %open = $self->sife_closed_scaf($track);

	if(scalar keys %open == 1){
		my ($k) = keys %open;
		return ($open{$k}{ort}, $k, $open{$k}{min});
	}else{
		my ($a, $b) = keys %open;
		my ($ort_a,  $ort_b ) = ($open{$a}{ort},  $open{$b}{ort} );
		my ($head_a, $tail_a) = ($open{$a}{min} == 1, $open{$a}{max} == $self->{scaf_degree}[$a]);
		my ($head_b, $tail_b) = ($open{$b}{min} == 1, $open{$b}{max} == $self->{scaf_degree}[$b]);
		
		my @left  = ($ort_a, $a, $open{$a}{min});
		my @right = ($ort_b, $b, $open{$b}{min});
		
		if($head_a and $head_b){					#both heads are in
			($ort_a eq "+") ? (return (@right, @left)) : (return (@left, @right));
		}elsif($tail_a and $tail_b){					#both tails are in
			($ort_a eq "+") ? (return (@left, @right)) : (return (@right, @left));
		}elsif($head_a and $tail_b){					#head, tail
			($ort_a eq "+") ? (return (@right, @left)) : (return (@left, @right));
		}elsif($tail_a and $head_b){					#tail, head
			($ort_a eq "+") ? (return (@left, @right)) : (return (@right, @left));
		}
print STDERR $track, " wrong track in get_track osc\n"; #@left; #TODO
	}
}
sub clean_up_scaf{
	my $self = shift;
	
	for(sort {$a <=> $b} keys %{$self->scaf}){
		for(my $i = $#{$self->scaf->{$_}}; $i >= 0; $i--){
			splice(@{$self->scaf->{$_}}, $i, 1) unless $self->scaf->{$_}[$i];
		}
		delete $self->scaf->{$_} if @{$self->scaf->{$_}} == 0;
	}
}


#################
# Scaffold Link
#################

sub set_scaf_link{
	my $self = shift;

	$self->prompt("Scaffold Link");

	my %visited;
	for my $os (keys %{$self->scaf_adj}){
		my $ros = $self->negative($os);
		next if defined $self->scaf_adj->{$ros};	# scape non-terminal scaf
		next if $visited{$os};
	
		my $link = $os;
		$self->traverse($os, \$link, \%visited);
		$self->{scaf_link}{$link} = 1;
	}
	
	open OUT, ">3.scaffold.link" or die "3.scaffold.link\n";	
	print OUT "$_\n" for keys %{$self->{scaf_link}};
}
sub traverse{
	my ($self, $os, $link_ref, $visited) = @_;

	return unless defined $self->scaf_adj->{$os};
	
	my $next_os = $self->scaf_adj->{$os};
	${$link_ref} .= $next_os;
	$visited->{$self->negative($next_os)} = 1;	# avoid reverse link

	$self->traverse($next_os, $link_ref, $visited);
}


#################
# Merge Scaffold
#################

sub merge_scaf{
	my $self = shift;

	$self->prompt("Merge Scaffold");

	for(keys %{$self->scaf_link}){
		my @os = $_ =~ /([+-]\d+)/g;
		my @contigs;
		for(@os){
			my ($o, $s) = $_ =~ /([+-])(\d+)/;
			next if @{$self->scaf->{$s}} == 0;
			if($o eq "+"){
				push @contigs, @{$self->scaf->{$s}};
			}else{
				push @contigs, $self->reverse_scaf($self->scaf->{$s});
			}
			delete $self->scaf->{$s};
		}

		my ($fo, $fs) = $os[0] =~ /([+-])(\d+)/;
		if($fo =~ /^\+/){
			@{$self->{scaf}{$fs}} = @contigs;
		}else{
			@{$self->{scaf}{$fs}} = $self->reverse_scaf(\@contigs);
		}
	}

	my $n_scaf = scalar keys %{$self->scaf};
	my $n_cntg = 0;
	$n_cntg += @{$self->scaf->{$_}} for keys %{$self->scaf};
	print STDERR "# scaffolds\t$n_scaf\n# contigs\t$n_cntg\n# average\t", $n_cntg/$n_scaf, "\n";

	open OUT, ">3.scaffold.merged_scaf" or die "$! 3.scaffold.merged_scaf\n";
	for( sort {abs $a <=> abs $b} keys %{$self->scaf}){
		print OUT "$_\n", join( "\n", @{$self->scaf->{$_}} ), "\n";
	}
	
}


#################
# Sort Scaffold
#################

sub sort_scaf{
	my $self = shift;
	
	$self->prompt("Sort Scaffold");
	$self->get_scaf_length;

	my $i = 0;
	for(reverse sort { $self->scaf_length->{$a} <=> $self->scaf_length->{$b} } keys %{$self->scaf}){
		for(@{$self->scaf->{$_}}){
			if(/^[+-][TQ]Contig\d+\.\d+$/){							# target contig
				push @{$self->{sorted_scaf}{$i}}, $_;
			}elsif(/^[+-]MContig\d+\.\d+$/){
				if(defined $self->track_domina->{$_}){
					my $t = $self->track_domina->{$_};
					push @{$self->{sorted_scaf}{$i}}, $t;
				}elsif(defined $self->track_domina->{$self->negative($_)}){
					my $t = $self->track_domina->{$self->negative($_)};
					push @{$self->{sorted_scaf}{$i}}, $self->reverse_track($t);
				}else{  print STDERR "Domina not found: $_\n"; }
			}else{ print STDERR "Unknown contig name:$_|sort_scaf\n"; }
		}
		$i++;
	}

	open  OUT, ">3floor_plan" or die "$! 3floor_plan\n";
	my $i = 0;
	for(sort {$a <=> $b} keys %{$self->sorted_scaf}){
		my $j = 1;
		for( @{$self->sorted_scaf->{$_}} ){ 
			my $n = "Contig$i.$j";
			$self->{plan}{$n} = $_;
			print OUT "$n\t$_\n";
			$j++; 
		}
		$i++;
	}
}
sub get_scaf_length{
	my $self = shift;
	
	for my $s (keys %{$self->scaf}){
		my $len = 0;
		for(@{$self->scaf->{$s}}){
			if(/^[+-]TContig\d+\.\d+$/){						# target contig
				my ($n) = $_ =~ /^[+-]T(Contig\d+\.\d+)$/;
				$len += length $self->target->{$n};
			}elsif(/^[+-]MContig\d+\.\d+$/){					# merged contig 
				my $t = $self->track_domina->{$_};
				$len += $self->get_track_length($t);
			}elsif(/^[+-]QContig\d+\.\d+$/){					#query contig
				my ($n) = $_ =~ /^[+-]Q(Contig\d+\.\d+)$/;
				 $len += length $self->query->{$n};
			}else{ print STDERR "Unknow contig name:$_|get_scaf_length\n"; }
		}
		$self->{scaf_length}{$s} = $len;
	}

	open OUT, ">3.scaffold.length" or die "$! 3.scaffold.length\n";
	print OUT $self->scaf_length->{$_}, "\n" for sort {$a <=> $b} keys %{$self->scaf_length};
}
sub add_unique_scaf{
	my $self = shift;
	
	$self->prompt("Add Unique Contigs");

	while( my ($key, $value) = each( %{$self->unique_scaf} )){
		$self->{scaf}{$key} = $value;
#        	print "$key => $value\n";
    	}

	my $n_scaf = scalar keys %{$self->unique_scaf};
	my $n_cntg = 0;
	$n_cntg += @{$self->unique_scaf->{$_}} for keys %{$self->unique_scaf};
	
	($n_scaf > 0) ? (print STDERR "# unique scaffolds\t$n_scaf\n# unique contigs\t$n_cntg\n# average\t", $n_cntg/$n_scaf, "\n") : (print STDERR "# unique scaffolds\t0\n# unique contigs\t$n_cntg\n");

	my $t_scaf = scalar keys %{$self->scaf},
	my $t_cntg = 0;
	$t_cntg += @{$self->scaf->{$_}} for keys %{$self->scaf};
	print STDERR "# total scaffolds:\t", scalar keys %{$self->scaf}, "\n";
	print STDERR "# total contigs:  \t$t_cntg\n# average:\t", $t_cntg/$t_scaf, "\n";
}


###########
# Util
###########

#Contig225.9	0	121376	Contig17.157	736	122147	121163	+	124466	122179
sub by_pcap_name{
	my $self = shift;

	my ($ai, $af) = $a =~ /^>?[+-]?[TQM]?Contig(\d+)\.(\d+)/;
	my ($bi, $bf) = $b =~ /^>?[+-]?[TQM]?Contig(\d+)\.(\d+)/;
	
	if($ai == $bi){
		$af <=> $bf;
	}else{
		$ai <=> $bi;
	}
}

sub add{
	my ($self, $a, $b) = @_;	
	return "+" if( ($a ne "-" and $b ne "-") or ($a eq "-" and $b eq "-") );
	return "-";
}
sub negative{
	my ($self, $a) = @_;	
	
	my ($sign, $body) = $a =~ /^([+-])(.*)$/;
	$self->add("-", $sign).$body;
}

sub reverse_scaf{
	my ($self, $ss) = @_;

	my @rss;
	for(reverse @$ss){
		push @rss, $self->negative($_);
	}
	@rss;
}

sub reverse_track{
	my ($self, $trk) = @_;

	$trk =~ tr/+-/-+/ ;
	@blks = $trk =~ /([+-][TQ]Contig\d+\.\d+\(\d+\:\d+\))/g;
	@blks = reverse @blks;
	join("", @blks);
}

sub get_track_length{
	my ($self, $track) = @_;

	my @len = $track =~ /[+-][TQ]Contig\d+\.\d+\(\d+\:(\d+)\)/g;	
	sum(0, @len);
}

sub prompt{
	my $self = shift;
	my $date = `date`;
	my ($time) = $date =~ /(\d\d:\d\d:\d\d)/;
	print STDERR join("", @_), "\n$time\n";
}



1;
