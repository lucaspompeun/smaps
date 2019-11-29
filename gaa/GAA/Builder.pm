# 	class GAA::Builder
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


package GAA::Builder;

# constructor
sub new{
	my $class = shift;
	my $self = {@_};
	bless $self, $class;

	print STDERR ">>>\n>>> 4. Builder\n>>>\n";
	die "Builder->new(
        plan_file =>,	
        target_file =>, 
        query_file  => , 
        target_qual_file =>, 
        query_qual_file  =>,)\n"
	unless (-e $self->plan_file or defined $self->plan) and 
           (-e $self->target_file and -e $self->query_file or 
            -e $self->target_qual_file and -e $self->query_qual_file);

	$self->set_plan unless defined $self->plan;
	$self->merged_prefix("merged") unless defined $self->merged_prefix;
	$self->merged_contig($self->merged_prefix . ".fa"); 
	$self->merged_qual($self->merged_prefix . ".qual"); 

	return $self; 
}

# accessor methods
sub plan_file{   $_[0]->{plan_file}   = $_[1] if defined $_[1]; $_[0]->{plan_file};}	# sorted name => track
sub target_file{ $_[0]->{target_file} = $_[1] if defined $_[1]; $_[0]->{target_file}; }
sub query_file{  $_[0]->{query_file}  = $_[1] if defined $_[1]; $_[0]->{query_file};  }
sub target_qual_file{ $_[0]->{target_qual_file} = $_[1] if defined $_[1]; $_[0]->{target_qual_file}; }
sub query_qual_file{  $_[0]->{query_qual_file}  = $_[1] if defined $_[1]; $_[0]->{query_qual_file};  }

sub plan{   $_[0]->{plan}   = $_[1] if defined $_[1]; $_[0]->{plan}; }
sub target{ $_[0]->{target} = $_[1] if defined $_[1]; $_[0]->{target}; }
sub query{  $_[0]->{query}  = $_[1] if defined $_[1]; $_[0]->{query}; }
sub target_qual{ $_[0]->{target_qual} = $_[1] if defined $_[1]; $_[0]->{target_qual}; }
sub query_qual{  $_[0]->{query_qual}  = $_[1] if defined $_[1]; $_[0]->{query_qual};  }

sub merged_prefix{ $_[0]->{merged_prefix} = $_[1] if defined $_[1]; $_[0]->{merged_prefix}; }
sub merged_contig{ $_[0]->{merged_contig} = $_[1] if defined $_[1]; $_[0]->{merged_contig}; }
sub merged_qual{   $_[0]->{merged_qual}   = $_[1] if defined $_[1]; $_[0]->{merged_qual}; }


sub build{
	my $self = shift;
	
	$self->build_contigs if $self->target_file and $self->query_file;
	$self->build_quals   if $self->target_qual_file and $self->query_qual_file;
	$self->prompt(">>> Done Builder <<<");
}


###############
# Contig
###############

sub build_contigs{
	my $self = shift;
	
	$self->target( $self->read_fa($self->target_file) );
	$self->query(  $self->read_fa($self->query_file) );
	
	$self->prompt("Build Contigs => $self->{merged_contig}");

	open OUT, ">$self->{merged_contig}" or die "$! out_contig\n";
	for(sort {$self->by_pcap_name} keys %{$self->plan}){
		my $str = $self->build_contig($self->{plan}{$_});
		print STDERR $self->{plan}{$_}, " 0-len plan, wrong track? $_\n" unless $str;
		my $len = length $str;

		print OUT ">$_\t$len\n";
		for(my $i = 0; $i < $len; $i += 60){
			print OUT substr($str, $i, 60), "\n";
		}
	}
}

#+TContig1.1(0:3753)+QContig8564.1(3682:232)

sub build_contig{
	my ($self, $track) = @_;
	
	if($track =~ /^[+-][TQ]Contig\d+\.\d+$/){
		my ($o, $c, $n) = $track =~ /^([+-])([TQ])(Contig\d+\.\d+)$/;
		if($c eq "T"){
			($o eq "+") ? (return $self->target->{$n}) : (return $self->rc($self->target->{$n}) );
		}else{
			($o eq "+") ? (return $self->query->{$n})  : (return $self->rc($self->query->{$n}) );
		}
	}else{
		my $m_cntg;
		my @blks = $track =~ /([+-][TQ]Contig\d+\.\d+\(\d+\:\d+\))/g;	
		for(@blks){
			my ($o, $c, $n, $s, $l) = /^([+-])([TQ])(Contig\d+\.\d+)\((\d+)\:(\d+)\)$/;
			
			next unless $l > 0; # 0-length is just for marking contig
			
			my $cntg;
			($c eq "T") ? ($cntg = $self->{target}{$n}) : ($cntg = $self->{query}{$n});
			$cntg = substr $cntg, $s, $l;
			$cntg = $self->rc($cntg) if $o eq "-";
			$m_cntg .= $cntg;			
		}
		return $m_cntg;
	}
}


sub rc{
	my ($self, $str) = @_;
	
	$str =~ tr/ACGTacgt/TGCAtgca/;
	reverse $str;
}


###############
# Qual
###############

sub build_quals{
	my $self = shift;
	
	$self->target_qual( $self->read_fa($self->target_qual_file) );
	$self->query_qual(  $self->read_fa($self->query_qual_file) );
	
	$self->prompt("Build Quals   => $self->{merged_qual}");

	open OUT, ">$self->{merged_qual}" or die "$! merged_qual";
	for(sort {$self->by_pcap_name} keys %{$self->{plan}}){
		my @qual = split /\s+/, $self->build_qual($self->{plan}{$_});
		my $len = @qual;
		print OUT ">$_\t$len\n";
		my $i;

		for($i = 0; $i+59 < $len; $i += 60){
			print OUT join(" ", @qual[$i..$i+59]), "\n";
		}
		#print OUT @qual[$i..$#qual], "\n";  ## modified by lye Dec 13, 2010
		print OUT join(" ", @qual[$i..$#qual]), "\n";
	}
}


sub build_qual{
	my ($self, $track) = @_;
	
	if($track =~ /^[+-][TQ]Contig\d+\.\d+$/){
		my ($o, $c, $n) = $track =~ /^([+-])([TQ])(Contig\d+\.\d+)$/;
		if($o eq "+"){
			($c eq "T") ? (return $self->target_qual->{$n}) : (return $self->query_qual->{$n});
		}else{
			my $qual;
			($c eq "T") ? ($qual = $self->target_qual->{$n}) : ($qual = $self->query_qual->{$n});
			my @quals = split / +/, $qual;
			return join(" ", reverse @quals);
		}
	}else{
		my $m_qual;
		my @blks = $track =~ /([+-][TQ]Contig\d+\.\d+\(\d+\:\d+\))/g;	
	
		for(@blks){
			my ($o, $c, $n, $s, $l) = /^([+-])([TQ])(Contig\d+\.\d+)\((\d+)\:(\d+)\)$/;
			
			next unless $l > 0;
			
			my $qual;
			($c eq "T") ? ($qual = $self->target_qual->{$n}) : ($qual = $self->query_qual->{$n});
			my @quals = split / +/, $qual;
			@quals = @quals[$s .. $s+$l-1];
			@quals = reverse @quals if $o eq "-";
			$m_qual .= join(" ", @quals) . " ";
		}

		$m_qual =~ s/ $//;
		$m_qual =~ s/^ //;
		return $m_qual;
	}
}


#####################
### Util
#####################

sub prompt{
	my $self = shift;
	my $date = `date`;
	my ($time) = $date =~ /(\d\d:\d\d:\d\d)/;
	print STDERR join("", @_), "\n$time\n";
}

sub set_plan{
	my $self = shift;
	
	$self->prompt("Read floor plan");

	open IN, $self->plan_file or die "$! plan_file\n";
	while(<IN>){
		chomp;
		my ($ord, $track, $len) = split /\t/, $_;
		$self->{plan}{$ord} = $track;
#		$self->{ord_length}{$ord} = $len;
	} 
}


sub read_fa{
	my ($self, $f) = @_;

	$self->prompt("Read Fasta:\t", $f);

	open IN, $f or die "$! fasta file\n";
	my $h;
	my $name = "";
	while(<IN>){
		chomp;

		if(/^>/){ 
			$_ =~ s/^>//; 
			($name) = split /\s/, $_;
		}else{ 
                        my $spacer = "";
                        $spacer = " " if $_ =~ /\d+/; # if quality values, spacer is added by lye @ Dec 13, 2010
			$h->{$name} .= $_ . $spacer;
		} 	
	}
        $h  =~ s/ $//;
	$h;
}

#for "[>+-]Contig0.1" format
sub by_pcap_name{
	my $self = shift;

	my ($ai, $af) = $a =~ /^>?[+-]?Contig(\d+)\.(\d+)/;
	my ($bi, $bf) = $b =~ /^>?[+-]?Contig(\d+)\.(\d+)/;
	
	if($ai == $bi){
		$af <=> $bf;
	}else{
		$ai <=> $bi;
	}
}

1;
