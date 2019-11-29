#       class GAA::Gap
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


package GAA::Gap;
use List::Util qw(sum);

sub new{
	my $class = shift;
	my $self  = {@_};
	bless $self, $class;
	
	print STDERR ">>>\n>>> Gap\n>>>\n";
	die "Gap->new(
        plan_file =>, 
        gap_t_file =>, 
        target_file =>, 
        merged_file =>, [qual_file =>])\n" 
	unless -e $self->plan_file and -e $self->gap_t_file and -e $self->target_file and -e $self->merged_file;

	return $self;
}

#accessor methods
sub plan_file{   $_[0]->{plan_file}   = $_[1] if defined $_[1]; $_[0]->{plan_file}; }
sub plan{        $_[0]->{plan}        = $_[1] if defined $_[1]; $_[0]->{plan}; }
sub gap_t_file{  $_[0]->{gap_t_file}  = $_[1] if defined $_[1]; $_[0]->{gap_t_file}; }
sub gap_t{       $_[0]->{gap_t}       = $_[1] if defined $_[1]; $_[0]->{gap_t}; }
sub target_file{ $_[0]->{target_file} = $_[1] if defined $_[1]; $_[0]->{target_file}; }
sub contig_len { $_[0]->{contig_len}  = $_[1] if defined $_[1]; $_[0]->{contig_len};  }
sub merged_file{ $_[0]->{merged_file} = $_[1] if defined $_[1]; $_[0]->{merged_file}; }
sub merged{      $_[0]->{merged}      = $_[1] if defined $_[1]; $_[0]->{merged}; }
sub qual_file{   $_[0]->{qual_file}   = $_[1] if defined $_[1]; $_[0]->{qual_file}; }

sub rplan{       $_[0]->{rplan}       = $_[1] if defined $_[1]; $_[0]->{rplan}; }
sub nplan{       $_[0]->{nplan}       = $_[1] if defined $_[1]; $_[0]->{nplan}; }
sub gap{         $_[0]->{gap}         = $_[1] if defined $_[1]; $_[0]->{gap}; }

sub generate_gap{
	my $self = shift;
	
	$self->read_plan;
#	$self->dump_plan;
	$self->read_gap_t;
#	$self->dump_gap_t;
	$self->read_contig_len;
#	$self->dump_contig_len;
	$self->read_merged;

	$self->set_rplan;
	$self->dump_rplan;

	$self->calc_shrink;
	$self->dump_shrink;
	$self->get_negative;
	$self->dump_negative;
	
	$self->align_neighbor;
	$self->get_neighbor_alignment;
	$self->dump_neighbor_alignment;

	$self->build_link;
	$self->dump_link;

	$self->build_ntrack;
	$self->dump_ntrack;

	$self->merge_neighbor;
	$self->rename_contig;
	$self->dump_gap;
	$self->dump_nplan;
	if(-e $self->{qual_file}){
		$self->read_qual;
		$self->merge_neighbor_qual;
		$self->rename_qual;
	}
	
	$self->prompt(">>> Done Gap <<<");
}


sub calc_shrink{
	my $self = shift;
	
	$self->prompt("Calculate shrink");

	for my $s (keys %{$self->gap_t}){
		for(my $c = 0; $c < @{$self->gap_t->{$s}}; $c++){

			my $aa = "TContig$s.".($c+1);
			my $bb = "TContig$s.".($c+2);
			my $ra = $self->rplan->{$aa};
			my $rb = $self->rplan->{$bb};
			next if($ra eq $rb);

			my ($o1, $s1, $c1) = $ra =~ /([+-])(\d+)\.(\d+)/;
			my ($o2, $s2, $c2) = $rb =~ /([+-])(\d+)\.(\d+)/;
	
# dump wrong plan
print STDERR "wrong rplan: $ra=p-osc,$rb=c-osc\n" 
if $o1 ne $o2 or $s1 ne $s2 or ($c2-$c1 != $o1."1");
		
			my $p1 = $self->plan->{$s1}[$c1-1];
			my $p2 = $self->plan->{$s2}[$c2-1];
					
			if ($o1 eq "-"){
				$p1 = $self->reverse_track($p1);
				$p2 = $self->reverse_track($p2);
			}

			my $shrink = 0;
			$shrink   += $self->get_tail("+$aa", $p1);
			$shrink   += $self->get_head("+$bb", $p2);
			
			my $gap = $self->{gap_t}->{$s}[$c] - $shrink;
			my $rep = "$s1.$c1";
			   $rep = "$s2.$c2" if $o1 eq "-";
			$self->{shrink}{"$rep"} = $shrink;
			$self->{gshrink}{"$rep"} = $gap;
		}	
	}
}


sub dump_shrink{
	my $self = shift;
	
	$self->prompt("Dump Shrink:\tg.shrink");
	
	open OUT, ">g.shrink" or die "$! g.shrink\n";
	for(sort {$self->by_pcap_name} keys %{$self->{shrink}}){
		print OUT "Contig$_\t", $self->{shrink}{$_}, "\t", $self->{gshrink}{$_}, "\n";
	}
}

sub get_negative{
	my $self = shift;
	for(keys %{$self->{gshrink}}){
		if($self->{gshrink}->{$_} < 0){
			$self->{negative}{$_} = $self->{gshrink}->{$_};
		}
	}
}

sub dump_negative{
	my $self = shift;
	open OUT, ">g.negative" or die "$! g.negative\n";
	for(keys %{$self->{negative}}){
		print OUT "$_\t", $self->{negative}->{$_}, "\n";
	}
}


sub align_neighbor{
	my $self = shift;
	
	$self->prompt("Align Neighbors");
	if(-e "g.align_neighbor"){
		unlink "g.align_neighbor" or warn "Could not unlink file g.align_neighbor: $!";
	}
	
	my $cnt = 0;
	for(keys %{$self->{negative}}){
		$self->align_neighbor_of($_);
		$cnt++;
		print "$cnt " if $cnt%10 == 0;
	}
	print "\n$cnt pairs of contigs are aligned.\n";
}

sub align_neighbor_of{
	my ($self, $a) = @_;
	
	my ($s, $c) = $a =~ /(\d+)\.(\d+)/;
	my $b = $s.".".($c+1);
	
	open OA, ">g.a" or die "$! g.a\n";
	print OA ">Contig$a\n", $self->{merged}->{"Contig".$a}, "\n";
	open OB, ">g.b" or die "$! g.b\n";
	print OB ">Contig$b\n", $self->{merged}->{"Contig".$b}, "\n";
	
	my $blat_out = `blat g.b g.a -fastMap -noHead g.tmp`;
	`cat g.tmp >> g.align_neighbor`;
}


sub get_neighbor_alignment{
	my $self = shift;

	my $tip_rate = 0.2;
	my (%hh, %h);

	open IN, "g.align_neighbor" or warn "$! g.align_neighbor\n";

	while (<IN>){
		my @ff = split /\t/, $_;
		next if $ff[0] < 50;
		my ($s1, $c1) = $ff[9]  =~ /Contig(\d+).(\d+)/;
		my ($s2, $c2) = $ff[13] =~ /Contig(\d+).(\d+)/;
		push @{$hh{$ff[9]}}, $_ if($s1 == $s2 and $c1 == $c2 - 1);
	}

	print STDERR "# overlaps:\t", scalar keys %hh, "\n";

	for my $a (keys %hh){
		for(@{$hh{$a}}){
			my @ff = split /\t/, $_;
			next if $ff[10] - $ff[12] > $tip_rate*$ff[0] or $ff[15] > $tip_rate*$ff[0];
			next if $ff[8] eq "-";
			if($self->{overlap}->{$a}){
				my @f2 = split /\t/, $self->{overlap}->{$a};
				$self->{overlap}->{$a} = $_ if $ff[0] > $f2[0];
			}else{
				$self->{overlap}->{$a} = $_;			
			}
		}
	}
}

sub dump_neighbor_alignment{
	my $self = shift;
		
	print STDERR "# proper overlaps:\t", scalar keys %{$self->{overlap}}, "\n";

	open OUT, ">g.neighbor.psl" or die "$! g.neighbor.psl\n";
	for(sort {$self->by_pcap_name} keys %{$self->{overlap}}){
		print OUT $self->{overlap}->{$_};
	}
}

sub build_link{
	my $self = shift;

	$self->prompt("Build Links");
	my ($ps, $pc);
	for(reverse sort {$self->by_pcap_name} keys %{$self->{overlap}}){
		push @{$self->{link}{$_}}, $self->{overlap}{$_};

		my ($s, $c) = $_ =~ /Contig(\d+)\.(\d+)/;
		if($s == $ps and $c == $pc - 1){
			push @{$self->{link}{$_}}, @{$self->{link}{"Contig$ps.$pc"}};
			delete $self->{link}{"Contig$ps.$pc"};
		}
		$ps = $s;
		$pc = $c;
	}
}

sub dump_link{
	my $self = shift;

	open OUT, ">g.link" or die "g.link: $!";
	for(sort {$self->by_pcap_name} keys %{$self->{link}}){
#		print OUT ">$_\n";
		for(@{$self->{link}{$_}}){
			print OUT;
		}
	}
}

sub build_ntrack{
	my $self = shift;

	$self->prompt("Build Neighbor Track");
	for(sort {$self->by_pcap_name} keys %{$self->{link}}){
		my @ll = @{$self->{link}{$_}};
		my @ff = split /\t/, $ll[0];
		$self->{ntrack}{$_} = "$ff[9](0:$ff[12])$ff[13]($ff[16]:";

		for(my $i = 1; $i < @ll; $i++){
			@ff = split /\t/, $ll[$i];
			$self->{ntrack}{$_} .= "$ff[12])$ff[13]($ff[16]:";
		}
		$self->{ntrack}{$_} .= "$ff[14])";
	}
}

sub dump_ntrack{
	my $self = shift;

	open OUT, ">g.ntrack" or die "g.ntrack: $!";
	for(sort {$self->by_pcap_name} keys %{$self->{ntrack}}){
		print OUT "$_\t", $self->{ntrack}{$_}, "\n";
	}
}

sub merge_neighbor{
	my $self = shift;
	$self->prompt("Merge Neighbors");

	for(sort {$self->by_pcap_name} keys %{$self->{ntrack}}){
		my $track = $self->{ntrack}{$_};
		my @ii = $track =~ /(Contig\d+\.\d+\(\d+\:\d+\))/g;
		my $str;
		for(@ii){
			my ($n, $s, $e) = $_ =~ /(Contig\d+\.\d+)\((\d+)\:(\d+)\)/;
			if($e >= $s){
				$str .= substr $self->{merged}{$n}, $s, $e-$s if $e>=$s;
			}else{
				$str = substr $str, 0, length $str - $s + $e;
				print STDERR "$track\n";
			}
		}
		
		my @aa = $track =~ /(Contig\d+\.\d+)\(\d+\:\d+\)/g;
		my $n = shift @aa;
		$self->{merged}{$n} = $str;

		$n =~ s/Contig//;
		my ($s, $c) = $aa[-1] =~ /Contig(\d+)\.(\d+)/;
		$self->{gshrink}{$n} = $self->{gshrink}{"$s.$c"};
		for(@aa){
			delete $self->{merged}{$_};
			$_ = s/Contig//;
			delete $self->{gshrink}{$_};
		}
	}
}

sub rename_contig{
	my $self = shift;

	$self->prompt("Rename Contigs");
	my %hh;
	for(sort {$self->by_pcap_name} keys %{$self->{merged}}){
		my ($s, $c) = $_ =~ /Contig(\d+)\.(\d+)/;
		push @{$hh{$s}}, $c;
	}

	open MRG, ">".$self->merged_file or die "$self->merged_file: $!";
	
	my $i = 0;
	for my $s (sort {$a <=> $b} keys %hh){
		my $j = 1;
		for my $c (@{$hh{$s}}){
			my $new_name = "Contig$i.$j";
			my $old_name = "Contig$s.$c";
			$self->dump_contig_multi_lines(MRG, $new_name, $self->{merged}{$old_name}, 60);
			
			$self->{gap}{$new_name} = $self->{gshrink}{"$s.$c"} if defined $self->{gshrink}{"$s.$c"};

			if(defined $self->{ntrack}{$old_name}){
				$self->{nplan}{$new_name} = $self->{ntrack}{$old_name};
			}else{
				$self->{nplan}{$new_name} = $old_name;
			}
			$j++;
		}
		$i++;
	}
}

sub dump_contig_multi_lines{
	my ($self, $OUT, $name, $str, $len) = @_;
	print $OUT ">$name\n";
	for(my $i = 0; $i < length $str; $i+=$len){
		print $OUT substr($str, $i, $len), "\n";
	}
}


sub dump_gap{
	my $self = shift;
	
	$self->prompt("Dump Gap");
	
	open OUT, ">g.gap" or die "g.gap: $!";

	for(sort {$self->by_pcap_name} keys %{$self->{gap}}){
		print OUT $_, "\t";
		if($self->{gap}{$_} < 10){ 
			print OUT "10\n";
		}else{
			print OUT $self->{gap}{$_}, "\n";
		}
	}
}

sub dump_nplan{
	my $self = shift;
	
	$self->prompt("Dump Neighbor Plan");

	open OUT, ">g.nplan" or die "g.nplan: $!";
	for(sort {$self->by_pcap_name} keys %{$self->{nplan}}){
		print OUT $_, "\t", $self->{nplan}{$_}, "\n";
	}
}


sub get_tail{
	my ($self, $t, $trk) = @_;
	
	@blks = $trk =~ /([+-][TQ]Contig\d+\.\d+)\((\d+\:\d+)\)/g;

	return 0 if @blks == 0;
	
	while($blks[0] ne $t){
		shift @blks;
		shift @blks;
	}
	
	my $contig = substr $blks[0], 2;
	my ($s, $l) = split /:/, $blks[1];
	my $shrink = $s + $l - $self->contig_len->{$contig};
		
	for(my $i = 3; $i < @blks; $i += 2){
		my ($s, $l) = split /:/, $blks[$i];
		$shrink += $l;
	}
	
	$shrink;
}

sub get_head{
	my ($self, $t, $trk) = @_;
	
	@blks = $trk =~ /([+-][TQ]Contig\d+\.\d+)\((\d+\:\d+)\)/g;

	return 0 if @blks == 0;
	
	while($blks[-2] ne $t){
		pop @blks;
		pop @blks;
	}
	
	my $contig = substr $blks[-2], 2;
	my ($s, $l) = split /:/, $blks[-1];
	my $shrink = 0 - $s;
		
	for(my $i = 1; $i < @blks -2; $i += 2){
		my ($s, $l) = split /:/, $blks[$i];
		$shrink += $l;
	}
	
	$shrink;
}


###############
# Read & Dump
###############

sub read_plan{
	my $self = shift;
	
	$self->prompt("Read Plan:\t", $self->plan_file);
	open IN, $self->plan_file or die "$! plan_file\n";
	while(<IN>){
		chomp;
		my ($name, $track) = split;
		my ($s, $c) = $name =~ /Contig(\d+)\.(\d+)/;
		$self->{plan}{$s}[$c-1] = $track;
	}

}

sub dump_plan{
	my $self = shift;

	$self->prompt("Dump Plan:\tg.plan");
	open OUT, ">g.plan" or die "$! g.plan\n";

	for my $s (sort {$a <=> $b} keys %{$self->plan}){
		my $i = 0;
		for my $t (@{$self->plan->{$s}}){
			$i++;
			print OUT "Contig$s.$i\t", $t, "\n";
		}
	}
}


# plan format
# +TContig0.66
# +TContig0.66(0:28055)+QContig0.67(0:868)
# {s}[c-1]=track

# rplan format
# {TContig1.1}=+22.1

sub set_rplan{
	my $self = shift;
	
	$self->prompt("set rPlan");
	
	for my $s (keys %{$self->plan}){
		my $c = 0;
		for my $track (@{$self->plan->{$s}}){
			$c++;
			my @ot = $track =~ /([+-])(TContig\d+\.\d+)/g;
		
			for(my $i = 0; $i < @ot; $i += 2){
 				my ($o, $t) = @ot[$i, $i+1];
				$self->{rplan}{$t} = "$o$s.$c";
			}
		}
	}
}

sub dump_rplan{
	my $self = shift;
	
	$self->prompt("Dump rPlan:\tg.rplan");

	open OUT, ">g.rplan" or die "$! g.rplan\n";
	for(sort {$self->by_pcap_name} keys %{$self->rplan}){
		print OUT "$_\t", $self->rplan->{$_}, "\n";
	}
}

# {s}[c-1]=gap
sub read_gap_t{
	my $self = shift;
	
	$self->prompt("Read Gap:\t", $self->gap_t_file);
	open IN, $self->gap_t_file or die "$! gap_t_file\n";
	while(<IN>){
		chomp;
		my ($name, $gap) = split;
		my ($s, $c) = $name =~ /^Contig(\d+)\.(\d+)$/;
		$self->{gap_t}{$s}[$c-1] = $gap;				# target gap
	}
}

sub dump_gap_t{
	my $self = shift;
	
	$self->prompt("Dump Gap:\tg.gap_t");
	open OUT, ">g.gap_t" or die "$! g.gap_t\n";
	for my $s (sort {$a <=> $b} keys %{$self->gap_t}){
		my $i = 0;
		for my $c (@{$self->gap_t->{$s}}){
			$i++;
			print OUT "Contig$s.$i $c\n" if defined $c;
		}
	}	
}

sub read_contig_len{
	my $self = shift;
	
	$self->prompt("Read Contigs' length:\t", $self->target_file);
	open IN, $self->target_file or die "$! target_file\n";
	
	my $name = "";
	while(<IN>){
		chomp;
		if(/^>/){
			$_ =~ s/^>//; 
			($name) = split /\s/, $_;			
		}else{
			$self->{contig_len}->{$name} += length $_;
		}
	}
}

sub dump_contig_len{
	my $self = shift;
	
	$self->prompt("Dump Contigs' lengths:\tg.contig_len");
	open OUT, ">g.contig_len" or die "$! g.contig_len\n";
	
	for(sort {$self->by_pcap_name} keys %{$self->contig_len}){
		print OUT $_, "\t", $self->contig_len->{$_}, "\n";
	}
}


sub read_merged{
	my $self = shift;
	
	$self->prompt("Read Merged Assembly:\t", $self->{merged_file});
	
	open IN, $self->{merged_file} or die "$! merged\n";
	my $name;
	while(<IN>){
		chomp;
		if(/^>/){ 
			$_ =~ s/^>//; 
			($name) = split /\s/, $_;
		}else{ 
			$self->{merged}{$name} .= $_;
		} 	
	}
	print STDERR "# merged contigs:     \t", scalar keys %{$self->{merged}}, "\n";
}


sub dump_merged{
	my $self = shift; 
	
	open OUT, ">g.merged" or die "$! merged\n";
	for(sort {$self->by_pcap_name} keys %{$self->{merged}}){
		print OUT ">$_\n";
		for(my $i = 0; $i < length $self->{merged}->{$_}; $i += 60){
			print OUT substr($self->{merged}->{$_}, $i, 60), "\n";
		}
	}
}


##########
# Quals
##########

sub read_qual{
	my $self = shift;
	
	$self->prompt("Read Quality File:\t", $self->{qual_file});
	
	open IN, $self->{qual_file} or die "$! qual\n";
	my $name;
	while(<IN>){
		chomp;
		if(/^>/){ 
			$_ =~ s/^>//; 
			($name) = split /\s/, $_;
		}else{ 
			$self->{qual}{$name} .= $_." ";
		} 	
	}
	print STDERR "# merged quals:     \t", scalar keys %{$self->{qual}}, "\n";
}


sub dump_qual{
	my $self = shift; 
	
	open OUT, ">g.qual" or die "$! qual\n";
	for(sort {$self->by_pcap_name} keys %{$self->{qual}}){
		my @qual = split /\s+/, $self->{qual}{$_};
		my $len = @qual;
		print OUT ">$_\t$len\n";
		my $i;
		for($i = 0; $i+59 < $len; $i += 60){
			print OUT join(" ", @qual[$i..$i+59]), "\n";
		}
		print OUT join(" ", @qual[$i..$#qual]), "\n";
	}
}

sub merge_neighbor_qual{
	my $self = shift;
	$self->prompt("Merge Neighbors Quals");

	for(sort {$self->by_pcap_name} keys %{$self->{ntrack}}){
		my $track = $self->{ntrack}{$_};
		my @ii = $track =~ /(Contig\d+\.\d+\(\d+\:\d+\))/g;
		my @newq;
		for(@ii){
			my ($n, $s, $e) = $_ =~ /(Contig\d+\.\d+)\((\d+)\:(\d+)\)/;
			my @qual = split /\s+/, $self->{qual}{$n};
			if($e > $s){
				push @newq, @qual[$s..$e-1];
			}else{
				pop @newq for(1..$s-$e);
			}
		}
		
		my @aa = $track =~ /(Contig\d+\.\d+)\(\d+\:\d+\)/g;
		my $n = shift @aa;
		$self->{qual}{$n} = join(" ", @newq);
		for(@aa){
			delete $self->{qual}{$_};
		}
	}
}

sub rename_qual{
	my $self = shift;

	$self->prompt("Rename Quals");
	open MRG, ">".$self->qual_file or die "$self->qual_file: $!";
	for(sort {$self->by_pcap_name} keys %{$self->{nplan}}){
		my $ori = $self->{nplan}{$_};
		 ($ori) = $ori =~ /(Contig\d+\.\d+)/;
		$self->dump_qual_multi_lines(MRG, $_, $self->{qual}{$ori}, 60);
	}
}

sub dump_qual_multi_lines{
	my ($self, $OUT, $name, $str, $l) = @_;
	my @qual = split /\s+/, $str;
	my $len = @qual;
	
	print $OUT ">$name\t$len\n";
	my $i;
	for($i = 0; $i+$l-1 < $len; $i+=$l){
		print $OUT join(" ", @qual[$i..($i+$l-1)]), "\n";
	}
	print $OUT join(" ", @qual[$i..$#qual]), "\n";
}


##########
# tools
##########

sub reverse_track{
	my ($self, $trk) = @_;

	$trk =~ tr/+-/-+/ ;
	@blks = $trk =~ /([+-][TQ]Contig\d+\.\d+\(\d+\:\d+\))/g;
	@blks = reverse @blks;
	join("", @blks);
}


sub by_pcap_name{
	my $self = shift;

	my ($ai, $af) = $a =~ /^>?[+-]?[TQM]?[Contig]*(\d+)\.(\d+)/;
	my ($bi, $bf) = $b =~ /^>?[+-]?[TQM]?[Contig]*(\d+)\.(\d+)/;
	
	if($ai == $bi){
		$af <=> $bf;
	}else{
		$ai <=> $bi;
	}
}


sub prompt{
	my $self = shift;
	my $date = `date`;
	my ($time) = $date =~ /(\d\d:\d\d:\d\d)/;
	print STDERR join("", @_), "\n$time\n";
}

1;
