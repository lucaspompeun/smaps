# 	class GAA::Validator
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


package GAA::Validator;

use List::Util qw[min max];

sub new{
	my $class = shift;
	my $self = {@_};
	bless $self, $class;
	
	print STDERR ">>>\n>>> Validator\n>>>\n";
	die "Validator->new(track_file =>, target_file =>)\n" 
	unless (-e $self->track_file or defined $self->track) and -e $self->target_file;

	return $self;
}

sub target_file{$_[0]->{target_file}= $_[1] if defined $_[1]; $_[0]->{target_file}; }
sub track_file{ $_[0]->{track_file} = $_[1] if defined $_[1]; $_[0]->{track_file}; }
sub scaf_degree{$_[0]->{scaf_degree}= $_[1] if defined $_[1]; $_[0]->{scaf_degree}; }
sub t_length{   $_[0]->{t_length}   = $_[1] if defined $_[1]; $_[0]->{t_length}; }
sub track{      $_[0]->{track}      = $_[1] if defined $_[1]; $_[0]->{track}; }
sub track_bad{  $_[0]->{track_bad}  = $_[1] if defined $_[1]; $_[0]->{track_bad}; }
sub break_point{$_[0]->{break_point}= $_[1] if defined $_[1]; $_[0]->{break_point}; }

sub validate{
	my $self = shift;

	$self->set_track unless defined $self->track;
	$self->set_scaf_degree_n_contig_length;

	$self->remove_track_bad;
	$self->rescue_track_bad;

	$self->prompt(">>> Done Validator <<<");		
	$self->track;
}


####################
# Break track
####################

sub remove_track_bad{
	my $self = shift;

	$self->prompt("# merged tracks:\t", scalar keys %{$self->track}, "\nRemedy bad tracks");

	for my $t (keys %{$self->track}){
		my @tt = $t =~ /([+-]TContig\d+\.\d+)/g;
		if(not $self->proper_track(\@tt)){
			$self->{track_bad}{$t} = 1;			#$self->get_track_length($t);
			delete $self->track->{$t};
		}
	}

	$self->prompt("# tracks:  \t", scalar keys %{$self->track}, 
	            "\n# to break:\t", scalar keys %{$self->{track_bad}});

	open  OUT, ">v.track.good" or die "$! v.track.good\n";
	print OUT "$_\t", $self->track->{$_}, "\n" for keys %{$self->track};
	open  OUT, ">v.track.bad" or die "$! v.track.bad\n";
	print OUT "$_\t", $self->track->{$_}, "\n"  for keys %{$self->track_bad};
}


sub rescue_track_bad{
	my $self = shift;
	
	$self->break_track_bad;

	for(keys %{$self->{track_broken}}){
		for(@{$self->{track_broken}{$_}}){
			$self->{track}{$_} = 1;
		}
	}

	my $n_broken_track = 0;
	$n_broken_track += @{$self->{track_broken}{$_}} for keys %{$self->{track_broken}};
	$self->prompt("# broken track:\t$n_broken_track\n# total tracks:\t", scalar keys %{$self->track});
	
	open OUT, ">v.track" or die "$! v.track\n";
	print OUT "$_\n" for(keys %{$self->track});
}
sub break_track_bad{
	my $self = shift;
	
	return unless defined $self->{track_bad};

	for(keys %{$self->track_bad}){
		$self->get_break_point_of($_);
		$self->break_track_bad_of($_);
	}
	
	open OUT, ">v.break_point" or die "$! v.break_point\n";
	for(keys %{$self->{break_point}}){
		print OUT "$_\t", $self->{break_point}{$_}, "\n";
	}
	open OUT, ">v.track.broken" or die "$! v.track.broken\n";
	for(keys %{$self->{track_broken}}){
		print OUT "$_\t", $self->{break_point}{$_}, "\n", join("\n", @{$self->{track_broken}{$_}}), "\n\n";
	}
}
sub get_break_point_of{
	my ($self, $bt) = @_;
	
	my @t = $bt =~ /([+-]TContig\d+\.\d+)/g;
	my @itvl;
	my ($s, $e) = (0, 1);
	for($e = 1; $e < @t; $e++){
		my @tt = @t[$s..$e];
		unless($self->proper_track(\@tt)){
			for(my $k = $e-1; $k >= $s; $k--){
				my @tt = @t[$k..$e];
				unless($self->proper_track(\@tt)){
					push @itvl, $s, $k;
					push @itvl, $k+1, $e-1 if $k<$e-1;	# inclusive
					$s = $e;
					last;
				}	
			}
		}
	}
	push @itvl, $s, $e-1 if $e == @t;					# last boundary is beyond array scope
	$self->{break_point}{$bt} = join(" ", @itvl);
}
sub repair_end{
	my ($self, $track, $ht) = @_;
	
	my ($o, $n, $s, $l) = $track =~ /^([+-])T(Contig\d+\.\d+)\((\d+)\:(\d+)\)/;
	
	if(($ht eq "T" and $o eq "+") or ($ht eq "H" and $o eq "-")){
		$l = $self->t_length->{$n} - $s;
		return $o."T$n($s:$l)";
	}else{
		$l += $s;
		return $o."T$n(0:$l)";
	}
}
sub break_track_bad_of{
	my ($self, $bt) = @_;

	my @itvl = split / /, $self->break_point->{$bt};
	my @t  = $bt =~ /([+-]TContig\d+\.\d+)/g;
	my @tq = $bt =~ /([+-][TQ]Contig\d+\.\d+\(\d+\:\d+\))/g;
	my @trks;
	if(@itvl > 0 and $itvl[0] == 0){					# first match of q
		my @trk;
		my $j = -1;
		for(my $k = -1; $k < $itvl[1];){
			$k++ if $tq[++$j] =~ /^[+-]T/;
		}
		if($j > 0){
			$tq[$j] = $self->repair_end($tq[$j], "T");
			push @trks, join("", @tq[0..$j]);
		}
		shift @itvl; shift @itvl;
	}

	if(@itvl > 0 and $itvl[-1] == $#t){					# last match of q
		my @trk;
		my $j = -1;
		for(my $k = -1; $k < $itvl[-2];){
			$k++ if $tq[++$j] =~ /^[+-]T/;
		}
		if($j < $#tq){
			$tq[$j] = $self->repair_end($tq[$j], "H");
			push @trks, join("", @tq[$j..$#tq]);
		}
		pop @itvl; pop @itvl;
	}
	
	if(@itvl > 0){
		my $k = -1;
		my @trk;
		for(my $j = 0; $j < @tq; $j++){
			$k++ if $tq[$j] =~ /^[+-]T/;
			if($k >= $itvl[0] and $k < $itvl[1]){
				push  @trk, $tq[$j];
			}elsif($k == $itvl[1]){
				push @trk, $tq[$j];
				if($itvl[1] > $itvl[0]){
					$trk[0] = $self->repair_end($trk[0], "H");
					$trk[-1]= $self->repair_end($trk[-1],"T");
					push @trks, join("", @trk);
				}
				shift @itvl;
				shift @itvl;
				@trk = ();
			}
		}
	}
	push @{$self->{track_broken}{$bt}}, @trks;
}


sub proper_track{
	my ($self, $t) = @_; # $t is array ref
	
	return 1 if $self->consecutive_contigs($t) and $self->proper_scaf_adj($t);
	return 0;
}
sub consecutive_contigs{
	my ($self, $t) = @_;

	my %sc;
	for(@$t){
		my ($o, $s, $c) = $_ =~ /([+-])TContig(\d+)\.(\d+)/;
		push @{$sc{$s}}, $c;
	}
	
	for(keys %sc){
		my $cc = $sc{$_};
		next if @$cc == 1;

		my $diff = $cc->[1] - $cc->[0];		# Allow reverse increasing (-5.5-5.6) connection
		return 0 unless abs $diff == 1;

		for(my $i = 1; $i < $#$cc; $i++){
			return 0 if $cc->[$i+1] - $cc->[$i] != $diff;
		}
	}
	1;
}
sub proper_scaf_adj{
	my ($self, $t) = @_;
	
# ort constistent
	my %sco;
	for(@$t){
		my ($o, $s, $c) = $_ =~ /([+-])TContig(\d+)\.(\d+)/;
		if(not defined $sco{$s}{ort}){
			$sco{$s}{ort} = $o;
		}elsif($sco{$s}{ort} ne $o){
			return 0;
		}
		push @{$sco{$s}{contig}}, $c;
	}
	
# sife closed scaf
	for(keys %sco){
		delete $sco{$_} if @{$sco{$_}{contig}} == $self->scaf_degree->[$_];
	}
	
# proper adj
	my %open = %sco;
	if(scalar keys %open < 2){
		return 1;
	}elsif(scalar keys %open > 2){
		return 0;
	}else{
		my ($s1, $s2) = keys %open;
		my ($d1, $d2) = ($self->{scaf_degree}[$s1], $self->{scaf_degree}[$s2]);
		
		my ($o1, $o2) = ($open{$s1}{ort},  $open{$s2}{ort} );
		my ($h1, $t1) = sort {$a <=> $b} ($open{$s1}{contig}[0], $open{$s1}{contig}[-1]);
		my ($h2, $t2) = sort {$a <=> $b} ($open{$s2}{contig}[0], $open{$s2}{contig}[-1]);
		$h1 = ($h1 == 1);
		$t1 = ($t1 == $d1);
		$h2 = ($h2 == 1);
		$t2 = ($t2 == $d2);			
		if( ($h1 and $h2) or ($t1 and $t2) ){		#both heads/tails are in
			return 1 if $o1 ne $o2;
		}elsif( ($h1 and $t2) or ($t1 and $h2) ){	#tail, head
			return 1 if $o1 eq $o2;
		}
		return 0;
	}
}


##################################

sub reverse_track{
	my ($self, $trk) = @_;

	$trk =~ tr/+-/-+/;
	@blks = $trk =~ /([+-][TQ]Contig\d+\.\d+\(\d+\:\d+\))/g;
	@blks = reverse @blks;
	
	join("", @blks);
}

sub get_track_of{
	my ($self, $q_name) = @_;
	
	my $qm = $self->q_match->{$q_name};
	my @itvl = $self->{q_interval}->{$q_name} =~ /(\d+)\-(\d+)\;/g;
	
	for(my $i = 0; $i < @itvl; $i += 2){
		my @m = @$qm[$itvl[$i] .. $itvl[$i+1]-1];

		$self->remove_covered_match(\@m);

		my $track = $self->log_track(\@m);

		if($itvl[$i] == 0){						# first match: keep query head
			my $rl = $self->get_q_relation($qm->[0]);
			if($rl ne "H"){
				my %ff = $self->match_to_fields($qm->[0]);
				$track =~ s/^[+-]TContig\d+\.\d+\(\d+\:\d+\)/\+Q$ff{q_name}\(0\:$ff{q_start}\)/;
			}
		}
		if($itvl[$i+1]-1 == $#$qm){					# last match: keep query tail
			my $rl = $self->get_q_relation($qm->[-1]);
			if($rl ne "T"){
				my %ff = $self->match_to_fields($qm->[-1]);
				my $len = $ff{q_length} - $ff{q_end};
				$track =~ s/[+-]TContig\d+\.\d+\(\d+\:\d+\)$/\+Q$ff{q_name}\($ff{q_end}\:$len\)/;
			}
		}

		$track = $self->tidy_track($track);
		$self->{track}{$track} = 1;
	}
} 
sub remove_covered_match{
	my ($self, $m) = @_;

	return if $#$m < 1;

	my (%c, %p);
	%p = $self->match_to_fields($m->[0]);
	for( my $i = 1; $i < $#$m; $i++){
		%c = $self->match_to_fields($m->[$i]);
		if($c{q_end} <= $p{q_end}){
			print STDERR $c{t_name}, " is fully covered in ", $c{q_name}, "\n";
			splice @$m, $i--, 1;
		}else{ %p = %c; }
	}
}
sub log_track{				# keep ref terminal tips, chop intermediate tips
	my ($self, $m) = @_;

	my $track;
	my ($first, $last) = @$m[0, -1];
	
	my %f = $self->match_to_fields($first);
	if($f{ort} eq "+"){   							# keep target head tip
		$track .= "+T$f{t_name}(0:$f{t_start})";
	}else{
		$track .= "-T$f{t_name}($f{t_end}:" . ($f{t_length}-$f{t_end}) . ")";
	}

	my %pre_f; $pre_f{q_end} = $f{q_start};
	for(@$m){
		%f = $self->match_to_fields($_);

		if($pre_f{q_end} <= $f{q_start}){
			my $gap_q = $f{q_start} - $pre_f{q_end};
			$track .= "+Q$f{q_name}($pre_f{q_end}:" . $gap_q . ")" if $gap_q;
			$track .= "$f{ort}T$f{t_name}($f{t_start}:" . ($f{t_end}-$f{t_start}) . ")"; 
			%pre_f = %f;
		}elsif($pre_f{q_end} < $f{q_end}){
			$track .= "+Q$f{q_name}($pre_f{q_end}:" . ($f{q_end}-$pre_f{q_end}) . ")";
			my $overlap = $self->fields_to_match(%pre_f)."\n".$self->fields_to_match(%f);
			$self->{overlap}{"$f{q_name}($f{q_start}"} = $overlap;
			$track .= "$f{ort}T$f{t_name}($f{t_end}:0)";					#log t_name
			%pre_f = %f;
		}else{ print STDERR $f{t_name}, " fully covered should've been removed ", $f{q_name}, "\n"; }
	}

	if($f{ort} eq "+"){							# keep target tail tip
		$track .= "+T$f{t_name}($f{t_end}:" . ($f{t_length}-$f{t_end}) . ")";
	}else{
		$track .= "-T$f{t_name}(0:$f{t_start})";
	}

	return $track;
}

sub get_q_rel_str{
	my $self = shift;

	$self->prompt("Query Relation String");
	
	for(keys %{$self->q_match}){
		$self->{q_rel_str}{$_} = $self->get_q_rel_str_of($_);
	}
}
sub get_q_rel_str_of{
	my ($self, $q_name) = @_;

	my $str;
	for(@{$self->q_match->{$q_name}}){
		my $rl = $self->get_q_relation($_);
		if($rl eq "<"){
			my %f = $self->match_to_fields($_);
			$rl = 1 if @{$self->t_match->{$f{t_name}}} > 1; 
		}elsif($rl eq "H" or $rl eq "T"){
			my %f = $self->match_to_fields($_);
			my $t_rls = $self->t_rel_str->{$f{t_name}};
			$rl = 2 unless $t_rls =~ /^T?[^HT]*H?$/;	
		}else{
			$rl = 0;
		}
		$str .= $rl;
	}
	$str;
}
sub get_t_rel_str_of{
	my ($self, $t_name) = @_;
	
	my $str;
	for(@{$self->t_match->{$t_name}}){
		$str .= $self->get_t_relation($_);
	}
	$str;
}
sub get_t_rel_str{
	my $self = shift;
	
	$self->prompt("Target Relation String");

	for(keys %{$self->t_match}){
		$self->{t_rel_str}{$_} = $self->get_t_rel_str_of($_);
	}
	
	open OUT, ">2.t_rel_str" or die $!, "2.t_rel_str\n";
	for(sort {$self->by_pcap_name} keys %{$self->t_rel_str}){
		print OUT "$_\t", $self->t_rel_str->{$_}, "\n";
	}
}

sub get_q_interval{
	my $self = shift;
	
	$self->prompt("Get q intervals");
	
	for(keys %{$self->q_match}){
		$self->get_q_interval_of($_);
	}
	
	open OUT, ">2.q_interval" or die $!, "2.q_interval\n";
	for(sort {$self->by_pcap_name} keys %{$self->q_match}){
		print OUT "$_\t", $self->q_rel_str->{$_}, "\t";
		print OUT $self->{q_interval}{$_} if $self->{q_interval}{$_};
		print OUT "\n";
	}
}
sub strip_noisy_ends{
	my ($self, $str) = @_;

	my $start = 0;
	my $end = length $str;

	if( my ($heads) = $str =~ m/^(HH+)/){		# strip multi-heads
		$start = length $heads;
	}
	if( my ($tails) = $str =~ m/(TT+)$/){		# strip multi-tails
		$end -= length $tails;
	}
	
	my @a = split //, $str;
	for(my $i = $start+1; $i < $end; $i++){		# start = last H
		$start = $i if $a[$i] eq "H";
	}
	for(my $i = $end-2; $i > $start; $i--){ 	# end = first T
		$end = $i if $a[$i] eq "T";
	}

	return ($start, $end) if $end - $start > 1;
}
sub get_q_interval_of{
	my ($self, $q_name) = @_;

	my $str = $self->q_rel_str->{$q_name};
	if($str eq "H" or $str eq "T" or $str eq "<"){
		$self->{q_interval}{$q_name} .= "0-1;";		# [0, 3)
		return;
	}

	my ($start, $end) = $self->strip_noisy_ends($str);
	my $qm = $self->q_match->{$q_name};
	my @rls = split //, $str;
	my $s = $start;
	my $e = $s;
	my @a = ();
	for(my $i = $start; $i < $end; $i++){
		push @a, $qm->[$i];
		if($rls[$i] =~ /[HT<]/ and $self->ort_consistent(\@a) and $self->consecutive(\@a) and $self->short_overlap(\@a) ){
			$e++;
		}else{
			unless($self->consecutive(\@a)){	# split at both ends of non-consecutive region
				my @b;
				push @b, pop @a, pop @a;
				my $e2 = $e;
				while($self->consecutive(\@b)){
					$e--;
					push @b, pop @a;
				}
				$self->{q_interval}{$q_name} .= "$e-$e2;" if $e2 - $e > 1;
			}

			$self->{q_interval}{$q_name} .= "$s-$e;" if $e - $s > 1;
			
			if($rls[$i] =~ /[HT<]/){
				$s = $e = $i--;
			}else{
				$s = $e = $i+1;
			}
			@a = ();
		}
	}
	$self->{q_interval}{$q_name} .= "$s-$e;" if $e - $s > 1;		# [0, 3)
}


#####################
# Check Consistency
#####################

sub group_contig_by_scaf{
	my ($self, $a) = @_;
	my %sc;
	for(@$a){
		my @parts = split /\s+/, $_;
		my $t_name = $parts[3];
		my ($s, $c) = $t_name =~ /(\d+)\.(\d+)/;
		push @{$sc{$s}}, $c;
	}
	%sc;
}
sub consecutive{
	my ($self, $a) = @_;
	
	my %sc = $self->group_contig_by_scaf($a);
	for(keys %sc){
		my $cc = $sc{$_};
		next if @$cc == 1;
		my $diff = $cc->[1] - $cc->[0];
		return 0 unless abs $diff == 1;

		for(my $i = 1; $i < $#$cc; $i++){
			return 0 if $cc->[$i+1] - $cc->[$i] != $diff;
		}
	}
	1;
}

sub ort_consistent{
	my ($self, $a) = @_;

	my %ort;
	for(@$a){
		my @parts = split /\s+/, $_;
		if(not defined $ort{$s}){
			$ort{$s} = $ps[7];
		}elsif($ort{$s} ne $ps[7]){
			return 0;
		}
	}
	1;
}

sub short_overlap{
	my ($self, $a) = @_;
	
	return 1 if @$a < 2;
	
	my %p = $self->match_to_fields($a->[-2]);
	my %c = $self->match_to_fields($a->[-1]);
	
	return $p{q_end} - $c{q_start} < $self->max_overlap;
}


##################
# Utils
##################

sub prompt{
	my $self = shift;
	my $date = `date`;
	my ($time) = $date =~ /(\d\d:\d\d:\d\d)/;
	print STDERR join("", @_), "\n$time\n";
}

#Contig225.9	0	121376	Contig17.157	736	122147	121163	+	124466	122179
sub match_to_fields{
	my ($self, $match) = @_;

	my @parts = split /\t/, $match;
	die "match format error:$match\n" if @parts < 10;
	
	my %h;
	($h{q_name}, $h{q_start}, $h{q_end}, $h{t_name}, $h{t_start}, $h{t_end}, $h{score}, $h{ort}, $h{q_length}, $h{t_length}) = @parts[0..9];
	
	%h;
}

sub fields_to_match{
	my ($self, %h) = @_;
	"$h{q_name}\t$h{q_start}\t$h{q_end}\t$h{t_name}\t$h{t_start}\t$h{t_end}\t$h{score}\t$h{ort}\t$h{q_length}\t$h{t_length}";
}

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

sub by_q_start{
	my $self = shift;

	my ($q_start_a) = $a =~ /^>?[+-]?Contig\d+\.\d+\s+(\d+)/;
	my ($q_start_b) = $b =~ /^>?[+-]?Contig\d+\.\d+\s+(\d+)/;
	
	$q_start_a <=> $q_start_b;
}

sub by_t_start{
	my $self = shift;

	my @aa = split /\t/, $a;
	my @bb = split /\t/, $b;
	my $t_start_a = $aa[4];
	my $t_start_b = $bb[4];
	
	$t_start_a <=> $t_start_b;
}

sub get_t_relation{
	my ($self, $match) = @_;

	#reverse q
	my $rl = $self->get_q_relation($match);

	my @fields = split /\t/, $match;
	my $ort = $fields[7];
	$rl =~ tr/HT/TH/ if $ort eq "-";

	return $rl;
}

#target rl query
sub get_q_relation{
	my ($self, $match) = @_;

#	my ($q_name, $q_start, $q_end, $t_name, $t_start, $t_end, $score, $ort, $q_length, $t_length)
	my @fields = split /\t/, $match;
	my ($q_start, $q_end, $t_start, $t_end, $ort, $q_length, $t_length) = @fields[1..2, 4..5, 7..9];
	
	#reverse t
	($t_start, $t_end) = ($t_length - $t_end, $t_length - $t_start) if $ort eq "-";

	my $t_tail = $t_length - $t_end;
	my $q_tail = $q_length - $q_end;
	
	my $tip = $self->tip_size;
	if(($t_start <= $tip or $q_start <= $tip) and ($t_tail <= $tip or $q_tail <= $tip) ){ #proper match
		if($t_start < $q_start){	# "<" is strict, ">" is actually >=
			($t_tail < $q_tail) ? (return "<") : (return "T");
		}else{
			($t_tail < $q_tail) ? (return "H") : (return ">");
		}
	}else{ return 0; } 
}


###########
# Read 
###########

sub set_track{
	my $self = shift;

	$self->prompt("Set Track");

	open IN, $self->track_file or die "$! track_file\n";
	while(<IN>){
		chomp;
		my ($track) = split /\s+/, $_;
		$self->{track}{$track} = 1;
	}
}
sub set_scaf_degree_n_contig_length{
	my $self = shift;

	$self->prompt("Set Scaf Degree and Contig Length");

	open IN, $self->target_file or die "$! target_file\n";
	my $name = "";
	while(<IN>){
		chomp;
		if(/^>/){
			$name = substr $_, 1;
			my ($s, $c) = $_ =~ /Contig(\d+)\.(\d+)/; 
			$self->{scaf}{$s}[$c-1] = "+T$_";	#distinguish old contig from merged contig by if have "Contig"
		}else{
			$self->{t_length}{$name} += length;
		}
	}

	for(keys %{$self->{scaf}}){
		$self->{scaf_degree}[$_] = @{$self->{scaf}{$_}};
	}
	open OUT, ">v.scaffold.degree" or die "$! v.scaffold.degree\n"; 
	print OUT "$_\n" for(@{$self->scaf_degree});
	open OUT, ">v.t_length" or die "$! v.t_length\n";
	for(sort {$self->by_pcap_name} keys %{$self->t_length}){
		print OUT "$_\t", $self->t_length->{$_}, "\n";
	}
}
sub set_qt_match{
	my $self = shift;
	
	$self->prompt("Read Match:\t", $self->match_file);
	
	open IN, $self->match_file or die "$!:match_file\n";
	while(<IN>){
		chomp;
		my @fields = split /\s+/, $_;
		my ($q_name, $t_name) = @fields[0, 3];
		$self->{qt_match}{$q_name}{$t_name} = $_;
	}

	open OUT, ">v.match.qt" or die $!, "v.match.qt\n";
	for(sort {$self->by_pcap_name} keys %{$self->qt_match}){
		my @m;
		for my $t (keys %{$self->qt_match->{$_}}){
			push @m, $self->qt_match->{$_}{$t};
		}
		@m = sort {$self->by_q_start} @m;
		print OUT join("\n", @m), "\n";
	}
}
sub set_overlap{
	my $self = shift;
	
	$self->prompt("Read Overlap:\t", $self->overlap_file);
	
	open IN, $self->overlap_file or die "$! overlap_file\n";
	while(<IN>){
		chomp;
		my $m1 = <IN>;
		my $m2 = <IN>;
		chomp $m1;
		chomp $m2;
		$self->{overlap}{$_}{m1} = $m1;
		$self->{overlap}{$_}{m2} = $m2;
	}
}

1;
