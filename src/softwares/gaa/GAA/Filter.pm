# 	class GAA::Filter
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


package GAA::Filter;
use List::Util qw[min max];

sub new{
	my $class = shift;
	my $self = {@_};
	bless $self, $class;
	
	print STDERR ">>>\n>>> 2. Filter\n>>>\n";
	die "Filter->new( 
        match_file =>, 
        tip_size => 90, 
        max_overlap => 1000)\n" 
    unless -e $self->match_file;
	
	$self->tip_size(90)      unless defined $self->{tip_size};
	$self->max_overlap(1000) unless defined $self->{max_overlap};

	return $self;
}

sub match_file{ $_[0]->{match_file} = $_[1] if defined $_[1]; $_[0]->{match_file}; }
sub tip_size{   $_[0]->{tip_size}   = $_[1] if defined $_[1]; $_[0]->{tip_size}; }
sub max_overlap{$_[0]->{max_overlap}= $_[1] if defined $_[1]; $_[0]->{max_overlap}; }
sub overlap{    $_[0]->{overlap}    = $_[1] if defined $_[1]; $_[0]->{overlap}; }

sub q_match{   $_[0]->{q_match}   = $_[1] if defined $_[1]; $_[0]->{q_match}; }
sub t_match{   $_[0]->{t_match}   = $_[1] if defined $_[1]; $_[0]->{t_match}; }
sub q_rel_str{ $_[0]->{q_rel_str} = $_[1] if defined $_[1]; $_[0]->{q_rel_str}; }
sub t_rel_str{ $_[0]->{t_rel_str} = $_[1] if defined $_[1]; $_[0]->{t_rel_str}; }
sub q_interval{$_[0]->{q_interval}= $_[1] if defined $_[1]; $_[0]->{q_interval}; }
sub track{     $_[0]->{track}     = $_[1] if defined $_[1]; $_[0]->{track}; }
sub track_bad{ $_[0]->{track_bad} = $_[1] if defined $_[1]; $_[0]->{track_bad}; }
sub track_adj{ $_[0]->{track_adj} = $_[1] if defined $_[1]; $_[0]->{track_adj}; }

sub filt{
	my $self = shift;

	$self->set_q_match;
	$self->set_t_match;
	$self->get_t_rel_str;
	$self->get_q_rel_str;
	$self->get_q_interval;	
	$self->get_track;
	$self->get_track_adj;
	$self->merge_tracks;
	$self->prompt(">>> Done Filter <<<");		
	$self->track;
}

sub merge_tracks{
	my $self = shift;
	
	$self->prompt("Linking");
	
	my $adj = $self->track_adj;
	for my $mid (keys %$adj){
		if( scalar keys %{$adj->{$mid}} == 2 and $mid =~ /^T/){ 	# Only T joint
			my ($pre, $post) = keys %{$adj->{$mid}};
			$self->merge_track($pre, $mid, $post);
			$links++;
		}
	}

	$self->prompt("# links:\t\t$links\n# track after merge:\t", scalar keys %{$self->track});
	open OUT, ">2track" or die "$! 2track\n";
	for(keys %{$self->track}){
		print OUT $_, "\n";
	}
}

sub merge_track{
	my ($self, $pre, $mid, $post) = @_;
	
	my $adj = $self->track_adj;
	my $t1 = $adj->{$pre}{$mid};
	my $t2 = $adj->{$mid}{$post};
	my @blks1 = $t1 =~ /([+-][TQ]Contig\d+\.\d+\(\d+:\d+\))/g;
	my @blks2 = $t2 =~ /([+-][TQ]Contig\d+\.\d+\(\d+:\d+\))/g;
	my $blk1 = pop @blks1;
	my $blk2 = shift @blks2;
	my ($o1, $n1, $s1, $l1) = $blk1 =~ /^([+-])([TQ]Contig\d+\.\d+)\((\d+):(\d+)\)$/g;
	my ($o2, $n2, $s2, $l2) = $blk2 =~ /^([+-])([TQ]Contig\d+\.\d+)\((\d+):(\d+)\)$/g;

	if($n1 eq $n2 and $o1 eq $o2){
		my ($o, $n, $s, $l) = ($o1, $n1);
		($o eq "+") ? ($s = $s1, $l = $l2 - $s1) : ($s = $s2, $l = $l1 - $s2);
		my $t = join("", @blks1)."$o$n($s:$l)".join("", @blks2);
		
		$adj->{$pre}{$post} = $t;
		$adj->{$post}{$pre} = $self->reverse_track($t);
		delete $adj->{$pre}{$mid};
		delete $adj->{$mid}{$pre};
		delete $adj->{$mid}{$post};
		delete $adj->{$post}{$mid};
		
		$self->{track}{$t} = 1;
		delete $self->track->{$t1};
		delete $self->track->{$t2};
		delete $self->track->{$self->reverse_track($t1)};
		delete $self->track->{$self->reverse_track($t2)};
	}else{ print STDERR "Can't merge\n$t1\n$t2\n";}
}

sub tidy_track{
	my ($self, $trk) = @_;

	@blks = $trk =~ /([+-][TQ]Contig\d+\.\d+)\((\d+):(\d+)\)/g;
	my ($pn, $ps, $pl) = @blks[0..2];
	for(my $i = 3; $i < @blks; $i += 3){
		my ($n, $s, $l) = @blks[$i..$i+2];
		if($n eq $pn){
			if($ps+$pl==$s or $s+$l==$ps){
				$ps = $blks[$i+1] = min($s, $ps);
				$pl = $blks[$i+2] = $l + $pl;
				splice @blks, $i-3, 3;
				$i -= 3;
=h			}elsif($pl == 0){ splice @blks, $i-3, 3; $i -= 3;	# consecutive 0-length 			
			}elsif($l  == 0){ splice @blks, $i, 3; $i -= 3;
=cut
			}else{ print STDERR "What's this! $trk\n"; }
		}else{
			($pn, $ps, $pl) = ($n, $s, $l);
		}
	}
	
	for(my $i = 0; $i < @blks; $i += 3){
		$blks[$i] .= "(";
		$blks[$i+1] .= ":";
		$blks[$i+2] .= ")";
	}
	
	join("", @blks);
}

sub get_track_adj{
	my $self = shift;

	$self->prompt("# track before merge:\t", scalar keys %{$self->track});

	for my $track (keys %{$self->track}){			# Only TContig has neighbor
		my ($first) = $track =~ /^[+-]([TQ]Contig\d+\.\d+)\(\d+:\d+\)/;
		my ($last)  = $track =~ /[+-]([TQ]Contig\d+\.\d+)\(\d+:\d+\)$/;
		
print STDERR "first=$first -> last=$last, track=$track already defined adj\n" 
if defined $self->{track_adj}{$first}{$last} or defined $self->{track_adj}{$last}{$first};
next unless $first and $last;

		$self->{track_adj}{$first}{$last} = $track;
		$self->{track_adj}{$last}{$first} = $self->reverse_track($track);
	}
	
	open OUT, ">2.track_adj" or die "$! 2.track_adj\n";
	my $adj = $self->track_adj;
	for my $mid (keys %$adj){
		if( scalar keys %{$adj->{$mid}} == 2 ){
			my ($pre, $post) = keys %{$adj->{$mid}};
			print OUT $adj->{$pre}{$mid}, "\n";
			print OUT $adj->{$mid}{$post}, "\n\n";
		}
	}
}

sub reverse_track{
	my ($self, $trk) = @_;

	$trk =~ tr/+-/-+/;
	@blks = $trk =~ /([+-][TQ]Contig\d+\.\d+\(\d+\:\d+\))/g;
	@blks = reverse @blks;
	
	join("", @blks);
}

sub get_track{
	my $self = shift;
	
	$self->prompt("Tracing");

	for(keys %{$self->q_match}){
		$self->get_track_of($_);
	}
	
	open OUT, ">2.track.by_q" or die $!, "2.track.by_q\n";
	for(keys %{$self->{track}}){
		print OUT "$_\t", $self->track->{$_}, "\n";
	}
	
	open OUT, ">2.overlap" or die "$! 2.overlap\n";
	for(sort {$self->by_pcap_name} keys %{$self->{overlap}}){
		print OUT $_, "\n", $self->overlap->{$_}, "\n";
	}
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
				my %ff = $self->match_to_fields($qm->[0]);
				my %ff = $self->match_to_fields($qm->[-1]);
				my $len = $ff{q_length} - $ff{q_end};
				$track =~ s/[+-]TContig\d+\.\d+\(\d+\:\d+\)$/\+Q$ff{q_name}\($ff{q_end}\:$len\)/;
			}
		}

		$track = $self->tidy_track($track);
		$self->{track}{$track} = "+Q$q_name";
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
	if($f{ort} eq "+"){  					# keep target head tip, +Q
		my $len = $f{t_start};
		$track .= "+T$f{t_name}(0:$len)";		# if $len > 0;
	}else{
		my $len = $f{t_length}-$f{t_end};
		$track .= "-T$f{t_name}($f{t_end}:$len)";	# if $len > 0;
	}

	my %pre_f;# $pre_f{q_end} = -1;
	for(@$m){
		%f = $self->match_to_fields($_);

		if(not %pre_f){		#first match of q; #$pre_f{q_end} == -1){
			my $len_t = $f{t_end} - $f{t_start};
			$track .= "$f{ort}T$f{t_name}($f{t_start}:$len_t)";
			 
			%pre_f = %f;
		}elsif($pre_f{q_end} <= $f{q_start}){
			my $len_q = $f{q_start} - $pre_f{q_end};
			$track .= "+Q$f{q_name}($pre_f{q_end}:$len_q)";		# if $len_q;
		
			my $len_t = $f{t_end} - $f{t_start};
			$track .= "$f{ort}T$f{t_name}($f{t_start}:$len_t)";
			 
			%pre_f = %f;
		}elsif($pre_f{q_end} < $f{q_end}){
			my $len_q = $f{q_end}-$pre_f{q_end};
			$track .= "+Q$f{q_name}($pre_f{q_end}:$len_q)";
			
			my $overlap = $self->fields_to_match(%pre_f)."\n".$self->fields_to_match(%f);
			$self->{overlap}{"$f{q_name}($pre_f{q_end}"} = $overlap;
			
			my $start;
			($f{ort} eq "+") ? ($start = $f{t_end}) : ($start = $f{t_start});
			$track .= "$f{ort}T$f{t_name}($start:0)";		#log t_name
			%pre_f = %f;
		}else{ print STDERR $f{t_name}, " fully covered should've been removed ", $f{q_name}, "\n"; }
	}

	if($f{ort} eq "+"){							# keep target tail tip
		my $len = $f{t_length}-$f{t_end};
		$track .= "+T$f{t_name}($f{t_end}:$len)"; 	# if $len > 0; # replace it with q_tip if rl != T
	}else{
		my $len = $f{t_start};
		$track .= "-T$f{t_name}(0:$len)";		# if $len > 0;
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
#print STDERR "$q_name\t$str\n" unless $str =~ /^H*[^HT]*T*$/;
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
	for(sort {$self->by_pcap_name} keys %{$self->{t_rel_str}}){
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
				unshift @b, pop @a;
				unshift @b, pop @a;
				my $e2 = $e;
				while($self->consecutive(\@b)){
					$e--;
					unshift @b, pop @a;
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

sub get_ort{
	my ($self, $a) = @_;

	my %ort;
	for(@$a){
		my @parts = split /\s+/, $_;
		my $t_name = $parts[3];
		my ($s, $c) = $t_name =~ /(\d+)\.(\d+)/;
		$ort{$s} = $parts[7] unless $ort{$s};
	}
	%ort;
}

sub consecutive{
	my ($self, $a) = @_;
	
	my %sc = $self->group_contig_by_scaf($a);
	my %so = $self->get_ort($a);
	for(keys %sc){
		my $cc = $sc{$_};
		next if @$cc == 1;
		my $diff = $cc->[1] - $cc->[0];
		return 0 unless (($diff == 1 and $so{$_} eq "+") or ($diff == -1 and $so{$_} eq "-"));

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


###################
# Read Match
###################

sub set_q_match{
	my $self = shift;
	
	$self->prompt("Query Match:\t", $self->match_file);
	
	open IN, $self->match_file or die "$!:match_file\n";
	while(<IN>){
		chomp;
		my @fields = split /\s+/, $_;
		my $q_name = $fields[0];
		push @{$self->{q_match}{$q_name}}, $_;
	}
	
	for my $q_name (keys %{$self->{q_match}}){
		for(@{$self->q_match->{$q_name}}){
			if($self->get_q_relation($_) eq ">"){
				delete $self->q_match->{$q_name};
				last;
			}
		}
	}
	
	my %qt_match;
	for my $q_name (keys %{$self->q_match}){
		for(@{$self->q_match->{$q_name}}){
			my @fields = split /\s+/, $_;
			my ($q_name, $t_name) = @fields[0, 3];
			if($qt_match{$q_name}{$t_name}){
				delete $self->q_match->{$q_name};
				last;
			}else{
				$qt_match{$q_name}{$t_name} = 1;
			}
		}
	}

	open OUT, ">2.match.qi" or die $!, "2.match.qi\n";
	for(sort {$self->by_pcap_name} keys %{$self->q_match}){
		@{$self->q_match->{$_}} = sort {$self->by_q_start} @{$self->q_match->{$_}};
		print OUT join("\n", @{$self->q_match->{$_}}), "\n";
	}
}

sub set_t_match{
	my $self = shift;

	for my $q_name (keys %{$self->q_match}){
		for(@{$self->q_match->{$q_name}}){
			my @fields = split /\s+/, $_;
			my $t_name = $fields[3];
			push @{$self->{t_match}{$t_name}}, $_;
		}
	}

	open OUT, ">2.match.ti" or die $!, "2.match.ti\n";
	for(sort {$self->by_pcap_name} keys %{$self->{t_match}}){
		@{$self->t_match->{$_}} = sort {$self->by_t_start} @{$self->t_match->{$_}};
		print OUT join("\n", @{$self->t_match->{$_}}), "\n";
	}
}


##############
# Other
##############

sub print_adj{
	my $self = shift;

	print "Track Adjacency >>>>\n";
	my $adj = $self->track_adj;
	for my $f (keys %$adj){
		for my $s (keys %{$adj->{$f}}){
			print "$f\t$s\t", $adj->{$f}{$s}, "\n";
		}
	}
	print "<<<<\n";
}

sub get_link_track{
	my $self = shift;
	
	for(keys %{$self->track}){
		my @blks = $_ =~ /([+-]TContig\d+\.\d+)\(\d+:\d+\)/g;
		my $t_link = join("", @blks);
		$t_link =~ s/T//g;
		$self->{link_track}{$t_link} = $_;
	}
	open OUT, ">2.t_link_track" or die $!, "2.t_link_track\n";
	print OUT "$_\t", $self->link_track->{$_}, "\n" for(keys %{$self->link_track});
}

1;
