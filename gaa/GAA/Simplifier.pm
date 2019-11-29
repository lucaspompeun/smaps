#	class GAA::Simplifier
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


package GAA::Simplifier;
#use strict;

# constructor
sub new{
	my $class = shift;
	my $self = {@_};
	bless $self, $class;

	print STDERR ">>>\n>>> 1. Simplifier\n>>>\n";
	die "Simplifie->new( match_file =>, score_cutoff => 100)\n" 
	unless -e $self->match_file;
	
	$self->score_cutoff(100) unless defined $self->score_cutoff;	#Yao: score_cutoff = 100
	$self->gap_cutoff(50)    unless defined $self->gap_cutoff;	#Yao: gap_cutoff   = 50
	$self->gap_rate(0.05)    unless defined $self->gap_rate;	#Yao: gap_rate     = 0.02
	$self->OVERLAPCF(0.80)   unless defined $self->OVERLAPCF; 	# overlap cutoff to drop a match

	return $self;
}

=h
my $USAGE = <<"END_USAGE";
Usage  :  $! query_vs_target.psl [score_cutoff=100, OVERLAPCF=0.80]
Options:
          blat target query [-ooc=11.ooc] -fastMap query_vs_target.psl
END_USAGE
=cut

#die $USAGE if @ARGV < 1;
#if (scalar(@ARGV) < 1) { print(STDERR, $USAGE); exit(1); }

# even ort is -, qStart and qEnd are the positions in original query, whereas, blockSizes and qStarts are
# from reverse compliment, where qStart = qSize - revQEnd

# finding repeats part can be excluded if use pslReps before this script

#accessor methods
sub match_file{   $_[0]->{match_file}   = $_[1] if defined $_[1]; $_[0]->{match_file}; }
sub OVERLAPCF{    $_[0]->{OVERLAPCF}    = $_[1] if defined $_[1]; $_[0]->{OVERLAPCF}; }
sub score_cutoff{ $_[0]->{score_cutoff} = $_[1] if defined $_[1]; $_[0]->{score_cutoff}; }
sub gap_cutoff{   $_[0]->{gap_cutoff}   = $_[1] if defined $_[1]; $_[0]->{gap_cutoff}; }
sub gap_rate{     $_[0]->{gap_rate}     = $_[1] if defined $_[1]; $_[0]->{gap_rate}; }



sub simplify{
	my $self = shift;
	
	$self->get_unique_match;
	$self->prompt(">>> Done Simplifier <<<");
}


my @hits;        # all hits
my @sortedHits;  # all hits sorted by score
my @oneCon;
my @sortCon;
my ($qName, $qStart, $qEnd, $tName, $tStart, $tEnd, $score, $ort);
my %ctgMarked;
my %bacMarked;
my %lengthOfBac;    # length of contigs and bacs
my %lengthOfCtg;


sub get_unique_match{
	my $self = shift;

	$self->prompt("Read match:\t", $self->match_file);

	open FH, $self->match_file or die $!, $self->match_file, "\n";
	open GAP, ">1.match.gap" or die $!, ">1.match.gap\n";

	<FH>;<FH>;<FH>;<FH>;<FH>;  # remove heading

	my $i = 0;
	while (my $line = <FH>) {
		chomp $line;
		next if (!$line);
		my @arr = split(' ', $line);

		#Yao: discard match with ave(gap_size) > 50=gap_cutoff
		#if( ($arr[6] and $arr[7]/$arr[6] > $self->gap_cutoff) or ($arr[4] and $arr[5]/$arr[4] > $self->gap_cutoff) ){ 
			if( $arr[7] > $self->gap_rate*$arr[0] or $arr[5] > $self->gap_rate*$arr[0]) {
				print GAP "$line\n";
				next;
			}
		#}

		my ($match,$misMatch,$repMatch,$qInsertNum,$tInsertNum) = @arr[0..2, 4, 6];
		($ort,$qName,$qStart,$qEnd,$tName,$tStart,$tEnd) = @arr[8..9, 11..13, 15..16];  
		$lengthOfBac{$qName} = $arr[10];
		$lengthOfCtg{$tName} = $arr[14];
		$score = &GetScore($match, $repMatch, $misMatch, $qInsertNum, $tInsertNum);

		#Yao: discard match with score < score_cutoff=100
		next if $score < $self->score_cutoff;

		push @{$hits[$i]}, ($qName, $qStart, $qEnd, $tName, $tStart, $tEnd, $score, $ort);
		$i++;    
	}
	close (FH);
	my $total = $i;  # total number of matches

# mark all repetitive matches. since repeats tend to collapse in contigs,
# we need to check them in finished seqs


	$self->prompt("Get unique match => 1match.unique");

	# sort hits by contigs
	my @sortedHitsByContig = sort {$a->[3] cmp $b->[3]} @hits;

	# go through each contig and find repeats and mask them
	# find only very similar repeats, for others inclusive or very different repeats
	# filtered out by going through scores

	my $i = 0;
	push (@{$oneCon[$i]}, @{$sortedHitsByContig[$i]});

	my $numCon = 1;
	my $lastName = $oneCon[$i][3];
	for $i (1..$total-1) {
		my $thisName = $sortedHitsByContig[$i][3];
		if ($thisName ne $lastName) {    # find repeat
			if ($numCon > 1) {       # find repeat
				my @sortCon = sort { $a->[4] <=> $b->[4] } @oneCon;     # by starting pos

				# MaskRep($numCon);
				$jj = 0;
				$tName = $sortCon[$jj][3]; # are all the same for contig names
				my ($ptStart, $ptEnd, $pScore) = @{$sortCon[$jj]}[4..6];
				my $pMask = 0;

				for $jj (1..$numCon-1) {
					($tStart, $tEnd, $score) = @{$sortCon[$jj]}[4..6];

					my $over = MIN($ptEnd, $tEnd) - MAX($ptStart, $tStart);
					my $p1 = $over / ($tEnd - $tStart);
					my $p2 = $over / ($ptEnd - $ptStart);
					my $ds = abs ($score - $pScore);

					# if overlap 98% for both matches, and small diff in score
					# dds2 one match score 2; blat 1; very stringent condition
					if ($p1 >= 0.98 && $p2 >= 0.98 && $ds < 10 && $score < 800) {
						if (!$pMask) {
							my ($qName, $qtStart, $qtEnd) = @{$sortCon[$jj-1]}[0..2];
							$bacMarked{$qName} .= " $qtStart-$qtEnd";
							$ctgMarked{$tName} .= " $ptStart-$ptEnd";
						}
						my ($qName, $qtStart, $qtEnd) = @{$sortCon[$jj]}[0..2];
						$bacMarked{$qName} .= " $qtStart-$qtEnd";
						$ctgMarked{$tName} .= " $tStart-$tEnd";
						$pMask = 1;
					} else {
						$pMask = 0;
					}

					($ptStart, $ptEnd, $pScore) = ($tStart, $tEnd, $score);
				}
			}
			@oneCon = ();
			$numCon = 1;
			push (@{$oneCon[$numCon-1]}, @{$sortedHitsByContig[$i]});
			$lastName = $thisName;
		} else {
			push (@{$oneCon[$numCon++]}, @{$sortedHitsByContig[$i]});
		}
	}
	# sort hits by score from high to low
	@sortedHits = sort { $b->[6] <=> $a->[6] } @hits;


	# go through matches from highest score to lowest, find unique matches
	open OUT, ">1match.unique" or die $!, "1match.unique\n";
	#my $index = 0;
	for my $y (0..$total-1) {
		($qName, $qStart, $qEnd, $tName, $tStart, $tEnd, $score, $ort) = @{$sortedHits[$y]};
		my $doneCtg = 0;
		my $doneBac = 0;
		if ($ctgMarked{$tName}) {
			$doneCtg = $self->CheckMatch($tName, $tStart, $tEnd, $ctgMarked{$tName});
		}
		if ($bacMarked{$qName}) {
			$doneBac = $self->CheckMatch($qName, $qStart, $qEnd, $bacMarked{$qName});
		}
		
		if ($doneCtg && $doneBac) {
			next;
		} elsif ($doneCtg) {                         # mark bac part if contig is marked already
			$bacMarked{$qName} .= " $qStart-$qEnd"; # so bac part won't have minor matches later
			next;
		} elsif ($doneBac) {                            # same for contig part
			$ctgMarked{$tName} .= " $tStart-$tEnd";
			next;
		}

		#push @{$uniHits[$index++]}, ($qName, $qStart, $qEnd, $tName, $tStart, $tEnd, $score, $ort);
		print OUT "$qName\t$qStart\t$qEnd\t$tName\t$tStart\t$tEnd\t$score\t$ort\t$lengthOfBac{$qName}\t$lengthOfCtg{$tName}\n";

		$ctgMarked{$tName} .= " $tStart-$tEnd";
		$bacMarked{$qName} .= " $qStart-$qEnd";
	}
	close(OUT);
}




##############
# Subroutines
##############


sub prompt{
	my $self = shift;
	my $date = `date`;
	my ($time) = $date =~ /(\d\d:\d\d:\d\d)/;
	print STDERR join("", @_), "\n$time\n";
}

sub GetScore {
	my ($m, $rep, $mis, $qGap, $tGap) = @_;
	return ($m + ($rep>>1) - $mis - $qGap - $tGap);
}

# check if better matches from the same position exist
sub CheckMatch {   
	my $self = shift;

	my ($name, $start, $end, $marked) = @_;
	my $len = $end - $start;
	my @oneMatch = split(' ', $marked);
	
	for(@oneMatch) {
		my @crr = split('-', $_);         
		my $overlap = MIN($crr[1], $end) - MAX($crr[0], $start);        

		return 1 if $overlap > $len * $self->OVERLAPCF;
	}
	return 0;  # $flag = 1 if already considered a better match
}

sub MIN{ ($_[0] <= $_[1]) ? (return $_[0]) : (return $_[1]); } 
sub MAX{ ($_[0] <= $_[1]) ? (return $_[1]) : (return $_[0]); }

1;


