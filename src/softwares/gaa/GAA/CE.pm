#	class GAA::CE
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


package GAA::CE;

=h
To run gaa.pl with CE option, these two paramenters need to be supplied:
 -pair_status_target <newb.454PairStatus.txt> -pair_status_query <cabog.454PairStatus.txt>

for example, 
gaa.pl -t Newbler.query_contigs.bases -q Cabog.reference_contigs.bases -T Newbler.query_contigs.quals -Q Cabog.reference_contigs.quals -g newbler.gap.txt -pair_status_target newb.454PairStatus.txt -pair_status_query cabog.454PairStatus.txt -o OUTPUT_DIR &


To produce 454PairStatus.txt file from newbler, run gsmapper or prepare the alignment file using other aligners in the format as below.

runMapping contigs.bases.newb *.sff 

Above will generate 454PairStatus.txt.

$ head 454PairStatus.txt
Template	Status	Distance	Left Accno	Left Pos	Left Dir	Right Accno	Right Pos	Right Dir	Left Distance	Right Distance
FC3FKPP01BUVA8	TruePair	3967	Contig18.13	24188	-	Contig18.13	20221	+		
FC3FKPP01BTJNG	TruePair	3820	Contig23.11	44530	+	Contig23.11	48350	-		
FC3FKPP01AFDGW	TruePair	2498	Contig41.19	5103	-	Contig41.19	2605	+		
FC3FKPP01CEJ64	TruePair	2498	Contig41.19	5103	-	Contig41.19	2605	+		
FC3FKPP01E3K7I	TruePair	3735	Contig117.3	98721	-	Contig117.3	94986	+		
FC3FKPP01A7AS8	FalsePair	-	Contig129.5	2438	+	Contig129.6	1986	-	276	1986
FC3FKPP01BJX6Q	FalsePair	-	Contig133.3	28972	+	Contig140.4	1428	-	5803	1428
FC3FKPP01AWH9A	TruePair	2636	Contig40.10	122575	+	Contig40.10	125211	-		
FC3FKPP01DCED4	MultiplyMapped	-	Contig18.9	1244	-	Repeat		

score_cutoff => minimum match cutoff of BLAT alignment.
gap_cutoff   => minimum gap size for CE evaluation.
num_cutoff   => minumum coverage of read pairs for CE evaluation.
=cut

# constructor
sub new {
	my $class = shift;
	my $self = {@_};
        bless $self, $class;

   	print STDERR ">>>\n>>> 0. CE\n>>>\n";
	die "CE->new(
		match_file =>, 
		pair_status_target =>, 
		pair_status_query =>, 
		score_cutoff => 100, 
		gap_cutoff => 50, 
		num_cutoff => 15
	)\n"
  	unless -e $self->match_file or -e $self->pair_status_target or -e $self->pair_status_query;

        $self->score_cutoff(100) unless defined $self->score_cutoff;    
        $self->gap_cutoff(50)    unless defined $self->gap_cutoff;      
        $self->num_cutoff(10)    unless defined $self->num_cutoff;

	return $self;
}

sub match_file{   $_[0]->{match_file}   = $_[1] if defined $_[1]; $_[0]->{match_file}; }
sub pair_status_target{   $_[0]->{pair_status_target}   = $_[1] if defined $_[1]; $_[0]->{pair_status_target}; }
sub pair_status_query{   $_[0]->{pair_status_query}   = $_[1] if defined $_[1]; $_[0]->{pair_status_query}; }
sub score_cutoff{ $_[0]->{score_cutoff} = $_[1] if defined $_[1]; $_[0]->{score_cutoff}; }
sub gap_cutoff{   $_[0]->{gap_cutoff}   = $_[1] if defined $_[1]; $_[0]->{gap_cutoff}; }
sub num_cutoff{   $_[0]->{num_cutoff}   = $_[1] if defined $_[1]; $_[0]->{num_cutoff}; }
sub lib_insertU{$_[0]->{lib_insertU}= $_[1] if defined $_[1]; $_[0]->{lib_insertU};}
sub lib_insert{ $_[0]->{lib_insert} = $_[1] if defined $_[1]; $_[0]->{lib_insert}; }
sub lib_sd{     $_[0]->{lib_sd}     = $_[1] if defined $_[1]; $_[0]->{lib_sd}; }
sub mate_pair_t{ $_[0]->{mate_pair_t} = $_[1] if defined $_[1]; $_[0]->{mate_pair_t}; }
sub mate_pair_q{ $_[0]->{mate_pair_q} = $_[1] if defined $_[1]; $_[0]->{mate_pair_q}; }
sub contig_len_t{ $_[0]->{contig_len_t} = $_[1] if defined $_[1]; $_[0]->{contig_len_t}; }
sub contig_len_q{ $_[0]->{contig_len_q} = $_[1] if defined $_[1]; $_[0]->{contig_len_q}; }

sub ce{
	my $self = shift;

 #	$self->lib_sd(300);
        $self->lib_insertU(8000);

	$self->mate_pair_t($self->pair_status($self->pair_status_target));
	$self->mate_pair_q($self->pair_status($self->pair_status_query));

    	$self->get_ces;
	$self->prompt(">>> Done CE <<<");
}

sub get_ces {
	my $self = shift;
	
        $self->prompt("Read match:\t", $self->match_file);

        open FH, $self->match_file or die $!, $self->match_file, "\n";
        open GAP, ">ce.psl" or die $!, ">ce.psl\n";

 	<FH>;<FH>;<FH>;<FH>;<FH>;  # remove heading

        my $i = 0;
        while (my $line = <FH>) {
                chomp $line;
                next if (!$line);
                my @arr = split(' ', $line);
  		my ($match,$misMatch,$repMatch,$qInsertNum,$qBaseInsert,$tInsertNum,$tBaseInsert) = @arr[0..2, 4..7];

		next if ($match < $self->score_cutoff || $qBaseInsert < $self->gap_cutoff || $tBaseInsert < $self->gap_cutoff);
                my ($ort,$qName,$qSize,$qStart,$qEnd,$tName,$tStart,$tEnd) = @arr[8..13, 15..16];
                my ($blockCount,$blockSizes,$qStarts,$tStarts) = @arr[17..20];	
		next if ($blockCount < 2);
		$self->{contig_len_t}{$tName} = $arr[14] unless defined $self->{contig_len_t}{$tName};
		$self->{contig_len_q}{$qName} = $arr[10] unless defined $self->{contig_len_q}{$qName};

		my @sizes = split(",", $blockSizes);
		my @qs = split(",", $qStarts);
		my @ts = split(",", $tStarts);
		
		my $ce = "\t";
		for my $k (0..$blockCount-2) {
			my $qe = $qs[$k] + $sizes[$k];
			my $te = $ts[$k] + $sizes[$k];
			my $q1 = $qe;
 			my $q2 = $qs[$k+1];
			if ($ort eq '-') {
				$q1 = $qSize - $qs[$k+1];
				$q2 = $qSize - $qe;
			}
			
			my $gapQ = $qs[$k+1] - $qe;
			my $gapT = $ts[$k+1] - $te;
			next if ($gapQ < $self->gap_cutoff && $gapT < $self->gap_cutoff);
	
			my $ce_valueq = $self->get_gap_ce($self->mate_pair_q, $qName, $q1, $q2, $self->{contig_len_q}{$qName});
			my $ce_valuet = $self->get_gap_ce($self->mate_pair_t, $tName, $te, $ts[$k+1], $self->{contig_len_t}{$tName});

			if ($ce_valueq != 0 && $ce_valuet != 0) {
				$ce_valueq = sprintf("%.2f", $ce_valueq);
				$ce_valuet = sprintf("%.2f", $ce_valuet);
				$ce .= ("$qName,".$q1.",".$q2.",".$ce_valueq);
				$ce .= (",$tName,".$te.",".$ts[$k+1].",".$ce_valuet.";");
			}
		}
		if ($ce =~ /Contig/) {
 			print GAP $line, $ce, "\n";
		}
	}
}

sub get_gap_ce {
	my ($self, $h, $name, $pos1, $pos2, $len) = @_;
	my $ce = 0;

	my ($sum, $n) = (0, 0);
	for my $i (0..$len-1) {
		next unless defined $h->{$name}{$i};
		my @ends = split(' ', $h->{$name}{$i});
		foreach my $end (@ends) {
			if ($pos1 >= $i && $pos1 <= $end && $pos2 >= $i && $pos2 <= $end) {
				$sum += $end - $i;	
				$n++;
			}
			last if ($i > $pos1 || $i > $pos2);
		}			
	}
	$ce = ($sum / $n - $self->lib_insert) * sqrt($n) / $self->lib_sd if ($n > $self->num_cutoff);
 	return sprintf("%.2f", $ce);
}

sub pair_status{
        my ($self, $file) = @_;

	my $h;

        $self->prompt("Read 454PairStatus.txt :\t", $file);

        open IN, $file or die "$! pair_status $file\n";
	<IN>;
        my ($sum, $sum2, $n) = (0, 0, 0);
        while(<IN>){
                chomp;
		my ($read_name, $flag, $distance, $l_na, $l_pos, $l_ori, $r_na, $r_pos, $r_ori) = split;
		next if (/MultiplyMapped|Unmapped/ || !defined $r_ori || $distance eq '-');
		#  pairs that have distance below a cutoff say 3*u
		if ($flag eq 'TruePair') {
			$sum += $distance;
     			$sum2 += $distance * $distance;
       			$n++;
		}
		next if ($l_na ne $r_na || $distance > $self->lib_insertU);
		if ($l_pos > $r_pos) {
			($l_pos, $r_pos) = ($r_pos, $l_pos);
		}
		
		$h->{$l_na}{$l_pos} .= " $r_pos";   # attach all right positions to the end, some PE's might have a same location
     	}
   	my $ave = $sum / $n;  # generate library mean and sd
  	my $sd = sqrt(($sum2 - $sum * $ave)/($n - 1));
	$self->lib_insert($ave) unless defined $self->lib_insert;
	$self->lib_sd($sd) unless defined $self->lib_sd;
	print STDERR "lib_insert ", $self->lib_insert, "\n";
	print STDERR "lib_sd ", $self->lib_sd, "\n";
	$h;
}	


sub prompt{
        my $self = shift;
        my $date = `date`;
        my ($time) = $date =~ /(\d\d:\d\d:\d\d)/;
        print STDERR join("", @_), "\n$time\n";
}

1;
