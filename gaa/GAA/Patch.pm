# 	class GAA::Patch
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


package GAA::Patch;

# constructor
sub new {
	my $class = shift;
	my $self = {@_};
        bless $self, $class;

   	print STDERR ">>>\n>>> Plan Patching\n>>>\n";
	die "Patch->new(plan_file =>, ce_file =>)\n"
  	unless -e $self->plan_file and -e $self->ce_file;
  	
	return $self;
}

# accessor methods
sub plan_file{ $_[0]->{plan_file} = $_[1] if defined $_[1]; $_[0]->{plan_file}; }
sub ce_file{   $_[0]->{ce_file}   = $_[1] if defined $_[1]; $_[0]->{ce_file}; }
sub ce_cutoff{ $_[0]->{ce_cutoff} = $_[1] if defined $_[1]; $_[0]->{ce_cutoff}; }

sub patch{
	my $self = shift;
	
	$self->ce_cutoff(3.3);
	
	$self->read_plan;
	$self->read_ce;
	$self->fix;
	rename "3floor_plan", "3floor_plan.befor_pair";
	$self->dump_plan;	
}

sub fix{
	my $self = shift;
	$self->prompt("Fixing");

	for my $t (keys %{$self->{ce}}){
		for my $q (keys %{$self->{ce}{$t}}){
			my $patch = $self->{ce}{$t}{$q};
			my @tmp   = split /\t/, $patch;
			   $patch = $tmp[-1];
			my @pp    = split /;/, $patch;
			for(@pp){
				$self->fix_a($_, $tmp[8], $tmp[10], $tmp[14]);	# Todo
				
			}
		}
	}
}

sub fix_a{
	my ($self, $patch, $ort, $q_length, $t_length) = @_;
	
	my ($q_name, $q_start, $q_end, $q_ce, $t_name, $t_start, $t_end, $t_ce) = split /,/, $patch;
	
	return if abs($t_ce) < $self->ce_cutoff or abs($t_ce) < abs($q_ce)+1;
	return if $q_ce > $t_ce and $q_end - $q_start < $t_end - $t_start;
	return if $q_ce < $t_ce and $q_end - $q_start > $t_end - $t_start;
	print "Fixing ($q_name, $q_start, $q_end, $q_ce, $t_name, $t_start, $t_end, $t_ce, $ort, $q_length, $t_length)\n";
	
	my $m_name = $self->{rplan}{$t_name};
	my $plan   = $self->{plan}{$m_name};
	if($plan =~ /^[+-]TContig\d+\.\d+$/){
		my ($t_ort, $t) = split /T/, $plan;
		my $q_ort = $self->add($t_ort, $ort);
		$plan  = $t_ort."T$t_name(0:$t_start)";
		$plan .= $q_ort."Q$q_name($q_start:". ($q_end-$q_start) . ")";
		$plan .= $t_ort."T$t_name($t_end:" . ($t_length-$t_end) . ")";
	}elsif($plan =~ /[+-]TContig\d+\.\d+\(\d/){
		my $i = $self->test_patch($plan, $t_name, $q_name, $t_start, $t_end);
		return if $i < 0;
		
		my @tq = $plan =~ /([+-][TQ]Contig\d+\.\d+\(\d+\:\d+\))/g;
		my @ti = $tq[$i] =~ /([+-])([TQ])(Contig\d+\.\d+)\((\d+)\:(\d+)\)/;
		my $t_ort = $ti[0];
		my $q_ort = $self->add($t_ort, $ort);
		
		$plan = "";
		$plan .= join("", @tq[0..$i-1]) if $i > 0;
		$plan .= $t_ort."T$t_name($ti[3]:" .($t_start-$ti[3]). ")";
		$plan .= $q_ort."Q$q_name($q_start:". ($q_end-$q_start) . ")";
		$plan .= $t_ort."T$t_name($t_end:" .($ti[4]-$t_end). ")";
		$plan .= join("", @tq[$i+1..$#tq]) if $i < $#tq;
	}
	
	$self->{plan}{$m_name} = $plan;						# make the change
}

sub test_patch{
	my ($self, $plan, $t_name, $q_name, $t_start, $t_end) = @_;
	
	my @tq = $plan =~ /([+-])([TQ])(Contig\d+\.\d+)\((\d+)\:(\d+)\)/g;
	for(my $i = 0; $i < @tq; $i+= 5){
		if($tq[$i+1] eq "T" and $tq[$i+2] eq $t_name){
			if($tq[$i+3] < $t_start and $t_end < $tq[$i+3]+$tq[$i+4]){
				if( ($i > 4 and $tq[$i-4] eq "Q" and $tq[$i-3] eq $q_name) or 
				    ($i < @tq and $tq[$i+6] eq "Q" and $tq[$i+7] eq $q_name) ){
					return $i/5;    
				}
			}
		}
	}
	-1;
}


sub read_plan{
	my $self = shift;
	$self->prompt("Reading ", $self->plan_file);
	
	open IN, $self->plan_file or die "plan_file: $!";
	while(<IN>){
		chomp;
		my ($c, $t) = split /\t/, $_;
		$self->{plan}{$c}  = $t;
		
		my @tt = $t =~ /T(Contig\d+\.\d+)/g;
		for(@tt){
			$self->{rplan}{$_} = $c;
		}
	}
}

sub dump_plan{
	my $self = shift;
	$self->prompt("dump new plan");
	
	open OUT, ">3floor_plan" or die $!;
	for(sort {$self->by_pcap_name} keys %{$self->{plan}}){
		print OUT $_, "\t", $self->{plan}{$_}, "\n";
	}
}

sub dump_rplan{
	my $self = shift;
	$self->prompt("dump rplan");
	
	open OUT, ">p.rplan" or die $!;
	for(sort {$self->by_pcap_name} keys %{$self->{rplan}}){
		print OUT $self->{rplan}{$_}, "\t", $_, "\n";
	}
}

sub read_ce{
	my $self = shift;
	$self->prompt("Parsing ", $self->ce_file);
	
	open IN, $self->ce_file or die "ce_file: $!";
	while(<IN>){
		chomp;
		next unless $self->reliable($_);
		
		my @ff = split /\t/, $_;
		my ($q_name, $t_name) = ($ff[9], $ff[13]);
		$self->{ce}{$q_name}{$t_name} = $_;
	}

}

sub dump_ce{
	my $self = shift;
	
	open OUT, ">p.reliable" or die "p.reliable: $!";
	for my $q (sort {$self->by_pcap_name} keys %{$self->{ce}}){
		for my $t (sort {$self->by_pcap_name} keys %{$self->{ce}{$q}}){
			print OUT $self->{ce}{$q}{$t}, "\n";
		}
	}	
}

sub reliable{
	my ($self, $m) = @_;
	my @ff = split /\t/, $m;
	
	return 0 if $ff[0] < 1000;
	return 0 if $ff[5] > 1000 or $ff[7] > 1000;
	return 0 if $ff[5] * 2 > $ff[0] or $ff[7] * 2 > $ff[0];
	return 0 unless $self->get_q_relation_psl($m, 100);
	
	1;
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

sub prompt{
        my $self = shift;
        my $date = `date`;
        my ($time) = $date =~ /(\d\d:\d\d:\d\d)/;
        print STDERR join("", @_), "\n$time\n";
}

sub get_q_relation_psl{
	my ($self, $match, $tip) = @_;

	my @fields = split /\t/, $match;
	my ($q_length, $q_start, $q_end, $t_length, $t_start, $t_end, $ort) = @fields[10..12, 14..16, 8];
	
	#reverse t
	($t_start, $t_end) = ($t_length - $t_end, $t_length - $t_start) if $ort eq "-";

	my $t_tail = $t_length - $t_end;
	my $q_tail = $q_length - $q_end;
	
	if(($t_start <= $tip or $q_start <= $tip) and ($t_tail <= $tip or $q_tail <= $tip) ){ #proper match
		if($t_start < $q_start){	# "<" is strict, ">" is actually >=
			($t_tail < $q_tail) ? (return "<") : (return "T");
		}else{
			($t_tail < $q_tail) ? (return "H") : (return ">");
		}
	}else{ return 0; } 
}

sub add{
	my ($self, $a, $b) = @_;	
	return "+" if( ($a ne "-" and $b ne "-") or ($a eq "-" and $b eq "-") );
	return "-";
}

1;

