#!/usr/bin/perl

use strict;
use warnings;

use Sort::Array qw(Sort_Table Discard_Duplicates);
use Bio::SeqIO;

###################################################################################################################################
#   G L O B A L    V A R I A B L E S
###################################################################################################################################

	my (
		$aln_file,				# file containing the input alignment
		$rate_file,				# file containing the rates
		$keep_frac,				# fraction of sites to keep
		$base_pairing_file,		# file containing which bases are paired (for rRNA stems)
		$base_pairing,			# 1|0: whether or not base pair information is present
		$partition_file,		# file containing the data partition definitions
		$partitioning,			# 1|0: whether or not partitions are present in the data
		$output_file,			# file to contain the output alignment
		$output_format,			# file format of the output alignment
		$random_sets,			# number of random alignments to be generated

		$aln_in,				# contains the input alignment
		$taxa,					# contains list of taxa in input order
		$input_seq_length,		# alignment length
		$aln_out,				# contains the output alignment
	);
	my $default_output_format = 'nexus';


###################################################################################################################################
#   P A R S E    C O M M A N D    L I N E    O P T I O N S
###################################################################################################################################

	unless (($ARGV[0]) and (substr($ARGV[0],0,1) eq "-") and
			($ARGV[1]) and
			($ARGV[2]) and (substr($ARGV[2],0,1) eq "-") and
			($ARGV[3]) and
			($ARGV[4]) and (substr($ARGV[4],0,1) eq "-") and
			($ARGV[5]) and
			($ARGV[6]) and (substr($ARGV[6],0,1) eq "-") and
			($ARGV[7]))
		{usage()}
	
	for (my $i=0; $i<scalar(@ARGV); $i+=2) {
		if ($ARGV[$i] eq "-a") {$aln_file=$ARGV[$i+1]}
		elsif ($ARGV[$i] eq "-r") {$rate_file=$ARGV[$i+1]}
		elsif ($ARGV[$i] eq "-f") {$keep_frac=$ARGV[$i+1]}
		elsif ($ARGV[$i] eq "-pa") {$base_pairing=1; $base_pairing_file=$ARGV[$i+1]}
		elsif ($ARGV[$i] eq "-pt") {$partitioning=1; $partition_file=$ARGV[$i+1]}
		elsif ($ARGV[$i] eq "-o") {$output_file=$ARGV[$i+1]}
		elsif ($ARGV[$i] eq "-of") {$output_format=$ARGV[$i+1]}
		elsif ($ARGV[$i] eq "-ra") {$random_sets=$ARGV[$i+1]}
		else {usage()}
	}
	unless ($aln_file and $keep_frac and $output_file) {usage()}
	unless (-e $aln_file) {die "error -- input file $aln_file does not exist\n"}
	unless ($rate_file or $random_sets) {usage()}
	if ($rate_file) {unless (-e $rate_file) {die "error -- input file $rate_file does not exist\n"}}
	if ($random_sets) {unless ($random_sets =~ /^\d+$/) {die "error -- number of random alignments (-ra) must be integer\n"}}
	unless (($keep_frac =~ /^[\d\.]+$/) and ($keep_frac <= 1) and ($keep_frac > 0)) {die "error -- fraction of sites to keep must be decimal number between 0 and 1\n"}
	unless ($base_pairing) {$base_pairing = 0}
	unless ($partitioning) {$partitioning = 0}
	unless ($output_format) {$output_format = $default_output_format}

	sub usage {			# prints warning about incorrect command line usage
		print "\nusage:\n"; 
		print "\nmandatory parameters\n";
		print "   -a    alignment_file  (must be in fasta format)\n";
		print "   -r    rates_file  (must be in the HyPhy output format)\n";
		print "           if you wish to remove random sites (-ra option), this is not mandatory\n";
		print "   -f    fraction of characters to keep\n";
		print "   -o    output alignment  (in nexus format)\n";
		print "\noptional parameters\n";
		print "   -pa   file with list of paired bases\n";
		print "   -pt   file with definitions of data partitions  (one per line)\n";
		print "   -of   output format: fasta|nexus  (default: ",$default_output_format,")\n";
		print "   -ra   number of alignments with random sites stripped to be generated\n";
		print "\n";
		exit;
	}


###################################################################################################################################
#   A C Q U I R E    I N P U T    D A T A
###################################################################################################################################

	$aln_in = parse_alignment($aln_file);
	$input_seq_length = scalar(keys(%$aln_in));
	print "alignment read from $aln_file\n  $input_seq_length sites\n  ",scalar(@$taxa)," taxa\n";
	if ($rate_file) {
		parse_rates_file($rate_file);
		my $rates_found = 0; foreach my $char (keys %$aln_in) {if (defined $aln_in->{$char}->{rate}) {++$rates_found}}
		print "rates read from $rate_file\n  $rates_found rate values read\n";
	}
	if ($base_pairing) {parse_base_pairings_file($base_pairing_file)}
	if ($partitioning) {parse_data_partitions_file($partition_file)}

###################################################################################################################################
#   O R G A N I Z E    S I T E S    F O R    O U T P U T    A N D    C A L L    R O U T I N E    T O    S A V E    O U T P U T
###################################################################################################################################

	print "saving output to $output_file\n";
	if ($rate_file) {
		my $list;
		foreach my $char (keys(%$aln_in)) {
			push @$list, $aln_in->{$char}->{'rate'}.','.$char;
		}
		sort_according_to_list($list);
		save_output($output_file);
	} else {
		for (my $rep = 1; $rep <= $random_sets; ++$rep) {
			my $list;
			foreach my $char (keys(%$aln_in)) {
				push @$list, rand(1).','.$char;    # this will lead to the sites being sorted according to random rates (i.e. random sites will be removed)
			}
			sort_according_to_list($list);
			save_output($output_file.'.'.get_zeros(length($random_sets)-length($rep)).$rep);
		}
	}


###################################################################################################################################
#   subroutine  ---  S O R T    A C C O R D I N G    T O    L I S T
###################################################################################################################################

	sub sort_according_to_list {
		my $in = shift;
		my $list; $list = sort_list_ascending($in);
		my $reduced_hash;
		my $continue; $continue=1;
		for (my $i=0; $i < scalar(@$list); ++$i) {
			if ($continue) {
				$list->[$i] =~ /.+?,(.+)/;
				my $char = $1;
				if ($aln_in->{$char}->{'linked'}) {
					if ( ( ( scalar(keys(%$reduced_hash)) + 2 ) / scalar(@$list) ) < $keep_frac ) {
						$reduced_hash->{$char} = 1;
						$reduced_hash->{$aln_in->{$char}->{'linked'}} = 1;
					} else {
						$continue=0;
					}
				} else {
					if ( ( ( scalar(keys(%$reduced_hash)) + 1 ) / scalar(@$list) ) < $keep_frac ) {
						$reduced_hash->{$char} = 1;
					} else {
						$continue=0;
					}
				}
			} else {
				$i = scalar(@$list)
			}
		}
		my $counter; $counter = 0;
		foreach my $char (sort(keys(%$reduced_hash))) {
			++$counter;
			$aln_out->[$counter] = $char;
			$aln_in->{$char}->{'output_pos'} = $counter;
		}	
	}

###################################################################################################################################
#   subroutine  ---  S A V E    O U T P U T
###################################################################################################################################

	sub save_output {
		my $file = shift;
		open(FH,'>'.$file) || die "error -- unable to write to $file\n";
		my $max_ID_length;
		if ($output_format eq 'nexus') {
			# print nexus header
				print FH "#NEXUS\n[TITLE: site stripping: keeping ",($keep_frac*100),"\% ";
				if ($rate_file) {print FH "slowest"} else {print FH "randomly picked"}
				print FH" sites of $aln_file]\n\n";
				print FH "begin data;\ndimensions ntax=",scalar(@$taxa)," nchar=",(scalar(@$aln_out)-1),";\nformat interleave datatype=DNA missing=N gap=-;\n\nmatrix\n";
			# calculate maximum sequence ID length
				$max_ID_length = 0;
				foreach my $taxon (@$taxa) {
					if (length($taxon) > $max_ID_length) {$max_ID_length = length($taxon);}
				}
		}
		for (my $taxon_nr = 0; $taxon_nr < scalar(@$taxa); ++$taxon_nr) {
			if ($output_format eq 'fasta') {print FH '>',$taxa->[$taxon_nr],"\n";}
			if ($output_format eq 'nexus') {print FH $taxa->[$taxon_nr],get_spaces(($max_ID_length - length($taxa->[$taxon_nr])) + 2);}
			for (my $pos = 1; $pos < scalar(@$aln_out); ++$pos) {
				my $char = $aln_out->[$pos];
				print FH substr($aln_in->{$char}->{string},$taxon_nr,1);
			}
			print FH "\n";
		}
		if ($output_format eq 'nexus') {
			print FH ";\nend;\n\n";
			if ($partitioning) {
				my $parts;
				for (my $i=1; $i < scalar(@$aln_out); ++$i) {
					my $char = $aln_out->[$i];
					my $part = $aln_in->{$char}->{'partition'};
					my $output_pos = $aln_in->{$char}->{'output_pos'};
					push (@{$parts->{$part}},$output_pos);
				}
				foreach my $part (keys(%$parts)) {
					print FH "charset $part = @{$parts->{$part}};\n"
				}
			}
			if ($base_pairing) {
				my ($done,@pairs);
				for (my $i=1; $i < scalar(@$aln_out); ++$i) {
					my $char1 = $aln_out->[$i];
					if ($aln_in->{$char1}->{'linked'}) {
						my $char2 = $aln_in->{$char1}->{'linked'};
						my $pos1 = $aln_in->{$char1}->{'output_pos'};
						my $pos2 = $aln_in->{$char2}->{'output_pos'};
						unless ($done->{$char1}) {
							push @pairs,"$pos1\:$pos2";
						}
						$done->{$char1} = 1;
						$done->{$char2} = 1;
					}
				}
				print FH "pairs ",join(', ',@pairs),";\n";
			}
		}
		close FH;
	}

###################################################################################################################################
#   subroutines  ---  C A L C U L A T I O N S
###################################################################################################################################

	sub sort_list_ascending {
		my $list = shift;
		my $out;
		@$out = Sort_Table(
			cols      => '2',
			field     => '1',
			sorting   => 'ascending',
			structure => 'csv',
			separator => "\,",
			data      => $list,
		);
		return $out;
	}
	
	sub get_spaces {
		my $nr = shift;
		my $out; $out = "";
		for (my $i=0; $i<$nr; ++$i) {$out .= ' ';}
		return $out;
	}
	
	sub get_zeros {
		my $nr = shift;
		my $out; $out = "";
		for (my $i=0; $i<$nr; ++$i) {$out .= '0';}
		return $out;
	}

###################################################################################################################################
#   subroutines  ---  F I L E    P A R S E R S
###################################################################################################################################

	sub file_check {		# checks whether file exists and dies if this is not the case
		my $file = shift;
		(-e $file) || die "error -- file not found: $file\n";
	}
	
	sub parse_alignment {		# parses the alignment using Bio::SeqIO
		my $file = shift;
		file_check($file);
		my $seqio = Bio::SeqIO->new('-file' => $file, '-format' => 'fasta');
		my $data;
		while (my $seq = $seqio->next_seq) {
			push @$taxa,$seq->id;
			my $seqa;
			@$seqa = split('',$seq->seq);
			for (my $i=0; $i < scalar(@$seqa); ++$i) {
				my $char = 'i'.sprintf("%20d",($i+1));
				$data->{$char}->{'string'} .= $seqa->[$i];
			}
		}
		return $data;
	}
	
	sub parse_rates_file {		# parses the file with rates from HyPhy
		my $file = shift;
		file_check($file);
		open(FH,$file); my @array_in = <FH>; close FH;
		for (my $i=1; $i < scalar(@array_in); ++$i) {
			my $char = 'i'.sprintf("%20d",($i));
			if ($array_in[$i] =~ /^(.+?)\,/) {
				$aln_in->{$char}->{'rate'} = $1;
			} elsif ($array_in[$i] =~ /^(.+?)\t/) {
				$aln_in->{$char}->{'rate'} = $1;
			} else {
				die "FATAL ERROR -- could not derive rate from position $char from $file -- make sure file is the two-column output from HyPhy and that it uses either comma or tab character to separate fields\n"
			}
		}
	}
	
	sub parse_base_pairings_file {		# parses the file with base-pair information
		my $file = shift;
		file_check($file);
		open(FH,$file); my @a = <FH>; close FH;
		my $counter; $counter = 0;
		foreach my $line (@a) {
			if ($line =~ /(\d+):(\d+)/) {
				my ($d1,$d2) = ($1,$2);
				my $char1 = 'i'.sprintf("%20d",$d1);
				my $char2 = 'i'.sprintf("%20d",$d2);
				$aln_in->{$char1}->{'linked'} = $char2;
				$aln_in->{$char2}->{'linked'} = $char1;
				my $pair_rate = ( ( $aln_in->{$char1}->{'rate'} + $aln_in->{$char2}->{'rate'} ) / 2 );
				$aln_in->{$char1}->{'rate'} = $pair_rate;
				$aln_in->{$char2}->{'rate'} = $pair_rate;
			}
		}
		foreach my $char (keys(%$aln_in)) {
			unless ($aln_in->{$char}->{'linked'}) {
				$aln_in->{$char}->{'linked'} = 0;
			}
		}
	}
	
	sub parse_data_partitions_file {
		my $file = shift;
		file_check($file);
		open(FH,$file); my @a = <FH>; close FH;
		foreach my $line (@a) {
			$line =~ /\s*(.+?)\s*=\s*(.+)$/;
			my $name = $1;
			my $def = $2;
			my @parts = split (' ',$def);
			foreach my $part (@parts) {
				if ($part =~ /\\/) {
					$part =~ /(\d+)-(.+?)\\(\d+)/;
					my ($start,$stop,$interval) = ($1,$2,$3);
					if ($stop eq '.') {$stop = $input_seq_length;}
					for (my $i = $start; $i <= $stop; $i += $interval) {
						my $char = 'i'.sprintf("%20d",$i);
						$aln_in->{$char}->{'partition'} = $name;
					}
				} elsif ($part =~ /-/) {
					$part =~ /(\d+)\-(\d+)/;
					my ($start,$stop) = ($1,$2);
					if ($stop eq '.') {$stop = $input_seq_length;}
					for (my $i = $start; $i <= $stop; ++$i) {
						my $char = 'i'.sprintf("%20d",$i);
						$aln_in->{$char}->{'partition'} = $name;
					}
				} else {
					$part =~ /(\d+)/;
					my $char = 'i'.sprintf("%20d",$1);
					$aln_in->{$char}->{'partition'} = $name;
				}
			}
		}
		foreach my $char (keys(%$aln_in)) {
			unless ($aln_in->{$char}->{'partition'}) {
				$aln_in->{$char}->{'partition'} = 0;
			}
		}
	}