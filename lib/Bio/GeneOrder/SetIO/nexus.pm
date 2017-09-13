#
# BioPerl module for Bio::GeneOrder::SetIO::nexus
#
#	Based on code written by Dennis Lavrov
#	Adapted for Bioperl by Walker Pett
#
# Copyright Dennis Lavrov, Walker Pett
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::GeneOrder::SetIO::nexus - nexus GeneOrder input/output stream

=head1 SYNOPSIS

    Do not use this module directly.  Use it via the L<Bio::GeneOrder::SetIO> class, as in:

    use Bio::GeneOrder::SetIO;

    $out  = Bio::GeneOrder::SetIO->new( -file     => "outputfilename" ,
                                        -format   => "nexus",
                                        -encoding => "MPME");
    $out->write_set($set);


=head1 DESCRIPTION

This object can write L<Bio::GeneOrder::Set> objects to a NEXUS matrix file 
formatted with a specified encoding.

A nice feature of this module is that it can be used to write distance matrices 
in a variety of encodings relevant to the phylogenetic analysis of genome arrangements.

Supported encodings:

 'MPBE'              For Maximum Parsimony on Binary Encodings
 'MPME'              For Maximum Parsimony on Multi-State Encodings
 'copies'            Encodes a character state for each gene as the number of copies of that gene
 'breakpoint'        Encodes the breakpoint distance matrix.
 'inversion'         Encodes the inversion distance matrix.
 'common_intervals'  Encodes the common interval distance matrix.
 'DCJ'               Encodes the Double-Cut and Join distance matrix.

For a detailed description of each of these encodings, see:

-MPBE

	Cosner, M.E., Jansen, R.K., Moret, B.M.E., Raubeson, L.A., Wang, L.-S., 
	Warnow, T., and Wyman, S.K. (2000). An empirical comparison of phylogenetic 
	methods on chloroplast gene order data in Campanulaceae. In Comparative 
	Genomics: Empirical and Analytical Approaches to Gene Order Dynamics, Map 
	Alignment, and the Evolution of Gene Families (ed. D. Sankoff and J. Nadeau), 
	pp. 99-121. Kluwer Academic Pub., Dordrecht, Netherlands.

	M.E. Cosner, R.K. Jansen, B.M.E. Moret, L.A. Raubeson, L. Wang, T. Warnow, and S.K. 
	Wyman. A new fast heuristic for computing the breakpoint phylogeny and experimental 
	phylogenetic analyses of real and synthetic data. In Proc. 8th Int'l Conf. on Intelligent 
	Systems for Mol. Biol. ISMB-2000, pages 104-115, 2000.

-MPME

	Wang, L.-S., Jansen, R.K., Moret, B.M.E., Raubeson, L.A., and Warnow, T. 
	(2002). Fast phylogenetic methods for genome rearrangement evolution: An em- 
	pirical study. In Proc. 7th Pacific Symp. on Biocomputing (PSB'02), pp. 524-535. 
	World Scientific Pub.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHORS

Walker Pett, Dennis Lavrov

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# 'Let the code begin...

package Bio::GeneOrder::SetIO::nexus;

use strict;
use Thread;
use Bio::GeneOrder::Distance;

use base qw(Bio::GeneOrder::SetIO);

#The set of 77 allowed character state symbols in PAUP except '?' for missing characters
my @ALL_SYMBOLS = qw(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z 
				 a b c d e f g h i j k l m n o p q r s t u v w x y z
				 0 1 2 3 4 5 6 7 8 9
				 ! @ # $ % ^ & * - + _ = < > . ' : | / \ `
				 );
				 
my ($NTAX,$NCHAR,$TAXLABELS,%STATELABELS,$STATELABELS,%SYMBOL_TABLE,%MATRIX,$CHARWIDTH,$MAX,
	$NOCHARSTATELABELS,$EXTENDED_SYMBOLS,$MARGIN,$SPACE,$WARNED,@SYMBOLS,$DISTANCE,$MISSING);

#We declare this encodings array for other programs that determine
#whether this module supports encodings
our @ENCODINGS = qw(MPBE MPME copies);
push @ENCODINGS, Bio::GeneOrder::Distance->supported_distances;

=head2 new

 Title   : new
 Usage   : $setio = new Bio::GeneOrder::SetIO( -format   => 'nexus',
                                               -file     => 'filename',
                                               -encoding => 'encoding');
 Function: returns a new Bio::GeneOrder::SetIO object to handle nexus files
 Returns : Bio::GeneOrder::SetIO::nexus object
 Args    : -file     => name of file to read in or to write, with ">"
           -fh       => alternative to -file param - provide a filehandle
                        to read from or write to
           -format   => gene order file format to process or produce
           -encoding => the character state encoding to use when writing
                        or interpreting gene order sets.
                        See DESCRIPTION for more info.

=cut

sub _initialize {
  my($self,@args) = @_;

  my ( $encoding) = $self->_rearrange( [qw(ENCODING)], @args );
  
  $self->{'encoding'} = defined $encoding ? $encoding : 'MPBE';

  $self->_initialize_io(@args);
  1;
}

=head2 next_set

 Title   : next_set
 Note    : reading from NEXUS files is not supported

=cut

sub next_set {
    my ($self) = @_;

	return $self->throw("Reading of NEXUS files is not supported.");
}

=head2 write_set

 Title   : write_set
 Usage   : $stream->write_set(@set)
 Function: writes the nexus-format object (.nex) into the stream
 Returns : 1 for success and 0 for failure
 Args    : Bio::GeneOrder::Set object

=cut

sub write_set {
    my ( $self, $set ) = @_;

	my @results = $self->get_matrix($set,$self->{'encoding'});
	my $matrix = ${ $results[0] }[0];

	#actually print everything to nexus file
	$self->_print("\#NEXUS\n");
	if($DISTANCE){
		$self->_print("begin distances;\ndimensions ntax = $NTAX;\n".
					  "format triangle=both labels;\n");
	}else{
		$self->_print("begin data;\ndimensions ntax = $NTAX  nchar = $NCHAR;\n".
					  "format missing=$MISSING transpose  respectcase  symbols= \"@SYMBOLS[0..$MAX]\";\n");
	}
	$self->_print("TAXLABELS\n$TAXLABELS;\n") or return ;
	$self->_print("charSTATELABELS\n$STATELABELS;\n") unless( $NOCHARSTATELABELS);
	$self->_print("matrix\n$matrix;\nend;") or return ;
	
	$self->flush if $self->_flush_on_write && defined $self->_fh;
    return 1;

}

=head2 get_matrix

 Title   : get_matrix
 Note    : get the matrix and TAXLABELS etc. programatically

=cut

sub get_matrix {
	 my ( $self, $set, $encoding, $width,$a,$b) = @_;
	 
	 ($NTAX,$NCHAR,$TAXLABELS,%STATELABELS,%SYMBOL_TABLE,%MATRIX,$CHARWIDTH,$MAX,$STATELABELS,
	$NOCHARSTATELABELS,$EXTENDED_SYMBOLS,$MARGIN,$SPACE,$WARNED,@SYMBOLS,$DISTANCE,$MISSING) = ();
	
	$SPACE = 3;
	my @sorted = $set->orders;
	$NTAX = scalar @sorted;
	
	#We need space to print the numbers at the top of the matrix
	$CHARWIDTH = length(sprintf("%d",scalar @sorted));

	#multi-state encoding
	if( $encoding eq 'MPME'){
		$EXTENDED_SYMBOLS=1;
	 	
		$MISSING = '?';

		foreach my $order ( @sorted){
			$TAXLABELS .= "\t'".$order->name."'\n";
			my $str = $order->order;
			while( $str =~ /(?=(^|\s|\-)([^\-\s]*)\s+(\-)?(\S*))/g){
				my ($strandA,$geneA,$strandB,$geneB) = ($1,$2,$3,$4);
				#print $order->name."\t$strandA'$geneA'\t$strandB'$geneB'\n";
				next if $geneA eq '' || $geneB eq '';
				#print "$strandA$geneA $strandB$geneB\n";
				
				my $character = $strandA eq '-' ? 
					"5_".$geneA : "3_".$geneA;
					
				my $state = $strandB eq '-' ? 
					"3_".$geneB : "5_".$geneB;
				
				$MARGIN = length($state) +$SPACE > $MARGIN ? length($state) +$SPACE : $MARGIN;

				push(@{ $STATELABELS{ $character}}, $state) unless grep($_ eq $state, @{ $STATELABELS{ $character}});
				$MATRIX{ $order->name}{ $character} = $state;
				
				$character = $strandB eq '-' ? 
					"3_".$geneB : "5_".$geneB;

				$state = $strandA eq '-' ? 
					"5_".$geneA : "3_".$geneA;

				$MARGIN = length($state) +$SPACE > $MARGIN ? length($state) +$SPACE : $MARGIN;

				push(@{ $STATELABELS{ $character}}, $state) unless grep($_ eq $state, @{ $STATELABELS{ $character}});
				$MATRIX{ $order->name}{ $character} = $state;
			}
			
			if($order->is_circular){
				$order->order =~ /^(\-)?(\S+).+\s(\-)?(\S+)$/;
				my ($strandB,$geneB,$strandA,$geneA) = ($1,$2,$3,$4);
				next if $geneA eq '' || $geneB eq '';
				my $character = $strandA eq '-' ? 
					"5_".$geneA : "3_".$geneA;
					
				my $state = $strandB eq '-' ? 
					"3_".$geneB : "5_".$geneB;
				
				$MARGIN = length($state) +$SPACE > $MARGIN ? length($state) +$SPACE : $MARGIN;

				push(@{ $STATELABELS{ $character}}, $state) unless grep($_ eq $state, @{ $STATELABELS{ $character}});
				$MATRIX{ $order->name}{ $character} = $state;
				
				$character = $strandB eq '-' ? 
					"3_".$geneB : "5_".$geneB;

				$state = $strandA eq '-' ? 
					"5	_".$geneA : "3_".$geneA;

				$MARGIN = length($state) +$SPACE > $MARGIN ? length($state) +$SPACE : $MARGIN;

				push(@{ $STATELABELS{ $character}}, $state) unless grep($_ eq $state, @{ $STATELABELS{ $character}});
				$MATRIX{ $order->name}{ $character} = $state;
			}
		}

	#binary encoding
	}elsif( $encoding eq 'MPBE'){

		$self->throw("MPBE encoding is not currently implemented");

	#copies encoding
	}elsif( $encoding eq 'copies'){

		$MISSING = '0';
		
		my $key = $set->_key();
		
		my @genes = keys %{ $key->{'index'} };
		
		#print join "\n", sort @genes , "\n";
		#return;
	
		foreach my $order ( @sorted){
			my $name = $order->name;
			$TAXLABELS .= "\t'".$name."'\n";

			foreach my $gene ( @genes ){
				if( my @count = grep($gene eq $_, $order->genes)){
					$MATRIX{ $name}{ $gene} = @count;

					for(my $i=0;$i<@count+1;$i++){
						push(@{ $STATELABELS{ $gene}},$i) unless grep($_ == $i, @{ $STATELABELS{ $gene}});
					}
					$MARGIN = length($gene) +$SPACE > $MARGIN ? length($gene) +$SPACE : $MARGIN;
				}else{
					$MATRIX{ $name}{ $gene} = 0;
				}
			}
		}
		
		@SYMBOLS = qw(0 1 2 3 4 5 6 7 8 9 a b c d e f g h i j k l m n o p q r s t u v w x y z);
		$NOCHARSTATELABELS = 1;

	#breakpoint encoding
	}elsif( grep($_ eq $encoding, Bio::GeneOrder::Distance->supported_distances) ){
		
		$MISSING = '?';
		my @threads;
		
		foreach my $order ( @sorted){
			$TAXLABELS .= "\t'".$order->name."'\n";
			$MARGIN = length($order->name) +$SPACE if length($order->name) +$SPACE > $MARGIN;

			foreach my $order2 ( @sorted){
				my $count = $order->distance->$encoding($order,$order2);

				$MAX = $count if $count > $MAX;
				$CHARWIDTH = (length($count)+1) if (length($count)+1) > $CHARWIDTH;
				$MATRIX{ $order->name}{ $order2->name} = $count;
				$MATRIX{ $order2->name}{ $order->name} = $count;
			}
		}
		
		map( $_->join, @threads);
		
		@SYMBOLS = (0..$MAX);
		map( @{ $STATELABELS{ $_->name}} = @SYMBOLS, @sorted);
		
		$DISTANCE = 1;
		$NOCHARSTATELABELS = 1;
		
	#shared encoding
	}else{

		$self->throw("encoding: ".$encoding." not recognized");
	}
	
	if($EXTENDED_SYMBOLS){ @SYMBOLS = @ALL_SYMBOLS;}
	
	foreach my $character (sort keys %STATELABELS){
		$NCHAR++;
		$STATELABELS .= sprintf("%-3s'%s'/",$NCHAR,$character);

		for(my $i=0;$i<@{ $STATELABELS{$character}};$i++){
			$STATELABELS{ $character}->[$i] =~ s/\s+//g;
			$STATELABELS .= "'".$STATELABELS{ $character}->[$i]."' ";
			#add a label for this particular character state
			#if we have more states than symbols, throw an exception
			$self->throw("Number of character states for character '$character' exceeds maximum number allowed: $i") 
				if($i == scalar @SYMBOLS);
			#PAUP only allows 32 character states on 32-bit machines
			if ($i == 32 && !$WARNED && !$DISTANCE){
				$self->warn("Number of character states exceeds maximum number allowed in 32-bit PAUP (32)");
				$WARNED++;
			}
			#PAUP only allows 64 character states on 64-bit machines
			if ($i == 64 && !$WARNED && !$DISTANCE){
				$self->warn("Number of character states exceeds maximum number allowed in 64-bit PAUP (32)");
				$WARNED++;
			}

			$SYMBOL_TABLE{ $character}{ $STATELABELS{ $character}->[$i]} = $SYMBOLS[$i];
			$MAX = $i if $i > $MAX;
		}
		$STATELABELS .= ",\n";
	}
	
	my @matrix;
	
	if($width){
		my $tax_width = int(($width - $MARGIN)/($CHARWIDTH+1));

		my $i=0;
		while($i < $NTAX){
			my $chunk = $self->_get_chunk($i,$tax_width)."\n";
			push @matrix, $chunk;
			$i += $tax_width;
		}
	}else{
		my $chunk = $self->_get_chunk(0,$NTAX);
		push @matrix, $chunk;
	}
	
	return (\@matrix,$CHARWIDTH,$MARGIN,$NTAX);
}

sub _get_chunk {
	my ($self,$start,$end) = @_;
	
	#build strings from the STATELABELS and matrix
	my $chunk;
	
	foreach my $character (sort keys %STATELABELS){
		$chunk .= sprintf("%-".$MARGIN."s","'".$character."'");

		my @orders = sort keys %MATRIX;
		for(my $u=$start;$u<$start+$end;$u++){
			last if $u >= @orders;
			
			$chunk .= defined $MATRIX{ $orders[$u]}{ $character} ? 
						sprintf("%-".($CHARWIDTH+1)."s",$SYMBOL_TABLE{ $character}{ $MATRIX{ $orders[$u]}{ $character}}) : 
						sprintf("%-".($CHARWIDTH+1)."s",$MISSING);
		}
		$chunk .= "\n";
	}
	
	return $chunk;
}

=head2 supported_encodings

 Title   : supported_encodings
 Note    : Get a list of encodings supported by this module

=cut

sub supported_encodings {
	my $self = shift;

	return @ENCODINGS;
}

1;
