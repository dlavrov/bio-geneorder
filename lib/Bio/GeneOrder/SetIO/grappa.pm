#
# BioPerl module for Bio::GeneOrder::SetIO::grappa
#
#	Based on code written by Dennis Lavrov
#	Adapted for Bioperl by Walker Pett
#
# Copyright Dennis Lavrov, Walker Pett
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::GeneOrder::SetIO::grappa - grappa GeneOrder input/output stream

=head1 SYNOPSIS

    Do not use this module directly.  Use it via the Bio::GeneOrder::SetIO class, as in:

    use Bio::GeneOrder::SetIO;

    $out  = Bio::GeneOrder::SetIO->new( -file     => "outputfilename" ,
                                        -format   => "grappa");
    $out->write_set($set);


=head1 DESCRIPTION

This object can write GeneOrder set objects to a GRAPPA formatted .so
file for use with GRAPPA of Moret et al.

For information about GRAPPA and its associated file format, see:

     http://www.cs.unm.edu/~moret/GRAPPA/

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

package Bio::GeneOrder::SetIO::grappa;

use strict;

use base qw(Bio::GeneOrder::SetIO);

sub _initialize {
  my($self,@args) = @_;

  $self->_initialize_io(@args);
  1;
}


=head2 next_set

 Title   : next_set
 Usage   : not implemented

=cut

sub next_set {
    my ($self) = @_;

	return $self->throw("Reading of GRAPPA files is not supported.");
}

=head2 write_set

 Title   : write_set
 Usage   : $stream->write_set(@set)
 Function: writes the grappa-format object (.so) into the stream
 Returns : 1 for success and 0 for failure
 Args    : Bio::GeneOrder::Set object

=cut

sub write_set {
    my ( $self, $set ) = @_;

	my (%numbers,$dupe);
	
	#$self->throw("Analysis of gene orders with unequal gene content is not supported in GRAPPA 2.0") unless $set->flush;
	
	my $i=1;
	foreach my $order ($set->orders){
		my %used;

		$self->_print(">".$order->name."\n");
		my @pi = $order->pi;
		next, $self->warn("Multichromosomal genomes not supported in GRAPPA 2.0\nskipping ",$order->name) if scalar(@pi) > 1;
		
		my @order;
		foreach my $pi (@pi){
			#If we've already seen a gene of the same name in the same order
			#we have duplicate genes and GRAPPA will complain.
			my @ppi = $pi->pi;
			foreach my $name (@ppi){
				$dupe = 1 if(defined $used{ abs($name)});
		
				#If we haven't seen this gene name for this order, count up
				if(!defined $numbers{ abs($name)}){
					$numbers{ abs($name)} = $i++;
				}
				push @order, $name/abs($name)*$numbers{ abs($name)};
	
				#We have used this gene name
				$used{ abs($name)}++;
			}
		}
		$self->_print("@order\n");
	}
	
	#$self->warn("Analysis of gene orders with unequal gene content is not supported in GRAPPA") 
	#	if($dupe);


	$self->flush if $self->_flush_on_write && defined $self->_fh;
    return 1;	
}

1;