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

Bio::GeneOrder::SetIO::fasta - geneorder raw fasta file input/output stream

=head1 SYNOPSIS

    Do not use this module directly.  Use it via the L<Bio::GeneOrder::SetIO> class, as in:

    use Bio::GeneOrder::SetIO;

    $out  = Bio::GeneOrder::SetIO->new( -file     => ">outputfilename" ,
                                        -format   => "go" );
    out->write_Set($set);

	$in  = Bio::GeneOrder::SetIO->new( -file     => "inputfilename" ,
                                        -format   => "go" );
    $set = $in->next_set();


=head1 DESCRIPTION

This object can write L<Bio::GeneOrder::Set> objects to a raw gene order file.
It can also read L<Bio::GeneOrder::Set> objects from GeneOrder formatted files.
Each gene order is contained in a GeneOrder formatted file with the name of the
gene order on a line starting with '>', followed by a line containing a space-delimted 
sequence of gene names prepended with transcriptional polarity symbols '+' or '-'.

For example:

>Lampsilis ornata|NC_005335
-nad4 +nad4L -atp8 +trnD -atp6 -cox3 -cox1 -cox2 ... etc.

;Lines can be commented with a semi-colon.  Blank lines are ignored.

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

package Bio::GeneOrder::SetIO::fasta;

use strict;
use Bio::GeneOrder;

use base qw(Bio::GeneOrder::SetIO);
use vars qw(%REV %FILTER %SWITCH $LINEAR);

BEGIN {
	$LINEAR = '~';
	%REV = ( '+' => '-',
			 '-' => '',
			 '' => '-');
	%SWITCH = ( 1	=> '', -1	=> '-', 0	  => undef,
			  ''	=> 1,   '-'	=> -1 , undef => 0,'+' => 1);
	%FILTER = ();
}


=head2 new

 Title   : new
 Usage   : $setio = new Bio::GeneOrder::SetIO( -format   => 'nexus',
                                               -file     => 'filename');
 Function: returns a new Bio::GeneOrder::SetIO object to handle nexus files
 Returns : Bio::GeneOrder::SetIO::nexus object
 Args    : -file     => name of file to read in or to write, with ">"
           -fh       => alternative to -file param - provide a filehandle
                        to read from or write to
           -format   => gene order file format to process or produce
=cut

sub _initialize {
  my($self,@args) = @_;

  $self->_initialize_io(@args);
  1;
}

=head2 next_set

 Title   : next_set
 Usage   : $set = $stream->next_set()
 Function: retrieves a GeneOrder::Set object from the stream
 Returns : Bio::GeneOrder::Set object

=cut

sub next_set {
    my ($self,$test,$verbose) = @_;

	return 1 if($test);

	my (@orders,@pi,$entry,$name,$order);
	
	while (defined ($entry = $self->_readline) ) {
		chomp $entry;
		if ( $entry =~ /^>(.*)/ ) {
			if(defined $name){
				#print "$name\n";
				push @orders, Bio::GeneOrder->new(@pi,-name => $name);
				@pi = ();
			}
			$name = $1;
		
		}elsif( $entry !~ /\S/ || $entry =~ /^;/){
		}elsif(defined $name){
			push @pi, $entry;
		}else{
			$self->throw("gene order file '". $self->file. "' not formatted correctly");
		}
	}
	
	if(defined $name){
		push @orders, Bio::GeneOrder->new(@pi,-name => $name);
		@pi = ();
	}
	
	return @orders;
}

=head2 write_set

 Title   : write_set
 Usage   : $stream->write_set(@set)
 Function: writes the geneorder-format object (.go) into the stream
 Returns : 1 for success and 0 for failure
 Args    : Bio::GeneOrder::Set object

=cut

sub write_set {
    my ( $self, $set ) = @_;
	
	foreach my $order ($set->orders){
		$self->_print(">".$order->name."\n");
		my @pi = $order->pi;
		
		foreach my $pi (@pi){
			my $string = join ' ', map $SWITCH{abs($_)/$_}.$order->{'key'}->{'name'}->{abs($_)}, grep( ! defined $order->{'key'}->{'filt'}->{$order->{'key'}->{'name'}->{abs($_)}} || $order->{'key'}->{'filt'}->{$order->{'key'}->{'name'}->{abs($_)}} == 0, $pi->pi );
			$string = "$LINEAR ".$string unless($pi->is_circular);
			$self->_print("$string\n");
		}
		
		$self->_print("\n");
	}		

	$self->flush if $self->_flush_on_write && defined $self->_fh;
    return 1;

}

1;
