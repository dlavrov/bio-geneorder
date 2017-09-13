#
# BioPerl module for Bio::GeneOrder::SetIO
#
#	Based on code written by Dennis Lavrov
#	Adapted for Bioperl by Walker Pett
#
# Copyright Dennis Lavrov, Walker Pett
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::GeneOrder::SetIO - Handler for GeneOrder::SetIO Formats

=head1 SYNOPSIS

    use Bio::GeneOrder::SetIO;

    $in  = Bio::GeneOrder::SetIO->new( -file => "inputfilename" ,
                                       -format => 'grappa');
    $out = Bio::GeneOrder::SetIO->new( -file => ">outputfilename" ,
                                       -format => 'nexus');

    while ( my $set = $in->next_set() ) {
	    $out->write_set($set);
    }

  # Now, to actually get at the set object and the order objects it
  # contains, use the Bio::GeneOrder::Set and Bio::GeneOrder methods

    use Bio::GeneOrder::SetIO;

    $in  = Bio::GeneOrder::SetIO->new( -file => "inputfilename" ,
                                       -format => 'grappa');

    while ( my $set = $in->next_set() ) {
       foreach my $order ( $set->orders){
           print $order->name,"\n";
       }
    }


  # The SetIO system does have a filehandle binding like SeqIO

    use Bio::GeneOrder::SetIO;

    $in  = Bio::GeneOrder::SetIO->newFh( -file => "inputfilename" ,
                                         -format => 'grappa');
    $out = Bio::GeneOrder::SetIO->newFh( -format   => 'nexus',
                                         -encoding => 'binary');

    # World's shortest Grappa<->NexusMPBE format converter:
    print $out $_ while <$in>;


=head1 DESCRIPTION

Bio::GeneOrder::SetIO is a handler module for the formats in the SetIO set
(eg, Bio::GeneOrder::SetIO::grappa).  It is based on the L<Bio::SeqIO> module 
and its methods can be used in almost exactly the same way. The GeneOrder::SetIO
system provides support for reading and writing Bio::GeneOrder::Set objects to 
file types that can be used with software like GRAPPA, PAUP, etc for performing 
analyses on gene order data.  It is important to note that SetIO will not attempt 
to guess file formats as SeqIO does. The '-format' argument is required for methods 
that use it.

See L<Bio::SeqIO> for a more in-depth explanation of BioPerl's I/O subsystems.

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

package Bio::GeneOrder::SetIO;

use strict;

use Symbol();

use base qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : $stream = Bio::GeneOrder::SetIO->new( -file => $filename,
                                                 -format => 'Format')
 Function: Returns a new gene order set stream
 Returns : A Bio::GeneOrder::SetIO stream initialised with the appropriate format
 Args    : Named parameters:
             -file     => $filename
             -fh       => filehandle to attach to
             -format   => format

=cut

sub new {
	my ($caller,@args) = @_;
	my $class = ref($caller) || $caller;

	# or do we want to call SUPER on an object if $caller is an
	# object?
	if( $class =~ /Bio::GeneOrder::SetIO(\S+)/ ) {
		my ($self) = $class->SUPER::new(@args);
		$self->_initialize(@args);
		return $self;
	} else {

		my %param = @args;
		@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	
		if (!defined($param{'-file'}) && !defined($param{'-fh'})) {
		  $class->throw("file argument provided, but with an undefined value") if exists($param{'-file'});
		  $class->throw("fh argument provided, but with an undefined value") if (exists($param{'-fh'}));
		}
	
		if( !defined $param{'-format'}){
		  $class->throw("format argument provided, but with an undefined value") if exists($param{'-format'});
		  $class->throw("format argument is required") if !exists($param{'-format'});
		}
	
		my $format = $param{'-format'} ||
			$class->_guess_format( $param{-file} || $ARGV[0] );
	
		$format = "\L$format";	# normalize capitalization to lower case
	
		return unless( $class->_load_format_module($format) );
		return "Bio::GeneOrder::SetIO::$format"->new(@args);
    }
}

=head2 newFh

 Title   : newFh
 Usage   : $fh = Bio::GeneOrder::SetIO->newFh(-file=>$filename,-format=>'Format')
 Function: does a new() followed by an fh()
 Example : $fh = Bio::GeneOrder::SetIO->newFh(-file=>$filename,-format=>'Format')
           $sequence = <$fh>;   # read a sequence object
           print $fh $sequence; # write a sequence object
 Returns : filehandle tied to the Bio::GeneOrder::SetIO::Fh class
 Args    :

=cut

sub newFh {
  my $class = shift;
  return unless my $self = $class->new(@_);
  return $self->fh;
}

=head2 fh

 Title   : fh
 Usage   : $obj->fh
 Function:
 Example : $fh = $obj->fh;      # make a tied filehandle
           $sequence = <$fh>;   # read a sequence object
           print $fh $sequence; # write a sequence object
 Returns : filehandle tied to Bio::GeneOrder::SetIO class
 Args    : none

=cut


sub fh {
  my $self = shift;
  my $class = ref($self) || $self;
  my $s = Symbol::gensym;
  tie $$s,$class,$self;
  return $s;
}

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  $self->_initialize_io(@args);
  1;
}

=head2 next_set

 Title   : next_set
 Usage   : $set = stream->next_set
 Function: Reads the next set object from the stream and returns it.

           Certain driver modules may encounter entries in the stream
           that are either misformatted or that use syntax not yet
           understood by the driver. If such an incident is
           recoverable, the driver will issue a warning. In the case of
           a non-recoverable situation an exception will be thrown.  Do
           not assume that you can resume parsing the same stream
           after catching the exception. Note that you can always turn
           recoverable errors into exceptions by calling
           $stream->verbose(2).

 Returns : a Bio::GeneOrder::Set sequence object
 Args    : none

See L<Bio::Root::RootI>, L<Bio::GeneOrder::Set>

=cut

sub next_set {
   my ($self, $seq) = @_;
   $self->throw("Sorry, you cannot read from a generic Bio::GeneOrder::SetIO object.");
}

=head2 write_set

 Title   : write_set
 Usage   : $stream->write_set($set)
 Function: writes the $set object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::GeneOrder::Set object

=cut

sub write_set {
    my ($self, $seq) = @_;
    $self->throw("Sorry, you cannot write to a generic Bio::GeneOrder::SetIO object.");
}

=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL SetIO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_format_module {
	my ($self, $format) = @_;
	my $module = "Bio::GeneOrder::SetIO::".$format;
	my $ok;

	eval {
		$ok = $self->_load_module($module);
	};
	if ( $@ ) {
		print STDERR <<END;
$self::$format cannot be found
Exception $@
For more information about the GeneOrder::SetIO system please see the 
GeneOrder::SetIO docs. This includes ways of checking for formats at 
compile time, not run time
END
		;
	}
	return $ok;
}

=head2 _guess_format

 Title   : _guess_format
 Usage   : $obj->_guess_format($filename)
 Function: guess format based on file suffix
 Example :
 Returns : guessed format of filename (lower case)
 Args    :

=cut

sub _guess_format {
   my $class = shift;
   return unless $_ = shift;
   return 'nexus'      if /\.(nex|nxs|nexus)$/i;
   return 'grappa'     if /\.so$/i;
}

sub DESTROY {
	my $self = shift;
	$self->close();
}

sub TIEHANDLE {
	my ($class,$val) = @_;
	return bless {'setio' => $val}, $class;
}

sub READLINE {
	my $self = shift;
	return $self->{'setio'}->next_set() unless wantarray;
	my (@list, $obj);
	push @list, $obj while $obj = $self->{'setio'}->next_set();
	return @list;
}

sub PRINT {
	my $self = shift;
	$self->{'setio'}->write_set(@_);
}

1;





