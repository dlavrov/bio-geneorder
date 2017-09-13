#
# BioPerl module for Bio::GeneOrder::gene
#
#	Based on code written by Dennis Lavrov
#	Adapted for Bioperl by Walker Pett
#
# Copyright Dennis Lavrov, Walker Pett
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::GeneOrder::permutation - Basic permutation structure

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::GeneOrder> 
class.  

=head1 DESCRIPTION

The object is a simple representation of the order of genes on a chromosome

=head1 FEEDBACK

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

# Let the code begin...

package Bio::GeneOrder::permutation;
use strict;

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::GeneOrder::permutation( -source	=>  $source,
                                                      -circular	=>  $i,
                                                      -pi       =>  $arrayRef);
 Function: Returns a new Bio::GeneOrder::permutation object 
 Returns : A Bio::GeneOrder::permutation object
 Args    : -source            => an accession number or other reference	
           -circular          => a value of 1 or 0 [default 0]
           -pi                => a reference to an array of signed integers

=cut

sub new {
	my ($caller, @args) = @_;
	my $self = $caller->SUPER::new(@args);
	bless $self, $caller;

	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys

	$caller->throw("source argument provided, but with an undefined value") 
		if( exists($param{'-source'}) && !defined($param{'-source'}) );
	$caller->throw("circular argument provided, but with an undefined value") 
		if( exists($param{'-circular'}) && !defined($param{'-circular'}) );
	$caller->throw("pi argument is required") 
		if( !defined($param{'-pi'}) );
	
	$self->source($param{'-source'}) if defined $param{'-source'};
	$self->is_circular(defined $param{'-circular'} ? $param{'-circular'} : 0);
	$self->pi(@{$param{'-pi'}});
	
	return $self;
}

=head2 reorder

 Title   : reorder
 Usage   : $pi->reorder($key);
 Function: Reorders genes so that the specified gene is first and in the same orientation.
 Returns : 1 if reordered, 0 if not

=cut

sub reorder {
	my ($self,$key) = @_;
	
	$self->throw("reorder requires an argument") unless defined $key;
	
	my $r = 0;
	
	my @pi = $self->pi;
	
	if( grep(abs($_) == $key, @pi) ){
		
		while(1){
			my $gene = shift @pi;
			push @pi, $gene;
			last if abs($pi[0]) == $key;
		}

		if($pi[0] < 0){
			@pi = reverse(@pi);
			map($_ *= -1, @pi);
			unshift @pi, pop @pi;
		}
		
		$self->pi(@pi);
		$r = 1;
	}
	
	return $r;
}

=head2 pi

 Title   : pi
 Usage   : $obj->pi( 'all' )
 Function: Get/set the permutation structure
 Returns : Integer


=cut

sub pi {
	my ($self,@pi) = @_;

	if( @pi){
		if($pi[0] eq 'all'){
			@pi = unpack("s*", $self->{'pi'});
			shift @pi;
			
			return @pi;
		}
		
		my @return = @pi;
		unshift @pi, $self->is_circular;
		$self->{'pi'} = pack("s*", @pi);
		
		return @return;
	}else{
		
		@pi = unpack("s*", $self->{'pi'});
		
		shift @pi;
		@pi = grep($self->{'key'}->{'filt'}->{ $self->{'key'}->{'name'}->{abs($_)} } == 0, @pi);
		return @pi;
	}
}

=head2 is_circular

 Title   : is_circular
 Usage   : my $circle = $pi->is_circular();
 Function: Returns true if the source sequence is circular
 Returns : A boolean value

=cut

sub is_circular {
	my ($self, $value) = @_;
	
	my $circular = 0;
	my @pi = ();
	
	if(defined $self->{'pi'}){
		@pi = unpack("s*", $self->{'pi'});
		$circular = shift @pi;
	}

	if( defined $value){
		$circular = $value;
		unshift @pi, $circular;
		$self->{'pi'} = pack("s*", @pi);
	}
	
	return $circular;
}

=head2 source

 Title   : source
 Usage   : $obj->source()
 Function: Get/set the source
 Returns : String


=cut

sub source {
	my ($self, $value) = @_;

	if( defined $value){
		$self->{'source'} = $value;
	}
	
	return $self->{'source'};
}

=head2 _key

 Title   : _key
 Usage   : $permutation->_key();
 Function: Get/set the gene name/number key.

=cut

sub _key {
	my ($pi,$key,$map) = @_;

	
	if(defined $key){
		if(defined $map){
			my @pi = $pi->pi('all');
			map($_ = (abs($_)/$_)*$key->{'index'}->{ $map->{ $pi->{'key'}->{'name'}->{abs($_)} } }, @pi);
			$pi->pi(@pi);
		}else{
			my @pi = $pi->pi('all');
			map($_ = (abs($_)/$_)*$key->{'index'}->{ $pi->{'key'}->{'name'}->{abs($_)} }, @pi);
			$pi->pi(@pi);
		}
		
		$pi->{'key'} = $key;
	}
	
	return $pi->{'key'};
}

1;

