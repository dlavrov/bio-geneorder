#
# BioPerl module for Bio::GeneOrder::Distance
#
# Copyright Dennis Lavrov, Walker Pett
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

our $VERSION = '0.1';

=head1 NAME

Bio::GeneOrder::Distance - An object for obtaining distance measurements between gene orders

=head1 SYNOPSIS

    DO THIS
    
=head1 DESCRIPTION

Bio::GeneOrder::Distance is an object that is used to calculate
distance measures between gene orders using their underlying
permutation structure. It uses an XS library of C routines
for computing distances based on:

breakpoints
adjacencies
inversions
translocations
common intervals
DCJ
TDRL

In addition, relevant correction estimators and median solvers are available
where appropriate.  Distances can be calculated under insertions, deletions
duplications and/or multichromosoal genomes where applicable/available.

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

package Bio::GeneOrder::Distance;

use warnings;
use strict;
use Bio::GeneOrder;

use base qw(Bio::Root::Root);
use vars qw(%REV %SWITCH);

our @DISTANCES = qw(adjacencies breakpoints inversions translocations DCJ common_intervals);

BEGIN {
	%REV = ( '+' => '-',
			 '' => '-',
			 '-' => '' );
	%SWITCH = ( 1	=> '', -1	=> '-', 0	  => undef,
			  ''	=> 1,   '-'	=> -1 , undef => 0,'+' => 1);
}

require XSLoader;
XSLoader::load('Bio::GeneOrder::Distance', $VERSION);

#Create a variable to hole our singleton instance
our $INSTANCE;

=head2 new

 Title   : new
 Usage   : $obj = new Bio::GeneOrder::Distance();
 Function: Returns a singleton Bio::GeneOrder::Distance object 
 Returns : A Bio::GeneOrder::Distance object

=cut

sub new {
	my ($caller,@args) = @_;
	
	unless(defined $INSTANCE){
		$INSTANCE = $caller->SUPER::new(@args);
		bless $INSTANCE, $caller;
		
		$INSTANCE->{'cache'} = ();
	}
	
	return $INSTANCE;
}

=head2 pack_order

 Title   : packed
 Usage   : $arrayRef = $distanceObj->packed($geneOrder);
 Function: Returns the permutation of the gene order as an array of packed arrays of shorts
 Returns : An array reference

=cut

sub pack_order {
	my ($self,$order) = @_;
	
	my (@packed,%key);
	
	my @pi = $order->pi;
	foreach my $pi (@pi){
		my $circular = $pi->is_circular;
		
		my @p = $pi->pi;
		
		unshift @p, $circular;
		push @packed, pack("s*",@p);
	}
	
	return \@packed;
}

=head2 pack_order_reduce

 Title   : packed
 Usage   : $arrayRef = $distanceObj->packed($geneOrder);
 Function: Returns the permutation of the gene order as an array of packed arrays of shorts
 Returns : An array reference

=cut

sub pack_order_reduce {
	my ($self,$order) = @_;
	
	my (@packed,%key);
	
	my @pi = $order->pi;
	
	map(@key{map(abs($_), $_->pi)} = (), @pi);
	
	my $i = 1;
	
	map($key{$_} = $i++, sort {$a <=> $b} keys %key);
	foreach my $pi (@pi){
		my $circular = $pi->is_circular;
		
		my @p = $pi->pi;
		
		map($_ = $_/abs($_)*$key{abs($_)}, @p);
		unshift @p, $circular;
		push @packed, pack("s*",@p);
	}
	
	return \@packed;
}

=head2 adjacencies

 Title   : adjacencies
 Usage   : $adjacencies = $distanceObj->adjacencies($geneOrderA,$geneOrderB);
 Function: Returns the number of adjacencies between two gene orders
 Returns : A scalar value

=cut

sub adjacencies {
	my ($self,$orderA,$orderB) = @_;
	
	my $adjacencies; 
	
	if( defined $self->_cache('adjacencies')->{"$orderA"}{"$orderB"} ){
		$adjacencies = $self->_cache('adjacencies')->{"$orderA"}{"$orderB"};
	}else{
		$adjacencies = adjacencies_xs($self->pack_order($orderA),$self->pack_order($orderB));
	
		$self->_cache('adjacencies')->{"$orderA"}{"$orderB"} = $adjacencies;
	}
	
	my @adjacencies = unpack("s*",$adjacencies) if defined $adjacencies;
	my @adj_list = ();
	
	while(@adjacencies){
		my @adj = ( $SWITCH{abs($adjacencies[0])/$adjacencies[0]}.$orderA->{'key'}->{'name'}->{abs($adjacencies[0])},
					$SWITCH{abs($adjacencies[1])/$adjacencies[1]}.$orderA->{'key'}->{'name'}->{abs($adjacencies[1])} );
		push @adj_list, \@adj;
		shift @adjacencies;
		shift @adjacencies;
	}
	
	return @adj_list;
}

=head2 breakpoints

 Title   : breakpoints
 Usage   : $breakpoints = $distanceObj->breakpoints($geneOrderA,$geneOrderB);
 Function: Returns the number of breakpoints between two GeneOrder objects.
 Returns : Scalar value

=cut

sub breakpoints {
	my ($self,$orderA,$orderB) = @_;

	my $breakpoints;
	
	if( defined $self->_cache('breakpoints')->{"$orderA"}{"$orderB"} ){
		$breakpoints = $self->_cache('breakpoints')->{"$orderA"}{"$orderB"};
	}else{
		$breakpoints = breakpoints_xs($self->pack_order($orderA),$self->pack_order($orderB));
			  	
		$self->_cache('breakpoints')->{"$orderA"}{"$orderB"} = $breakpoints;
	}
	
	return $breakpoints;
}

=head2 inversions

 Title   : inversions
 Usage   : $inversions = $distanceObj->inversions($geneOrderA,$geneOrderB);
 Function: Returns the number of inversions between two GeneOrder objects.
 Returns : Scalar value

=cut

sub inversions {
	my ($self,$orderA,$orderB) = @_;

	my $inversions;
	
	if( defined $self->_cache('inversions')->{"$orderA"}{"$orderB"} ){
		$inversions = $self->_cache('inversions')->{"$orderA"}{"$orderB"};
	}else{
		$inversions = inversions_xs($self->pack_order_reduce($orderA),$self->pack_order_reduce($orderB));
				
		$self->_cache('inversions')->{"$orderA"}{"$orderB"} = $inversions;
	}
	
	return $inversions;
}

=head2 DCJ

 Title   : DCJ
 Usage   : $inversions = $distanceObj->DCJ($geneOrderA,$geneOrderB);
 Function: Returns the Double-Cut and Join distance between two GeneOrder objects.
 Returns : Scalar value

=cut

sub DCJ {
	my ($self,$orderA,$orderB) = @_;

	my $DCJ;
	
	if( defined $self->_cache('DCJ')->{"$orderA"}{"$orderB"} ){
		$DCJ = $self->_cache('DCJ')->{"$orderA"}{"$orderB"};
	}else{
		$DCJ = DCJ_xs($self->pack_order_reduce($orderA),$self->pack_order_reduce($orderB));
		$self->_cache('DCJ')->{"$orderA"}{"$orderB"} = $DCJ;
	}
	
	return $DCJ;
}

=head2 _cache

 Title   : _cache
 Usage   : $breakpoints = $distanceObj->_cache('adjacencies');
 Function: Returns a reference to the cached hash of distance metrics 
           for gene orders that have been passed to this object.
 Returns : Hash reference

=cut

sub _cache {
	my ($self,$arg) = @_;

	if($arg eq 'clear'){
		$self->{'cache'} = {};
		return;
	}
	unless(defined $self->{'cache'}->{$arg}){
		$self->{'cache'}->{$arg} = {};
	}
	
	return $self->{'cache'}->{$arg};
}

=head2 supported_distances

 Title   : supported_distances
 Note    : Get a list of distances supported by this module

=cut

sub supported_distances {
	my $self = shift;

	return @DISTANCES;
}

1;
