#
# BioPerl module for Bio::GeneOrder::Set
#
#	Based on code written by Dennis Lavrov
#	Adapted for Bioperl by Walker Pett
#
# Copyright Dennis Lavrov, Walker Pett
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::GeneOrder::Set - An object for reading/writing/manipulating sets of gene orders

=head1 SYNOPSIS

    use Bio::GeneOrder;
	use Bio::GeneOrder::Set;
    use Bio::Seq;

	#Gene orders can be created using a Bio::SeqI compliant object 
	#containing at least one feature of type CDS,tRNA,mRNA,or rRNA.
    
    my $go  = Bio::GeneOrder->new($seqobj);
	my $go2  = Bio::GeneOrder->new($seqobj2);
	my $go3  = Bio::GeneOrder->new($seqobj3);
	...
    
    #Sets can be created from multiple gene orders
    
    my $goSet = Bio::GeneOrder::Set->($go,$go2,$go3 ...);

	#You can add/remove orders to a set individually
	#via the traditional array functions

	$goSet->push( $go4);
	$goSet->shift();

	#You can also remove orders matching certain criteria
	#Like the number of boundaries they share with other orders
	#in the set, or the number of genes they have, etc.
	
	$goSet->filter_orders( -unique => 1);
	$goSet->filter_orders( -min_shared => 3);
    
=head1 DESCRIPTION

Bio::GeneOrder::Set is an object that is used to represent the order of genes
in a sequence.  It stores information about the transcriptional polarity
for each gene as well as gene boundaries within the gene order.  Use this
object in conjunction with Bio::GeneOrderSet and Bio::GeneOrderSetIO 
to construct sets of gene orders that can be prepared for analysis
with software like PAUP and GRAPPA.

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

package Bio::GeneOrder::Set;

use strict;
use Bio::GeneOrder;
use Bio::GeneOrder::Distance;
use Storable;

use base qw(Bio::Root::Root);
use vars qw(%OFILTER %GFILTER);

BEGIN {
	%OFILTER = ();
	%GFILTER = ();
}
    

=head2 new

 Title   : new
 Usage   : $obj = new Bio::GeneOrder::Set( $order1, $order2, $order3 ... );
 Function: Returns a new Bio::GeneOrder::Set object 
 Returns : A Bio::GeneOrder::Set object
 Args    : Takes a list of 0 or more Bio::GeneOrder objects
           OR
           -file        => retrieves Set object stored in a file using the save method.

=cut

sub new {
	my ($caller, @args) = @_;
	
	my $self;
	
	if($args[0] eq '-file'){
		$caller->throw("file argument provided but with an undefined value") 
			unless $args[1];
		$self = retrieve($args[1]) || $caller->throw("Set object could not be opened from $args[1]");
		
		my @orders = $self->orders;
		for(my $i = 0; $i<scalar @orders;$i++){
			$self->{'indices'}{ $orders[$i]->name } = $i;
		}
		
		return $self;
	}
		

	$self = $caller->SUPER::new(@args);
	bless $self, $caller;

	if($args[0] eq '-verbose'){
		shift @args;
		$self->{'verbose'}=1;
	}

	$self->{'orders'} = ();
	map( $_->filtered(0), @{$self->{'orders'}});
	
	$self->throw("at least one GeneOrder object is required to initialize a Set object") 
		unless( @args);
	
	$self->push(@args);
	
	$self->{'distance'} = Bio::GeneOrder::Distance->new();

	return $self;
}

=head2 orders

 Title   : orders
 Usage   : my $order = $geneOrderSet->orders('name');	   
 Function: Returns the gene order object with the requested name
           or an array of gene orders if no name is provided.
           With no arguments returns unfiltered gene orders.
 Returns : A Bio::GeneOrder object or an array of Bio::GeneOrder objects
 Args    : -all             => Returns all orders, filtered or unfiltered.
           -filtered        => Returns only filtered gene orders
           -name            => Returns the gene order with the specified name.

=cut

sub orders {
	my ($self,@args) = @_;
	
	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	
	$self->throw("name argument provided, but with an undefined value") 
		if( !defined $param{'-name'} && exists $param{'-name'});
	
	if(defined $param{'-name'} ){
		return $self->{ $param{'-name'} };
	}
	
	my @return;
	if(defined $param{'-all'} ){
		@return = @{ $self->{'orders'} };
	}elsif(defined $param{'-filtered'}){
		@return = grep( $_->filtered, @{ $self->{'orders'} } );
	}else{
		@return = grep( $_->filtered == 0, @{ $self->{'orders'} } );
	}
	
	return sort {$a->name cmp $b->name} @return;
}

=head2 index

 Title   : index
 Usage   : my $index = $geneOrderSet->index('name');	   
 Function: Returns the index of the gene order object with the requested name
 Returns : A Bio::GeneOrder object or an array of Bio::GeneOrder objects

=cut

sub index {
	my ($self,$name) = @_;
	
	$self->throw("name argument not provided") 
		if( !defined $name );
		
	return $self->{'indices'}{ $name };
		
}

=head2 genes

 Title   : genes
 Usage   : my @genes = $geneOrderSet->genes;	   
 Function: Returns the unique gene names in this set.  
           With no arguments, returns unfiltered gene names.
 Returns : An array of scalars
 Args    : -all             => Returns all genes, filtered or unfiltered.
           -filtered        => Returns only filtered genes.
=cut

sub genes {
	my ($self,@args) = @_;
	
	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	
	$self->throw("name argument provided, but with an undefined value") 
		if( !defined $param{'-name'} && exists $param{'-name'});
		
	my %genes;
	map( @genes{$_->genes(-all=>1)} = (),$self->orders);
	
	my @genes = keys %genes;
	
	if(defined $param{'-name'}){
		return sort grep($_->name eq $param{'-name'}, @genes);
	}
	
	if(defined $param{'-all'} ){
		return sort @genes;
	}elsif(defined $param{'-filtered'}){
		return sort grep($self->{'key'}->{'filt'}->{$_}, @genes);
	}else{
		return sort grep($self->{'key'}->{'filt'}->{$_} == 0, @genes);
	}
}

=head2 flush

 Title   : flush
 Usage   : my $bool = $geneOrderSet->flush();	   
 Function: Returns true if all gene orders in the set are of the same length
 Returns : Scalar value

=cut

sub flush {
	
	my @orders = shift->orders;
	
	my $no = $orders[0]->no_genes;
	
	return grep($_->no_genes != $no, @orders) ? 0 : 1;
}

=head2 no_orders

 Title   : no_orders
 Usage   : my $no_orders = $geneOrderSet->no_orders();	   
 Function: Returns the number of unfiltered gene order objects in the set
 Returns : Scalar value

=cut

sub no_orders {
	
	my @count = CORE::shift->orders;
	return scalar @count;
}

=head2 push

 Title   : push
 Usage   : $orderSet->push( $order1);
 Function: Adds a GeneOrder object to the GeneOrder set
 Returns : Boolean value if successful.
 Args    : Takes a list of 0 or more Bio::GeneOrder objects.

=cut

sub push {
	my ($self,@orders) = @_;

	my $i = 1;
	foreach my $order (@orders){
		$self->throw("arguments must be of class Bio::GeneOrder") 
			if ref($order) ne 'Bio::GeneOrder';

		my $nameA = $order->name;
		
		if( defined $self->{ $nameA }){
			$nameA .= "_".$self->{'no_orders'};
			$order->name($nameA);
		}

		CORE::push @{ $self->{'orders'} }, $order;
		
		$self->{ $nameA } = $order;
		$self->{'no_orders'}++;
	}
	
	#update key
	my (@genes,%types);
	map( eval{ push @genes, $_->genes(-all => 1) }, @orders);
	map( eval{ @types{ keys %{ $_->_key()->{'type'} } } = values %{ $_->_key()->{'type'} } }, @orders);
	
	my %saw;
	@saw{@genes} = ();
	@genes = keys %saw;
	
	map( $self->{'key'}->{'filt'}->{$_} = 0, @genes);
	
	$i = 1;
	map($self->{'key'}->{'index'}->{$_} = $i++, @genes);
	%{ $self->{'key'}->{'name'} } = reverse(%{ $self->{'key'}->{'index'} });
	$self->{'key'}->{'type'} = \%types;
	
	$self->_key($self->{'key'});
	
	#update filters
	$self->filter_genes(%GFILTER) if %GFILTER;
	$self->filter_orders(%OFILTER) if %OFILTER;
	
	#update index values
	my @orders = $self->orders;
	$i = 0;
	map( $self->{'indices'}{ $_->name } = $i++, @orders);
	
	return 1;
}

=head2 filter_orders

 Title   : filter_orders
 Usage   : $orderSet->filter_orders( -min_shared    =>  3,
                                     -min_neighbors => 4  );
 Function: Markes gene orders in the set as filtered based on the specified criteria.
 Returns : An array of Bio::GeneOrder objects
 Args    : -invert          => Filters those orders to which the parameters do not apply.
                               if this is the only option provided, inverts the current filter.
           -unfilter		=> Unfilters those orders to which the parameters apply.
           -name            => Removes orders whose names match a given name or regular expression.
           -unique          => Removes all but one gene order of those that are identical.
           					   if used with the -name argument, removes those gene orders that
           					   are identical to the gene order specified by -name
           -flush			=> Ensures all orders have the same number of genes.
           -min_genes		=> Removes any order that has less than a certain number of genes.
           -max_genes		=> Removes any order that has more than a certain number of genes.
		   -max_copies      => Removes any orders that have more than a specified number of copies
                               of any gene. If used with -name, applies only to the gene with the
                               indicated name.
           -min_copies      => Removes any orders that have less than a specified number of copies
                               of any gene. If used with -name, applies only to the gene with the
                               indicated name.
           -max_{distance}	=> Removes all orders whose {distance} is greater than a certain 
                               number of rearrangements. If used with -name, removes those gene orders that
           					   are greater than the specified distance from the gene order specified by -name.
           -min_{distance}  => Removes all orders whose {distance} is less than a certain 
                               number of rearrangements. If used with -name, removes those gene orders that
           					   are less than the specified distance from the gene order specified by -name.
           -cluster_size    => Specifies the minimum number of orders that constitute a viable
                               group when applying upper and lower limits on distances.
                               For example, if -cluster_size is set to 2, and -min_inversions 
                               is set to 4, then any orders that are less than 4 inversions 
                               from less than 2 other orders will be filtered.

=cut

sub filter_orders {
	my ($self,@args) = @_;
	
	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	
	$self->throw("name argument provided, but with an undefined value") 
		if( !defined $param{'-name'} && exists $param{'-name'});
		
	$self->throw("invert argument provided, but with an undefined value") 
		if( !defined $param{'-invert'} && exists $param{'-invert'});

	$self->throw("min_shared argument provided, but with an undefined value") 
		if( !defined $param{'-min_shared'} && exists $param{'-min_shared'});
	
	$self->throw("max_shared argument provided, but with an undefined value") 
		if( !defined $param{'-max_shared'} && exists $param{'-max_shared'});

	$self->throw("max_genes argument provided, but with an undefined value") 
		if( !defined $param{'-max_genes'} && exists $param{'-max_genes'});

	$self->throw("min_genes argument provided, but with an undefined value") 
		if( !defined $param{'-min_genes'} && exists $param{'-min_genes'});

	$self->throw("max_copies argument provided, but with an undefined value") 
		if( !defined $param{'-max_copies'} && exists $param{'-max_copies'});

	$self->throw("min_copies argument provided, but with an undefined value") 
		if( !defined $param{'-min_copies'} && exists $param{'-min_copies'});
	
	$self->throw("max_breakpoints argument provided, but with an undefined value") 
		if( !defined $param{'-max_breakpoints'} && exists $param{'-max_breakpoints'});
	
	$self->throw("min_breakpoints argument provided, but with an undefined value") 
		if( !defined $param{'-min_breakpoints'} && exists $param{'-min_breakpoints'});

	$self->throw("cluster_size argument provided, but with an undefined value") 
		if( !defined $param{'-cluster_size'} && exists $param{'-cluster_size'});

	$self->throw("cluster_size argument provided, but without any additional limit arguments") 
		if( defined $param{'-cluster_size'} && !defined $param{'-min_genes'} && 
			!defined $param{'-max_genes'} && !defined $param{'-min_copies'} && 
			!defined $param{'-max_copies'} && 
				!grep(defined $param{"-max_$_"}, Bio::GeneOrder::Distance->supported_distances) &&
				!grep(defined $param{"-min_$_"}, Bio::GeneOrder::Distance->supported_distances) );

	$self->throw("max_shared must be an integer") 
		if( defined $param{'-max_shared'} && $param{'-max_shared'} =~ /D/);
	$self->throw("min_shared must be an integer") 
		if( defined $param{'-min_shared'} && $param{'-min_shared'} =~ /D/);
	$self->throw("max_genes must be an integer") 
		if( defined $param{'-max_genes'} && $param{'-max_genes'} =~ /D/);
	$self->throw("min_genes must be an integer") 
		if( defined $param{'-min_genes'} && $param{'-min_genes'} =~ /D/);
	$self->throw("max_breakpoints must be an integer") 
		if( defined $param{'-max_breakpoints'} && $param{'-max_breakpoints'} =~ /D/);
	$self->throw("min_breakpoints must be an integer") 
		if( defined $param{'-min_breakpoints'} && $param{'-min_breakpoints'} =~ /D/);
	$self->throw("max_copies must be an integer") 
		if( defined $param{'-max_copies'} && $param{'-max_copies'} =~ /D/);
	$self->throw("min_copies must be an integer") 
		if( defined $param{'-min_copies'} && $param{'-min_copies'} =~ /D/);
	
	$self->throw("cluster_size must be a positive integer") 
		if( defined $param{'-cluster_size'} && ($param{'-cluster_size'} <= 0 || $param{'-cluster_size'} =~ /D/));
	
	#Add our options to the filter
	@OFILTER{ keys %param } = values %param;
	
	#But we don't include invert in the stored filter
	delete $OFILTER{'-invert'};
	
	if($param{'-unfilter'}){
		delete @OFILTER{ keys %param };
	}
	
	my @filtered = ();
	
	if(scalar(keys %param) == 1 && (keys %param)[0] eq '-invert'){
		#If invert is our only option provided, invert the filter
		map $_->filtered( $_->filtered ? 0 : 1 ), @{ $self->{'orders'} };
		
	}elsif(%param){
		#-name
		if(defined $param{'-name'} && !defined $param{'-max_copies'} && !defined $param{'-min_copies'}){
			if($param{'-name'} =~ /^\/.*\/(.*)/){
				#Delete genes that match the regexp passed
				my $tags = $1;
				my $name =~ s/^\/([^\/]*)\/.*/$1/;
					
				my $regexp;
				if($tags =~ /i/){
					$regexp = qr/$name/i
				}else{
					$regexp = qr/$name/;
				}
				
				my @matched = grep($_->name =~ $regexp, @{ $self->{'orders'} });
				foreach my $match (@matched){
					if($param{'-unique'}){
						my @unique = grep( $match->breakpoints($_) == 0, 
												  $self->orders);
						foreach(@unique){
							CORE::push @filtered, $self->filter_orders( '-name' => $_->name );
						}
					}else{
						CORE::push @filtered, $match;
						if($param{'-unfilter'}){
							$match->filtered(0);
						}else{
							$match->filtered(1);
						}
					}
				}
			}elsif($param{'-unique'}){
				my @unique = grep( $self->{ $param{'-name'} }->breakpoints($_) == 0, 
										  $self->orders);
				foreach(@unique){
					CORE::push @filtered, $self->filter_orders( '-name' => $_->name );
				}
			}else{
				my @f;
				if(my @d = grep(defined $param{"-min_$_"}, Bio::GeneOrder::Distance->supported_distances)){
					foreach my $d (@d){
						CORE::push @f, grep( $self->distance->$d($self->{ $param{'-name'} },$_) < $param{"-min_$d"}, $self->orders);
					}
				}
				if(my @d = grep(defined $param{"-max_$_"}, Bio::GeneOrder::Distance->supported_distances)){
					foreach my $d (@d){
						CORE::push @f, grep( $self->distance->$d($self->{ $param{'-name'} },$_) > $param{"-max_$d"}, $self->orders);
					}
				}
				foreach(@f){
					CORE::push @filtered, $self->filter_orders( '-name' => $_->name );
				}
				
				#If we haven't performed any distance filtering, then simply filter the orders that match the name
				unless(@f){
					CORE::push @filtered, $self->{ $param{'-name'} };
					foreach my $match (@filtered){
						if($param{'-unfilter'}){
							$match->filtered(0);
						}else{
							$match->filtered(1);
						}
					}
				}
			}
		}else{
			#my %temp = %OFILTER;
			#If we set a value for -min_neighbors, use it
			my $cluster_size = defined $param{'-cluster_size'} ? $param{'-cluster_size'} : 1;
			
			my %bound_count;
			($bound_count{$_}++) for map ($_->bounds, @{ $self->{'orders'} });
			my $max_bound_count = (sort {$bound_count{$b} <=> $bound_count{$a}} keys %bound_count)[0];
			
			my $matched =0;
			foreach my $order (@{ $self->{'orders'} }){
				#If this order has already been filtered, we need not filter it again
				unless( $order->filtered ){
					#-unique
					if( $param{'-unique'}){
						my $num_identical = grep( $order->breakpoints($_) == 0, 
												  $self->orders);
		
						if( $num_identical > 1){
								$matched++;
								CORE::push @filtered, $self->filter_orders( '-name' => $order->name, );
						}
					}
			
					#-min_genes
					if( defined $param{'-min_genes'}){
						if($order->no_genes < $param{'-min_genes'}){
							$matched++;
							CORE::push @filtered, $self->filter_orders( '-name' => $order->name, )
								unless($param{'-invert'});
						}
					}
					
					#-max_genes
					if( defined $param{'-max_genes'}){
						if($order->no_genes > $param{'-max_genes'}){
							$matched++;
							CORE::push @filtered, $self->filter_orders( '-name' => $order->name, )
								unless($param{'-invert'});
						}
					}
					
					#-max_copies
					if(defined $param{'-max_copies'}){
						if(defined $param{'-name'}){
							if($param{'-max_copies'} < $order->genes(-name => $param{'-name'}) ){
								$matched++;
								CORE::push @filtered, $self->filter_orders( '-name' => $order->name, )
									unless($param{'-invert'});
							}
						}else{
							my $num_copied;
							my @genes = $order->genes;
							foreach my $gene ($self->genes){
								my @named = grep($_ eq $gene, @genes);
								$num_copied++ if $param{'-max_copies'} < scalar(@named);
							}
							
							if($num_copied){
								$matched++;
								CORE::push @filtered, $self->filter_orders( '-name' => $order->name, )
									unless($param{'-invert'});
							}
						}
					}
			
					#-min_copies
					if(defined $param{'-min_copies'}){
						if(defined $param{'-name'}){
							if($param{'-min_copies'} > $order->genes(-name => $param{'-name'}) ){
								$matched++;
								CORE::push @filtered, $self->filter_orders( '-name' => $order->name, )
									unless($param{'-invert'});
							}
						}else{
							my $num_copied;
							my @genes = $order->genes;
							foreach my $gene ($self->genes){
								my @named = grep($_ eq $gene, @genes);
								$num_copied++ if $param{'-min_copies'} > scalar(@named);
							}
							
							if($num_copied){
								$matched++;
								CORE::push @filtered, $self->filter_orders( '-name' => $order->name, )
									unless($param{'-invert'});
							}
						}
					}
			
					#-flush
					if( $param{'-flush'}){
						if($max_bound_count != $order->bounds){
							$matched++;
							CORE::push @filtered, $self->filter_orders( '-name' => $order->name, )
								unless($param{'-invert'});
						}
					}
					
					#-min_{distance}
					if(my @d = grep(defined $param{"-min_$_"}, Bio::GeneOrder::Distance->supported_distances)){
						#For each gene order, count the number of other gene orders 
						#for which the specified arguments apply
						
						foreach my $d (@d){
							my $num_neighbors = grep( $self->distance->$d($order,$_) < $param{"-min_$d"}, $self->orders);
							
							if( ($num_neighbors-1) > $cluster_size){
								$matched++;
								CORE::push @filtered, $self->filter_orders( '-name' => $order->name, )
									unless($param{'-invert'});
							}
						}
					}
			
					#-max_{distance}
					if(my @d = grep(defined $param{"-max_$_"}, Bio::GeneOrder::Distance->supported_distances)){
						#For each gene order, count the number of other gene orders 
						#for which the specified arguments apply
						
						foreach my $d (@d){
							my $num_neighbors = grep( $self->distance->$d($_,$order) > $param{"-max_$d"}, $self->orders);
						
							if( $num_neighbors > $cluster_size){
								$matched++;
								CORE::push @filtered, $self->filter_orders( '-name' => $order->name, )
									unless($param{'-invert'});
							}
						}
					}
					
					if($param{'-invert'} && $matched == 0){
						CORE::push @filtered, $self->filter_orders( '-name' => $order->name, );
					}
				}
			}
			
			#%OFILTER = %temp;
		}
	}else{
	#If there are no options provided
	#clear the filter
		map( $_->filtered(0), @{$self->{'orders'}});
		%OFILTER = ();
	}
	
	#update index values
	my $i = 0;
	$self->{'indices'} = {};
	map $self->{'indices'}{$_->name} = $i++ , $self->orders;
		
	return @filtered;
}

=head2 filter_genes

 Title   : filter_genes
 Usage   : my @filtered = $geneOrderSet->filter_genes( -type => 'gene type',
                                                       -name => /regexp/    );
 Function: Filters genes of specified name/type from gene order set and returns the number filtered
 Returns : Scalar value
 Args	 : -type    =>  A sequence feature primary  tag, ie. 'gene', 'tRNA', 'mRNA', 'rRNA' etc.
           -name    =>  A string or regular expression that identifies gene names to filter
                        ie. To match an exact string use a bareword: 'name'.
                        To match an inexact string, supply a regular expression: /'name'/i
           -local   =>  Removes genes from the set that are not represented in every gene order.
                        

=cut

sub filter_genes {
	my ($self,@args) = @_;

	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	
	$self->throw("name argument provided, but with an undefined value") 
		if( !defined $param{'-name'} && exists $param{'-name'});
	
	$self->throw("type argument provided, but with an undefined value") 
		if( !defined $param{'-type'} && exists $param{'-type'});
	
	@GFILTER{ keys %param } = values %param;
	
	my @matched;
	my @genes = keys %{$self->{'key'}->{'index'}};
	my $filtered;
	
	if(%param){
		if( defined $param{'-name'} && $param{'-name'} =~ /^\/.*\/(.*)/){
			#Delete genes that match the regexp passed
			my $tags = $1;
			my $name = $param{'-name'};
			$name =~ s/^\/([^\/]*)\/.*/$1/;
				
			my $regexp;
			if($tags =~ /i/){
				$regexp = qr/$name/i
			}else{
				$regexp = qr/$name/;
			}
			
			if($param{'-invert'}){
				@matched = grep( $_ !~ $regexp, @genes);
			}else{
				@matched = grep( $_ =~ $regexp, @genes);
			}
			
		}elsif(defined $param{'-name'}){
			#Delete genes that equal the name passed
			if($param{'-invert'}){
				@matched = grep( $_ ne $param{'-name'}, @genes);
			}else{
				@matched = grep( $_ eq $param{'-name'}, @genes);
			}
		}elsif( defined $param{'-type'}){
			#Delete genes that equal the type passed
			if($param{'-invert'}){
				@matched = grep( $self->{'key'}->{'type'}->{$_} ne $param{'-type'}, @genes);
			}else{
				#map( eval{ print "$_: ".$self->{'key'}->{'type'}->{$_}."\n"}, @genes);
				@matched = grep( $self->{'key'}->{'type'}->{$_} eq $param{'-type'}, @genes);
			}
		}
		
		if(@matched){
			if($param{'-unfilter'}){
				map( $self->{'key'}->{'filt'}->{$_} = 0, @matched);
			}else{
				map( $self->{'key'}->{'filt'}->{$_} = 1, @matched);
			}
			$filtered = 1;
		}
		
	}elsif( keys %{ $self->{'key'}->{'filt'} } ){
		map( $self->{'key'}->{'filt'}->{$_} = 0, @genes);
		$filtered = 1;
	}
	
	if($filtered){
		$self->_key($self->{'key'});
		$self->distance->_cache('clear');
	}
	
	return @matched;
}

=head2 ofilter

 Title   : ofilter
 Usage   : my %filter = $geneOrderSet->ofilter();
 Function: Returns the filter currently set on gene orders in the set
 Returns : A hash of parameters for the filter_orders method
                        

=cut

sub ofilter {

	return %OFILTER;
}

=head2 gfilter

 Title   : gfilter
 Usage   : my %filter = $geneOrderSet->gfilter();
 Function: Returns the filter currently set on genes in the set
 Returns : A hash of parameters for the filter_genes method
                        

=cut

sub gfilter {

	return %GFILTER;
}

=head2 rename_genes

 Title   : rename_genes
 Usage   : my $genes = $geneOrderSet->rename_genes( -list =>  @list_of_names,
                                                    -list =>  @list_of_names2,
                                                    -list =>  ...      );
         OR
           my $genes = $geneOrderSet->rename_genes( -table => $filename   );

 Function: Renames genes from names matched in a list to the first name in the list.
		   Any number of lists may be supplied to this function.
		   Alternatively, you may rename genes using a synonym table read from a file.
		   Returns an the number of renamed genes.  
 Returns : A scalar value.
 Args	 : -list    =>  A list of names to match to genes that will be renamed.
           -table   =>  The name of a file containing a synonymous gene name table.
                        The first line of this file should read ">Bio::GeneOrder::synonymous_genes".
                        Subsequent lines should contain pipe "|" delimited lines of synonymous 
                        gene names.  The first name in each line will be used as the primary name
                        for that gene.

=cut

sub rename_genes {
	my ($self,@args) = @_;

	#our list to return
	my @renamed = ();
	#our list of argument lists
	my @lists =();

	if( grep $_ eq '-table' || $_ eq 'table', @args){
		$self->throw("-table argument requires exactly one table name") 
			unless( scalar @args == 2);
		$self->throw("-table argument supplied incorrectly") 
			unless( $args[0] eq '-table');

		open TABLE, $args[1] or die "could not open synonym table file $args[1]: $!";
		
		my $entry = <TABLE>;
		$self->throw("synonym table not formatted correctly") 
			unless( $entry =~ />Bio::GeneOrder::synonyms/);

		while(<TABLE>){
			#ignore blank lines
			next if !($_ =~ /\S/);
			#push each line as an array of names onto our list of lists
			chomp $_;
			my @list = split ';',$_;

			$self->throw("lists require at least two synonymous gene names") 
				unless( scalar @list >= 2);
			
			push @lists, \@list;
		}
		
		close TABLE or die "could not close synonym table file $args[1]: $!";
			
	#if lists are provided
	}elsif( grep $_ eq '-list' || $_ eq 'list', @args){
		$self->throw("-list argument supplied incorrectly") 
			unless( $args[0] eq '-list');
		
		#initialize some vars
		my $push_list = 0;
		my $list_count = 0;
		my @list = ();

		#go through each argument
		for(my $i=0;$i< scalar @args;$i++){
			#at each '-list' argument, if we've already counted some list items
			#mark them for adding to our list of lists
			#allow user to use either '-list' or 'list'
			if($args[$i] eq '-list' || $args[$i] eq 'list'){
				$push_list = $list_count ? 1 : 0;
			}else{
				#push each synonym into our list
				push @list, $args[$i];
				$list_count++;
			}

			#push our list argument onto our list of lists
			if( $push_list || $i == scalar @args -1){
				$self->throw("-list argument requires at least two synonymous names") 
					if $list_count<2;
				
				my @pusher = @list;
				push @lists, \@pusher;

				$push_list = 0;
				$list_count = 0;
				@list = ();
			}
		}
	}else{
		$self->throw("rename_genes requires one argument of type list or table") ;
	}
	
	my @genes = keys %{$self->{'key'}->{'index'}};
	
	my $map = {};

	#rename each gene from our list of argument lists
	foreach my $list (@lists){
		my $rename = shift @$list;
		$rename =~ s/\t//g;
		
		last unless @genes;

		foreach my $name ( @$list){
		
			my (@genesA,@genesB);
			
			if($name =~ /^\/.*\/(.*)/){
				my $tags = $1;
				$name =~ s/^\/([^\/]*)\/.*/$1/;
				
				my $regexp;
				if($tags =~ /i/){
					$regexp = qr/$name/i;
				}else{
					$regexp = qr/$name/;
				}

				foreach(@genes){
					if($_ =~ $regexp){
						push @genesA, $_;
						$map->{$_} = $rename;
					}else{
						push @genesB, $_;
						$map->{$_} = $_;
					}
				}
			}else{
				foreach(@genes){
					if($_ eq $name){
						push @genesA, $_;
						$map->{$_} = $rename;
					}else{
						push @genesB, $_;
						$map->{$_} = $_;
					}
				}
			}
			
			if(@genesA){
				@genes = @genesB;
				push @renamed, $rename;
			}
		}
	}
	
	if(@renamed){
		my @new = (@renamed,@genes);
		
		my $new_key = {};
		
		my $i=1;
		foreach(keys %{$self->{'key'}->{'index'}}){
			$new_key->{'index'}->{$map->{$_}} = $i++ unless defined $new_key->{'index'}->{$map->{$_}};
			$new_key->{'name'}->{$new_key->{'index'}->{$map->{$_}}} = $map->{$_};
			$new_key->{'type'}->{$map->{$_}} = $self->{'type'}->{$_};
			$new_key->{'filt'}->{$map->{$_}} = 0;
		}
		
		$self->_key($new_key,$map);
		$self->distance->_cache('clear');
		
		$self->filter_genes(%GFILTER) if %GFILTER;
	}
	
	return @renamed;
}

=head2 reorder

 Title   : reorder
 Usage   : $geneOrderSet->reorder('name');
 Function: Reorders gene orders in the set that contain the specified gene so that the
           specified gene is first and in the same orientation, and returns the number of reordered
           gene orders.
 Returns : Scalar value

=cut

sub reorder {
	my ($self,$name) = @_;
	
	my $num_reordered = 0;
	
	foreach( @{ $self->{'orders'}} ){
		if( my @m = $_->genes(-name => $name) ){
			$num_reordered += $_->reorder($name);
		}
	}
	
	return $num_reordered;
}

=head2 unique

 Title   : unique
 Usage   : $geneOrderSet->unique();
 Function: Filters duplicate gene orders, and returns the list of orders filtered out.
 Returns : Scalar value.

=cut

sub unique {
	my $self = CORE::shift;

	$OFILTER{'-unique'} = 1;
	return $self->filter_orders(%OFILTER);
}

=head2 distance

 Title   : distance
 Usage   : $geneOrderA->distance();
 Function: Get a distance object
 Returns : Bio::GeneOrder::Distance object

=cut

sub distance {
	return shift->{'distance'};
}

=head2 save

 Title   : save
 Usage   : $geneOrderSet->save('filename');
 Function: Save the GeneOrder set object to a file for later recovery.
           Depending on the size of your set, this could be a large file.
 Returns : 1 for success, 0 for failure.

=cut

sub save {
	my ($self,$file) = @_;

	store $self, $file || $self->throw("GeneOrder::Set could not be saved to $file");

	return 1;
}

sub max_genes {
	my $self = CORE::shift;

	my $max_genes = 0;
	foreach my $order ($self->orders){
		$max_genes = $order->no_genes > $max_genes ? $order->no_genes : $max_genes;
	}
	
	return $max_genes;
}

=head2 _sort_orders

 Title   : _sort_orders
 Usage   : $orderSet->_sort_orders();
 Function: Sorts the internal array of gene orders alphabetically.

=cut

sub _sort_orders {
	my $self = shift;

	@{ $self->{'orders'}} = sort {$a->name cmp $b->name} @{ $self->{'orders'}};

}

=head2 _key

 Title   : _key
 Usage   : $orderSet->_key();
 Function: Get/set the gene name/number key.

=cut

sub _key {
	my ($self,$key,$map) = @_;

	
	if(defined $key){
		map( $_->_key($key,$map), $self->orders);
		
		$self->{'key'} = $key;
	}
	
	return $self->{'key'};
}

1;
	
	

	

	
