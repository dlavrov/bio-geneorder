#
# BioPerl module for Bio::GeneOrder
#
#	Based on code written by Dennis Lavrov
#	Adapted for Bioperl by Walker Pett
#
# Copyright Dennis Lavrov, Walker Pett
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

our $VERSION = '0.4';

=head1 NAME

Bio::GeneOrder - An object for obtaining and manipulating gene orders

=head1 SYNOPSIS

    use Bio::GeneOrder;
    use Bio::Seq;

    #Gene orders can be created using a Bio::SeqI compliant object 
    #containing at least one feature of type "gene".
    
    my $go  = Bio::GeneOrder->new($seqobj);
    
    #All genes of a given type can be purged from the gene order
    
    $go->purge(-type => "tRNA");
    
    #The source sequence object can be retrieved from a GeneOrder object
    #for use with other Bioperl objects.  If the GeneOrder object was created
    #from a NEXUS matrix or from a GRAPPA file, it will likely not contain
    #a source sequence object, and this method will raise an exception.
    
    my $source = $go->seq();

    #Retrieve the number of genes or gene boundaries in the order

    my $no_genes = $go->no_genes;
    my $no_bounds = $go->bounds;

    #Retrieve a gene by name or by number or simply retrieve an
    #array of all the genes

    $geneObj = $go->genes('cox2');   #Access genes directly
    my @genes = $go->genes;     #Returns an array

    #All genes of a particular type can be removed from the order
    #where the type is one of CDS, tRNA, mRNA or rRNA
	
    my @purged_genes = $go->purge( -type  =>  'tRNA' );

    #Gene orders can be randomized for performing sampling
    #of random orders with equal gene content
	
    my $random_order = $go->shuffle();

    #Genes can be renamed to account for synonymous gene names
    #The first name in each list provided is used as the
    #primary name.

    $go->rename_genes( -list => 'cox1','COI','COX1',
                       -list => 's-rRNA','12S ribosomal RNA','12S rRNA' );

    #A file with a table of synonymous names can also be used
    #Refer to the method description for how this file should be
    #formatted.

    $go->rename_genes( -table => 'synonyms.txt' );
	
    #Two GeneOrder objects can be compared to determine which
    #gene boundaries are shared between them, and how many.
    
    my @shared_bounds = $go->shared_bounds($go2);
    
    my $number_shared = scalar @shared_bounds;
    
=head1 DESCRIPTION

Bio::GeneOrder is an object that is used to represent the order of genes
in a sequence.  It stores information about the transcriptional direction
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

package Bio::GeneOrder;

use warnings;
use strict;
use Bio::Seq;
use Bio::GeneOrder::permutation;
use Bio::GeneOrder::Distance;
use Storable;

use base qw(Bio::Root::Root);
use vars qw(%REV %FILTER %SWITCH $LINEAR);

BEGIN {
	$LINEAR = '~';
	%REV = ( '+' => '-',
			 '' => '-',
			 '-' => '' );
	%SWITCH = ( 1	=> '', -1	=> '-', 0	  => undef,
			  ''	=> 1,   '-'	=> -1 , undef => 0,'+' => 1);
	%FILTER = ();
}

=head2 new

 Title   : new
 Usage   : $obj = new Bio::GeneOrder( $seqObj1, $seqObj2, ...
                                      -name   =>    'optional name' );
           $obj = new Bio::GeneOrder( $orderString1, $orderString2, ...
                                      -name   =>    'required name' );
           $obj = new Bio::GeneOrder( $permutation1, $permutation2, ...
                                      -name   =>    'required name' );
           $obj = new Bio::GeneOrder( $orderString, $seqObj, $permutation ...
                                      -name   =>    'required name' );
		   $obj = new Bio::GeneOrder( -file   =>    'geneOrder.obj' );
 Function: Returns a new Bio::GeneOrder object 
 Returns : A Bio::GeneOrder object initialized with a Bio::SeqI compliant object
 Args    : Accepts an array of Bio::SeqI compliant objects, gene order strings, or
           Bio::GeneOrder::permutation objects with the additional optional arguments.
           -name              => a string representing the name of this gene order
                                 required if -order is provided
                                 [default is species from source sequence]
           -file              => A file name from which to load a GeneOrder object.
                                 These files are created using the save method.

=cut

sub new {
	my ($caller, @args) = @_;
	my $self = $caller->SUPER::new(@args);
	bless $self, $caller;
	
	my (@seqs,@perms,@orders);
	while(defined $args[0] && $args[0] ne '-name' && $args[0] ne '-file'){
		if(ref($args[0]) =~ /Bio::Seq/){
			push @seqs, shift @args;
		}elsif(ref($args[0]) =~ /Bio::GeneOrder::permutation/){
			push @perms, shift @args;
		}else{
			push @orders, shift @args;
		}
	}
	
	$caller->throw("invalid number of arguments or invalid argument order") 
		if scalar @args % 2 != 0;

	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys

	$caller->throw("file argument provided, but with an undefined value") 
		if( !defined($param{'-file'}) && exists($param{'-file'}) );

	if( defined $param{'-file'}){
		return retrieve($param{'-file'}) || $self->throw("GeneOrder object could not be opened from ".$param{'-file'});
	}

	$caller->throw("no sequences provided") 	
		unless(@seqs || @orders);
	$caller->throw("-name argument is required when not providing sequence objects") 	
		unless(@seqs || defined $param{'-name'});
	$caller->throw("name argument provided, but with an undefined value")
		if( exists $param{'-name'} && !defined($param{'-name'}));
	$caller->throw("circular argument provided, but with an undefined value")
		if( exists $param{'-circular'} && !defined($param{'-circular'}));
	
	$self->{'key'} = {};
	$self->{'name'} = $param{'-name'} if defined $param{'-name'};
	$self->is_circular($param{'-circular'}) if defined $param{'-circular'};
	
	my @permutations;
	
	my $i=1;
	if( @seqs){
		my @accessions;
		foreach my $seq (@seqs){
			my @pi;
			
			my $accession = $seq->accession_number();
			my $circular = $seq->is_circular ? 1 : 0;
			
			my @sources = grep( ($_->primary_tag() eq 'CDS' || $_->primary_tag() =~ /(t|r)RNA/), $seq->get_SeqFeatures());
			
			#Sort the genes in the order by their midpoints to adjust for genes that may not appear
			#in order in the sequence object
			my @new_sources;
			
			foreach my $s (@sources){
				foreach( $s->location->each_Location ){
					($s->{'start'},$s->{'end'}) = ($_->start,$_->end);
					my $ss = {};
					%{$ss} = %{$s};
					bless $ss, ref($s);
					push @new_sources, $ss;
				}
			}
			
			@new_sources = sort { ($a->{'start'} + $a->{'end'})/2 <=> ($b->{'start'} + $b->{'end'})/2 } @new_sources;
		
			if(!defined $self->{'name'}){
				if(defined $seq->species()){
					my @species = $seq->species()->classification();
					$self->{'classification'} = \@species;
					$self->{'name'} = $species[0];
					$self->{'name'} =~ s/mitochondrion\s?//i;
				}else{
					$self->{'name'} = ($seq->desc() =~ /^(\S+)\s+(\S+)/);
				}
			}
				
			foreach my $source (@new_sources){
				my $gene_name;
				
				if( $source->has_tag('gene')){
					$gene_name = ($source->get_tag_values('gene'))[0];
				}elsif( $source->has_tag('standard_name')){
					$gene_name = ($source->get_tag_values('standard_name'))[0];
				}elsif( $source->has_tag('product')){
					$gene_name = ($source->get_tag_values('product'))[0];
				}elsif( $source->has_tag('locus_tag')){
					$gene_name = ($source->get_tag_values('locus_tag'))[0];
				}elsif( $source->has_tag('db_xref')){
					$gene_name = ($source->get_tag_values('db_xref'))[0];
				}else{
					$self->warn($source->type." at ".$source->start.'..'.$source->end." in ".$self->{'name'}." has an unknown name");
					$gene_name = 'unk';
				}
		
				unless($gene_name =~ /\(...\)/){
					if( $source->has_tag('codon_recognized')){
						my $codon = ($source->get_tag_values('codon_recognized'))[0];
						unless( $gene_name =~ /$codon/i){
							$gene_name .= "($codon)";
						}
					}elsif( $source->has_tag('note')){
						my $note = ($source->get_tag_values('note'))[0];
						if($note =~ /codon[s]?[\s|_]recognized:\s?(...)/){
							my $codon = $1;
							unless( $gene_name =~ /$codon/i){
								$gene_name .= "($codon)";
							}
						}
					}
				}
				
				$gene_name =~ s/\s/_/g;
				
				unless(defined $self->{'key'}->{'index'}->{$gene_name}){
					$self->{'key'}->{'name'}->{$i} = $gene_name;
					$self->{'key'}->{'index'}->{$gene_name} = $i++;
					$self->{'key'}->{'type'}->{$gene_name} = $source->primary_tag;
				}
				
				push @pi, !defined $source->strand || $source->strand == 0 ? $self->{'key'}->{'index'}->{$gene_name} : $self->{'key'}->{'index'}->{$gene_name} * $source->strand;
			}
			
			push @permutations, Bio::GeneOrder::permutation->new( -circular => $circular,
																  -source   => $accession,
																  -pi		=> \@pi );
			push @accessions, $accession;
		}
	
		$self->{'source'} = join ';' , @accessions;
		$self->{'name'} .= '|'.$self->{'source'};
	}
	
	if(@orders){	
		
		foreach my $order (@orders){
			my @pi;
			my $circular = 1;
			$circular = 0 if($order =~ s/^$LINEAR\s*//);
			
			my @names = map( /^[\-+]?(\S+)/, split / /, $order);
			my @strands = map( /^([\-+])?\S+/, split / /, $order);
			map(eval{$_ = '+' unless $_},@strands);

			for(my $u=0;$u<@names;$u++){
				$names[$u] =~ s/\s/_/g;
				
				unless(defined $self->{'key'}->{'index'}->{$names[$u]}){
					$self->{'key'}->{'name'}->{$i} = $names[$u];
					$self->{'key'}->{'index'}->{$names[$u]} = $i++;
					
					if($names[$u] =~ /trn/i){
						$self->{'key'}->{'type'}->{$names[$u]} = 'tRNA';
					}elsif($names[$u] =~ /rrn/i){
						$self->{'key'}->{'type'}->{$names[$u]} = 'rRNA';
					}else{
						$self->{'key'}->{'type'}->{$names[$u]} = 'CDS';
					}
				}
				
				push @pi, $SWITCH{$strands[$u]} == 0 ? $self->{'key'}->{'index'}->{$names[$u]} : $self->{'key'}->{'index'}->{$names[$u]} * $SWITCH{$strands[$u]};
			}
			
			push @permutations, Bio::GeneOrder::permutation->new( -circular => $circular,-pi => \@pi );
		}
		
	}
	
	if(@perms){
	
		foreach my $permutation (@perms){
		
			my $key = $permutation->_key;
			my %names = keys %{$key->{'index'}};
			
			
		}
	}

	$self->{'pi'} = \@permutations;
	$self->{'name'} =~ s/\s/_/g;
	$self->{'circular'} = ${$self->{'pi'}}[0]->is_circular && @{$self->{'pi'}} == 1 ? 1 : -1;
	map( $self->{'key'}->{'filt'}->{$_} = 0, keys %{$self->{'key'}->{'index'}});
	map( $_->{'key'} = $self->{'key'} , $self->pi );
	$self->filtered(0);
	
	$self->{'distance'} = Bio::GeneOrder::Distance->new();
	
	return $self;
}

=head2 name

 Title   : name
 Usage   : my $name = $geneOrder->name();
 Function: Returns the name of this gene order, usually species name from source sequence
 Returns : A string representing the name of the gene order	

=cut

sub name {
	my ($self, $value) = @_;

	if( defined $value){
		$self->{'name'} = $value;
	}
	
	return $self->{'name'};
}

=head2 classification

 Title   : classification
 Usage   : my $classification = $geneOrder->classification();
 Function: Returns the classification array of the organism represented by this gene order.
 Returns : An array

=cut

sub classification {
	my ($self, $value) = @_;
	
	if(defined $value){
		$self->{'classification'} = $value;
	}
	
	return defined $self->{'classification'} ? @{ $self->{'classification'} } : undef;
}

=head2 source

 Title   : source
 Usage   : my $acc = $geneOrder->source();
 Function: Returns the accession number from the source sequence, or nothing if undefined
 Returns : A string representing the accession number of the source sequence

=cut

sub source {
	
	return shift->{'source'};

}

=head2 is_circular

 Title   : is_circular
 Usage   : my $circle = $geneOrder->is_circular();
 Function: Returns true if the source sequence is circular
           returns false if gene order consists of more than one permutation
 Returns : A boolean value

=cut

sub is_circular {
	my ($self, $value) = @_;

	if( defined $value){
		$self->{'circular'} = $value;
	}
	return $self->{'circular'};
}

=head2 filtered

 Title   : filtered
 Usage   : $obj->filtered()
 Function: Get/set the filtered state
 Returns : Integer


=cut

sub filtered {
	my ($self,$state) = @_;

	if( defined $state){
		$self->{'filtered'} = $state;
	}
	return $self->{'filtered'};
}

=head2 pi

 Title   : pi
 Usage   : my $order = $geneOrder->pi();
 Function: Returns the permutations representing the gene order
 Returns : An array of Bio::GeneOrder::permutation objects

=cut

sub pi {
	return @{shift->{'pi'}};
}

=head2 order

 Title   : order
 Usage   : my $order = $geneOrder->order();
 Function: Returns gene order as an array of strings
 Returns : A scalar value

=cut

sub order {
	my $self = shift;
	
	my @order;
	
	foreach my $pi (@{$self->{'pi'}}){
		push @order, join ' ', map $SWITCH{abs($_)/$_}.$self->{'key'}->{'name'}->{abs($_)}, grep( ! defined $self->{'key'}->{'filt'}->{$self->{'key'}->{'name'}->{abs($_)}} || $self->{'key'}->{'filt'}->{$self->{'key'}->{'name'}->{abs($_)}} == 0, $pi->pi );
	}
	
	my $order = join "\n", @order;
	
	return $order;
}

=head2 genes

 Title   : genes
 Usage   : my $gene = $geneOrder->genes();	   
 Function: Returns an array of the names of the genes in this gene order
 Returns : An array of strings
 Args    : -all             => Returns all genes, filtered or unfiltered.
           -filtered        => Returns only filtered genes.

=cut

sub genes {
	my ($self,@args) = @_;
	
	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	
	$self->throw("name argument provided, but with an undefined value") 
		if( !defined $param{'-name'} && exists $param{'-name'});
		
	my @genes;
	foreach my $pi (@{$self->{'pi'}}){
		push @genes, map( $self->{'key'}->{'name'}->{abs($_)}, $pi->pi('all') );
	}
	if(defined $param{'-name'}){
		return sort grep($_ eq $param{'-name'}, @genes);
	}
	
	if(defined $param{'-all'} ){
		return sort @genes;
	}elsif(defined $param{'-filtered'}){
		return sort grep($self->{'key'}->{'filt'}->{$_}, @genes);
	}else{
		return sort grep($self->{'key'}->{'filt'}->{$_} == 0, @genes);
	}
}

=head2 bounds

 Title   : bounds
 Usage   : my @bounds = $geneOrder->bounds();
 Function: Returns the gene boundaries in the geneOrder
 Returns : An array of array references, each containing two genes

=cut

sub bounds {
	my $self = shift;
	my @bounds;
	
	foreach my $pi (@{$self->{'pi'}}){
		my @genes;
	
		if(%FILTER){
			@genes = map( $SWITCH{abs($_)/$_}.$self->{'key'}->{'name'}->{abs($_)}, grep( $self->{'key'}->{'filt'}->{$self->{'key'}->{'name'}->{abs($_)}} == 0, $pi->pi ) );
		}else{
			@genes = map( $SWITCH{abs($_)/$_}.$self->{'key'}->{'name'}->{abs($_)}, $pi->pi );
		}
		
		push @genes, $genes[0] if $pi->is_circular;
		
		for(my $i=0;$i<scalar @genes-1;$i++){
			my @bound = ($genes[$i],$genes[$i+1]);
			push @bounds, \@bound
		}
	}
	
	return @bounds;
}

=head2 no_genes

 Title   : no_genes
 Usage   : my $no_genes = $geneOrder->no_genes();
 Function: Returns total number of genes in geneOrder
 Returns : A scalar value

=cut

sub no_genes {
	my $self = shift;
	
	my $no_genes;
	
	foreach my $pi (@{$self->{'pi'}}){
		$no_genes += $pi->pi;
	}
		
	return $no_genes;
}

=head2 no_bounds

 Title   : no_bounds
 Usage   : my $no_bounds = $geneOrder->no_bounds();
 Function: Returns total number of bounds in geneOrder
 Returns : A scalar value

=cut

sub no_bounds {
	my $self = shift;
	
	my $no_bounds;
	
	foreach my $pi (@{$self->{'pi'}}){
		$no_bounds += $pi->pi;
		$no_bounds-- unless $pi->is_circular;
	}
		
	return $no_bounds;
}

=head2 adjacencies

 Title   : adjacencies
 Usage   : my $adjacencies = $geneOrderA->adjacencies( $geneOrderB);
 Function: Find the boundaries that two gene orders share and return them in an array.
 Returns : An array of shared boundaries.  Each shared boundary is an array reference
           with two elements, each of which is a reference to a gene object contained in that boundary.

=cut

sub shared_bounds {
	my ($orderA,$orderB) = @_;

	return $orderA->distance->adjacencies($orderA,$orderB);

}

=head2 breakpoints

 Title   : breakpoints
 Usage   : my $breakpoints = $geneOrderA->breakpoints( $geneOrderB);
 Function: Returns the number of breakpoints between two GeneOrder objects.
 Returns : Scalar value

=cut

sub breakpoints {
	my ($orderA,$orderB) = @_;

	return $orderA->distance->breakpoints($orderA,$orderB);
}

=head2 filter_genes

 Title   : filter_genes
 Usage   : my @filtered = $geneOrder->filter_genes( -type => 'gene type',
                                                    -name => 'regexp'    );
 Function: Purges genes of specified name/type from gene order and returns them
 Returns : An array of Bio::GeneOrder::gene objects. If no arguments are specified
           or the search terms match no genes, no genes are purged and this array is empty.
 Args	 : -type      =>  One of 'gene', 'tRNA', 'mRNA', or 'rRNA' etc.
           -name      =>  A string or regular expression that identifies gene names to purge
                          ie. To match an exact string use a bareword: 'name'.
                          To match an inexact string, supply a regular expression: /'name'/i
           -invert    =>  Filters those genes that do not match the given criteria.
           -unfilter  =>  Unfilters those genes that match the given criteria.
                        

=cut

sub filter_genes {
	my ($self,@args) = @_;

	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	
	$self->throw("name argument provided, but with an undefined value") 
		if( !defined $param{'-name'} && exists $param{'-name'});
	
	$self->throw("type argument provided, but with an undefined value") 
		if( !defined $param{'-type'} && exists $param{'-type'});
	
	@FILTER{ keys %param } = values %param;
	
	my @matched;
	my @genes = keys %{$self->{'key'}->{'index'}};
	
	
	
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
				@matched = grep( $self->{'key'}->{'type'}->{$_} eq $param{'-type'}, @genes);
			}
		}
		
		if(@matched){
			if($param{'-unfilter'}){
				map( $self->{'key'}->{'filt'}->{$_} = 0, @matched);
			}else{
				map( $self->{'key'}->{'filt'}->{$_} = 1, @matched);
			}
			
			$self->{'distance'}->_cache('clear');
		}
		
	}elsif( keys %{ $self->{'key'}->{'filt'} } ){
		map( $self->{'key'}->{'filt'}->{$_} = 0, @genes);
		$self->{'distance'}->_cache('clear');
		$self->_key($self->{'key'});
	}
	
	return @matched;
}

=head2 rename_genes

 Title   : rename_genes
 Usage   : my @genes = $geneOrderA->rename_genes( -list =>  @list_of_names,
                                                  -list =>  @list_of_names2,
                                                  -list =>  ...      );
         OR
           my @genes = $geneOrderA->rename_genes( -table => $filename   );

 Function: Renames genes from names matched in a list to the first name in the list.
		   Any number of lists may be supplied to this function.
		   Alternatively, you may rename genes using a synonym table read from a file.
		   Returns an array of the renamed genes.  
 Returns : An array of Bio::GeneOrder::gene objects.
 Args	 : -list    =>  A list of names to match to genes that will be renamed.
           -table   =>  The name of a file containing a synonymous gene name table.
                        The first line of this file should read ">Bio::GeneOrder::synonyms".
                        Subsequent lines should contain semicolon-delimited lines of synonymous 
                        gene names.  The first name in each line will be used as the primary name
                        for that gene.  In place of a literal name, a synonymous gene name may
                        also be a regular expression.  rename_genes will recognize a synonymous
                        gene name quoted with forward slashes '/' as a regular expression.

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
			$new_key->{'type'}->{$map->{$_}} = $self->{'key'}->{'type'}->{$_};
			$new_key->{'filt'}->{$map->{$_}} = 0;
		}
		
		$self->_key($new_key,$map);
		$self->filter_genes(%FILTER) if %FILTER;
		$self->{'distance'}->_cache('clear');
	}
	
	return @renamed;
}

=head2 reorder

 Title   : reorder
 Usage   : $geneOrder->reorder('name');
 Function: Reorders genes so that the specified gene is first and in the same orientation.
 Returns : 1 if reordered, 0 if not

=cut

sub reorder {
	my ($self,$name) = @_;
	
	$self->throw("reorder requires an argument") unless defined $name;
	
	my $r = 0;
	my $key = $self->{'key'}->{'index'}->{$name};
	
	foreach my $pi (@{$self->{'pi'}}){
		$r += $pi->reorder($key) if $pi->is_circular;
	}
	
	return $r;
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
 Usage   : $geneOrderA->save('filename');
 Function: Save the GeneOrder object to a file for later recovery.
 Returns : 1 for success, 0 for failure.

=cut

sub save {
	my ($self,$file) = @_;

	store $self, $file || $self->throw("GeneOrder could not be saved to $file");

	return 1;
}

=head2 _key

 Title   : _key
 Usage   : $orderSet->_key();
 Function: Get/set the gene name/number key.

=cut

sub _key {
	my ($self,$key,$map) = @_;

	
	if(defined $key){
		map($_->_key($key,$map), $self->pi);	
		$self->{'key'} = $key;
	}
	
	return $self->{'key'};
}


1;
