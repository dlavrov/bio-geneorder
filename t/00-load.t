#!perl -T

use Test::More tests => 1;

BEGIN {
	use_ok( 'Bio::GeneOrder' );
}

diag( "Testing Bio::GeneOrder $Bio::GeneOrder::VERSION, Perl $], $^X" );
