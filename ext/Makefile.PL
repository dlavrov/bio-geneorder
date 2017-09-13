use ExtUtils::MakeMaker;

$CC = 'g++';
$ENV{'MACOSX_DEPLOYMENT_TARGET'} = 10.3;

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'		=> 'Bio::GeneOrder::Distance',
    'VERSION_FROM'	=> 'Distance.pm',
    'CC'		=> $CC,
    'CCFLAGS'		=> "$$CFLAGS -fPIC",
    'INC'		=> '-I./libd',
    'LD'		=> 'env MACOSX_DEPLOYMENT_TARGET=10.3 $(CC)',
    'XSOPT'		=> '-C++',
    'MYEXTLIB'		=> 'libd/libsw$(LIB_EXT)',
    'clean'		=> { 'FILES' => 'libd/*.o libd/*.a' }
);

sub MY::postamble{
'
$(MYEXTLIB): 
	DEFINE=\'$(DEFINE)\'; CC=\'$(PERLMAINCC)\'; CFLAGS=\'$(CCFLAGS)\'; export DEFINE INC CC CFLAGS; \
		cd libd && $(MAKE) CC=\'$(PERLMAINCC)\' CFLAGS=\'$(CCFLAGS) $(DEFINE)\' DEFINE=\'$(DEFINE)\' libsw$(LIB_EXT) -e
			
';
}
