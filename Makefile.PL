use ExtUtils::MakeMaker;
WriteMakefile(
    'NAME'		=> 'Bio::GeneOrder',
    'VERSION_FROM'	=> 'lib/Bio/GeneOrder.pm',
    'DISTNAME'		=> 'Bio-GeneOrder',
    'PREREQ_PM'		=> {'Test::More'	=> 0,
			    'Bio::Perl'		=> 1.006,
			   },
    'EXE_FILES'		=> ['bin/gogo'],
);
