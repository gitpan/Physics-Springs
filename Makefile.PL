use 5.006;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'		=> 'Physics::Springs',
    'VERSION_FROM'	=> 'lib/Physics/Springs.pm', # finds $VERSION
    'PREREQ_PM'		=> {
                             Physics::Particles => 1.00,
                             Test::More => 0,
                             Data::Dumper => 0,
			     Sub::Assert => 0,
                           }, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?    ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM => 'lib/Physics/Springs.pm', # retrieve abstract from module
       AUTHOR     => 'Steffen Mueller <springs-module at steffen-mueller dot net>') : ()),
);
