#!/usr/bin/perl

$NRGSUITE_GIT_PATH = "/Users/francisgaudreault/Development/NRGsuite/NRGsuite";
$EXECUTABLES_GIT_PATH = "/Users/francisgaudreault/Development/Executables";

for($i=0; $i<@ARGV; $i++){
    if($ARGV[$i] eq '-n'){
	$NRGSUITE_GIT_PATH = $ARGV[++$i];
    }elsif($ARGV[$i] eq '-e'){
	$EXECUTABLES_GIT_PATH = $ARGV[++$i];
    }elsif($ARGV[$i] eq '-p'){
	$PLATFORM = $ARGV[++$i];
    }
}

sub Clean {
    `rm -rf $NRGSUITE_GIT_PATH/FlexAID/WRK`;
    `rm -rf $NRGSUITE_GIT_PATH/GetCleft/WRK`;
    `rm -rf $NRGSUITE_GIT_PATH/*.pyc`;
    `rm -rf $NRGSUITE_GIT_PATH/*.pyo`;
    `rm -rf $NRGSUITE_GIT_PATH/FlexAID/*.pyc`;
    `rm -rf $NRGSUITE_GIT_PATH/FlexAID/*.pyo`;
    `rm -rf $NRGSUITE_GIT_PATH/GetCleft/*.pyc`;
    `rm -rf $NRGSUITE_GIT_PATH/GetCleft/*.pyo`;
}

sub Build {
    `cp -r $EXECUTABLES_GIT_PATH/$PLATFORM/* $NRGSUITE_GIT_PATH/`;
}

if(!$PLATFORM){ die "no platform entered. exit.\n"; }

if($PLATFORM eq 'MacOSX64'){
	&Clean;
	&Build;

}elsif($PLATFORM eq 'Linux32'){
	&Clean;
	&Build;
	
}elsif($PLATFORM eq 'Linux64'){
	&Clean;
	&Build;
	
}else{
	die "unknown platform: $PLATFORM. exit.\n";
}

print "Done!\n";
