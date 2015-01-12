$NRGSUITE_GIT_PATH = "/Users/lmorency/Documents/PhD/Programs/NRGsuite";
$INSTALLATION_DIR = "/Applications/NRGsuite";
$VERSION = '2.47d';

for($i=0; $i<@ARGV; $i++){
    if($ARGV[$i] eq '-v'){
	$VERSION = $ARGV[++$i];
    }elsif($ARGV[$i] eq '-n'){
	$NRGSUITE_GIT_PATH = $ARGV[++$i];
    }
}

$cmd = "pkgbuild --identifier \"Najmanovich Research Group\" --install-location $INSTALLATION_DIR --root $NRGSUITE_GIT_PATH --filter .DS_Store --filter .git --filter .pyc --filter .gitignore --filter dist --filter NRGsuite_VS2010 --filter NRGsuite_Xcode --scripts $NRGSUITE_GIT_PATH/scripts NRGsuite_$VERSION"."_MacOSX64.pkg";

print "$cmd\n";
`$cmd`;
