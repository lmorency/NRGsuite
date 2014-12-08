#!/usr/bin/perl

@PLATFORMS = qw/MacOSX64 Linux32 Linux64 Win32 Win64/;
@LINUX_EXCLUDES = qw/.git .gitignore scripts dist .DS_Store NRGsuite_Xcode NRGsuite_VS2012/;
$EXECUTABLES_GIT_PATH = "/Users/lmorency/Documents/PhD/Programs/Executables";
$NRGSUITE_GIT_PATH = "/Users/lmorency/Documents/PhD/Programs/NRGsuite";
$PLATFORM = "MacOSX64";

$ENV{'COPYFILE_DISABLE'}=true;

sub build_excludes_string {
    $excludes = "";
    foreach $e (@LINUX_EXCLUDES){
	$excludes .= " --exclude \"$e\"";
    }
    return $excludes;
}

for($i=0; $i<@ARGV; $i++){
    if($ARGV[$i] eq '-v' || $ARGV[$i] eq '--version'){
	$VERSION = $ARGV[++$i];
    }elsif($ARGV[$i] eq '-p' || $ARGV[$i] eq '--platform'){
	$PLATFORM = $ARGV[++$i];
    }
}

if("@PLATFORMS"!~/\b$PLATFORM\b/){ die "Invalid platform. Possible platforms are @PLATFORMS\n"; }
if(!$VERSION){ die "No version entered. Please enter a version number with the -v or --version flag." }

print "Platform is $PLATFORM\n";

if($PLATFORM eq 'MacOSX64'){
    print "installing executables...\n";
    `perl install.pl -e "$EXECUTABLES_GIT_PATH" -n "$NRGSUITE_GIT_PATH" -p $PLATFORM`;
    print "building pkg...\n";
    `perl pkgbuild.pl -n "$NRGSUITE_GIT_PATH" -v $VERSION`;

}elsif($PLATFORM =~ /Linux/){
    print "installing executables...\n";
    `perl install.pl -e "$EXECUTABLES_GIT_PATH" -n "$NRGSUITE_GIT_PATH" -p $PLATFORM`;

    chdir("$NRGSUITE_GIT_PATH/..");
    print "now in ".`pwd`;
    $archive1 = "NRGsuite.tar.gz";
    $tar1 = "tar zcf $archive1 ".&build_excludes_string." NRGsuite";
    print "tarring archives...\n";
    `$tar1`;
    `mv $archive1 $NRGSUITE_GIT_PATH/dist`;

    chdir("$NRGSUITE_GIT_PATH/dist");
    print "now in ".`pwd`;
    $archive2 = "NRGsuite_$VERSION"."_$PLATFORM.tar";
    $tar2 = "tar cf $archive2 $archive1 LICENSE install_linux.sh";
    `$tar2`;
    `rm $archive1`;
    
}elsif($PLATFORM =~ /Win/){
    print "The installer for this platform can only be installed using the NSIS application on Windows.\n";
}
