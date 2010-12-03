#!/usr/bin/perl -w
# $Id: slicer.pl,v 1.2 2003/08/08 18:39:02 oliver Exp $
# first: source the gromacs environment, 
# eg source ~/bin/GMXRC_new.bash


$ANA_D = ".";
$MAKEFILE = $ANA_D."/Makefile";  # this is Makefile.grid

$Z0 = 0;
$BZ = 8;      # box 0 .. BZ (nm)
$DZ = 0.6;    # slices (nm)

for ($Z1 = $Z0; $Z1 <= $BZ-$DZ; $Z1+=$DZ) {   
    $Z2 = $Z1 + $DZ;
    $XYP_NAME = sprintf("%.2f",$Z1) . "_xyp";
    $XYP_NAME =~ tr/./_/;
    $XYP = "Grid.d/".$XYP_NAME.".eps";
    print "slicer: [$Z1 -> $Z2] $XYP\n";
    system "make -f $MAKEFILE  Z1=$Z1 Z2=$Z2 Grid.d/xyp.eps";
    rename "Grid.d/xyp.eps", $XYP;
}

print "Show eps files:\n";
system "ls Grid.d/*_xyp.eps | xargs -i gv {}";
