#!/usr/bin/perl -w
# $Id: a_flux.pl,v 1.2 2002/06/05 19:52:12 oliver Exp $

use strict;
use lib ".";            # put Statistics.pm and Messages.pm into that dir 
use Messages; 
use Statistics;

use FileHandle;
use Getopt::Long;

my $Version_ID = q($Id: a_flux.pl,v 1.2 2002/06/05 19:52:12 oliver Exp $ );

# imported from Messages
# DEBUG level says how chatty we are:
#  <0  dead silent
#   0  no diagnostic messages, only output
#   1  a few additional important informations
#   3  verbose
#   5  extra verbose
#  10  start with real debugging stuff
#  15  many intermediate results
#  16  print input files line-by-line
#  20  everything   
$DEBUG = 1;

my $DataFn = "flux.dat";
my $OutputFn = "flux_wsm.dat";
my $TauDistFn = "tau_cross.dat";
my $WINSIZE  = 100;
my $BINWIDTH = 4;
my ($opt_verbose, $opt_version, $opt_help);

# how to indicate a comment line in the data file (regex)
my $commentchar = '[;#@]';



sub print_usage {
    my $Prog = Progname;
    print <<EOF;
usage: $Prog [OPTIONS] [FILE]
Analyse and smooth data form FILE [$DataFn], produced by g_flux:
INPUT:    t PHI(t) PHI+(t) PHI-(t) |PHI(t)| tau
OUTPUT:   t ... window average ... (without tau)
  distribution:
          tau  N   N/N_tot  N/(N_tot*Bin_width)
           
1. partition time in windows of size SIZE and average all values
   (data smoothing: increases resolution in y and decreases it in x)
2. histogram/distribution of residency/crossing times tau    



The OPTIONS can be abbreviated; defaults are given in [].        
--help              help
--verbose           (alias for DEBUGLEVEL=3)
--version           print version id
--debug=DEBUGLEVEL  0 (almost) quiet, <0 silent,
                    1..3 verbose, 
                    >5 really debugging stuff (20 max) [$DEBUG]

--window=SIZE       window size in ps [$WINSIZE]
--binwidth=SIZE     bin width for histogram in ps [$BINWIDTH]

--output=FILE       [$OutputFn]
--dist=FILE         distribution [$TauDistFn]
EOF
    exit;
}

use Carp;
sub vec_add ($$) {
    # in  \$a \$b
    # out \(a+b)
    my ($a, $b) = @_;
    my (@this, $i);
    croak "vec_add() failed: unlike vector length" unless @$a == @$b;
    for ($i=0; $i < @$a; $i++) { push @this, @$a->[$i] + @$b->[$i] };
    return \@this;
};

sub vec_sub ($$) {
    # in  \$a \$b
    # out \(a-b)
    my ($a, $b) = @_;
    my (@this, $i);
    croak "vec_sub() failed: unlike vector length" unless @$a == @$b;
    for ($i=0; $i < @$a; $i++) { push @this, @$a->[$i] - @$b->[$i] };
    return \@this;
};

sub vec_sprod ($$) {
    # in  $r \$a
    # out \(r*a)
    my ($r, $a) = @_;
    my (@this, $i);
    for ($i=0; $i < @$a; $i++) { push @this, $r*@$a->[$i] };
    return \@this;
};

sub xvec_sprod ($$$$) {
    # indexed vec_sprod
    # in  $r \$a
    # out \(r*a) (but only act on a[m..n])
    my ($r, $a, $m, $n) = @_;
    my (@this, $i);
    for ($i=0; $i < @$a; $i++) { 
        push @this, ($i >= $m and $i <= $n) ? $r*@$a->[$i] : @$a->[$i];
    };
    return \@this;
};
    



######################################
#
# MAIN
#
######################################


# options
&GetOptions( "output=s" => \$OutputFn,
             "dist=s"   => \$TauDistFn,
             "window=f" => \$WINSIZE,
             "binwidth=f" => \$BINWIDTH,
             "debug=i" => \$DEBUG,
             "verbose" => \$opt_verbose,
             "version" => \$opt_version,
             "help" => \$opt_help) or &print_usage;
if ($opt_help) {&print_usage};
if ($opt_version) { msg (0, "%s\n", $Version_ID); exit 0 };
if ($opt_verbose) { $DEBUG = 3 };


open_log or die "Cannot open logfile for writing, stopped ";

if ($ARGV[0]) {
    $DataFn = $ARGV[0];
    -e $DataFn or 
        die "The  file $DataFn to read from does not exist, stopped ";
};

open (IN,  "< $DataFn") or die "Error opening input '$DataFn', stopped";
open (OUT, "> $OutputFn") or die "Error opening output '$OutputFn', stopped";

printf OUT "# window-averaged flux data, window = %g ps\n", $WINSIZE;

my @data;

# <>: read from all files on the commandline or STDIN:
SLURP: while (<IN>) {
    chomp;
    next SLURP if /^ *$commentchar/ || /^ *$/;
    msg (20, "[$.] [$_]\n");
    push @data, [ split ]; 
};

my ($t, $tau, @tau);
my %phi = (net  => 0,
           up   => 0,
           dwn  => 0,
           tot  => 0);
my $i;
my $window = [0, 0, 0, 0, 0, 0];   # pointer to @window
my $N_dat = 0;
my $t_start = $data[0]->[0];
my $t_end   = $t_start + $WINSIZE;

LINE: foreach (@data) {
    ($t, $phi{net}, $phi{up}, $phi{dwn}, $phi{tot}, $tau) =  @$_;
    
    if ($phi{tot} > 0) {
        $N_dat++;
        $window =  vec_add ($window,  $_);

        msg (16, "t=%8g phi{net}=%g phi{up}=%g phi{dwn}=%g phi{tot}=%g " .
             "tau=%g\n",  
             ($t, $phi{net}, $phi{up}, $phi{dwn}, $phi{tot}, $tau));
        msg (16, "%-4s%6g phi{net}=%g phi{up}=%g phi{dwn}=%g phi{tot}=%g " .
             "tau=%g N_dat=%g\n", "win", @$window, $N_dat);

        # collect tau with the correct weight
        # (g_flux already averages over all crossings that 
        #  finish at the same t)
        for ($i=0; $i<$phi{tot}; $i++) {push @tau, $tau; };
    };

    if ($t >= $t_end or $_ == $data[-1]) {      
        $window->[0] = ($N_dat) ? $window->[0]/$N_dat : ($t_start + $t)/2;
        $window->[5] = ($N_dat) ? $window->[5]/$N_dat : 0;
        $window = xvec_sprod (1/($t - $t_start), $window, 1, 4);

        printf OUT "%6g %g %g %g %g    %g\n", @$window;

        msg (15, "window [%g -> %g]\n", $t_start, $t);
        msg (15, "t=%6g phi{net}=%g phi{up}=%g phi{dwn}=%g " . 
             "phi{tot}=%g tau=%g\n\n", @$window);

        $t_end   += $WINSIZE;
        $t_start  = $t;
        $N_dat  = 0;
        $window = [0, 0, 0, 0, 0, 0];
    };
};

{
    my $HIST = new FileHandle;
    my (%hist, @avg_tau);

    open ($HIST,"> $TauDistFn") 
        or die "Error opening output '$TauDistFn', stopped";

    %hist = histogram ($BINWIDTH, 0, \@tau);
    @avg_tau = w_avg (\%hist);

    logger (0, "<tau> = %g±%g\n", @avg_tau);
    printf $HIST "# <tau> = %g±%g\n", @avg_tau;

    write_N_histogram ($HIST, \%hist);
    ascii_plot_histogram (\%hist);
};


exit(0);	
