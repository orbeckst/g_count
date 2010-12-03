#!/usr/bin/perl -w
# $Id: Statistics.pm,v 1.1 2002/06/05 19:51:25 oliver Exp $
# $Log: Statistics.pm,v $
# Revision 1.1  2002/06/05 19:51:25  oliver
# Auxiliary perl libraries (which I tend to use quite often). Put them in a directory of your choice and put that directory in turn in the 'use lib DIR' statement
#
# Revision 1.6  2002/06/03 13:03:46  oliver
# - corrected histogram(): now the center of a bin determines the bin's
#   value
# - use  Carp/croak for error messages
# - moved obsolete plot_histogram() to end
#
# Revision 1.5  2002/03/08 16:17:27  oliver
# new functions: weighted averages etc:
# input: \%data = \(x1 => p(x1), x2 => p(x2), ...)
# 	w_average             expectation value E(X) of distrib. %data
# 	w_standarddeviation   sigma = sqrt(Var(X))
# 	w_avg                 (mu, sigma)
#
# Revision 1.4  2002/03/04 14:06:52  oliver
# * split plot_histogram() in more focused routines:
#   - ascii_plot_histogram() : plots the graph on the terminal
#   - write_N_histogram()    : writes histogram to file
#   - write_P_histogram()    : writes probability distribution to file
#     CAVEAT: It actaully writes (x, N(x)/N_tot) but for a proper pdf I
#             should also divide by the binwidth.
#   The old plot_histogram() is retained in EXPORT_OK for backward
#   compatibility.
# * Decreased PLOTHISTMAX to 55
# * Added a 'return 1' (to make 'use Statistics' happy)
#
# Revision 1.3  2002/02/28 17:52:22  oliver
# plot_histogram() accepts an optional hash with the same keys as the
# histogram but with a second set of values which are printed after the
# values as an annotation
#
# Revision 1.2  2002/02/25 20:26:04  oliver
# * added plot_histogram() and PLOTHISTMAX (max number of chars in the
#   ascii hist plot)
# * use POSIX::ceil
#
# Revision 1.1.1.1  2002/02/25 13:01:38  oliver
# Perl packages; currently living in /sansom/gfij/oliver/Gromacs/Scripts/Perl
#
# use it as:
#
# use lib "/sansom/kir/oliver/local/lib/perl";
# use Statistics; 


package      Statistics;
require      Exporter;
@ISA       = qw(Exporter);
@EXPORT    = qw(avg poisson average   standarddeviation poissonerror 
		w_avg     w_average w_standarddeviation
		histogram ascii_plot_histogram 
		write_N_histogram write_P_histogram);
@EXPORT_OK = qw($PLOTHISTMAX plot_histogram);

use Carp;
use Messages;
use POSIX    qw(ceil);

# maximum number of characters in the histogram ascii graph
$PLOTHISTMAX = 55;

# make 'use Statistics' happy
return 1;

sub avg {
    # calculate average and standard deviation of a _ref_ to a list of values
    # and return the array ($avg, $stddev)
    my ($data) = @_;
    return (&average($data), &standarddeviation($data));
};

sub w_avg {
    # calculate weighted average (expectation value) and standard 
    # deviation of a _ref_ to a hash %data = (x1 => p1, x2 => p2, ...)
    # and return the array ($w_avg, $w_stddev)
    my ($data) = @_;
    return (&w_average($data), &w_standarddeviation($data));
};


sub poisson {
    # calculate average and error of a _ref_ to a list of values
    # assuming a Poisson distribution of N observations
    # and return the array ($avg, $error)
    my ($data) = @_;
    return (&average($data), &poissonerror($data));
};


sub average {
    my ($data) = @_;
    my ($sum, $x) = (0, 0);
    foreach $x (@$data) {
	$sum += $x;
    };

    return @$data > 0 ? $sum/@$data : 0;
};

sub poissonerror {
    # for Poisson-distribution f(x)=1/mu exp(-x/mu)
    # with mean mu =>  s=mu/sqrt(N)
    my ($data) = @_;
    return @$data > 0 ?  &average($data)/sqrt(scalar @$data) : 0;
};

sub standarddeviation {
    my ($data) = @_;
    my ($sum, $avg, $x) = (0, 0, 0);

    $avg = &average($data);

    foreach $x (@$data) {
	$sum += ($x - $avg)**2;
    };

    return @$data > 0 ? sqrt($sum/@$data) : 0;
};
	

sub w_average {
    # weighted average 
    # %data = (x1 => p1, x2 => p2)
    # <x> = Sum p_i x_i
    my ($data) = @_;
    my ($sum, $x, $norm) = (0, 0, 0);
    foreach $x (keys %$data) {
	$sum  += $data->{$x} * $x;
	$norm += $data->{$x};
#	printf "w_average: x=%g p(x)=%g sum=%g norm=%g\n", 
#	       $x,$data->{$x},$sum,$norm;
    };
    return $norm > 0 ? $sum/$norm : 0;
};


sub w_standarddeviation {
    my ($data) = @_;
    my ($sum, $avg, $x, $norm) = (0, 0, 0, 0);

    $avg = &w_average($data);

    foreach $x (keys %$data) {
	$sum += $data->{$x} * ($x - $avg)**2;
	$norm += $data->{$x};
    };

    return $norm > 0 ? sqrt($sum/$norm) : 0;
};


sub histogram  {
    # Create a histogram from a list of numbers
    # in: ($bwidth  $nbins  \@array)
    # (set bwidth or nbins)
    # returns the hash (x1 => N1, x2 => N2, ...) in LIST context
    # and the binwidth in SCALAR
    # (a hash is easier to sort in order to determine 
    # N_min and N-max for example)
    my ($bwidth, $nbins, $array) = @_;
    my ($max, $min);
    my (@sorted, $i, $j, $x, $x_centered);
    my (%hist, $time);
    
    @sorted = sort { $a <=> $b } @$array;
    
    ($max, $min) = ($sorted[-1], $sorted[0]);

    # For non neg. numbers: set explicitly min=0
    $min = 0 if ($min > 0);

    if ($bwidth > 0) {
	$nbins =  ceil( ($max - $min)/(1.0*$bwidth) );
    } 
    elsif ($nbins) {
	$bwidth = ($max - $min)/(1.0*$nbins);
    } 
    else {
	croak "histogram(): Provide either BINWIDTH or NBINS!"
    };
    
    return $bwidth unless wantarray;

    logger (3, "histogram: N=%g  x_max=%g x_min=%g (set to x_min=%g) ".
	    "N_bins=%g  binwidth=%g\n",
	 scalar @sorted, $max, $sorted[0], $min, $nbins, $bwidth);
    
    # define bins so that $x is the maximum value contained in the bin
    # except for the min value:
    #      1             2        N_bin-1     N_bin
    # |min  min+bw|   min+2*bw| ..... max-bw|    max|   
    #
    # NB: using the upper bound of the bin to be 'the' x-value of the
    # bin biases the calculation of averages towards high x values
    # with increasing binwidth. It would be a better idea to use the
    # center of the bin.

    for ($j=1, $i=0; $j<=$nbins; $j++)
    {
	$x = $min + $j * $bwidth;
	$x_centered = $x - $bwidth/2.0;

	$hist{$x_centered} = 0;
	msg (11, "i=$i x=$x x_centered=$x_centered sorted[i]="
	     .$sorted[$i]."\n") 
	    unless ($i >= @sorted);  
	# comparison is not exact, therefore epsilon=0.0001
	while($i < @sorted &&  ($sorted[$i] - $x) < 0.0001  )
	{
	    $hist{$x_centered}++;
	    $i++;
	};
    }

    # we get here only if we did not leave before in SCALAR context
    return %hist;
};



##########################################
# hacking plot_histogram() into pieces
# (perhaps make a class out of histogram ?)
##########################################

sub write_N_histogram {
    # plot_histogram FILEHANDLE \%histogram (\%annotation)
    # annotation is a hash with the same keys as histogram; its
    # entries are printed after the keys (user responsibility to keep them
    # short enough)
    my $HISTFILE = shift;
    my %hist  = %{ shift @_};
    my %annot = %{ shift @_} if (@_);
    my ($x);
    
    foreach $x (sort {$a <=> $b} (keys %hist)) {
	my $buf;	
	if (%annot) {
	    $buf = sprintf "%g %g   # %8.3f\n", 
	    $x, $hist{$x}, $annot{$x};
	} else {
	    $buf = sprintf "%g %g\n", 
	    $x, $hist{$x};
	};
	printf $HISTFILE $buf;
    };
};

sub write_P_histogram {
    # plot_histogram FILEHANDLE \%histogram (\%annotation)
    # print the probabilities (N/N_tot)
    # annotation is a hash with the same keys as histogram; its
    # entries are printed after the keys (user responsibility to keep them
    # short enough)
    my $PDFILE = shift;
    my %hist  = %{ shift @_};
    my %annot = %{ shift @_} if (@_);
    my ($Ntotal, $N, $x);
    
    # recalculate total number of observations
    $Ntotal = 0;
    foreach $N (values %hist) {
	$Ntotal += $N;
    };

    foreach $x (sort {$a <=> $b} (keys %hist)) {
	my $buf;	
	if (%annot) {
	    $buf = sprintf "%g %g   # %8.3f\n", 
	    $x, $hist{$x}/$Ntotal, $annot{$x};
	} else {
	    $buf = sprintf "%g %g\n", 
	    $x, $hist{$x}/$Ntotal;
	};
	printf $PDFILE $buf;
    };
};


sub ascii_plot_histogram {
    # ascii_plot_histogram \%histogram (\%annotation)
    # annotation is a hash with the same keys as histogram; its
    # entries are printed after the keys (user responsibility to keep them
    # short enough)
    my %hist  = %{ shift @_};
    my %annot = %{ shift @_} if (@_);
    my ($Ntotal, $N, @histsorted, $histmax, $histmin, $norm, $x, $k);
    
    @histsorted = sort {$a <=> $b} (values %hist);
    ($histmax, $histmin) = ($histsorted[-1], $histsorted[0]);

    # recalculate total number of observations
    $Ntotal = 0;
    foreach $N (values %hist) {
	$Ntotal += $N;
    };
    
    $norm = $PLOTHISTMAX/$histmax;
    logger (3, "histogram: count_max=$histmax count_min=$histmin\n");
    foreach $x (sort {$a <=> $b} (keys %hist)) {
	my $buf;
	if (%annot) {
	    $buf = sprintf "%5g %8g %5.2f %s", 
	           $x, $hist{$x}, $annot{$x}, $hist{$x}>0 ? "|" : ".";
	} else {
	    $buf = sprintf "%5g %8g %s", 
	           $x, $hist{$x}, $hist{$x}>0 ? "|" : ".";
	};

	logger (1, $buf);
	for ($k=1; $k <= $hist{$x} * $norm; $k++) {
	    logger (1, "#");
	};
	logger (1, "\n");
    };
    logger (1, "---- N_tot = $Ntotal ----\n\n");
};



####################################################################
# OBSOLETE -- retained for backward compatibility. Use the other 
# histogram functions! 
####################################################################

sub plot_histogram {
    # plot_histogram FILEHANDLE \%histogram (\%annotation)
    # annotation is a hash with the same keys as histogram; its
    # entries are printed after the keys (user responsibility to keep them
    # short enough)
    my $HISTFILE = shift;
    my %hist  = %{ shift @_};
    my %annot = %{ shift @_} if (@_);
    my ($Ntotal, $N, @histsorted, $histmax, $histmin, $norm, $x, $k);
    
    @histsorted = sort {$a <=> $b} (values %hist);
    ($histmax, $histmin) = ($histsorted[-1], $histsorted[0]);

    # recalculate total number of observations
    $Ntotal = 0;
    foreach $N (values %hist) {
	$Ntotal += $N;
    };

    

    $norm = $PLOTHISTMAX/$histmax;
    logger (3, "histogram: count_max=$histmax count_min=$histmin\n");
#    logger (3, "histogram: count_max=$histmax count_min=$histmin " .
#	    "bw=$bwidth\n");
    foreach $x (sort {$a <=> $b} (keys %hist)) {
	my $buf;
	printf $HISTFILE "%8g %8g\n", $x, $hist{$x};
#	printf PD   "%8g %6g\n",  $x, $hist{$x}/($bwidth * $Ntotal);
	
	if (%annot) {
	    $buf = sprintf "%5g %8g %5.2f %s", 
	           $x, $hist{$x}, $annot{$x}, $hist{$x}>0 ? "|" : ".";
	} else {
	    $buf = sprintf "%5g %8g %s", 
	           $x, $hist{$x}, $hist{$x}>0 ? "|" : ".";
	};

	logger (1, $buf);
	for ($k=1; $k <= $hist{$x} * $norm; $k++) {
	    logger (1, "#");
	};
	logger (1, "\n");
    };
    logger (1, "---- N_tot = $Ntotal ----\n\n");
};
