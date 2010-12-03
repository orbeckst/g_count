#!/usr/bin/perl -w
#
# $Id: Messages.pm,v 1.1 2002/06/05 19:51:25 oliver Exp $
# $Log: Messages.pm,v $
# Revision 1.1  2002/06/05 19:51:25  oliver
# Auxiliary perl libraries (which I tend to use quite often). Put them in a directory of your choice and put that directory in turn in the 'use lib DIR' statement
#
# Revision 1.4  2002/03/20 20:55:04  oliver
# new: $MSG_STDERR; if set  logger et al will write thir messages to
# stderr not stdout
#
# Revision 1.3  2002/02/25 20:30:19  oliver
# * new function (behaves like build ins)
#   Progname() -- return program's name as string; replace similar
#                 construct in the main program
#   date()     -- string with current date and time
# * open_log now chooses the filename and prints a first status line
#
# Revision 1.2  2002/02/25 14:05:03  oliver
# * appears to work
# * EXPORT_OK: if one wants to write from main to the logfile
#   (deprectated) or set an explicit default name for the logfile
# * added open_log: Logfile is now opened within the package (cleaner!!)
#
# global: $DEBUG
# use it as:
#
# use lib "/sansom/kir/oliver/local/lib/perl";
# use Messages; 

package      Messages;
require      Exporter;
@ISA       = qw(Exporter);
@EXPORT    = qw(Progname warning logger log_only msg open_log $DEBUG);
@EXPORT_OK = qw(LOG $MSG_STDERR date xe $SHELL_QUIET);

use POSIX qw(strftime);
use Cwd;
use File::Basename qw(basename);

$DEBUG = 0;
$MESG_STDERR = 0;   # normally print to stdout
$SHELL_QUIET = ">/dev/null 2>&1";

# somehow main is only happy if the package returns true 
# -- and $DEBUG=0 returns false... weird Perl.
return 1;

END {
    close (LOG);
};

sub Progname () {
    return basename $0;
};

sub date () {
    return POSIX::strftime("%c", localtime(time));
};

sub open_log () {
    my ($error, $LogFn);
    ($LogFn = Progname) =~ s(\.pl$)();  
    $LogFn .= ".log";
    if ($error = open (LOG, "> " . $LogFn)) {
	print LOG Progname . " [Logfile] --- ",
	"opened on ", date(), "\n";
	print LOG "in ", cwd(), "\n";
	print LOG 
	    "-----------------------------------------------------------\n\n";
    };
    return $error;
};


sub warning {
    my ($format, @args) = @_;
    $format = "Warning: " . $format;
    logger ( 0, $format, @args);
    return;
};

sub logger {
    my ($level, @args) = @_;
    do {
	if ($MSG_STDERR) {
	    printf STDERR @args;
	} else {
	    printf  @args;
	};
	printf LOG @args;
    } unless $DEBUG < $level;
    return;
};

sub log_only {
    my ($level, @args) = @_;
    do {
	printf LOG @args;
    } unless $DEBUG < $level;
    return;
};


sub msg {
    my ($level, @args) = @_;
    printf @args unless $DEBUG < $level;
    return;
};


# not really logging but...
sub xe ($) {
    my ($CMD) = @_;
    if (LOG) {
	logger (1, ">>> $CMD\n");
    } else {
	print STDERR ">>> $CMD\n";
    };
    system $CMD;
    return $!;
};
