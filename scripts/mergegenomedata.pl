#!/usr/bin/perl
# -*-Perl-*-
# Last changed Time-stamp: <2009-01-12 11:42:29 raim>
# $Id: mergegenomedata.pl,v 1.3 2009/01/13 01:25:28 raim Exp $

use strict;
use File::Basename;
use Getopt::Long;

$|=1; # flush print to stdout

my $write   = 0;
my $outfile = "genome.csv";
my $inDir   = "";
my $infiles = "";
my $matchcol= "0,1";
my $replace = "";
my $progress = 1; # report progress based on first column (ie. chromosome #)
# get command line parameters
usage() unless GetOptions("w"   => \$write,
			  "o=s" => \$outfile,
			  "i=s" => \$inDir,
			  "f=s" => \$infiles,
			  "m=s" => \$matchcol,
			  "p"   => \$progress,
			  "r=s" => \$replace,
			  "h"   => \&usage);
#usage() unless -e $inDir;

sub usage
{
    printf STDERR
	"\n\n Reads several tab-delimited chromosomal coordinate files";
    printf STDERR "\n from an input directory (PATH/*.csv)";
    printf STDERR "\n and either just checks whether coordinates ";
    printf STDERR "\n in each line of the files correspond to each other,";
    printf STDERR "\n or optionally also merges all value columns into one";
    printf STDERR "\n big coordinate value file\n\n";
    
    printf STDERR "\tUsage: @{[basename($0)]}  [options] \n\n";
    
    printf STDERR "-i=s\tdirectory with input *.csv files, current: $inDir\n";
    printf STDERR "\tOR\n";
    printf STDERR "-f=s\tlist of input files, ','-separated, alternative to -i\n, current: $infiles\n";
    printf STDERR "-m=s\tlist of columns (column numbers) which must match between input files (should be equal) and will not be duplicated, current: $matchcol\n";
    printf STDERR "-r=s\tlist of columns (IDs or column numbers) that will be skipped from the first input file only, current: $replace\n";
    printf STDERR "-p\tprint progress message, current: $progress\n";
    printf STDERR "-w\twrite output file, default name 'genome.csv', current: $write\n";
    printf STDERR "-o=s\talternative name for output file, current: $outfile\n\n";
    exit;
}

# flush message print to stdout
$| = 1;

# all files from calculateChromatinExperiments.pl and
# calculateChromatinStructures.pl
my @filelist;
if ( $infiles eq "" )
{
    my $dir = $inDir."/*.csv";
    @filelist = glob($dir);
}
else
{
    my @list = split /,/, $infiles; #/
    foreach ( @list )
    {
	$filelist[@filelist] = $_;
	print "added file $_\n";
    }
}

# open all csv files
my @io;
my $filecnt = 0;
foreach ( @filelist )
{
    open($io[$filecnt], "<$filelist[$filecnt]") or
	die("can't open csv file $filelist[$filecnt]\n");
    print "opened $filelist[$filecnt]\n";
    $filecnt++;
}

# which columns should match?
my @list = split /,/, $matchcol; #/
my @match;    
print "matching column(s)\t";
foreach ( @list )
{
    $match[@match] = $_;
    print "$_ ";
}
print "\n";

# which columns from first file to leave out?
my @list = split /,/, $replace; #/
my @repl;    
print "replacing column(s)\t";
foreach ( @list )
{
    $repl[@repl] = $_;
    print "$_ ";
}
print "\n";

# write all data to one big genome file or ...
if ( $write )
{
    open(outf, ">$outfile") or die "can't open outfile $outfile\n";
    print "writing $outfile ...\n";
}
# ... just check whether coordiantes are in sync between files
else
{
    print "checking correspondence of match columns ...\n";
}


# use first file for coordinate reference
my $fh = $io[0];
my $chrom = 0;
my @colnum = ();
while( <$fh> )
{
    my @values = ();
    @{$values[0]} = split /\t/;
    chomp(@{$values[0]});

    $colnum[0] = @{$values[0]};

    # progress message
    if ( $chrom != $values[0][0] && $progress )
    {
    	$chrom = $values[0][0];
    	print " ... finished\n" if $chrom > 1;
    	print "\tchromosome $chrom ... ";
    }

    # PRINT COORDINATE AND VALUE columns from first file
    if ( $write )
    {
	my $cnt = 0;
	for ( my $j=0; $j<@{$values[0]}; $j++ )
	{
	    my $printcol = 1;
	    for ( my $k=0; $k<@repl; $k++ )
	    {
		# replace is column number
		if ( $repl[$k] =~ /^\d+$/ )
		{
		    $printcol = 0 if $j == $repl[$k];
		}
		# replace is column ID, only for header line
		else
		{
		    if ( $values[0][$j] eq $repl[$k] )
		    {
			$printcol = 0;
			print "found replace value $repl[$k] at column $j\n";
			$repl[$k] = $j;
		    }
		}
	    }
	    next unless $printcol;
	    print outf "\t" if $cnt > 0;
	    print outf $values[0][$j];
	    $cnt++;
	}
    }
    
    # go through other files
    for ( my $i=1; $i<$filecnt; $i++ )
    {
	# move one line forward in file #i
	my $fhx = $io[$i];
	my $line = <$fhx>;
	
	@{$values[$i]} = split /\t/, $line;
	chomp(@{$values[$i]});
	
	# ERROR message
	# check whether match columns are equal!
	my $matches = 0;
	foreach ( @match )
	{
	    $matches++ if $values[$i][$_] eq $values[0][$_]	    
	}
	unless ( $matches == @match )
	{
	    print "\n\nERROR: $filelist[$i] coordinates differ from $filelist[0]\n";
	    foreach ( @match )
	    {
		print "VALUES: $values[0][$_]\t$values[$i][$_]\n";
	    }
	    exit;
	}
	# print value columns from other files
	if ( $write )
	{
	    for (my $j=0; $j < @{$values[$i]}; $j++ )
	    {
		my $printcol = 1;
		foreach ( @match )
		{
		    $printcol = 0 if $_ == $j;
		}
		print outf "\t".$values[$i][$j] if $printcol;
	    }
	}
    }
    print outf "\n" if $write;
}
print " ... finished\n";

print "FINISHED MERGING DATA INTO FILE $outfile\n" if $write;
# write chr coor and all values


__END__
