#!/usr/bin/perl

$infile = $ARGV[0];

open(IN, $infile) || die("Failed to open $infile");
$row = <IN>;
@parts = split(/\s+/, $row);
$header = shift(@parts);
for $i (0 .. $#parts) {
    print "$parts[$i]\n";
}
