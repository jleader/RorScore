#!/usr/bin/perl

use strict;
use warnings;

use Test::More;
use Path::Class qw(dir);;

require 'bin/rorschach.pl';

my $tdir = dir('t', 'data');
my $count = 0;

while (my $file = $tdir->next) {
    if ($file =~ m( / (\d\d) input \. txt $)x) {
        my $num = $1;
        my $outname = $num . 'output-expected.txt';
        my $expected = $tdir->file($outname)->slurp;
        my $result = main($file);
        is($result, $expected, "$num output should match $outname");
        ++$count;
    }
}

done_testing($count);
