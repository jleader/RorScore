#!/usr/bin/perl
use Module::Build;
my $build = Module::Build->new(
    dist_name=>'Rorshach::Score',
    module_name=>'Rorshach::Score',
    dist_version=>'0.1',
    dist_author=>'Jeremy Leader <jleader@alumni.caltech.edu>',
    license=>'perl',
    requires=>{
        perl=>'5.8.0',
    },
    build_requires=>{
        'Test::More'=>0,
    },
);
