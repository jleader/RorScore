#!/usr/bin/perl
use Module::Build;
Module::Build->new(
    dist_name           =>'Rorschach::Score',
    dist_abstract       =>
        'Script to create a structural summary from a Rorschach score sheet.',
    module_name         =>'Rorschach::Score',
    dist_version        =>'0.1',
    dist_author         =>'Jeremy Leader <jleader@alumni.caltech.edu>',
    license             =>'perl',
    requires            =>{
        perl                    =>'5.8.0',
        'Path::Class'           =>'0.25',
    },
    build_requires      =>{
        'Test::More'            =>'0.98',
    },
)->create_build_script;
