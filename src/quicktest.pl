#!/usr/local/bin/perl

if (!system("make")) {
    chdir("../test");
    system("make testcore");
}
