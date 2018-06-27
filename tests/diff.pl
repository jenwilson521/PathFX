#! /usr/bin/perl

use lib qw( lib );
use PathFX::Test;

my $test  = new PathFX::Test();
my $path  = "../results/pathfx_mike_demo";
my $truth = $path . '.truth';

foreach my $query (qw( DB00331 Montelukast )) {
	$test->compare_results( "$path/$query", "$truth/$query" );
}
