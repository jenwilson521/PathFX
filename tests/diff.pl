#! /usr/bin/perl
use Statistics::Descriptive;

my $target = shift;

my $truth = read_pathfx_results( 'Jenimab' );
my $test  = read_pathfx_results( 'Jenimab.macos.2018-06-13' );
my $stats = new Statistics::Descriptive::Full();

my $neighborhood = $truth->{ neighborhood };
foreach my $file (sort keys %$neighborhood) {
	if( ! exists( $test->{ neighborhood }{ $file } )) {
		warn "ERROR: Missing neighborhood file '$file'\n";
		next;
	}
	foreach my $a (sort keys %{ $neighborhood->{ $file }}) {
		if( ! exists( $test->{ neighborhood }{ $file }{ $a } )) {
			warn "ERROR: Missing neighborhood node '$a'\n";
			next;
		}
		foreach my $b (sort keys %{ $neighborhood->{ $file }{ $a }}) {
			if( ! exists( $test->{ neighborhood }{ $file }{ $a }{ $b } )) {
				warn "ERROR: Missing neighborhood edge '$a' -> '$b'\n";
				next;
			}
			my $i = $truth->{ neighborhood }{ $file }{ $a }{ $b };
			my $j = $test->{ neighborhood }{ $file }{ $a }{ $b };

			# print STDERR "$file\t$a\t$b\t$i\t$j\n";
			my $diff = abs( $i - $j );
			$stats->add_data( $diff );
		}
	}
}

printf "Neighborhood Difference\nmean=%.4f\tsdev=%.4f\n", $stats->mean(), $stats->standard_deviation();

$stats->clear();
my $associations = $truth->{ association };
foreach my $file ( sort keys %$associations ) {
	if( ! exists( $test->{ association }{ $file } )) {
		warn "ERROR: Missing association file '$file'\n";
		next;
	}
	my $phenotypes = $truth->{ association }{ $file }{ order };
	foreach my $phenotype (@$phenotypes) {
		if( ! exists( $test->{ association }{ $file }{ entries }{ $phenotype } )) {
			warn "ERROR: Missing association phenotype '$phenotype'\n";
			next;
		}

		my $i = $truth->{ association }{ $file }{ entries }{ $phenotype }{ probability };
		my $j = $test->{ association }{ $file }{ entries }{ $phenotype }{ probability };

		print STDERR "$file\t$phenotype\t$i\t$j\n";
		my $diff = abs( $i - $j );
		$stats->add_data( $diff );
	}
}

printf "Association Probability Difference\nmean=%.4f\tsdev=%.4f\n", $stats->mean(), $stats->standard_deviation();

# ============================================================
sub read_pathfx_results {
# ============================================================
	my $results = {};
	my $path = shift;
	opendir DIR, $path or die "Can't read '$path' $!";
	my @files = grep { /\.txt$/ } readdir DIR;
	closedir DIR;

	foreach my $file (@files) {
		if   ( $file =~ /neighborhood_\.txt$/ ) { $results->{ neighborhood }{ $file } = read_neighborhood( "$path/$file" );      } 
		elsif( $file =~ /assoc_table_\.txt$/ )  { $results->{ association }{ $file } = read_association_table( "$path/$file" ); }
	}
	return $results;
}

# ============================================================
sub read_neighborhood {
# ============================================================
	my $file = shift;
	my $neighborhood = {};
	open FILE, $file or die "Can't read '$file' $!";
	while( <FILE> ) {
		chomp;
		my ($a, $b, $score) = split /\t/;
		die "Redundant association" if exists $neighborhood->{ $a }{ $b };
		$neighborhood->{ $a }{ $b } = $score;
	}
	close FILE;

	return $neighborhood;
}

# ============================================================
sub read_association_table {
# ============================================================
	my $file  = shift;
	my $table = { entries => {}, order => [] };
	open FILE, $file or die "Can't read '$file' $!";
	my $header = <FILE>; chomp $header;
	my $headers = [ split /\t/, $header ];
	while( <FILE> ) {
		chomp;
		my @values = split /\t/;
		my $entry  = {};
		foreach my $field (@$headers) { $entry->{ $field } = shift @values; }
		my $phenotype = $entry->{ phenotype };
		die "Redundant phenotype found" if exists $table->{ entries }{ $phenotype };
		$table->{ entries }{ $phenotype } = $entry;
		push @{$table->{ order }}, $phenotype;
	}
	close FILE;

	return $table;
}
