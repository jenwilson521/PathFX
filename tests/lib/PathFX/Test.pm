package PathFX::Test;

use base PathFX;
use Statistics::Descriptive;

# ============================================================
sub new {
# ============================================================
	my ($class) = map { ref || $_ } shift;
	my $self = bless {}, $class;
	$self->init( @_ );
	return $self;
}

# ============================================================
sub init {
# ============================================================
	my $self = shift;
}

# ============================================================
sub compare_results {
# ============================================================
	my $self  = shift;
	my $test  = shift;
	my $truth = shift;
	my $stats = new Statistics::Descriptive::Full();

	$truth = $self->read_results( $truth );
	$test  = $self->read_results( $test );

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
				my $x = [ sort { $a <=> $b } @{$truth->{ neighborhood }{ $file }{ $a }{ $b }}];
				my $y = [ sort { $a <=> $b } @{$test->{ neighborhood }{ $file }{ $a }{ $b }}];

				# print STDERR "$file\t$a\t$b\t$i\t$j\n";
				foreach my $i ( 0 .. $#$x ) { 
					my $diff = abs( $x->[ $i ] - $y->[ $i ] );
					$stats->add_data( $diff );
				}
			}
		}
	}

	printf "Neighborhood Difference\n  mean=%.4f\tsdev=%.4f\n", $stats->mean(), $stats->standard_deviation();

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

			# print STDERR "$file\t$phenotype\t$i\t$j\n";
			my $diff = abs( $i - $j );
			$stats->add_data( $diff );
		}
	}

	printf "Association Probability Difference\n  mean=%.4f\tsdev=%.4f\n", $stats->mean(), $stats->standard_deviation();
}

1;
