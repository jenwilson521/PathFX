package PathFX;

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
sub read_results {
# ============================================================
	my $self    = shift;
	my $results = {};
	my $path = shift;
	opendir DIR, $path or die "Can't read '$path' $!";
	my @files = grep { /\.txt$/ } readdir DIR;
	closedir DIR;

	foreach my $file (@files) {
		if   ( $file =~ /neighborhood_\.txt$/ ) { $results->{ neighborhood }{ $file } = $self->read_neighborhood( "$path/$file" );      } 
		elsif( $file =~ /assoc_table_\.txt$/ )  { $results->{ association }{ $file }  = $self->read_association_table( "$path/$file" ); }
	}
	return $results;
}

# ============================================================
sub read_neighborhood {
# ============================================================
	my $self = shift;
	my $file = shift;
	my $neighborhood = {};
	open FILE, $file or die "Can't read '$file' $!";
	while( <FILE> ) {
		chomp;
		my ($a, $b, $score) = split /\t/;
		push @{$neighborhood->{ $a }{ $b }}, $score;
	}
	close FILE;

	return $neighborhood;
}

# ============================================================
sub read_association_table {
# ============================================================
	my $self  = shift;
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

1;
