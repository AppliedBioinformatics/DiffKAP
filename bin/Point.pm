package Point;

sub new {
	my $class = shift;
	my $self = {
		'x' => shift,
		'y' => shift,
	};
	bless $self, $class;
	return $self;
}

sub setX {
	my ($self, $x) = @_;
	$self->{'x'}=$x;
}

sub setY {
	my ($self, $y) = @_;
	$self->{'y'}=$y;
}

sub setXY {
	my ($self, $x, $y) = @_;
	if(ref($x) eq 'ARRAY'){
		$self->{'x'}=$x->[0];
		$self->{'y'}=$x->[1];
	} else {
		$self->{'x'}=$x;
		$self->{'y'}=$y;
	}
}

sub x {
	my ($self) = shift;
	return $self->{'x'};
}

sub y {
	my ($self) = shift;
	return $self->{'y'};
}

sub xy {
	my ($self) = shift;
	return [$self->{'x'}, $self->{'y'}];
}


sub print {
	my ($self) = shift;
	print "$self->{'x'} $self->{'y'}\n";
}

sub add {
	my ($self, $inPoint) = @_;
	$self->{'x'} += $inPoint->{'x'};
	$self->{'y'} += $inPoint->{'y'};
}

sub subtract {
	my ($self, $inPoint) = @_;
	$self->{'x'} -= $inPoint->{'x'};
	$self->{'y'} -= $inPoint->{'y'};
}

1;
