package CumulativePlot;

use Point;

sub new {
	my ($class, $pointArrayRef) = @_;
	my $self = {
		'data' => $pointArrayRef
	};
	bless $self, $class;
	return $self;
}

sub setData {
	my ($self, $pointArrayRef) = @_;
	$self->{'data'} = $pointArrayRef;
}

sub getY {
	my ($self, $x) = @_;
	my $y;
	foreach my $point (@{$self->{'data'}}){
		return $point->y if ($x == $point->x());
	}
	
	return undef;
}

sub getX {
	my ($self, $y) = @_;
	my $x;
	foreach my $point (@{$self->{'data'}}){
		return $point->x if ($y == $point->y);
	}
	
	return undef;
}

## returns a reference to an array of Point
sub getData {
	my ($self) = shift;
	return $self->{'data'};
}

sub getAllX {
	my ($self) = shift;
	my @x;
	foreach my $point (@{$self->{'data'}}){
		push(@x, $point->x);
	}	
	return \@x;
}

sub getAllY {
	my ($self) = shift;
	my @y;
	foreach my $point (@{$self->{'data'}}){
		push(@y, $point->y);
	}	
	return \@y;
}

sub add {
	my ($self, $point) = @_;
	push(@{$self->{'data'}}, $point);
}

sub getFirstNoneZeroPoint {
	my ($self) = shift;
	foreach my $point (@{$self->{'data'}}){
		return $point if ($point->y > 0);
	}	
	return;
}

sub getLastPoint {
	my ($self) = shift;
	my $size=$#{$self->{'data'}};
	return $self->{'data'}->[$size];
}

sub getMidYBtw2Points {
	my ($self, $x1, $x2, $xMid) = @_;
	my $y2 = $self->getY($x2);
	my $y1 = $self->getY($x1);
	my $yMid = ($xMid - $x1) / ($x2 - $x1) * ($y2 - $y1) + $y1;
	return $yMid;
}

sub getKneePoint {
	my ($self) = shift;
	$self->sortX;
	
	my $point1 = $self->getFirstNoneZeroPoint;
	my $point2 = $self->getLastPoint;
	my $maxDiff=0;
	my $maxDiffLocX;
	foreach my $curPoint (@{$self->{'data'}}) {
		if ($curPoint->x >= $point1->x && $curPoint->x <= $point2->x) {
			my $yMid = $self->getMidYBtw2Points($point1->x, $point2->x, $curPoint->x);
			my $curYDiff = $self->getY($curPoint->x) - $yMid;
			if ($maxDiff < $curYDiff) {
				$maxDiff = $curYDiff;
				$maxDiffLocX = $curPoint->x;
			}
		}
	}
	return new Point($maxDiffLocX, $self->getY($maxDiffLocX));
}


sub sortX {
	my ($self) = shift;
	my @sortedData = sort { $a->x <=> $b->x } @{$self->{'data'}}; 
	$self->{'data'} = \@sortedData;
}

sub sortY {
	my ($self) = shift;
	my @sortedData = sort { $a->y <=> $b->y } @{$self->{'data'}}; 
	$self->{'data'} = \@sortedData;
}


1;