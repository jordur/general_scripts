#!/usr/bin/perl
# 2011/10/14 - arbol sw.

package SOLiD_Reads;
sub new
{
    my $class = shift; # the first parameter is the class name (Reads_SOLiD)
	my $self = {
		_TAG_panel => shift,
		_TAG_x => shift,
		_TAG_y => shift,
		_TAG_type => shift,
		_TAG_additional => shift,
		_sequence  => shift
	};

    # Create and return object
    bless $self, $class;
    return $self;
}


sub create_from_tag_seq
{
    my ( $tag_id, $sequence )=@_; # assignation of parameters to the function
    chomp $tag_id;
    chomp $sequence;
    my @tag_aux=(split /,/, $tag_id);
    $tag_aux[0]=~ s/>//g;
    my @tag=(split /_/, $tag_aux[0]);
	if ($#tag_aux == 0) {
		$read = new SOLiD_Reads ($tag[0], $tag[1], $tag[2], $tag[3], "", $sequence);
	}
	else {
		my $additional="";
		foreach (@tag_aux[1..$#tag_aux]) {
			$additional=$additional . "," . $_;
		}
		$read = new SOLiD_Reads ($tag[0], $tag[1], $tag[2], $tag[3], $additional, $sequence);
	}
    return $read;
}


sub reads_compare 
{
# returns -1 if the read1 is considered minor than read2, +1 if read2 is considered minor than read1, and 0 if both reads are considered the same (without checking their respective F3, R3 or F5 type. For this aim, the bead TAG_ID will be checked (panel,X_coordinate,Y_coordinate)
	my ( $read1, $read2 ) = @_;

	if ( $read1->{_TAG_panel} < $read2->{_TAG_panel} ) {
		return -1;
	}
	elsif ( $read1->{_TAG_panel} > $read2->{_TAG_panel} ) {
		return +1;
	}
	else{
		if ( $read1->{_TAG_x} < $read2->{_TAG_x} ) {
			return -1;
		}
		elsif ( $read1->{_TAG_x} > $read2->{_TAG_x} ) {
			return +1;
		}
		else {
			if ( $read1->{_TAG_y} < $read2->{_TAG_y} ) {
				return -1;
			}
			elsif ( $read1->{_TAG_y} > $read2->{_TAG_y} ) {
				return +1;
			}
			else {
				return 0;
			}
		}
	}
}


sub return_file_lines
{
# this function returns the TAG_ID and sequence file lines that can be directly appended to the desired fasta/qual file
	my ($self)=@_;
	my $TAG_ID = ">" . $self->{_TAG_panel} . "_" . $self->{_TAG_x} . "_" . $self->{_TAG_y} . "_" . $self->{_TAG_type};
	if ( $self->{_TAG_additional} ne "" ) {
		$TAG_ID=$TAG_ID . $self->{_TAG_additional} . "\n";
	}
	else {
		$TAG_ID=$TAG_ID . "\n";
	}
	my $seq = $self->{_sequence} . "\n";
	return ( $TAG_ID, $seq);
}

1;
