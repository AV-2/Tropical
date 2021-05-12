package GraphTools;
use v5.20;
use strict;
use warnings;

use Data::Dumper;

=head1 Package Parameters
-\$ChipFiring::MAX\_ATTEMPTS = 13;

=head1 Attributes

Make a graph:

    GraphTools->new(
        v => $v,
        e => $e,
        t => $t,
        h => $h);

=cut

our $VERSION = '0.2';

use Class::Tiny qw/
v
e
h
t
linkString
incidentEdges
val
g
/;
    
#===================
sub BUILD {
    my ($self, $args) = @_;

    my @v = @{$args->{'v'}};
    my @e = @{$args->{'e'}};
    
    $self->v( {map {$_ => 1} @v} );
    $self->e( {map {$_ => 1} @e} );
    
    #make link string
    my $string;
    my %h = %{$self->h()};
    my %t = %{$self->t()};
    
    $string .= $h{$_} . ',' . $t{$_} . ' ' foreach (@e);
    $self->linkString($string);
    
    #make list of incident edges
    my %incident = (map {$_ => []} @v);
    map {push @{$incident{$h{$_}}}, $_ } keys %h;
    map {push @{$incident{$t{$_}}}, $_ } keys %t;
    $self->incidentEdges(\%incident);
    
    #calculate valencies
    my %valencies = (map {$_ => scalar @{$incident{$_}}} keys %incident);
    $self->val(\%valencies);    
    
    #calculate genus
    $self->g((scalar @e - scalar @v) + 1);
    
    return $self;
}

#give edge and one endpoint, return the other endpoint
sub otherEndpoint {
    my ($self, $edge, $endpoint) = @_;
    
    return ($self->h()->{$edge} != $endpoint) ? $self->h()->{$edge} : $self->t()->{$edge};
}

#=====================Graph processing===================


1;

__END__


=Changelog

0.1 -- Initial release
0.2 -- Implements edge contraction
