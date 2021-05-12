package GlueingDatum;
use v5.20;
use strict;
use warnings;
use GraphTools;

use Data::Dumper;
use Time::HiRes qw(gettimeofday tv_interval);
use Math::MatrixReal;
use Math::Round qw/round/;
use Try::Tiny;  

=head1 Package Parameters
-\$ChipFiring::MAX\_ATTEMPTS = 13;

=head1 Attributes
bv   vertices of the base tree

=cut

use Class::Tiny qw/
bv
be   
bh
bt
bval
basetree
d
gv
ge
Gv
Gvlabel
reverseGvlabel
Ge   
reverseGelabel
Gh
Gt
G
Gval
Gvdan
Gvnondan
Gedan
g
Hv
He   
Hh
Ht
Hval
GinH
H
phi
Aphi
rphi
mphi
dangling
signDet
signMult
isGD
/;
    
#===================
sub BUILD {
    my ($self, $args) = @_;
    
    #parameters
    my @bv = @{$args->{'bv'}};
    my @be = @{$args->{'be'}};
    my %bh = %{$args->{'bh'}};
    my %bt = %{$args->{'bt'}};
    my $d = $args->{'d'};
    my %Gvlabel;
    my %reverseGvlabel;
    my %Gelabel;
    my %reverseGelabel;
    my %Gh;
    my %Gt;
    
    
    #if no edge glueings set to empty hash
    $self->ge({}) if (not $self->ge());
    
    ################################Make base tree###############################
    $self->basetree(GraphTools->new(
        v => $args->{'bv'},
        e => $args->{'be'},
        t => $args->{'bt'},
        h => $args->{'bh'},
    ));
    
    #say Dumper($self->basetree());
    
    ################################Check it is indeed GlueingDatum###############################
    
    #------------Connectedness------------
    #say Dumper( scalar finestCommonCoarsening([map {@{$_}} values % { $self->gv() }], $d)   );
    
    if (scalar finestCommonCoarsening([(map {@{$_}} values %{ $self->gv() }), (map {[$_]} (0..$d-1))], $d-1) > 1) {
    #say Dumper([map {@{$_}} values % { $self->gv() }]);
    #say Dumper($d);
    #say Dumper(  finestCommonCoarsening([map {@{$_}} values % { $self->gv() }], $d)   );
    #say Dumper(finestCommonCoarsening([map {@{$_}} values %{ $self->gv() }], $d));
        $self->isGD(0);
        return $self;
    }
    
    
    ################################Make G###############################
    #the code below is very inefficient
    #in a sense I was expecting being able to use pairs as hash keys...
    my %phi = (v=>{}, e=>{});
    my %mphi = (v=>{}, e=>{});
    
    
    #make vertices of G
    my $i = 0;
    foreach my $v (@bv) {
       $Gvlabel{$v} = {};
       if (defined $self->gv()->{$v}) {
            foreach my $gluerelation (@{$self->gv()->{$v}}) {
                $Gvlabel{$v}->{$_} = $i foreach @$gluerelation;
                $reverseGvlabel{$i} = $gluerelation;
                $phi{v}->{$i} = $v;
                $mphi{v}->{$i} = scalar @$gluerelation;
                $i++;
            }
       }
       
    
       
       foreach (0..$d-1) {
            if (not defined $Gvlabel{$v}->{$_}) {
                $Gvlabel{$v}->{$_} = $i;
                $reverseGvlabel{$i} = [$_];
                $phi{v}->{$i} = $v;
                $i++;
            }
       }
       
    }
    
    my @Gv = (0..$i-1);
    $self->Gv(\@Gv);
    
    #make edges of G
    $i = 0;
    foreach my $e (@be) {
       $Gelabel{$e} = {};
       if (defined $self->ge()->{$e}) {
            foreach my $gluerelation (@{$self->ge()->{$e}}) {
                $Gelabel{$e}->{$_} = $i foreach @$gluerelation;
                $reverseGvlabel{$i} = $gluerelation;
                $phi{e}->{$i} = $e;
                $mphi{e}->{$i} = scalar @$gluerelation;
                $Gh{$i} = $Gvlabel{ $bh{$e} }->{$gluerelation->[0]};
                $Gt{$i} = $Gvlabel{ $bt{$e} }->{$gluerelation->[0]};
                $i++;
            }
       }
       
       
       
       foreach (0..$d-1) {
            if (not defined $Gelabel{$e}->{$_}) {
                $Gelabel{$e}->{$_} = $i;
                $reverseGelabel{$i} = [$_];
                $phi{e}->{$i} = $e;
                $mphi{e}->{$i} = 1;
                $Gh{$i} = $Gvlabel{ $bh{$e} }->{$_};
                $Gt{$i} = $Gvlabel{ $bt{$e} }->{$_};
                $i++;
            }
       }
    }
    
    
    $self->Ge([0..$i-1]);
    $self->Gh(\%Gh);
    $self->Gt(\%Gt);
    $self->g( (scalar @{ $self->Ge() })  - (scalar @{ $self->Gv() }) + 1);
    
    #make G
    $self->G(GraphTools->new(
        v => $self->Gv(),
        e => $self->Ge(),
        t => $self->Gt(),
        h => $self->Gh()) );
    #say '========== edges:', Dumper($self->Ge()), ' =============';
    #say '========== vertices:', Dumper($self->Gv()), ' =============';
    
    #say Dumper($self->G());
    
    #########################Mark dangling##########################
    my %Gval = %{ $self->G()->val() };
    my %Gvdan;
    my %Gedan;
    my @Gvnondan = @Gv;
    my $nnondan = 0;
    
    while ($nnondan - scalar @Gvnondan) { #while this changes from iteration to iteration
        $nnondan = scalar @Gvnondan;
        @Gvnondan = grep { my $v = $_;
                     my @incident = @{  $self->G()->incidentEdges()->{$v}  };
                             if ($Gval{$v} - (scalar grep {$Gedan{$_}} @incident) < 2) {
                                #mark $v as dangling, and also all its incident edges
                                $Gvdan{$v} = 1;
                                $Gedan{$_} = 1 foreach @incident;
                                0;
                            } else {1};
                    } @Gvnondan;
    }
    
    $self->Gvdan(\%Gvdan);
    $self->Gedan(\%Gedan);
    $self->Gval(\%Gval);
    
    
    
    #########################Make H##########################
    $i = 0;
    my %Hval = (map {$_ => 
                $Gval{$_} - (scalar grep {$Gedan{$_}}
                                    @{  $self->G()->incidentEdges()->{$_}  }
                            )} @Gvnondan);
    my @Hv = grep {$Hval{$_} > 2} @Gvnondan;
    
    my @He = (0.. (scalar @Hv) + $self->g() - 2);
    my %Ht;
    my %Hh;
    
    
    
    my $Aphi = {map {$_ => {map {$_ => 0} @be}} @He}; #edge length map
    my %GinH;
    
    
    foreach my $A (@Hv) {
        #pick incident edges that are not dangling, and not yet in an edge of H
        my @incident = grep {(not $Gedan{$_}) and (not defined($GinH{$_}))} 
                        @{ $self->G()->incidentEdges()->{$A} };
        
        foreach my $e (@incident) {
           #follow along
           $Ht{$i} = $A;
           $GinH{$e} = $i;
           my $B = $self->G()->otherEndpoint($e,$A);
           $Aphi->{$i}{$phi{e}->{$e}} += 1/$mphi{e}->{$e};
           
           
           while ($Hval{$B} == 2) { 
                $e = [grep {(not $Gedan{$_}) and (not ($_ == $e) )} 
                                @{ $self->G()->incidentEdges()->{$B} }]->[0];
                $GinH{$e} = $i;
                $B = $self->G()->otherEndpoint($e, $B);
                $Aphi->{$i}{$phi{e}->{$e}} += 1/$mphi{e}->{$e};
           }
           
           $Hh{$i} = $B;
           $i++;
        }
    }
    
    $self->Hv(\@Hv);
    $self->He(\@He);
    $self->Ht(\%Ht);
    $self->Hh(\%Hh);
    $self->Hval(\%Hval);
    $self->Aphi($Aphi);
    
    #make H
    $self->H(GraphTools->new(
        v => $self->Hv(),
        e => $self->He(),
        t => $self->Ht(),
        h => $self->Hh()) );
        
    #say '============== genus:', $self->g(), '===============';
    #say Dumper($self->H());
    
    #########################Make edge-length matrix##########################
    return $self if ((scalar @He) != scalar @{ $args->{'be'}} );
    
    my $Am = [map {my $h = $_; [map {$Aphi->{$h}{$_}} @be]} @He];
    
    try {
    $self->signDet(
        round(abs(Math::MatrixReal->new_from_cols($Am)->det())) 
    );
    
    
    my $numberLeaves = scalar grep {$self->basetree()->val()->{$_} == 1} @bv;
    $self->signMult($self->signDet()/2**$numberLeaves);
    };
    $self->isGD(1);
    return $self;
}

#=====================Graph processing===================

#=====================Small subroutines==================
sub finestCommonCoarsening {
    my @gd = @{ shift() };
    my $d = shift;
    
    for my $i (0..$d) {
        my %head;
        my @tail;
        for (@gd) {
                    if (elementIn($i,@$_)) {@head{@$_} = ()}
                    else {push @tail, $_}
                }
                    
        @gd = ([keys %head], @tail);
    }
    
    return @gd;
}

sub elementIn {
    my $a = shift;
    for (@_) {
        return 1 if ($a == $_);
    }
    return 0;
}

1;
