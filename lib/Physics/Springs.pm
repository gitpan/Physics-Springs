package Physics::Springs;

use 5.006;
use strict;
use warnings;

our $VERSION = '0.02';

use Carp;
use Physics::Particles;

use base 'Physics::Particles';

sub new {
   my $proto = shift;
   my $class = ref($proto) || $proto;

   my $self = $class->SUPER::new();
   $self->{_PhSprings_springs} = [];
   return $self;
}

sub add_spring {
   my $self   = shift;
   my %args   = @_;
   my $spring = {};
   my $k      = $args{k};
   my $p1     = $args{p1};
   my $p2     = $args{p2};
   my $l      = $args{l};

   defined($k) && defined($p1) && defined($p2)
     or croak("You need to supply several named arguments.");

   # default to being relaxed
   if (not defined $l) {
      my $dist = sqrt(
                  ( ($self->{p}[$p1]{x} - $self->{p}[$p2]{x}) )**2 +
                  ( ($self->{p}[$p1]{y} - $self->{p}[$p2]{y}) )**2 +
                  ( ($self->{p}[$p1]{z} - $self->{p}[$p2]{z}) )**2
                 );
      $l = abs($dist);
   }

   $spring->{k}   = $k;
   $spring->{p1}  = $p1;
   $spring->{p2}  = $p2;
   $spring->{len} = $l;

   push @{$self->{_PhSprings_springs}}, $spring;
   return $#{$self->{_PhSprings_springs}};
}

sub iterate_step {
   my $self      = shift;
   my @params    = @_;
   my $time_diff = $params[0];
   foreach my $spring (@{$self->{_PhSprings_springs}}) {
      my $p1 = $self->{p}[$spring->{p1}];
      my $p2 = $self->{p}[$spring->{p2}];
      my $l  = $spring->{len};
      my $k  = $spring->{k};

      my $dist = sqrt(
                  ( ($p1->{x} - $p2->{x}) )**2 +
                  ( ($p1->{y} - $p2->{y}) )**2 +
                  ( ($p1->{z} - $p2->{z}) )**2
                 );

      my $force1 = $k * ($dist-$l);
      my $force2 = -$force1;

      my $acc1 = $force1 / $p1->{m};
      my $acc2 = $force2 / $p2->{m};
      $acc1 = [
                map { $acc1 * ($p2->{$_} - $p1->{$_}) }
                    qw/x y z/
              ];
      $acc2 = [
                map { $acc2 * ($p2->{$_} - $p1->{$_}) }
                    qw/x y z/
              ];
  
      $p1->{x}  += $p1->{vx} * $time_diff +
                   $acc1->[0]*0.5*$time_diff**2;
      $p1->{y}  += $p1->{vy} * $time_diff +
                   $acc1->[1]*0.5*$time_diff**2;
      $p1->{z}  += $p1->{vz} * $time_diff +
                   $acc1->[2]*0.5*$time_diff**2;

      $p2->{x}  += $p2->{vx} * $time_diff +
                   $acc2->[0]*0.5*$time_diff**2;
      $p2->{y}  += $p2->{vy} * $time_diff +
                   $acc2->[1]*0.5*$time_diff**2;
      $p2->{z}  += $p2->{vz} * $time_diff +
                   $acc2->[2]*0.5*$time_diff**2;
  
      $p1->{vx} += $acc1->[0]*$time_diff;
      $p1->{vy} += $acc1->[1]*$time_diff;
      $p1->{vz} += $acc1->[2]*$time_diff;

      $p2->{vx} += $acc2->[0]*$time_diff;
      $p2->{vy} += $acc2->[1]*$time_diff;
      $p2->{vz} += $acc2->[2]*$time_diff;
   }

   $self->SUPER::iterate_step(@params);
}

1;
__END__

=head1 NAME

Physics::Springs - Simulate Particle Dynamics with Springs

=head1 SYNOPSIS

  use Physics::Springs;

=head1 ABSTRACT

  Simulate particle dynamics with springs.

=head1 DESCRIPTION

This module is intended as an add-on to the Physics::Particles module
(version 0.10) and may be used to simulate particle dynamics
including spring-like forces between any two particles you specify.

The module extends the API of Physics::Particle by one method which is
documented below. Please see the documentation to Physics::Particle for
more information about the API.

There are several particle properties required by Physics::Springs in
order to work: These are the x/y/z coordinates, the vx/vy/vz
velocity vector components, and a non-zero mass 'm'.

=head2 Methods

=over 2

=item add_spring

You may use the add_spring method to add springs to the system of particles.
Each spring has a starting and end particle, a relaxed length, and a spring
constant 'k'.

add_spring expects several named arguments. Required arguments are:
The spring constant k, the starting (p1) and end (p2) points.
Optional arguments are: The length of the relaxed spring 'len'.
If len is not specified, the current distance between p1 and p2 will be
used as the length of the relaxed spring.

=back

=head1 SEE ALSO

L<Physics::Particles>

L<Math::Project3D>, L<Math::Project3D::Plot> for a reasonably
simple way to visualize your data.

=head1 AUTHOR

Steffen Mueller, E<lt>springs-module at steffen-mueller dot netE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2003 by Steffen Mueller

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. 

=cut
