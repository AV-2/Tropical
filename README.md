# Tropical Perl
Several perl programs to manipulate graphs in the context of tropical geometry, specifically chip-firing and tropical morphisms. See the following articles and master thesis:

* (Var21)
* (DV20)
* (Var16)(Chapter 3)

## Aims 
* Graph Theory
    * Implement a graph object with basic functionality
    * Implement graph morphisms, including contraction morphisms
    * Walk around the moduli space of graphs by contracting one edge
* ChipFiring
    * Implement Dhar burning algorithm
    * Calculate reduced divisors
    * Certify whether two divisors are equivalent
    * Calculate rank of a divisor
    * Find all rank-1 divisors with a given degree
* Tropical morphisms
    * Calculate the length-map associated to a tropical morphism
    * Calculate the multiplicity of a tropical morphism
    * Certify whether two given tropical morphisms are isomorphic
    * Walk around the moduli space of tropical morphisms
* Glueing datums
    * Obtain the tropical morphism associated to a glueing datum
    * Deform glueing datums

## Installing
The following dependencies must be installed via cpan.
(!!to-do)

## Quickstart
The following code shows how to create Graph, Divisor, Morphism, and Glueing Datum objects.

### Creating a graph

### Querying information about a graph

### Graph contractions

### Creating a morphism

### Querying information about a morphism

### Creating a divisor

### Querying information about a divisor

### Pulling back a divisor

### Creating a glueing datum

### Deforming a glueing datum
## License
GPLv3

## Changelog
* v0.1 -- Initial release
    * Basic graph functionality, a lot of routines for chip-firing, basic tropical morphism and glueing datum.
