## Sage code and maple worksheets

[V. Neiger](https://www.unilim.fr/pages_perso/vincent.neiger), 
[B. Salvy](http://perso.ens-lyon.fr/bruno.salvy), and 
[G. Villard](http://perso.ens-lyon.fr/gilles.villard) 


Software material companion to a draft on modular polynomial composition, inverse composition, and modular power projection. 

*Page under construction.*

### Maple worksheets

### Sage code

The folder ``sage`` contains [SageMath](http://www.sagemath.org/) code. Some of
the functions require a recent version of SageMath such as SageMath 8.4.

Files in this folder:
* Four files with implementations of the algorithms for balanced basis
  computation, modular composition, inverse modular composition, and power
  projections. The code contains comments and offers to turn on a verbose mode
  which, when the function is called, prints some text explaining what are its
  main steps. It also gives timings for each step, although they should only be
  considered as indications for such prototype implementations.
* One file ``examples.sage`` containing simple example uses of these functions.
* One file ``test-suite.sage`` containing tests for each of the defined
  functions; simply running this file in Sage will run these functions
  on many instances and verify that they meet their specification.
* One file ``utils.sage``, containing a few basic helper functions.

Files in ``subroutines`` folder:
* code for computing _approximant bases_, used for reconstructing the
  denominator of polynomial matrix fraction representing linearly recurrent
  matrix sequences.
* code for performing _matrix division with remainder_, based on linear system
  solving over the univariate polynomials; this is used to find remainders
  modulo a balanced basis.
