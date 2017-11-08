[![DOI](https://zenodo.org/badge/89211508.svg)](https://zenodo.org/badge/latestdoi/89211508)

# Forcer
Forcer is a FORM program for the integration of massless four-loop propagator diagrams.

Forcer reduces integrals using parametrically solved IBP identities, similar to the three-loop Mincer program. Forcer has been used for the computation of the finite pieces of four-loop propagators and vertices, Mellin moments of splitting functions and coefficient functions, and for the computation of the five-loop beta function for Yang-Mills theory with fermions.

More information about the Forcer program can be found in our [research paper](http://arxiv.org/abs/1704.06650). If you use this software, please cite our work.

### Authors
Forcer was written by Ben Ruijl, Takahiro Ueda, and Jos Vermaseren. This work has been supported by the ERC Advanced Grant no. 320651, HEPGAME.

We are much indebted to Andreas Vogt for helping with the testing of the programs and making many useful suggestions about the efficiency.

### Usage

Forcer requires the latest development version of [FORM](https://github.com/vermaseren/form).
The generation scripts -- which are not needed for running the program -- require Python 2.7/Python 3.0, igraph, numpy, and sympy.

An example of Forcer input is the following program:

```
#-
#include forcer.h

L F =
 +1/<p1.p1>/.../<p6.p6>*Q.p3*Q.p4*vx(Q,p1,p5,p6)*vx(-p1,p2,p3)
                       *vx(-p5,-p6,p4)*vx(-Q,-p2,-p3,-p4)
 +1/<p1.p1>/.../<p5.p5>*vx(-Q,p2,p3)*vx(p1,-p2,p5)*vx(-p1,p4,Q)
                       *vx(-p3,-p4,-p5)*ex(p1,p4)
;

#call Forcer(msbarexpand=4)
B ep;
P +s;
.end
```

##### Input specification

After ```forcer.h``` is included, the input integral can be defined.  This is done by specifying the vertex structure using ```vx```. The external momentum should be called ```Q```. The propagators and momenta can simply be multiplied in, as shown in the example above. Insertions on lines can be specified using the ```ex``` function. In the second integral above, ```ex(p1,p4)``` means that there is a single ε on the propagator associated with momentum p1 and one on p4. The provided topologies will automatically be matched to Forcer's internal topologies. The dot products will also automatically be rewritten.

##### Output options

By calling ```Forcer()```, the integrals are computed. The optional argument ```msbarexpand``` can give the (unrenormalized) answer expanded in MSbar. Otherwise, the result will be given exactly in terms of the master integrals and rational coefficients. Other options include ```polyratfunexpand=div``` and ```polyratfunexpand=maxpow```, which enable the expansions of rational coefficients in ε at intermediate steps..

