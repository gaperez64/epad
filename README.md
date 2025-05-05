# Tools for the existential fragment of Presburger arithmetic with divisibility
The first goal of this repository is to have a library to manipulate systems
of divisibility constraints over linear (diophantine) polynomials conjoined
with equations and inequalities.
1. Normalization of such formulas (based on ["On the Complexity of Linear
   Arithmetic with Divisibility", by Joël Ouaknine, Antonia Lechner and Ben
   Worrell.](https://www.cs.ox.ac.uk/people/james.worrell/LICS-main.pdf))
2. Solvers for normalized formulas (based on the same work mentioned above,
   and an alternative one based on ["Integer Programming with GCD
   Constraints", by Rémy Défossez, Christoph Haase, Alessio Mansutti, and Guillermo
   A. Pérez.](https://epubs.siam.org/doi/10.1137/1.9781611977912.128))

## (Wrong) design choices
- We are only encoding systems, so no disjunctions.
- Variables are quantified over the integers for now. This is a bit off since
  after an affine change of variables we go to natural numbers instead. To
  remediate this, we are asking that a system of linear constraints is always
  given to represent that domain of the variables explicitly. (That is,
  working on a system of divisibilities with no further constraints implicitly
  tells the library that we have changed from integers to naturals as domain.)
- Everything is tuples (of tuples).
- We don't encode Ax = b but rather Ax - b = 0 and the last "coefficient" is
  the constant in the polynomial. Similar weirdness applies to inequalities.

## WIP: Lipshitz' normalization
Given a system of divisibilities and linear constraints, we can implement
Lipshitz' transformation into increasing normal form if we have:
- [X] A way of obtaining a (hybrid) semilinear representation of the solutions
  of a system of linear equations---which can be used to update the system of
  divisibilities via an affine change of variables while not increasing the
  number of variables
- [X] A way of enumerating the sign of the left-hand sides of divisibility
  constraints, towards obtaining left-hand sides with all positive
  coefficients, and adding them as inequalities to then remove them via the
  method described above
- [ ] A way of obtaining a finite representation of the module of the
  primitive polynomials of the left-hand sides
- [ ] A way of checking whether a given system is increasing for a given
  (semantic) order over the variables, based on the representation above.
  (NOTE: This should be as easy as checking whether all elements of the
  spanning set of the module have 0 in the entries larger than the leading
  variable of the LHS?)
- [ ] ...

## Dependencies
For now, Z3 and flint
