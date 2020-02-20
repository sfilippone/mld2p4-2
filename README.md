
                         MLD2P4  
MultiLevel Domain Decomposition Parallel Preconditioners Package
           based on PSBLAS (Parallel Sparse BLAS version 3.5)
    
Salvatore Filippone    (Cranfield University, UK)
Pasqua D'Ambra         (IAC-CNR, Naples, IT)
Daniela di Serafino    (Univ. of Campania "L. Vanvitelli", Caserta, IT)

---------------------------------------------------------------------

MLD2P4 (MultiLevel Domain Decomposition Parallel Preconditioners
Package based on PSBLAS) provides parallel Algebraic MultiGrid (AMG)
and Domain Decomposition preconditioners, to be used in the
iterative solution of linear systems.

The name of the package comes from its original implementation,
containing multilevel additive and hybrid Schwarz preconditioners,
as well as one-level additive Schwarz preconditioners. The current
version extends the original plan by including multilevel cycles
and smoothers widely used in multigrid methods. A purely algebraic
approach is applied to generate coarse-level corrections, so that
no geometric background is needed concerning the matrix to be
preconditioned.

MLD2P4 has been designed to provide scalable and easy-to-use 
preconditioners in the context of the PSBLAS (Parallel Sparse Basic
Linear Algebra Subprograms) computational framework and is used
in conjuction with the Krylov solvers available from PSBLAS. The
package employs object-oriented design techniques in Fortran 2003,
with interfaces to additional third party libraries such as MUMPS,
UMFPACK, SuperLU, and SuperLU_Dist, which can be exploited in building
multilevel preconditioners. The parallel implementation is based on
a Single Program Multiple Data (SPMD) paradigm; the inter-process
communication is based on MPI and is managed mainly through PSBLAS.


MAIN REFERENCE:

P. D'Ambra, D. di Serafino, S. Filippone,
MLD2P4: a Package of Parallel Algebraic Multilevel Domain Decomposition
Preconditioners in Fortran 95,
ACM Transactions on Mathematical Software, 37 (3), 2010, art. 30,
doi: 10.1145/1824801.1824808.


TO COMPILE

0. Unpack the tar file in a directory of your choice (preferrably
   outside the main PSBLAS directory).
1. run configure --with-psblas=<ABSOLUTE path of the PSBLAS install directory>
   adding the options for MUMPS, SuperLU, SuperLU_Dist, UMFPACK as desired.
   See MLD2P4 User's and Reference Guide (Section 3) for details.
2. Tweak Make.inc if you are not satisfied.
3. make; 
4. Go into the test subdirectory and build the examples of your choice.
5. (if desired): make install 


NOTES

- The single precision version is supported only by MUMPS and SuperLU;
  thus, even if you specify at configure time to use UMFPACK or SuperLU_Dist, 
  the corresponding preconditioner options will be available only from
  the double precision version.

- The preconditioners in MLD2P4 extend those of PSBLAS and are meant
  to be used with the PSBLAS Krylov solvers; so in an existing program
  you need to modify the type of the preconditioner object and its
  settings, but the rest of the application needs not be changed. 
   

The MLD2P4 team. 
---------------
Project lead:
Salvatore Filippone

Contributors:
Pasqua     D'Ambra
Daniela    di Serafino
Ambra	   Abdullahi Hassan
Alfredo    Buttari
