*deck %W%  %G%
      Subroutine Options
C***Begin prologue     Options
C***Date written       850601   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C***Keywords           M1, Link 1, $Route, Input, Options
C***Author             Martin, Richard (LANL)
C***Source             %W%   %G%
C***Purpose            Describes the available route card options.
C***Description
C     The sequence of steps MESA goes through in determining
C     a molecular wavefunction is determined by the options specified
C     in the $Route section.  A complete list of the available options
C     follows:
C
C     Overlay 1:
C       Coord     The coordinates are to be read in Cartesian form.
C       RdChk     Read the geometry from the checkpoint file.
C       InAU      The distances are in atomic units.
C       InRad     The angles are in radians.
C       Timing    Print timing statistics for this link.
C
C     Overlay 2:
C       NoPrtDis     Don't print the distance matrix.
C       NoPrtAng     Don't print the interatomic angles.
C       NoAngByZ     Use internal cutoffs instead of the Z-Matrix to
C                    determine which angles to print.  This is the default,
C                    and only choice, for direct Cartesian coordinate input.
C       NoTetra      Don't set angles within 0.001 degree of 109.471
C                    to ACos(-1/3).
C       NoPrtZ       Don't print the Z-Matrix and resulting coordinates.
C       NoCrowdT     Don't abort the run if two atoms are less than
C                    0.5 Angstroms apart.
C       NoSymm       Unconditionally turn symmetry off.
C                    Note that Symm is still called, and will determine the
C                    framework group, but the molecule is not oriented.
C       Coord        The coordinates were read in Cartesian form.
C       Timing       Print timing statistics for this link.
C
C     Overlay 3:
C       Print_Basis       Print the basis set information with normalized
C                         primitives.
C       Print_ECP         Print the effective core potential information.
C       Basis=String      The default basis set is denoted by the string.
C       Print_S           Print the overlap integrals.
C       Print_T           Print the kinetic energy integrals.
C       Print_V           Print the potential energy integrals.
C       Print_LP          Print the effective core potential integrals.
C       PreExponential=N  Cutoff to use on preexponential factors is
C                         10**(-N).  (Default=10**-15).
C       CutOff=N          Cutoff to use on magnitude of integrals actually
C                         stored is 10**(-N)).  (Default=10**-10).
C       Timing            Print timing statistics for this link.
C
C     Overlay 4:
C       Guess_Core    In this case the guess comes from diagonalizing the
C                     core Hamiltonian (H=T+V), where T is the kinetic
C                     energy and V is the nuclear attraction operator.
C                     This type of guess has the advantage that it will
C                     work for any basis set, but is generally a fairly
C                     poor starting point.
C       Guess_Huckel  An extended Huckel calculation is performed to obtain
C                     starting vectors.  The Huckel calculation is done
C                     in a mimimum(STO3G) basis set, and the occupied
C                     orbitals are projected onto the current basis using
C                     the corresponding orbital transformation.  This type
C                     of guess is also available for any basis set, and is
C                     usually fairly decent.
C       Guess_RdChk   The initial vectors are read from the checkpoint file
C                     and projected onto the current basis.  This is
C                     usually an excellent guess if you have performed a
C                     previous calculation in a  different basis, on a
C                     different electronic state, or at a nearby geometry.
C       Guess_RdRWF   The initial vectors are read from the read/write
C                     file, and projected onto the current basis.
C       Guess_RdInp   The initial guess is read from the input file.
C                     The input is signaled through a $Vectors section,
C                     followed by each of the occupied orbitals, read in
C                     a list directed format.
C       Alter         Alter the initial guess vectors.
C       Print_Guess   Print the initial orbitals.
C       Timing        Print timing statistics for this link.
C
C     Overlay 5:
C       PRINT_VECTOR           Print the SCF vector
C       CONVERGENCE=n          DIIS error for convergence, as 10**(-n)
C                               (default 10**-8)
C       SCFCYC=n               Maximum number of iterations
C       PSEUDOCANONICAL=n      1 to use pseudocanonical orbitals (default)
C                              0 not
C       EXTRAPOLATE=n          1 to extrapolate fock matrices
C                              0 not (default)
C       DIAGONALS=n            1 to use F(i,i)
C                              0 to use i (default)
C       STARTDIIS=n            Convergence at which to initiate DIIS, as
C                              10**(-n)    (default 10**4)
C       TIMING                 Collect and print timing statistics
C
C
C***References
C***Routines called
C***End prologue       Options
      Implicit Integer(A-Z)
      Return
      End
