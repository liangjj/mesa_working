*deck %W%  %G%
      Subroutine Intro
C***Begin prologue     Intro
C***Date written       850601   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C***Keywords           Intro
C***Author             Martin, Richard  (LANL)
C***Source             %W%   %G%
C***Purpose            General introduction to the MESA system.
C***Description
C     Introduction
C
C       Suppose you wanted to perform a Hartree-Fock calculation on the
C     water molecule.  An acceptable form of the input to MESA would
C     look like:
C
C
C       $Route
C       HF, 6-31G, Opt
C          $Title
C          Double-zeta geometry optimization for H2O.
C       $Geom
C          0 1
C       O
C       H1 O R
C          H2 O R    H1 Theta
C
C       R=0.96
C       Theta=104.5
C       $
C
C
C       This is an energy optimization at the Hartree-Fock level for the
C     water molecule.  The OH bond distance is initially 0.96
C     Angstroms, and the bond angle is 104.5 degrees.
C
C
C
C         Note that the deck  consists of a series of input sections,
C     the parts that begin with $...  Each input section ends when you
C     hit another $.  The end of one section is the beginning of
C     another section, except of course for the last one.  The input
C     sections may occur in any order within the deck. Within each
C     section are directives telling it what you want.  In the
C     example, the $Route section has a directive specifying a
C     Hartree-Fock calculation,  a directive that says to put a 6-31G
C     basis set on all the atoms, and a directive that says to
C     optimize the geometry.  This section is always required, and
C     essentially governs the type of calculation to be done. The
C     $Title section simply describes the job. You can put anything
C     you want to here.  It's a completely optional section.  The
C     $Geom section is required in some form in all calculations.  The
C     first directive card specifies the molecular charge and the spin
C     multiplicity of the state.  We're doing neutral water and the
C     singlet spin state.  The next three cards describe atoms in the
C     molecule.  Note that it's usually convenient to  distinguish
C     them with numerical identifiers like H1 and H2.  The H1-O bond
C     distance is specified by a parameter, R.  The H2-O bond distance
C     is specified by the same parameter, while the H2-O-H1 bond angle
C     is called Theta.  After all the atoms are specified, the
C     parameters are given values, in this case R=0.96, and
C     Theta=104.5.  The units are Angstroms for distances and degrees.
C
C
C         There are many input sections available so you can (we hope)
C     do what you want to do and put the programs through the
C     particular hoops you think are necessary.  We've tried to make
C     things flexible.  On the other hand, no one wants to prepare a
C     100 line input deck to do a calculation on the H2O molecule, so
C     most of the input sections default to reasonable procedures.
C     Usually you will be preparing $Route, $Title, and $Geom sections
C     and that's all.
C
C***References
C***Routines called
C***End prologue       Intro
      Implicit Integer(A-Z)
      Return
      End
