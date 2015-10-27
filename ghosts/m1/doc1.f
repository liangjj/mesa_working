*deck %W%  %G%
      subroutine doc()
C***Begin prologue     $Route
C***Date written       850601   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C***Keywords           M1, Link 1, Input, $Route
C***Author             Martin, Richard (LANL)
C***Source             %W%   %G%
C***Purpose            This section presents the basis route information.
C***Description
C       This section specifies the calculation MESA is to
C     perform.  It is always required.  The route is specified in free
C     format by keyword commands which specify the quantum mechanical
C     procedure, the basis set and other control features involving
C     optimization or other aspects of the execution.  Some keywords
C     are associated with options.  Three formats are available for
C     specifying a keyword:
C                            KeyWord
C                            KeyWord=Option
C                            KeyWord="Option1, Option2, ..."
C         If only one option is desired, either of the last two formats
C     is acceptable, for two or more only the last is acceptable.
C
C         The most important and widely used keywords are described
C     briefly below.  A complete list of available keywords and
C     options is given in the appendix.
C
C
C     Quantum mechanical procedure:
C         This keyword specifies the quantum mechanical method to be
C     used. Possibilities are:
C         HF      This requests that only an RHF calculation be
C     performed.  It is not available for states that are not the
C     lowest of their multiplicity, i.e. open-shell singlets can't be
C     done.
C
C         GVB     This requests that a GVB perfect pairing calculation be
C     performed.  It is usually accompanied by an option specifying
C     the number of pairs, etc.  This is described in more detail
C     later.  If the number of pairs is zero, this keyword can be used
C     to perform Hartree-Fock  calculations on states not the lowest
C     of their spin multiplicity, i.e. open shell singlets can be done.
C
C         CI      This requests that a configuration-interaction
C     calculation be performed.  It is usually specified in
C     conjunction with one of the previous keywords, which determines
C     where the one-electron orbitals used in the CI expansion are
C     obtained.  If this is the only method keyword, HF is assumed.
C
C         LCC     This requests that a linear coupled-cluster calculation
C     be performed.  It is usually specified in conjunction with one
C     of the previous keywords, which determines where the
C     one-electron orbitals used in the LCC expansion are obtained.
C     If this is the only method keyword, HF is assumed.
C
C         BASIS   This keyword specifies the default basis set to be
C     used.  A number of basis sets are stored on the Dat file. A
C     complete list can be found in the appendix along with
C     instructions on how to add new ones.  The most commonly used
C     bases are DZ, DZP, 6-31G, 6-31G**.  Basis sets can be mixed and
C     matched on various centers through the $Geom section if desired.
C
C
C     Optimization control:
C         These keywords control the calculation of energy derivatives
C     and the subsequent search for a stationary point on the potential
C     surface.
C
C         Opt      requests that a geometry optimization be performed
C     using analytically computed gradients.  The parameters defined
C     in the Z-Matrix are varied until the minimum on the
C     potential surface is found.  An equivalent way to specify this
C     option is Opt=Grad.
C
C         Opt=FP   requests that a geometry optimization be performed
C     using numerically computed gradients.  The parameters defined in
C     the Z-Matrix are varied in turn.  This option is available for all
C     types of wavefunctions, although it quickly becomes expensive
C     for systems with a large number of degrees of freedom.
C
C***References
C
C***Routines called
C***End prologue       $Route
C***Begin prologue     $Title
C***Date written       850601   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C***Keywords           M101, Link 101, Input, $Title
C***Author             Martin, Richard  (LANL)
C                      Binkley, Steven  (Gauss82)
C***Source             %W%   %G%
C***Purpose            Reads a descriptive title for the run.
C***Description
C                      The title section be up to six cards long.  It is
C                      echoed in the output, but serves no other purpose.
C
C***References
C***Routines called
C***End prologue       $Title
      return
      end 
