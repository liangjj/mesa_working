*deck @(#)doc101.f	5.1  11/6/94
      subroutine doc101()
c***begin prologue     $geom
c***date written       850601   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords
c                      m101, link 101, input,
c                      geometry, basis set, $geom,
c                      coord, rdchk, inau, inrad
c***author             martin, richard  (lanl)
c                      binkley, steven  (gauss82)
c***source             @(#)doc101.f	5.1   11/6/94
c***purpose
c                      reads the molecular geometry, and basis set.
c***description
c       the $geom section is read by m101 and is always required,
c     except in the case of a restart or a geom=rdchk or geom=chk directive.
c     it reads the molecular geometry and basis sets.
c     the input is free-field: the several items on a card may be
c     separated by one or more blanks.
c         the first set set of cards specifies the nuclei.  there are three
c     formats which you may want to use at one time or another.  they
c     are discussed in turn below.
c
c
c
c     z-matrix format:
c       element(basis,nuclear charge),n1,length,n2,alpha,n3,beta,j
c
c         the z-matrix format is the default in mesa.  a series of
c     cards, one for each center, is read in the format above.
c
c
c         element specifies the chemical nature of the nucleus.     it
c     may consist of simply the chemical symbol, such as 'h' or 'pt',
c     or it may begin with the chemical symbol followed immediately by
c     a numeric identifier, e.g. 'pt4' to designate the fourth
c     platinum atom.  dummy nuclei are denoted by the symbols 'x'
c     or'-'. this item is required for every center.  the remainder of
c     this field is optional and specifies the basis set to be used on
c     the center.  if unspecified, the default basis from the route
c     card is supplied.  if you wish to specify the basis here, enter it's
c     name (16 characters maximum) equated to basis
c     in the parentheses; e.g. (basis=sto-3g).  a list of the
c     basis sets known to mesa is available in mesadat. the nuclear charge may
c     also be entered here -- this option is useful for overwriting the
c     standard nuclear charge associated with the atomic symbol. for example,
c     (basis=nobasis,z=-1.0) can be used to introduce a -1.0 point charge
c     with no basis set being added to the center.
c
c
c          n1 specifies the (previously defined) nucleus for which the
c     internuclear distance from center n1 to center n will be given.
c     this item may be either an integer (with n1<n) giving the number
c     of the card defining the other center or the alphanumeric
c     label of the other center.
c
c
c         length is the internuclear distance r(n1,n).  it may be either
c     a positive floating point number in angstroms(unless modified by
c     the inau directive in $route) or an alphanumeric string of up to
c     16 characters.  in the latter case, the distance is represented
c     by a parameter, for which a value will be supplied in a later
c     section of $geom.  if a length is to be varied in an
c     optimization run, it must be specified as a parameter.  note
c     that a given parameter may be used to specify two distinct
c     lengths, in which case they are constrained to be equal.  the
c     items n1 and length are required for all nuclei after the first.
c     for the second nucleus, only element, n1, and length are
c     required.  the nucleus is placed on the z-axis.
c
c
c         n2 specifies the center for which the angle   alpha(n,n1,n2) is
c     given.  again, it may be either an integer or an alphanumeric
c     string which matches a previous element entry.  note that n1 and
c     n2 must represent different nuclei.
c
c         alpha is the internuclear angle.     this may be a floating
c     point number giving the angle in degrees (unless modified by the
c     inrad directive in $route) or an alphanumeric string representing
c     a parameter.  the numerical value of alpha must lie in the range
c     0<alpha<180 degrees.  n2 and alpha are required for all centers after
c     the second.  for the third center, only element,n1,length,n2,and
c     alpha are required.  this center is placed in the xz-plane.
c
c
c
c         n3 differs in meaning,as does beta, depending on the value of
c     the last item, j.  if j=0 or is omitted, n3 specifies the center
c     for which the internuclear dihedral angle beta(n,n1,n2,n3) will
c     be given.  as with n1 and n2, this may be either an integer or
c     an alphanumeric string matching a previous element entry.  this
c     center must be distinct from those specified by n1 and n2.
c
c         beta is the internuclear dihedral angle (if j=0,see below).
c     this may be a floating point number giving the angle in degrees
c     (unless modified by $route), or a(signed) alphanumeric string
c     representing a parameter.  the dihedral angle is defined as the
c     angle between the planes defined by (n,n1,n2) and (n1,n2,n3).
c     beta must lie in the range 0<beta<180 degrees.  the sign is positive if
c     the movement of the directed vector n1-->n2 towards the directed
c     vector n2-->n3 involves a right-handed screw motion.  figure 1
c     may be helpful to those of you who either flunked vector
c     analysis or are not too handy around the house.
c
c
c
c         j allows you some freedom over the specification of the
c     dihedral angle.  although it is always possible to specify the
c     center by a bond length, a bond angle and a dihedral angle(which
c     is what you must do if j=0 or is absent), it is sometimes more
c     convenient to use a bond length and two bond angles.  this is
c     called for by using j=11. in this case, the internuclear angle
c     beta(n,n1,n3) is specified.  as before, this may be either a
c     floating point value in degrees or an alphanumeric string
c     representing a parameter.  in the event that you specify two
c     bond angles, there will be two possible positions for the center
c     n.  this is fixed by the sign of j.  thus j=+1 if the triple
c     vector product n1-n(n1-n2 x n1-n3) is positive and j=-1 if the
c     product is negative.
c
c
c
c         if you have not specified any parameters, this completes the
c     $geom section.  if you have then you must specify one last set
c     of information, the initial values of the parameters.  this set
c     of data must be separated from the z-matrix by a blank card,
c     after which one card for each parameter is read in the form:
c
c         parameter=value,  d2e=value, type
c
c         parameter is the alphanumeric string defined by the preceeding
c     section and value is a floating point number representing the
c     initial value( in angstroms and degrees unless modified by
c     $route directives.
c
c         d2e is an optional field.  if  specified, it denotes how mesa
c     obtains initial guesses at the force constants to be used in the
c     optimization.  if omitted, internally stored guesses are
c     provided. the d2e value may be a floating point number, which
c     then specifies the initial diagonal second derivative associated
c     with this parameter( in hartrees/sq-bohr, or hartrees/sq-radian),
c     or the string 'numer', which causes the program to do a few
c     points numerically before doing an optimization step. 'numer' is
c     normally specified for several parameters, in which case the
c     program can get guesses at both the diagonal and off-diagonal
c     force constants involving them.  the d2e specification is
c     ignored for opt=ms and opt=fp runs.
c
c
c
c         type denotes  the status of the parameter.  it can be either
c     'constant', in which case the parameter is not varied in an
c     optimization, or 'variable', the default type.  this is included
c     for convenience only.  a 'constant' could have been just as
c     easily entered directly in the z-matrix as a floating point
c     number.
c
c
c
c
c
c
c     coord format:
c
c       element(basis), x, y, z
c
c
c         this case is specified by the coord directive in $route.  the
c     input section consists of a number of cards, one for each
c     center, which specify the element, basis set, and cartesian
c     coordinates of the element.
c
c
c         element is either a string (as above),or an integer specifying
c     the atomic number of the center.
c         x is a floating point number giving the x-coordinate of the
c     center in angstroms (unless modified by a directive in $route).
c         y is a floating point number giving the y-coordinate of the
c     center in angstroms (unless modified by a directive in $route).
c         z is a floating point number giving the z-coordinate of the
c     center in angstroms (unless modified by a directive in $route).
c
c
c
c
c
c         $geom
c         0 1
c         o
c         h1  o  0.96
c         h2  o  0.96      h1 104.5
c
c         this specifies the geometry for the water molecule.  the
c        bond distance is 0.96 angstroms and the angle is 104.5 degrees.
c
c
c         $geom
c         1 2
c         o(basis=6-31g)
c         h1(basis=sto-3g)   1  roh
c         h2 (basis=sto-3g)  1  roh      h1 theta
c
c         roh=0.96      constant
c         theta=104.5
c
c
c
c         this specifies a geometry for the water cation.  the hydrogens
c     are 0.96 angstroms from the oxygen and the distance will not be
c     varied in the optimization. the bond angle, theta, is initially
c     104.5 degrees, and will be optimized.  a double zeta basis (6-31g) will
c     be used on the oxygen, the hydrogens will use a minimal (sto-3g)
c     basis.
c
c
c
c
c
c         $geom
c         0 1
c         o    0.0       0.0    0.0
c         h    0.0       0.0   -0.96
c         h    0.0       0.0    0.96
c
c
c         this specifies yet another geometry for the water molecule.
c     the hydrogens are 0.96 angstroms from the oxygen along the
c     z-axis, and the bond angle is 180.0 degrees.
c
c
c
c
c
c         this section specifies a geometry suitable for optimizing the
c     hydrogen cyanide molecule.  if force or opt runs are  carried
c     out, the program is unable to handle the bond angles of 180 degrees
c     also be encountered in nearly linear situations such as ethynyl
c     groups in unsymmetrical molecules.  these situations can be
c     avoided by introducing dummy atoms along the angle bisector and
c     using the half-angle as the variable or constant.  this is
c     illustrated here and in the following example.
c
c
c
c
c
c          $geom
c          0 1
c          n
c          c  n rcn
c          x  c 1.0        n  half
c          o  c rco        x  half           n 180.0
c          h  o roh        c  coh            x  0.0
c
c          rcn=1.20
c          rco=1.3
c          roh=1.0
c          coh=105.0
c
c
c          in this section, half represents half of the nco angle which
c     is expected to be nearly linear.  note that a value of half less
c     the 90 degrees corresponds to a cis arrangement.
c
c
c
c
c***references
c
c***routines called
c
c***end prologue       $geom
      return
      end
