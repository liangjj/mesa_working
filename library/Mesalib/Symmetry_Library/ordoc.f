*deck @(#)ordoc.f	5.1  11/6/94
      subroutine ordoc
      implicit real*8(a-h,o-z)
c
c1orient
c
c                conventions for orienting molecules in cartesian space
c
c
c
c  a.  goals
c
c           1.  the molecule is oriented  so  as  to  simplify  the  3x3
c               transformation matrices.
c
c           2.  two z-matrices differing  only  in  the  values  of  the
c               internal  coordinates but which are identical as regards
c               the integer quantites  (such  as  occurs  on  subsequent
c               points  of  a  geometry  optimization) shall produce the
c               same orientation of the molecule.
c
c           3.  two different z-matrices for  the  same  molecule  shall
c               produce  the  same  coordinates,  save  for  a  possible
c               renumbering of the atoms.
c
c           4.  the number of molecular  orbital  coefficients  zero  by
c               symmetry should be maximized.
c
c  b.  general considerations
c
c           1.  a  right  handed  coordinate   system   will   be   used
c               throughout.
c
c           2.  the molecule will be translated so that  its  center  of
c               charge is at the origin.
c
c           3.  atoms will not be reordered relative to their order upon
c               input.
c
c           4.  the cartesian axes will be  considered  to  increase  in
c               priority in the order x < y < z.
c
c  c.  rules for the positioning of an axis
c
c           1.  an axis of rotation or a principal axis of charge may be
c               aligned with a cartesian axis in one of two ways.  which
c               way  will  be  decided  by  successively  applying   the
c               following tests until a definitive result is achieved:
c               1.  the third moment of charge should be positive.
c               2.  the sum of the projections of the atomic coordinates
c                   onto the reference axis should be positive.
c               3.  the first atom with a  non-zero  projection  on  the
c                   reference  axis should have a positive projection on
c                   that axis.
c
c           2.  if a rotation is neccessary in order to meet one of  the
c               above  criteria  it shall be a 180 degree rotation about
c               the axis defined below:
c                          reference         axis
c                            axis         of rotation
c
c                             x                y
c                             y                z
c                             z                x
c
c  d.  rules for  postioning  the  principal  axes  of  charge:  in  the
c      abscence  of any other rules, the principal axis corresponding to
c      the largest principal moment of charge will be aligned  with  the
c      highest priority cartesian axis available.
c
c  e.  special rules for asymmetric top molecules
c
c           1.  cs:  the molecular plane shall be made  coincident  with
c               the  xy plane.  note that this convention conflicts with
c               mulliken's suggestion, but that it  is  consistent  with
c               the  character  tables  of  cotton and of herzberg.  the
c               moleule is then rotated about the z-axis  in  accordance
c               with the rules given below for cn molecules.
c
c           2.  c2v:  following mulliken's recommendation for planar c2v
c               molecules,  the  molecular  plane  is  placed  in the yz
c               plane.  for non-planar molecules the following tests are
c               successively applied:
c               1.  the mirror plane with the most atoms is put  in  the
c                   yz plane.
c               2.  the mirror plane with the most non-hydrogen atoms is
c                   put in the yz plane.
c               3.  the mirror plane with the lowest  numbered  atom  is
c                   made coincident with yz.
c               4.  apply axes of charge rules as per section d.
c
c           3.  planar     d2h     molecules:  following      mulliken's
c               recommendation,  the molecular plane is placed in the yz
c               plane.  the molecule is rotated about the x axis so that
c               the  z  axis  shall  pass  through the greater number of
c               atoms, or, if this is not decisive, the  greater  number
c               of bonds.
c
c           4.  c2 molecules follow the general rules given below for cn
c               molecules.
c
c           5.  c1 molecules are not re-oriented.
c
c           6.  ci molecules are translated but not rotated.
c
c  f.  special rules for symmetric top molecules
c
c           1.  the unique axis is aligned with the z axis.
c
c           2.  definition:  a "circular-set" of atoms  is  composed  of
c               those  atoms lying in a plane which have the same atomic
c               number and which are equidistant from a  reference  axis
c               perpindicular to the plane.  atoms on the reference axis
c               are not included in any circular-sets.   a  circular-set
c               of atoms is generated by a proper rotation axis.
c
c           3.  definition:  the "key-atom" in a symmetric top  molecule
c               is  the  lowest  numbered  atom in the key circular-set.
c               the following tests are carried out successively to find
c               the key circular-set:
c               1.  which set is nearest the xy plane?
c               2.  which set has a positive projection on the z axis?
c               3.  which set is nearest to the z axis?
c               4.  which set is comprised  of  atoms  with  the  lowest
c                   atomic number?
c
c           4.  cn, cnh, sn:  the molecule is rotated about  the  z-axis
c               so  as  to  maximize  the number of pairs of heavy atoms
c               which  are  parallel  to  the  y  axis.   if   no   such
c               arrangement  is satisfactory then the key atom is placed
c               in  the  yz  plane  so  as  to  give  it  a  positive  y
c               coordinate.
c
c           5.  dn, dnh:  one of the c2 axes is made coincident with the
c               y  cartesian  axis.   the  tests  in 7 below are used to
c               decide which c2 axis is so positioned.
c
c           6.  dnd, cnv:  one of the vertical planes is made coincident
c               with  the  yz cartesian plane.  the tests in 7 below are
c               used to decide which plane is so positioned.
c
c           7.  test for orienting dn, dnh, dnd, and cnv molecules:
c               1.  maximize the projection of the key  atom  on  the  y
c                   axis/
c               2.  if two orientations give the maximum projection on y
c                   choose the one with the maximum projection on x.
c
c           8.  for molecules contained in the  xy  plane  the  standard
c               axis orientation rule (see section c) are applied to the
c               x axis to complete the  specification  of  a  molecule's
c               orientation.
c
c  g.  special rules for spherical top molecules
c
c           1.  definitions:  a "spherical-set" of atoms is composed  of
c               those  atoms  which  are equidistant from the origin and
c               which have the same atomic number.  spherical-sets shall
c               be  ordered  in  terms  of  increasing distance from the
c               origin and in terms of increasing atomic number  at  any
c               one  distance.   the  "key  atom" is the lowest numbered
c               atom in the first spherical-set.
c
c           2.  t, td, th
c               1.  align the three c2 axes with the cartesian  axes  so
c                   as to maximize the z-coordinate of the key atom.
c               2.  treat the molecule as a symmetric top.
c
c           3.  o, oh
c               1.  align the three c4 axes with the cartesian  axes  so
c                   as to maximize the z-coordinate of the key atom.
c               2.  treat the molecule as a symmetric top.
c
c           4.  i, ih
c               1.  align one of the c5 axes with the z axis  so  as  to
c                   maximize the z-coordinate of the key atom.
c               2.  treat the molecule as a symmetric top.
c
c
c  h.  references
c
c               cotton, f.  a.,  "the  chemical  applications  of  group
c               theory", 2nd ed., wiley-interscience, new york, 1971.
c
c               herzberg, g., "infrared and raman spectra of  polyatomic
c               molecules", d. van nostrand, princeton, 1945.
c
c               mulliken, r. s., "report on the notation for the spectra
c               of   polyatomic   molecules",   j.   chem.   phys.,  23,
c               1997(1955).
c
c?
c
c    lest the complier complain:
c
c     call rtrace(6hordoc ,1)
      return
      end
