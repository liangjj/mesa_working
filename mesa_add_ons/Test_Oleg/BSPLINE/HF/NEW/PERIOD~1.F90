
!=======================================================================
  MODULE periodic_table 
!=======================================================================
!   This module lists the labels of the elements of the periodic table 
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    SAVE

    CHARACTER(LEN=2), DIMENSION(103) ::  atom_label = &
     (/' H', 'He', 'Li', 'Be', ' B', ' C', ' N', ' O', ' F', 'Ne', &
       'Na', 'Mg', 'Al', 'Si', ' P', ' S', 'Cl', 'Ar', ' K', 'Ca', &
       'Sc', 'Ti', ' V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
       'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', ' Y', 'Zr', &
       'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
       'Sb', 'Te', ' I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
       'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
       'Lu', 'Hf', 'Ta', ' W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
       'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
       'Pa', ' U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', &
       'Md', 'No', 'Lw', '  ', '  ', '  ', '  ', '  ', '  ', '  '/)

    CHARACTER(LEN=2), DIMENSION(17) :: shell_order = &
     (/'1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '5s', &
       '5p', '4f', '5d'  '6s', '6p', '5f', '7s', '7p'/
 
    INTEGER, DIMENSION(17) :: shell_number = &
     (/   2,    6,  10,   12,   18,   28,   30,   36,   46,   48,  &
         54,   68,  78,   80,   86,  100,  102,  108/)

    CHARACTER(LEN=*), DIMENSION(13) :: shells = &
     (/'1s(2)', '.. 2s(2)2p(6)', '.. 3s(2)3p(6)',                  &
       '.. 3d(10)', ' .. 4s(2)4p(6)', '4d(10)', '.. 5s(2)5p(6)',   &
       '.. 4f(14)', ' '.. 5d(10)', '.. 6s(2)6p(6)' &
       '.. 5f(14)', ' .. 7s(2)7p(6)'/)

    INTEGER, DIMENSION(13) :: number_of_electrons = &
     (/  2,   10,   18,   28,   36,   46,   54,    ` &
        68,   78,   86,  100,  108/)

  END MODULE periodic_table



