!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! MODULES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 MODULE gen_data
!
!
   IMPLICIT NONE
!
!
   REAL*8, SAVE :: zero=0.e0, one=1.e0, two=2.e0
   REAL*8, SAVE :: convfs=2.419e-02, convh=27.2114e0
   REAL*8, SAVE :: bohr=0.529177e0, conv1=2.19473e+5
   REAL*8, SAVE :: pi=3.141592653589793
!
!
   REAL*8, DIMENSION(:), ALLOCATABLE, SAVE :: xn,yn,zn,wtnx,wtny,wtnz
   REAL*8, SAVE :: delrx,delry,delrz
   REAL*8, SAVE :: amp_opt,width_opt,t_int,t_fin,f_time
   INTEGER, SAVE :: n_shift

END MODULE gen_data
