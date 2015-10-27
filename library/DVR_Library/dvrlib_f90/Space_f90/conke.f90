!deck conke.f
!***begin prologue     conke
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            ke elements of DVR/FEM basis.
!***
!***references
!***routines called
!***end prologue       conke
!
  SUBROUTINE conke(kemat,norm,mat,tr,nglobal)
  USE dvr_global,   ONLY  : nreg, npt, iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nglobal
  INTEGER, DIMENSION(nreg)               :: tr
  REAL*8,  DIMENSION(*)                  :: mat
  REAL*8,  DIMENSION(nglobal,nglobal)    :: kemat
  REAL*8,  DIMENSION(nglobal)            :: norm
  CHARACTER (LEN=80)                     :: title
  LOGICAL                                :: bridge
  INTEGER                                :: row, reg, kei, start, end
  INTEGER                                :: keii, keij, kej, kejj
  INTEGER                                :: nfun, last, i, j
!
  WRITE(iout,1) nglobal
  kemat=0.d0
!
!     calculate the needed matrix elements
!
  row=1
  DO  reg=1,nreg
!
!        locate the matrix elements for this region and the previous region
!  
      kei=tr(reg)
      bridge=.true.
      start=2
      END=npt(reg)
!  
!     since all regions except the first and last have a bridge
!     function, offset the row and column indices so they start
!     at the second function.  the first function is included as
!     the bridge function from the previous interval.
!  
      keii=kei + END + 1
      keij=kei + 1
      IF(reg /= nreg) THEN
        kej=tr(reg+1)
        kejj=kej
      END IF
      IF(reg == 1) THEN
         start=1
         keii=kei
      END IF
      IF(reg == nreg) THEN
         bridge=.false.
      END IF
      nfun=END-start+1
!
!     get the global starting and ending value
!
      last=row+nfun-1
      WRITE(iout,2) reg, nfun, start, END, row, last, bridge
!
!     if we consider only the lower triangle of the full matrix
!     and recognize that use of the bloch operator ensures
!     hermiticity, then a function in region i connects with
!     functions in the i - 1 region and the i region.  in the
!     i - 1 region it can only connect to a bridge function, if
!     it exists which has a piece in the i region.
!     if the function in region i is a bridge function, it connects
!     to all lower functions in that region and with the first function
!     in the i + 1 region.
!
      IF(reg > 1) THEN
!  
!        do integrals between functions in region i and i - 1
!
!        ok we have a bridge function in region i - 1
!        it is a linear combination of F(n|i-1) + F(1|i).
!        we only need integrals from region i, suitably
!        weighted by the norms.
!  
         CALL makod_1(kemat(row,row-1),mat(keij),norm(row-1),npt(reg),nfun,prn(7))
      END IF
!
!        do integrals between functions all in region i
!        we already have the normalizations
!        if we have a bridge function in region i
!        it is a linear combination of F(n|i) + F(1|i+1).
!        we only need integrals from region i and i + 1, suitably
!        weighted by the norms.
!
         CALL makd_1(kemat(row,row),mat(keii),mat(kejj),norm(row),  &
                   bridge,npt(reg),nfun,nglobal,prn(7))
         row=row+nfun
  END DO
  DO  i=1,nglobal
      DO  j=1,i
          kemat(j,i)=kemat(i,j)
      END DO
  END DO
  IF(prn(8)) THEN
      title='kinetic energy hamiltonian'
      CALL prntrm(title,kemat,nglobal,nglobal,nglobal,nglobal,iout)
   END IF
1    FORMAT(/,5X,'size of global basis set = ',i4)
2    FORMAT(/,5X,'region = ',i3,1X,'number of functions = ',i3,/,5X,  &
    'starting function = ',i3,1X,'ending function = ',i3,  &
    /,5X,'global starting function = ',i3,1X, 'global ending function = ',i3,  &
    /,5X,'bridge function = ',l1)
4    FORMAT(10X,i4,12X,i4,1X,a38)
END SUBROUTINE conke



