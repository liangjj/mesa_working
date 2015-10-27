! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Final DVR Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck conham.f
!***begin prologue     conham
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            matrix elements of DVR/FEM basis.
!***
!***references
!***routines called
!***end prologue       conham
!
  SUBROUTINE conham(hmat,vmat,kmat,pmat,norm,mat,tr,pr,nglobal)
  USE dvr_global,   ONLY  : nreg, npt, iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nglobal
  INTEGER, DIMENSION(nreg)               :: tr, pr
  REAL*8, DIMENSION(*)                   :: mat
  REAL*8, DIMENSION(nglobal,nglobal)     :: hmat, kmat, pmat
  REAL*8, DIMENSION(nglobal)             :: vmat, norm
  CHARACTER (LEN=80)                     :: title
  LOGICAL                                :: bridge
  INTEGER                                :: i, j, row, reg, start, end
  INTEGER                                :: kei, kej, keii, keij, kejj
  INTEGER                                :: pri, prj, prii, prij, prjj
  INTEGER                                :: nfun, last
  hmat = 0.d0
  kmat = 0.d0
  pmat = 0.d0
!
!     calculate the needed matrix elements
!
  row=1
  DO  reg=1,nreg
!
!        locate the matrix elements for this region and the previous region
!  
      kei=tr(reg)
      pri=pr(reg)
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
      prii=pri + END + 1
      keij=kei + 1
      prij=pri + 1
      IF(reg /= nreg) THEN
        kej=tr(reg+1)
        prj=pr(reg+1)
        kejj=kej
        prjj=prj
      END IF
      IF(reg == 1) THEN
         start=1
         keii=kei
         prii=pri
      END IF
      IF(reg == nreg) THEN
         bridge=.false.
      END IF
      nfun=END-start+1
!
!     get the global starting and ending value
!
      last=row+nfun-1
      IF(prn(12)) then
         WRITE(iout,2) reg, nfun, start, END, row, last, bridge
      END IF
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
         CALL makod(kmat(row,row-1),mat(keij),pmat(row,row-1),mat(prij), &
                    norm(row-1),npt(reg),nfun,prn(7))
      END IF
!
!        do integrals between functions all in region i
!        we already have the normalizations
!        if we have a bridge function in region i
!        it is a linear combination of F(n|i) + F(1|i+1).
!        we only need integrals from region i and i + 1, suitably
!        weighted by the norms.
!
         CALL makd(kmat(row,row),mat(keii),mat(kejj),pmat(row,row),mat(prii), &
                   mat(prjj),norm(row),bridge,npt(reg),nfun,nglobal,prn(7))
         row=row+nfun
  END DO
  DO  i=1,nglobal
      DO  j=1,i
          kmat(j,i) = kmat(i,j)
          pmat(j,i) = pmat(i,j)
      END DO
  END DO
  hmat = kmat
!
  IF(prn(8)) THEN
     title='full kinetic energy'
     CALL prntrm(title,kmat,nglobal,nglobal,nglobal,nglobal,iout)
  END IF  
!
!     add the potential
!
   DO  i=1,nglobal
       hmat(i,i) = hmat(i,i) + vmat(i)
   END DO
   IF(prn(8)) THEN
      title='full hamiltonian'
      CALL prntrm(title,hmat,nglobal,nglobal,nglobal,nglobal,iout)
      title='full first derivative matrix'
      CALL prntrm(title,pmat,nglobal,nglobal,nglobal,nglobal,iout)
   END IF
1    FORMAT(/,5X,'size of global basis set = ',i5)
2    FORMAT(/,5X,'region = ',i3,1X,'number of functions = ',i5,/,5X,  &
    'starting function = ',i5,1X,'ending function = ',i5,  &
    /,5X,'global starting function = ',i5,1X, 'global ending function = ',i5,  &
    /,5X,'bridge function = ',l1)
4    FORMAT(10X,i4,12X,i4,1X,a38)
END SUBROUTINE conham



