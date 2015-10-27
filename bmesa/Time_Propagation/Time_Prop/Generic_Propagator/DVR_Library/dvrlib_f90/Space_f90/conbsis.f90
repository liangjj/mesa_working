! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Construct Final DVR Basis}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck conbsis.f
!***begin prologue     conbsis
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            DVR/FEM basis.
!***
!***references

!***routines called
!***end prologue       conbsis

  SUBROUTINE conbsis(x,xwt,f,df,ddf,norm,grid,q,wt,p,dp,ddp,nglobal)
!
  USE dvr_global,    ONLY   : nreg, npt, iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                            :: nglobal
  INTEGER, DIMENSION(nreg)           :: q, wt, p, dp, ddp                
  REAL*8, DIMENSION(nglobal)         :: x, xwt, norm
  REAL*8, DIMENSION(nglobal,nglobal) :: f, df, ddf
  REAL*8, DIMENSION(*)               :: grid
  REAL*8                             :: renorm
  CHARACTER (LEN=80)                 :: title
  LOGICAL                            :: bridge
  INTEGER                            :: reg, start, end, nfun, last, row 
!
  f=0.d0
  df=0.d0
  ddf=0.d0
!     calculate the needed functions
  row=1
  DO  reg=1,nreg
      IF(reg /= nreg) then
         renorm= grid(wt(reg+1))
      ELSE
         renorm=0.d0
      END IF
      bridge=.true.
      start=2
      END=npt(reg)
      IF(reg == 1) THEN
         start=1
      END IF
      IF(reg == nreg) THEN
         bridge=.false.
      END IF
      nfun = END - start + 1
      last = row + nfun - 1
      IF (prn(12)) then
          WRITE(iout,1) reg, nfun, start, END, &
                        row, last, bridge
      END IF


!        calculate the functions and derivatives in this region
!        region i, including the bridge function if present.
!        if there is a bridge function, it contains pieces
!        from region i and region i + 1.


      CALL filfun(x(row),xwt(row),f(row,row),   &
                  df(row,row),ddf(row,row),     &
                  norm(row),                    &
                  grid(q(reg)),grid(wt(reg)),   &
                  grid(p(reg)),grid(dp(reg)),   &
                  grid(ddp(reg)),renorm,bridge, &
                  npt(reg),nfun,nglobal,start)
      row = row + nfun
  END DO
  IF(prn(1)) THEN
     title='global coordinates'
     CALL prntrm(title,x,nglobal,1,nglobal,1,iout)
     title='global weights'
     CALL prntrm(title,xwt,nglobal,1,nglobal,1,iout)
     title='global norm'
     CALL prntrm(title,norm,nglobal,1,nglobal,1,iout)
  END IF
  IF(prn(2)) THEN
     title='global basis'
     CALL prntrm(title,f,nglobal,nglobal,nglobal,nglobal,iout)
     title='1 derivative of global basis'
     CALL prntrm(title,df,nglobal,nglobal,nglobal,nglobal,iout)
     title='2 derivative of global basis'
    CALL prntrm(title,ddf,nglobal,nglobal,nglobal,nglobal,iout)
  END IF
1    FORMAT(/,5X,'region = ',i5,1X,'number of functions = ',i5,/,5X,  &
    'starting function = ',i5,1X,'ending function = ',i5,  &
    /,5X,'global starting function = ',i5,1X, 'global ending function = ',i5,  &
    /,5X,'bridge function = ',l1)
END SUBROUTINE conbsis




