! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Fillfun}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle

!deck filfun.f
!***begin prologue     filfun
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            global function and derivatives on grid
!***description        fill up the global arrays from the regional arrays.
!***                   
!***references

!***routines called
!***end prologue       filfun

   SUBROUTINE filfun(q,wt,p,dp,ddp,norm,qi,wti,pi,dpi,  &
                  ddpi,renorm,bridge,ni,nfun,nglobal,start)
!
  USE input_output
!
  IMPLICIT NONE
  INTEGER                                   :: ni, nfun, nglobal, start
  REAL*8, DIMENSION(nglobal)                :: q, wt, norm
  REAL*8, DIMENSION(nglobal,nglobal)        :: p, dp, ddp
  REAL*8, DIMENSION(ni)                     :: qi, wti
  REAL*8, DIMENSION(ni,ni)                  :: pi, dpi, ddpi
  LOGICAL                                   :: bridge
  REAL*8                                    :: renorm, fac1, fac2, fac3
  INTEGER                                   :: cnti, i, cntj, j
!
  CHARACTER (LEN=80) :: title
!
!        get the normalizations and functions in
!        region i, including the bridge function if present.
!        if there is a bridge function, it contains a renormalization
!        factor involving the weight from regon i+1.
!

  cnti=start
  DO  i=1,nfun-1
      q(i)=qi(cnti)
      wt(i)=wti(cnti)
      norm(i) = SQRT ( 1.d0/ wt(i) )
      cntj=start
      DO  j=1,nfun
          p(j,i)=pi(cntj,cnti)*norm(i)
          dp(j,i)=dpi(cntj,cnti)*norm(i)
          ddp(j,i)=ddpi(cntj,cnti)*norm(i)
          cntj=cntj+1
      END DO
      cnti=cnti+1
  END DO
  q(nfun) = qi(cnti)
!
!     check if the last function in region i is a bridge function
!
   IF(bridge) THEN
!
!        yes it is.
!  
      wt(nfun) = wti(cnti) + renorm
      norm(nfun) = SQRT ( 1.d0 / wt(nfun) )
      cntj=start
      DO  j=1,nfun
          p(j,nfun)=pi(cntj,cnti)*norm(nfun)
          cntj=cntj+1
      END DO
   ELSE
!
!        no it is not.
!  
      wt(nfun)=wti(cnti)
      norm(nfun)= SQRT ( 1.d0/ wt(nfun) )
      cntj=start
      DO  j=1,nfun
          p(j,nfun)=pi(cntj,cnti)*norm(nfun)
          dp(j,nfun)=dpi(cntj,cnti)*norm(nfun)
          ddp(j,nfun)=ddpi(cntj,cnti)*norm(nfun)
          cntj=cntj+1
      END DO
   END IF
   RETURN
END SUBROUTINE filfun



