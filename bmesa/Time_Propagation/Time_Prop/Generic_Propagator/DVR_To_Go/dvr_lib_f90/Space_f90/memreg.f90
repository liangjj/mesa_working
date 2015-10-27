!*memreg
  SUBROUTINE MEMREG(q,wt,p,dp,ddp,ov,t,p_mom,v,words)
  USE dvr_global, ONLY  : nreg, npt  
  IMPLICIT NONE
  INTEGER, DIMENSION(nreg) :: q, wt, p, dp, ddp, ov, t, p_mom, v
  INTEGER, dimension(2)   :: words
  INTEGER                 :: i
!
  words(1)=0
  words(2)=0
  DO  i=1,nreg
      q(i) = words(1) + 1
      wt(i) = q(i)+npt(i)
      p(i) = wt(i)+npt(i)
      dp(i) = p(i)+npt(i)*npt(i)
      ddp(i) = dp(i)+npt(i)*npt(i)
      words(1)   = ddp(i) + npt(i)*npt(i)
      ov(i) = words(2) + 1
      t(i) = ov(i) + npt(i)*npt(i) 
      p_mom(i) = t(i) + npt(i)*npt(i)
      v(i) = p_mom(i) + npt(i)*npt(i)
      words(2)  = v(i) + npt(i)
  END DO
END SUBROUTINE memreg
