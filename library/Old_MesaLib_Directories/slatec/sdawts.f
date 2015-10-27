*deck sdawts
      subroutine sdawts (neq, iwt, rtol, atol, y, wt, rpar, ipar)
c***begin prologue  sdawts
c***subsidiary
c***purpose  set error weight vector for sdassl.
c***library   slatec (dassl)
c***type      single precision (sdawts-s, ddawts-d)
c***author  petzold, linda r., (llnl)
c***description
c-----------------------------------------------------------------------
c     this subroutine sets the error weight vector
c     wt according to wt(i)=rtol(i)*abs(y(i))+atol(i),
c     i=1,-,n.
c     rtol and atol are scalars if iwt = 0,
c     and vectors if iwt = 1.
c-----------------------------------------------------------------------
c***routines called  (none)
c***revision history  (yymmdd)
c   830315  date written
c   901009  finished conversion to slatec 4.0 format (f.n.fritsch)
c   901019  merged changes made by c. ulrich with slatec 4.0 format.
c   901026  added explicit declarations for all variables and minor
c           cosmetic changes to prologue.  (fnf)
c***end prologue  sdawts
c
      integer  neq, iwt, ipar(*)
      real  rtol(*), atol(*), y(*), wt(*), rpar(*)
c
      integer  i
      real  atoli, rtoli
c
c***first executable statement  sdawts
      rtoli=rtol(1)
      atoli=atol(1)
      do 20 i=1,neq
         if (iwt .eq.0) go to 10
           rtoli=rtol(i)
           atoli=atol(i)
10         wt(i)=rtoli*abs(y(i))+atoli
20         continue
      return
c-----------end of subroutine sdawts------------------------------------
      end
