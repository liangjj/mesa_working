*deck sdanrm
      real function sdanrm (neq, v, wt, rpar, ipar)
c***begin prologue  sdanrm
c***subsidiary
c***purpose  compute vector norm for sdassl.
c***library   slatec (dassl)
c***type      single precision (sdanrm-s, ddanrm-d)
c***author  petzold, linda r., (llnl)
c***description
c-----------------------------------------------------------------------
c     this function routine computes the weighted
c     root-mean-square norm of the vector of length
c     neq contained in the array v,with weights
c     contained in the array wt of length neq.
c        sdanrm=sqrt((1/neq)*sum(v(i)/wt(i))**2)
c-----------------------------------------------------------------------
c***routines called  (none)
c***revision history  (yymmdd)
c   830315  date written
c   901009  finished conversion to slatec 4.0 format (f.n.fritsch)
c   901019  merged changes made by c. ulrich with slatec 4.0 format.
c   901026  added explicit declarations for all variables and minor
c           cosmetic changes to prologue.  (fnf)
c***end prologue  sdanrm
c
      integer  neq, ipar(*)
      real  v(neq), wt(neq), rpar(*)
c
      integer  i
      real  sum, vmax
c
c***first executable statement  sdanrm
      sdanrm = 0.0e0
      vmax = 0.0e0
      do 10 i = 1,neq
        if(abs(v(i)/wt(i)) .gt. vmax) vmax = abs(v(i)/wt(i))
10      continue
      if(vmax .le. 0.0e0) go to 30
      sum = 0.0e0
      do 20 i = 1,neq
20      sum = sum + ((v(i)/wt(i))/vmax)**2
      sdanrm = vmax*sqrt(sum/neq)
30    continue
      return
c------end of function sdanrm------
      end
