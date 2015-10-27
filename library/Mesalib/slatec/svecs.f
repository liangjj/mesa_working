*deck svecs
      subroutine svecs (ncomp, lnfc, yhp, work, iwork, inhomo, iflag)
c***begin prologue  svecs
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (svecs-s, dvecs-d)
c***author  watts, h. a., (snla)
c***description
c
c  this subroutine is used for the special structure of complex valued
c  problems. mgsbv is called upon to obtain lnfc vectors from an
c  original set of 2*lnfc independent vectors so that the resulting
c  lnfc vectors together with their imaginary product or mate vectors
c  form an independent set.
c
c***see also  bvsup
c***routines called  mgsbv
c***common blocks    ml18jr
c***revision history  (yymmdd)
c   750601  date written
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  svecs
c
      dimension yhp(ncomp,*),work(*),iwork(*)
      common /ml18jr/ ae,re,tol,nxpts,nic,nopg,mxnon,ndisk,ntape,neq,
     1                indpvt,integ,nps,ntp,neqivp,numort,lnfcc,
     2                icoco
c***first executable statement  svecs
      if (lnfc .eq. 1) go to 5
      niv=lnfc
      lnfc=2*lnfc
      lnfcc=2*lnfcc
      kp=lnfc+2+lnfcc
      idp=indpvt
      indpvt=0
      call mgsbv(ncomp,lnfc,yhp,ncomp,niv,iflag,work(1),work(kp),
     1         iwork(1),inhomo,yhp(1,lnfc+1),work(lnfc+2),dum)
      lnfc=lnfc/2
      lnfcc=lnfcc/2
      indpvt=idp
      if (iflag .eq. 0  .and.  niv .eq. lnfc) go to 5
      iflag=99
      return
    5 do 6 k=1,ncomp
    6 yhp(k,lnfc+1)=yhp(k,lnfcc+1)
      iflag=1
      return
      end
