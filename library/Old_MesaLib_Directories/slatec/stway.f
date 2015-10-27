*deck stway
      subroutine stway (u, v, yhp, inout, stowa)
c***begin prologue  stway
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (stway-s, dstway-d)
c***author  watts, h. a., (snla)
c***description
c
c  this subroutine stores (recalls) integration data in the event
c  that a restart is needed (the homogeneous solution vectors become
c  too dependent to continue)
c
c***see also  bvsup
c***routines called  stor1
c***common blocks    ml15to, ml18jr, ml8sz
c***revision history  (yymmdd)
c   750601  date written
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  stway
c
      dimension u(*),v(*),yhp(*),stowa(*)
c
      common /ml8sz/ c,xsav,igofx,inhomo,ivp,ncomp,nfc
      common /ml15to/ px,pwcnd,tnd,x,xbeg,xend,xot,xop,info(15),istkop,
     1                knswot,kop,lotjp,mnswot,nswot
      common /ml18jr/ ae,re,tol,nxpts,nic,nopg,mxnon,ndisk,ntape,neq,
     1                indpvt,integ,nps,ntp,neqivp,numort,nfcc,
     2                icoco
c
c***first executable statement  stway
      if (inout .eq. 1) go to 100
c
c     save in stowa array and istkop
c
      ks=nfc*ncomp
      call stor1(stowa,u,stowa(ks+1),v,1,0,0)
      ks=ks+ncomp
      if (neqivp .eq. 0) go to 50
      do 25 j=1,neqivp
      ksj=ks+j
   25 stowa(ksj)=yhp(ksj)
   50 ks=ks+neqivp
      stowa(ks+1)=x
      istkop=kop
      if (xop .eq. x) istkop=kop+1
      return
c
c     recall from stowa array and istkop
c
  100 ks=nfc*ncomp
      call stor1(yhp,stowa,yhp(ks+1),stowa(ks+1),1,0,0)
      ks=ks+ncomp
      if (neqivp .eq. 0) go to 150
      do 125 j=1,neqivp
      ksj=ks+j
  125 yhp(ksj)=stowa(ksj)
  150 ks=ks+neqivp
      x=stowa(ks+1)
      info(1)=0
      ko=kop-istkop
      kop=istkop
      if (ndisk .eq. 0  .or.  ko .eq. 0) return
      do 175 k=1,ko
  175 backspace ntape
      return
      end
