*deck dstway
      subroutine dstway (u, v, yhp, inout, stowa)
c***begin prologue  dstway
c***subsidiary
c***purpose  subsidiary to dbvsup
c***library   slatec
c***type      double precision (stway-s, dstway-d)
c***author  watts, h. a., (snla)
c***description
c
c  this subroutine stores (recalls) integration data in the event
c  that a restart is needed (the homogeneous solution vectors become
c  too dependent to continue).
c
c***see also  dbvsup
c***routines called  dstor1
c***common blocks    dml15t, dml18j, dml8sz
c***revision history  (yymmdd)
c   750601  date written
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dstway
c
      integer icoco, igofx, indpvt, info, inhomo, inout, integ, istkop,
     1     ivp, j, k, knswot, ko, kop, ks, ksj, lotjp, mnswot, mxnon,
     2     ncomp, ndisk, neq, neqivp, nfc, nfcc, nic, nopg, nps, nswot,
     3     ntape, ntp, numort, nxpts
      double precision ae, c, pwcnd, px, re, stowa(*), tnd, tol, u(*),
     1     v(*), x, xbeg, xend, xop, xot, xsav, yhp(*)
c
      common /dml8sz/ c,xsav,igofx,inhomo,ivp,ncomp,nfc
      common /dml15t/ px,pwcnd,tnd,x,xbeg,xend,xot,xop,info(15),istkop,
     1                knswot,kop,lotjp,mnswot,nswot
      common /dml18j/ ae,re,tol,nxpts,nic,nopg,mxnon,ndisk,ntape,neq,
     1                indpvt,integ,nps,ntp,neqivp,numort,nfcc,
     2                icoco
c
c***first executable statement  dstway
      if (inout .eq. 1) go to 30
c
c        save in stowa array and istkop
c
         ks = nfc*ncomp
         call dstor1(stowa,u,stowa(ks+1),v,1,0,0)
         ks = ks + ncomp
         if (neqivp .lt. 1) go to 20
         do 10 j = 1, neqivp
            ksj = ks + j
            stowa(ksj) = yhp(ksj)
   10    continue
   20    continue
         ks = ks + neqivp
         stowa(ks+1) = x
         istkop = kop
         if (xop .eq. x) istkop = kop + 1
      go to 80
   30 continue
c
c        recall from stowa array and istkop
c
         ks = nfc*ncomp
         call dstor1(yhp,stowa,yhp(ks+1),stowa(ks+1),1,0,0)
         ks = ks + ncomp
         if (neqivp .lt. 1) go to 50
         do 40 j = 1, neqivp
            ksj = ks + j
            yhp(ksj) = stowa(ksj)
   40    continue
   50    continue
         ks = ks + neqivp
         x = stowa(ks+1)
         info(1) = 0
         ko = kop - istkop
         kop = istkop
         if (ndisk .eq. 0 .or. ko .eq. 0) go to 70
            do 60 k = 1, ko
               backspace ntape
   60       continue
   70    continue
   80 continue
      return
      end
