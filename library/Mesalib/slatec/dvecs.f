*deck dvecs
      subroutine dvecs (ncomp, lnfc, yhp, work, iwork, inhomo, iflag)
c***begin prologue  dvecs
c***subsidiary
c***purpose  subsidiary to dbvsup
c***library   slatec
c***type      double precision (svecs-s, dvecs-d)
c***author  watts, h. a., (snla)
c***description
c
c  this subroutine is used for the special structure of complex*16
c  valued problems. dmgsbv is called upon to obtain lnfc vectors from an
c  original set of 2*lnfc independent vectors so that the resulting
c  lnfc vectors together with their imaginary product or mate vectors
c  form an independent set.
c
c***see also  dbvsup
c***routines called  dmgsbv
c***common blocks    dml18j
c***revision history  (yymmdd)
c   750601  date written
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891009  removed unreferenced statement label.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dvecs
c
      integer icoco, idp, iflag, indpvt, inhomo, integ, iwork(*), k,
     1     kp, lnfc, lnfcc, mxnon, ncomp, ndisk, neq, neqivp, nic, niv,
     2     nopg, nps, ntape, ntp, numort, nxpts
      double precision ae, dum, re, tol, work(*), yhp(ncomp,*)
      common /dml18j/ ae,re,tol,nxpts,nic,nopg,mxnon,ndisk,ntape,neq,
     1                indpvt,integ,nps,ntp,neqivp,numort,lnfcc,
     2                icoco
c***first executable statement  dvecs
         if (lnfc .ne. 1) go to 20
            do 10 k = 1, ncomp
               yhp(k,lnfc+1) = yhp(k,lnfcc+1)
   10       continue
            iflag = 1
         go to 60
   20    continue
            niv = lnfc
            lnfc = 2*lnfc
            lnfcc = 2*lnfcc
            kp = lnfc + 2 + lnfcc
            idp = indpvt
            indpvt = 0
            call dmgsbv(ncomp,lnfc,yhp,ncomp,niv,iflag,work(1),work(kp),
     1                  iwork(1),inhomo,yhp(1,lnfc+1),work(lnfc+2),dum)
            lnfc = lnfc/2
            lnfcc = lnfcc/2
            indpvt = idp
            if (iflag .ne. 0 .or. niv .ne. lnfc) go to 40
               do 30 k = 1, ncomp
                  yhp(k,lnfc+1) = yhp(k,lnfcc+1)
   30          continue
               iflag = 1
            go to 50
   40       continue
               iflag = 99
   50       continue
   60    continue
      continue
      return
      end
