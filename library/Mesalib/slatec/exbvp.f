*deck exbvp
      subroutine exbvp (y, nrowy, xpts, a, nrowa, alpha, b, nrowb, beta,
     +   iflag, work, iwork)
c***begin prologue  exbvp
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (exbvp-s, dexbvp-d)
c***author  watts, h. a., (snla)
c***description
c
c  this subroutine is used to execute the basic technique for solving
c  the two-point boundary value problem
c
c***see also  bvsup
c***routines called  bvpor, xermsg
c***common blocks    ml15to, ml17bw, ml18jr, ml5mco, ml8sz
c***revision history  (yymmdd)
c   750601  date written
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c   910722  updated author section.  (als)
c***end prologue  exbvp
c
      dimension y(nrowy,*),a(nrowa,*),alpha(*),b(nrowb,*),beta(*),
     1         work(*),iwork(*),xpts(*)
      character*8 xern1, xern2
c
c     ****************************************************************
c
      common /ml8sz/ c,xsav,igofx,inhomo,ivp,ncomp,nfc
      common /ml18jr/ ae,re,tol,nxpts,nic,nopg,mxnon,ndisk,ntape,neq,
     1                indpvt,integ,nps,ntp,neqivp,numort,nfcc,
     2                icoco
      common /ml15to/ px,pwcnd,tnd,x,xbeg,xend,xot,xop,info(15),istkop,
     1                knswot,kop,lotjp,mnswot,nswot
      common /ml17bw/ kkkzpw,needw,neediw,k1,k2,k3,k4,k5,k6,k7,k8,k9,
     1                k10,k11,l1,l2,kkkint,lllint
c
      common /ml5mco/ uro,sru,eps,sqovfl,twou,fouru,lpar
c
c***first executable statement  exbvp
      kotc = 1
      iexp = 0
      if (iwork(7) .eq. -1) iexp = iwork(8)
c
c     compute orthonormalization tolerances.
c
   10 tol = 10.0**((-lpar-iexp)*2)
c
      iwork(8) = iexp
      mxnon = iwork(2)
c
c **********************************************************************
c **********************************************************************
c
      call bvpor(y,nrowy,ncomp,xpts,nxpts,a,nrowa,alpha,nic,b,
     1           nrowb,beta,nfc,iflag,work(1),mxnon,work(k1),ntp,
     2           iwork(18),work(k2),iwork(16),work(k3),work(k4),
     3           work(k5),work(k6),work(k7),work(k8),work(k9),
     4           work(k10),iwork(l1),nfcc)
c
c **********************************************************************
c **********************************************************************
c     if mgsbv returns with message of dependent vectors, we reduce
c     orthonormalization tolerance and try again. this is done
c     a maximum of 2 times.
c
      if (iflag .ne. 30) go to 20
      if (kotc .eq. 3  .or.  nopg .eq. 1) go to 30
      kotc = kotc + 1
      iexp = iexp - 2
      go to 10
c
c **********************************************************************
c     if bvpor returns message that the maximum number of
c     orthonormalizations has been attained and we cannot continue, then
c     we estimate the new storage requirements in order to solve problem
c
   20 if (iflag .ne. 13) go to 30
      xl = abs(xend-xbeg)
      zquit = abs(x-xbeg)
      inc = 1.5 * xl/zquit * (mxnon+1)
      if (ndisk .ne. 1) then
         nsafw = inc*kkkzpw + needw
         nsafiw = inc*nfcc + neediw
      else
         nsafw = needw + inc
         nsafiw = neediw
      endif
c
      write (xern1, '(i8)') nsafw
      write (xern2, '(i8)') nsafiw
      call xermsg ('slatec', 'exbvp',
     *   'in bvsup, predicted storage allocation for work array is ' //
     *   xern1 // ', predicted storage allocation for iwork array is '
     *   // xern2, 1, 0)
c
   30 iwork(1) = mxnon
      return
      end
