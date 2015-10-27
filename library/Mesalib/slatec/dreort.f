*deck dreort
      subroutine dreort (ncomp, y, yp, yhp, niv, w, s, p, ip, stowa,
     +   iflag)
c***begin prologue  dreort
c***subsidiary
c***purpose  subsidiary to dbvsup
c***library   slatec
c***type      double precision (reort-s, dreort-d)
c***author  watts, h. a., (snla)
c***description
c
c **********************************************************************
c   input
c *********
c     y, yp and yhp = homogeneous solution matrix and particular
c                     solution vector to be orthonormalized.
c     iflag = 1 --  store yhp into y and yp, test for
c                   reorthonormalization, orthonormalize if needed,
c                   save restart data.
c             2 --  store yhp into y and yp, reorthonormalization,
c                   no restarts.
c                   (preset orthonormalization mode)
c             3 --  store yhp into y and yp, reorthonormalization
c                   (when inhomo=3 and x=xend).
c **********************************************************************
c   output
c *********
c     y, yp = orthonormalized solutions.
c     niv = number of independent vectors returned from dmgsbv.
c     iflag = 0 --  reorthonormalization was performed.
c            10 --  solution process must be restarted at the last
c                   orthonormalization point.
c            30 --  solutions are linearly dependent, problem must
c                   be restarted from the beginning.
c     w, p, ip = orthonormalization information.
c **********************************************************************
c
c***see also  dbvsup
c***routines called  ddot, dmgsbv, dstor1, dstway
c***common blocks    dml15t, dml18j, dml8sz
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dreort
c
      double precision ddot
      integer icoco, iflag, igofx, ijk, indpvt, info, inhomo, integ,
     1     ip(*), istkop, ivp, j, k, kk, knswot, kop, l, lotjp, mflag,
     2     mnswot, mxnon, ncomp, ncompd, ndisk, neq, neqivp, nfc,
     3     nfcc, nfcp, nic, niv, nopg, nps, nswot, ntape, ntp, numort,
     4     nxpts
      double precision ae, c, dnd, dndt, dx, p(*), pwcnd, px, re, s(*),
     1     srp, stowa(*), tnd, tol, vnorm, w(*), wcnd, x, xbeg, xend,
     2     xop, xot, xsav, y(ncomp,*), yhp(ncomp,*), yp(*), ypnm
c
c     ******************************************************************
c
      common /dml8sz/ c,xsav,igofx,inhomo,ivp,ncompd,nfc
      common /dml15t/ px,pwcnd,tnd,x,xbeg,xend,xot,xop,info(15),istkop,
     1                knswot,kop,lotjp,mnswot,nswot
      common /dml18j/ ae,re,tol,nxpts,nic,nopg,mxnon,ndisk,ntape,neq,
     1                indpvt,integ,nps,ntp,neqivp,numort,nfcc,
     2                icoco
c
c **********************************************************************
c     begin block permitting ...exits to 210
c        begin block permitting ...exits to 10
c***first executable statement  dreort
            nfcp = nfc + 1
c
c           check to see if orthonormalization test is to be performed
c
c        ...exit
            if (iflag .ne. 1) go to 10
            knswot = knswot + 1
c        ...exit
            if (knswot .ge. nswot) go to 10
c     ......exit
            if ((xend - x)*(x - xot) .lt. 0.0d0) go to 210
   10    continue
         call dstor1(y,yhp,yp,yhp(1,nfcp),1,0,0)
c
c        ***************************************************************
c
c        orthogonalize the homogeneous solutions y
c        and particular solution yp.
c
         niv = nfc
         call dmgsbv(ncomp,nfc,y,ncomp,niv,mflag,s,p,ip,inhomo,yp,w,
     1               wcnd)
c
c           ************************************************************
c
c        check for linear dependence of the solutions.
c
         if (mflag .eq. 0) go to 50
c           begin block permitting ...exits to 40
               if (iflag .eq. 2) go to 30
                  if (nswot .le. 1 .and. lotjp .ne. 0) go to 20
c
c                    retrieve data for a restart at last
c                    orthonormalization point
c
                     call dstway(y,yp,yhp,1,stowa)
                     lotjp = 1
                     nswot = 1
                     knswot = 0
                     mnswot = mnswot/2
                     tnd = tnd + 1.0d0
                     iflag = 10
c           .........exit
                     go to 40
   20             continue
   30          continue
               iflag = 30
   40       continue
         go to 200
   50    continue
c           begin block permitting ...exits to 190
c              begin block permitting ...exits to 110
c
c                 ******************************************************
c
c              ...exit
                  if (iflag .ne. 1) go to 110
c
c                 test for orthonormalization
c
c              ...exit
                  if (wcnd .lt. 50.0d0*tol) go to 110
                  do 60 ijk = 1, nfcp
c              ......exit
                     if (s(ijk) .gt. 1.0d20) go to 110
   60             continue
c
c                 use linear extrapolation on logarithmic values of the
c                 norm decrements to determine next orthonormalization
c                 checkpoint.  other controls on the number of steps to
c                 the next checkpoint are added for safety purposes.
c
                  nswot = knswot
                  knswot = 0
                  lotjp = 0
                  wcnd = log10(wcnd)
                  if (wcnd .gt. tnd + 3.0d0) nswot = 2*nswot
                  if (wcnd .lt. pwcnd) go to 70
                     xot = xend
                     nswot = min(mnswot,nswot)
                     pwcnd = wcnd
                     px = x
                  go to 100
   70             continue
                     dx = x - px
                     dnd = pwcnd - wcnd
                     if (dnd .ge. 4) nswot = nswot/2
                     dndt = wcnd - tnd
                     if (abs(dx*dndt) .le. dnd*abs(xend-x)) go to 80
                        xot = xend
                        nswot = min(mnswot,nswot)
                        pwcnd = wcnd
                        px = x
                     go to 90
   80                continue
                        xot = x + dx*dndt/dnd
                        nswot = min(mnswot,nswot)
                        pwcnd = wcnd
                        px = x
   90                continue
  100             continue
c           ......exit
                  go to 190
  110          continue
c
c              *********************************************************
c
c              orthonormalization necessary so we normalize the
c              homogeneous solution vectors and change w accordingly.
c
               nswot = 1
               knswot = 0
               lotjp = 1
               kk = 1
               l = 1
               do 150 k = 1, nfcc
c                 begin block permitting ...exits to 140
                     srp = sqrt(p(kk))
                     if (inhomo .eq. 1) w(k) = srp*w(k)
                     vnorm = 1.0d0/srp
                     p(kk) = vnorm
                     kk = kk + nfcc + 1 - k
                     if (nfc .eq. nfcc) go to 120
c                 ......exit
                        if (l .ne. k/2) go to 140
  120                continue
                     do 130 j = 1, ncomp
                        y(j,l) = y(j,l)*vnorm
  130                continue
                     l = l + 1
  140             continue
  150          continue
c
               if (inhomo .ne. 1 .or. nps .eq. 1) go to 180
c
c                 normalize the particular solution
c
                  ypnm = ddot(ncomp,yp,1,yp,1)
                  if (ypnm .eq. 0.0d0) ypnm = 1.0d0
                  ypnm = sqrt(ypnm)
                  s(nfcp) = ypnm
                  do 160 j = 1, ncomp
                     yp(j) = yp(j)/ypnm
  160             continue
                  do 170 j = 1, nfcc
                     w(j) = c*w(j)
  170             continue
  180          continue
c
               if (iflag .eq. 1) call dstway(y,yp,yhp,0,stowa)
               iflag = 0
  190       continue
  200    continue
  210 continue
      return
      end
