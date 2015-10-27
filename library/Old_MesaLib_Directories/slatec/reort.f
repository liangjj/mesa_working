*deck reort
      subroutine reort (ncomp, y, yp, yhp, niv, w, s, p, ip, stowa,
     +   iflag)
c***begin prologue  reort
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (reort-s, dreort-d)
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
c***see also  bvsup
c***routines called  mgsbv, sdot, stor1, stway
c***common blocks    ml15to, ml18jr, ml8sz
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  reort
c
      dimension y(ncomp,*),yp(*),w(*),s(*),p(*),ip(*),
     1          stowa(*),yhp(ncomp,*)
c
c **********************************************************************
c
      common /ml8sz/ c,xsav,igofx,inhomo,ivp,ncompd,nfc
      common /ml15to/ px,pwcnd,tnd,x,xbeg,xend,xot,xop,info(15),istkop,
     1                knswot,kop,lotjp,mnswot,nswot
      common /ml18jr/ ae,re,tol,nxpts,nic,nopg,mxnon,ndisk,ntape,neq,
     1                indpvt,integ,nps,ntp,neqivp,numort,nfcc,
     2                icoco
c
c **********************************************************************
c***first executable statement  reort
      nfcp=nfc+1
c
c     check to see if orthonormalization test is to be performed
c
      if (iflag .ne. 1) go to 5
      knswot=knswot+1
      if (knswot .ge. nswot) go to 5
      if ((xend-x)*(x-xot) .lt. 0.) return
    5 call stor1(y,yhp,yp,yhp(1,nfcp),1,0,0)
c
c     ****************************************
c
c     orthogonalize the homogeneous solutions y
c     and particular solution yp.
c
      niv=nfc
      call mgsbv(ncomp,nfc,y,ncomp,niv,mflag,s,p,ip,inhomo,yp,w,wcnd)
c
c     ****************************************
c
c  check for linear dependence of the solutions.
c
      if (mflag .eq. 0)  go to 25
      if (iflag .eq. 2) go to 15
      if (nswot .gt. 1  .or.  lotjp .eq. 0) go to 20
   15 iflag=30
      return
c
c     retrieve data for a restart at last orthonormalization point
c
   20 call stway(y,yp,yhp,1,stowa)
      lotjp=1
      nswot=1
      knswot=0
      mnswot=mnswot/2
      tnd=tnd+1.
      iflag=10
      return
c
c     ****************************************
c
   25 if (iflag .ne. 1) go to 60
c
c     test for orthonormalization
c
      if (wcnd .lt. 50.*tol) go to 60
      do 30 ijk=1,nfcp
      if (s(ijk) .gt. 1.0e+20) go to 60
   30 continue
c
c     use linear extrapolation on logarithmic values of the norm
c     decrements to determine next orthonormalization checkpoint.
c     other controls on the number of steps to the next checkpoint
c     are added for safety purposes.
c
      nswot=knswot
      knswot=0
      lotjp=0
      wcnd=log10(wcnd)
      if (wcnd .gt. tnd+3.) nswot=2*nswot
      if (wcnd .ge. pwcnd) go to 40
      dx=x-px
      dnd=pwcnd-wcnd
      if (dnd .ge. 4) nswot=nswot/2
      dndt=wcnd-tnd
      if (abs(dx*dndt) .gt. dnd*abs(xend-x)) go to 40
      xot=x+dx*dndt/dnd
      go to 50
   40 xot=xend
   50 nswot=min(mnswot,nswot)
      pwcnd=wcnd
      px=x
      return
c
c     ****************************************
c
c     orthonormalization necessary so we normalize the homogeneous
c     solution vectors and change w accordingly.
c
   60 nswot=1
      knswot=0
      lotjp=1
      kk = 1
      l=1
      do 70 k = 1,nfcc
      srp=sqrt(p(kk))
      if (inhomo .eq. 1) w(k)=srp*w(k)
      vnorm=1./srp
      p(kk)=vnorm
      kk = kk + nfcc + 1 - k
      if (nfc .eq. nfcc) go to 63
      if (l .ne. k/2) go to 70
   63 do 65 j = 1,ncomp
   65 y(j,l) = y(j,l)*vnorm
      l=l+1
   70 continue
c
      if (inhomo .ne. 1  .or.  nps .eq. 1)  go to 100
c
c     normalize the particular solution
c
      ypnm=sdot(ncomp,yp,1,yp,1)
      if (ypnm .eq. 0.0)  ypnm = 1.0
      ypnm = sqrt(ypnm)
      s(nfcp) = ypnm
      do 80 j = 1,ncomp
   80 yp(j) = yp(j) / ypnm
      do 90 j = 1,nfcc
   90 w(j) = c * w(j)
c
  100 if (iflag .eq. 1) call stway(y,yp,yhp,0,stowa)
      iflag=0
      return
      end
