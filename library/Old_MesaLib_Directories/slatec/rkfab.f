*deck rkfab
      subroutine rkfab (ncomp, xpts, nxpts, nfc, iflag, z, mxnon, p,
     +   ntp, ip, yhp, niv, u, v, w, s, stowa, g, work, iwork, nfcc)
c***begin prologue  rkfab
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (rkfab-s, drkfab-d)
c***author  watts, h. a., (snla)
c***description
c
c **********************************************************************
c
c     subroutine rkfab integrates the initial value equations using
c     the variable-step runge-kutta-fehlberg integration scheme or
c     the variable-order adams method and orthonormalization
c     determined by a linear dependence test.
c
c **********************************************************************
c
c***see also  bvsup
c***routines called  bvder, deabm, derkf, reort, stor1
c***common blocks    ml15to, ml17bw, ml18jr, ml8sz
c***revision history  (yymmdd)
c   750601  date written
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  rkfab
c
      dimension p(ntp,*),ip(nfcc,*),u(ncomp,nfc,*),
     1          v(ncomp,*),w(nfcc,*),z(*),yhp(ncomp,*),
     2          xpts(*),s(*),stowa(*),work(*),iwork(*),
     3          g(*)
c
c **********************************************************************
c
      common /ml8sz/ c,xsav,igofx,inhomo,ivp,ncompd,nfcd
      common /ml15to/ px,pwcnd,tnd,x,xbeg,xend,xot,xop,info(15),istkop,
     1                knswot,kop,lotjp,mnswot,nswot
      common /ml18jr/ ae,re,tol,nxptsd,nic,nopg,mxnond,ndisk,ntape,neq,
     1                indpvt,integ,nps,ntpd,neqivp,numort,nfccd,
     2                icoco
      common /ml17bw/ kkkzpw,needw,neediw,k1,k2,k3,k4,k5,k6,k7,k8,k9,
     1                k10,k11,l1,l2,kkkint,lllint
c
      external bvder
c
c **********************************************************************
c  initialization of counters and variables.
c
c***first executable statement  rkfab
      kod = 1
      non = 1
      x = xbeg
      jon = 1
      info(1) = 0
      info(2) = 0
      info(3) = 1
      info(4) = 1
      work(1) = xend
      if (nopg .eq. 0)  go to 1
      info(3) = 0
      if (x .eq. z(1))  jon = 2
    1 nfcp1 = nfc + 1
c
c **********************************************************************
c *****beginning of integration loop at output points.******************
c **********************************************************************
c
      do 110 kopp = 2,nxpts
      kop=kopp
c
    5 xop = xpts(kop)
      if (ndisk .eq. 0)  kod = kop
c
c     step by step integration loop between output points.
c
   10 xxop = xop
      if (nopg .eq. 0)   go to 15
      if (xend.gt.xbeg.and.xop.gt.z(jon)) xxop=z(jon)
      if (xend.lt.xbeg.and.xop.lt.z(jon)) xxop=z(jon)
c
c **********************************************************************
   15 go to (20,25),integ
c     derkf integrator
c
   20 call derkf(bvder,neq,x,yhp,xxop,info,re,ae,idid,work,kkkint,
     1           iwork,lllint,g,ipar)
      go to 28
c     deabm integrator
c
   25 call deabm(bvder,neq,x,yhp,xxop,info,re,ae,idid,work,kkkint,
     1           iwork,lllint,g,ipar)
   28 if(idid .ge. 1) go to 30
      info(1) = 1
      if(idid .eq. -1) go to 15
      iflag = 20 - idid
      return
c
c **********************************************************************
c     gram-schmidt orthogonalization test for orthonormalization
c     (temporarily using u and v in the test)
c
   30 if (nopg .eq. 0)  go to 35
      if (xxop .ne. z(jon))  go to 100
      jflag=2
      go to 40
   35 jflag=1
      if (inhomo .eq. 3  .and.  x .eq. xend) jflag=3
c
   40 if (ndisk .eq. 0) non=numort+1
      call reort(ncomp,u(1,1,kod),v(1,kod),yhp,niv,
     1           w(1,non),s,p(1,non),ip(1,non),stowa,jflag)
c
      if (jflag .ne. 30) go to 45
      iflag=30
      return
c
   45 if (jflag .eq. 10) go to 5
c
      if (jflag .ne. 0)  go to 100
c
c **********************************************************************
c     store orthonormalized vectors into solution vectors.
c
      if (numort .lt. mxnon)  go to 65
      if (x .eq. xend) go to 65
      iflag = 13
      return
c
   65 numort = numort + 1
      call stor1(yhp,u(1,1,kod),yhp(1,nfcp1),v(1,kod),1,
     1           ndisk,ntape)
c
c **********************************************************************
c     store orthonormalization information, initialize
c     integration flag, and continue integration to the next
c     orthonormalization point or output point.
c
      z(numort) = x
      if (inhomo .eq. 1  .and.  nps .eq. 0)  c = s(nfcp1) * c
      if (ndisk .eq. 0)  go to 90
      if (inhomo .eq. 1)  write (ntape) (w(j,1), j = 1,nfcc)
      write(ntape) (ip(j,1), j = 1,nfcc),(p(j,1), j = 1,ntp)
   90 info(1) = 0
      jon = jon + 1
      if (nopg .eq. 1  .and.  x .ne. xop)  go to 10
c
c **********************************************************************
c     continue integration if we are not at an output point.
c
  100 if (idid .eq. 1)  go to 15
c
c     storage of homogeneous solutions in u and the particular
c     solution in v at the output points.
c
      call stor1(u(1,1,kod),yhp,v(1,kod),yhp(1,nfcp1),0,ndisk,ntape)
  110 continue
c **********************************************************************
c **********************************************************************
c
      iflag = 0
      return
      end
