*deck drkfab
      subroutine drkfab (ncomp, xpts, nxpts, nfc, iflag, z, mxnon, p,
     +   ntp, ip, yhp, niv, u, v, w, s, stowa, g, work, iwork, nfcc)
c***begin prologue  drkfab
c***subsidiary
c***purpose  subsidiary to dbvsup
c***library   slatec
c***type      double precision (rkfab-s, drkfab-d)
c***author  watts, h. a., (snla)
c***description
c
c **********************************************************************
c
c     subroutine drkfab integrates the initial value equations using
c     the variable-step runge-kutta-fehlberg integration scheme or
c     the variable-order adams method and orthonormalization
c     determined by a linear dependence test.
c
c **********************************************************************
c
c***see also  dbvsup
c***routines called  dbvder, ddeabm, dderkf, dreort, dstor1
c***common blocks    dml15t, dml17b, dml18j, dml8sz
c***revision history  (yymmdd)
c   750601  date written
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  drkfab
c
      integer icoco, idid, iflag, igofx, indpvt, info, inhomo, integ,
     1     ipar, istkop, ivp, j, jflag, jon,
     2     k1, k10, k11, k2, k3, k4, k5, k6, k7, k8, k9, kkkint,
     3     kkkzpw, knswot, kod, kop, kopp, l1, l2, lllint, lotjp,
     4     mnswot, mxnon, mxnond, ncomp, ncompd, ndisk, neediw, needw,
     5     neq, neqivp, nfc, nfcc, nfccd, nfcd, nfcp1, nic, niv, non,
     6     nopg, nps, nswot, ntape, ntp, ntpd, numort, nxpts, nxptsd,
     7     ip(nfcc,*), iwork(*)
      double precision ae, c, g(*), p(ntp,*), pwcnd, px, re,
     1     s(*), stowa(*), tnd, tol, u(ncomp,nfc,*),
     2     v(ncomp,*), w(nfcc,*), work(*), x, xbeg, xend, xop,
     3     xot, xpts(*), xsav, xxop, yhp(ncomp,*), z(*)
c
c     ******************************************************************
c
      common /dml8sz/ c,xsav,igofx,inhomo,ivp,ncompd,nfcd
      common /dml15t/ px,pwcnd,tnd,x,xbeg,xend,xot,xop,info(15),istkop,
     1                knswot,kop,lotjp,mnswot,nswot
      common /dml18j/ ae,re,tol,nxptsd,nic,nopg,mxnond,ndisk,ntape,neq,
     1                indpvt,integ,nps,ntpd,neqivp,numort,nfccd,
     2                icoco
      common /dml17b/ kkkzpw,needw,neediw,k1,k2,k3,k4,k5,k6,k7,k8,k9,
     1                k10,k11,l1,l2,kkkint,lllint
c
      external dbvder
c
c      *****************************************************************
c       initialization of counters and variables.
c
c     begin block permitting ...exits to 220
c        begin block permitting ...exits to 10
c***first executable statement  drkfab
            kod = 1
            non = 1
            x = xbeg
            jon = 1
            info(1) = 0
            info(2) = 0
            info(3) = 1
            info(4) = 1
            work(1) = xend
c        ...exit
            if (nopg .eq. 0) go to 10
            info(3) = 0
            if (x .eq. z(1)) jon = 2
   10    continue
         nfcp1 = nfc + 1
c
c        ***************************************************************
c        *****beginning of integration loop at output
c        points.******************
c        ***************************************************************
c
         do 210 kopp = 2, nxpts
            kop = kopp
            xop = xpts(kop)
            if (ndisk .eq. 0) kod = kop
c
   20       continue
c
c              step by step integration loop between output points.
c
c              begin block permitting ...exits to 190
c                 begin block permitting ...exits to 30
                     xxop = xop
c                 ...exit
                     if (nopg .eq. 0) go to 30
                     if (xend .gt. xbeg .and. xop .gt. z(jon))
     1                  xxop = z(jon)
                     if (xend .lt. xbeg .and. xop .lt. z(jon))
     1                  xxop = z(jon)
   30             continue
c
c                 ******************************************************
   40             continue
c                    begin block permitting ...exits to 170
                        go to (50,60), integ
c                       dderkf integrator
c
   50                   continue
                           call dderkf(dbvder,neq,x,yhp,xxop,info,re,ae,
     1                                 idid,work,kkkint,iwork,lllint,g,
     2                                 ipar)
                        go to 70
c                       ddeabm integrator
c
   60                   continue
                           call ddeabm(dbvder,neq,x,yhp,xxop,info,re,ae,
     1                                 idid,work,kkkint,iwork,lllint,g,
     2                                 ipar)
   70                   continue
                        if (idid .ge. 1) go to 80
                           info(1) = 1
c                    ......exit
                           if (idid .eq. -1) go to 170
                           iflag = 20 - idid
c     .....................exit
                           go to 220
   80                   continue
c
c                       ************************************************
c                           gram-schmidt orthogonalization test for
c                           orthonormalization (temporarily using u and
c                           v in the test)
c
                        if (nopg .eq. 0) go to 100
                           if (xxop .eq. z(jon)) go to 90
c
c                             ******************************************
c                                 continue integration if we are not at
c                                 an output point.
c
c           ..................exit
                              if (idid .ne. 1) go to 200
c                    .........exit
                              go to 170
   90                      continue
                           jflag = 2
                        go to 110
  100                   continue
                           jflag = 1
                           if (inhomo .eq. 3 .and. x .eq. xend)
     1                        jflag = 3
  110                   continue
c
                        if (ndisk .eq. 0) non = numort + 1
                        call dreort(ncomp,u(1,1,kod),v(1,kod),yhp,niv,
     1                              w(1,non),s,p(1,non),ip(1,non),stowa,
     2                              jflag)
c
                        if (jflag .ne. 30) go to 120
                           iflag = 30
c     .....................exit
                           go to 220
  120                   continue
c
                        if (jflag .ne. 10) go to 130
                           xop = xpts(kop)
                           if (ndisk .eq. 0) kod = kop
c              ............exit
                           go to 190
  130                   continue
c
                        if (jflag .eq. 0) go to 140
c
c                          *********************************************
c                              continue integration if we are not at an
c                              output point.
c
c           ...............exit
                           if (idid .ne. 1) go to 200
c                    ......exit
                           go to 170
  140                   continue
c
c                       ************************************************
c                           store orthonormalized vectors into solution
c                           vectors.
c
                        if (numort .lt. mxnon) go to 150
                        if (x .eq. xend) go to 150
                           iflag = 13
c     .....................exit
                           go to 220
  150                   continue
c
                        numort = numort + 1
                        call dstor1(yhp,u(1,1,kod),yhp(1,nfcp1),
     1                              v(1,kod),1,ndisk,ntape)
c
c                       ************************************************
c                           store orthonormalization information,
c                           initialize integration flag, and continue
c                           integration to the next orthonormalization
c                           point or output point.
c
                        z(numort) = x
                        if (inhomo .eq. 1 .and. nps .eq. 0)
     1                     c = s(nfcp1)*c
                        if (ndisk .eq. 0) go to 160
                           if (inhomo .eq. 1)
     1                        write (ntape) (w(j,1), j = 1, nfcc)
                           write (ntape)
     1                           (ip(j,1), j = 1, nfcc),
     2                           (p(j,1), j = 1, ntp)
  160                   continue
                        info(1) = 0
                        jon = jon + 1
c                 ......exit
                        if (nopg .eq. 1 .and. x .ne. xop) go to 180
c
c                       ************************************************
c                           continue integration if we are not at an
c                           output point.
c
c           ............exit
                        if (idid .ne. 1) go to 200
  170                continue
                  go to 40
  180             continue
  190          continue
            go to 20
  200       continue
c
c           storage of homogeneous solutions in u and the particular
c           solution in v at the output points.
c
            call dstor1(u(1,1,kod),yhp,v(1,kod),yhp(1,nfcp1),0,ndisk,
     1                  ntape)
  210    continue
c        ***************************************************************
c        ***************************************************************
c
         iflag = 0
  220 continue
      return
      end
