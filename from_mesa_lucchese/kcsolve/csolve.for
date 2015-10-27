      PROGRAM solver
c note IF npvec .NE. nbtot, q must be zero
      IMPLICIT REAL*8(a-h,o-z)
c       PARAMETER (nbig=800, nsmall=150,nchnl=20,maxnqv=50)
      PARAMETER (nbig=#maxbig, nsmall=#maxsmall,nchnl=#maxchan,maxnqv=50)
      COMPLEX*16  hdenom(nbig,nbig),smat(nsmall,nsmall),work(nbig)
      COMPLEX*16 tmat(nsmall,nsmall),zdotu,cov(nbig,nsmall)
      COMPLEX*16 htop(nbig,nsmall),htopp(nbig,nsmall),det
      COMPLEX*16 store(nsmall*nbig) 
      REAL*8 kchan(nchnl),hpp(nbig,nbig),hqq(nbig,nbig),
     1    hpq(nbig,nsmall),hfix(nbig),xsecmat(nchnl,nchnl)
c      COMMON /hold/ hdenom, hqq
      REAL*8 tq(nbig*maxnqv)
      INTEGER nlm(nchnl),kpvt(nbig),iclosed(nchnl),jkill(nbig),
     # kkill(nbig)
      EQUIVALENCE (hqq,hpq,hpp,store),(hdenom,htop)
      CHARACTER*8 istat,istatic,idump,ikill
      CHARACTER*9 iconq,icon
      CHARACTER*128 namkohn
c      DATA istatic/6hstatic/,idump/4hdump/,ikill/4hkill/,
c     #icon/9hcontractq/
      istatic='static  '
      idump='dump    '
      ikill='kill    '
      icon='contractq'
      OPEN(5,file='insolve')
      OPEN(6,file='outsolve')
c      CALL second(time)
      WRITE(6,*)'starting csolve'
      DO i=1,nbig
         DO j=1,nsmall
            cov(i,j)=0.
         ENDDO
      ENDDO
      READ(5,*)iprint,iopen,iestart
      READ(5,*)iconq
      IF(iopen.EQ.1)THEN
      READ(5,*)dele,coef
      WRITE(6,*)' energy defect and ci coefficient',dele,coef
      ENDIF
      READ(5,55)istat
 55   FORMAT(a8)
      WRITE (6,*) 'iconq', iconq
      IF(iconq.EQ.icon)THEN
         WRITE (6,*) 'calling drum'
         CALL drum
         WRITE (6,*) 'calling iosys'
         CALL iosys('read character "kohn filename" from rwf',
     #       0,0,0,namkohn)
         CALL iosys('open kohn as old',0,0,0,namkohn)
         CALL iosys('read integer "q-space nq" from kohn',1,nq,0,' ')
         CALL iosys('read integer "q-space nqvecs" from kohn',
     $        1,nqvecs,0,' ')
         WRITE(6,*)'unfolding contracted q-space vectors'
         WRITE(6,*)'nq from kohn=',nq
         WRITE(6,*)'nqvecs from kohn=',nqvecs
      ENDIF
      OPEN(7,file='bmesa',form='unformatted')
      OPEN(9,file='bcoef',form='unformatted')
      OPEN(77,file='btmat',form='unformatted')
      OPEN(8,file='bsolve',form='unformatted')
      OPEN(78,file='cmesa',form='unformatted')
c!!!!!!!!!!
      IF(istat.EQ.ikill)THEN
         READ(5,*)nkill
         READ(5,*)(kkill(i),i=1,nkill)
         WRITE(6,*)'The following q-space terms will be killed:'
         WRITE(6,*)(kkill(i),i=1,nkill)
      ENDIF
      IF(istat.EQ.istatic)THEN
      OPEN(99,file='bstatic',form='unformatted')
      ENDIF
      READ(7)nchan,nbtot,nfree,(nlm(i),i=1,nchan)
      READ(8)npvec,npdim
      WRITE(6,*)' npdim ',npdim
      IF(iconq.EQ.icon)THEN
         IF(nqvecs.NE.npdim)THEN
            WRITE(6,*)'Error stopping because nqvecs ', 
     $           'and npdim are different'
            STOP
         ENDIF
         IF((nq*nqvecs).GT.(maxnqv*nbig))THEN
            WRITE(6,*)'Error stopping because ', 
     $           '(nq*nqvecs).gt.(maxnqv*nbig)'
            STOP
         ENDIF
         CALL iosys('read real "q-space contraction vectors" from kohn',
     $        nq*nqvecs,tq,0,' ')
c      WRITE(6,*)'q-space vectors'
c      CALL matout(tq,nq,nqvecs,nq,nqvecs,6)
      ENDIF
      IF(npvec.NE.nbtot)THEN
      WRITE(6,*)' npvec and nbtot are unequal!', nbtot, npvec
      ENDIF
      WRITE(6,112)nchan,nbtot,nfree,(nlm(i),i=1,nchan)
 112  FORMAT(' number of channels=',i3/
     ^ ' dimension of P-space is ',i5/
     ^ ' number of free (lm) channels is =',i4/
     ^ ' number of lm terms per channel:'/(20i3))
      ntot=nbtot+nfree
      IF(nchan.GT.nchnl.OR.nfree.GT.nsmall.OR.
     ^ ntot.GT.nbig .OR. npdim .GT. nbig)THEN
      WRITE(6,*)' Error check dimensions'
      WRITE(6,*) 'nchan', nchan
      WRITE(6,*) 'nchnl', nchnl
      WRITE(6,*) 'nfree', nfree
      WRITE(6,*) 'nsmall', nsmall
      WRITE(6,*) 'ntot', ntot
      WRITE(6,*) 'nbig', nbig
      WRITE(6,*) 'npdim', npdim
      STOP
      ENDIF
      nall=ntot+npdim
      IF(iconq.NE.icon)THEN
         WRITE(9)nchan,nall,nfree
      ELSE
         nalll=ntot+nq
         WRITE(9)nchan,nalll,nfree
      ENDIF
      iw=nbtot+1
      READ(7)nener
c loop on energies
      DO 1 ie=1,nener
      DO 100 i=1,nbig
      DO 101 j=1,nbig
 101     hdenom(i,j)=0.
      DO 102 j=1,nsmall
 102  htopp(i,j)=0.
 100     CONTINUE
      READ(7)(kchan(ic),ic=1,nchan)
c*tnr* dropping the WRITE of kchan to unit 9******
c      WRITE(9)(kchan(ic),ic=1,nchan)
      WRITE(6,111) (kchan(i),i=1,nchan)
111   FORMAT(///,' ********************************',
     # ' new energy ********************************',
     # //,' channel momenta:',6e12.5)
c
c determine number of OPEN channels
c
      READ(78)(iclosed(ic),ic=1,nchan)
      WRITE(6,*)'iclosed',(iclosed(ic),ic=1,nchan)
      nopen=nchan
      DO i=1,nchan
         nopen=nopen-iclosed(i)
      ENDDO
      nofree=0
      notot=nbtot
      DO i=1,nopen
         notot=notot+nlm(i)
         nofree=nofree+nlm(i)
      ENDDO
      noall=notot+npdim
c!!!!!!!!!!!!!!!!!!!!!!!!
      IF(istat.EQ.ikill)THEN
c offset jkill by notot
         DO i=1,nkill
            jkill(i)=notot+kkill(i)
         ENDDO
      ENDIF
c ***tnr*** adding the WRITE to unit 9 of nofree and noall
      WRITE(6,*)'nofree, noall to bcoef=',nofree, noall
      IF(iconq.NE.icon)THEN
         WRITE(9)nofree,noall
      ELSE
         noalll=notot+nq
         WRITE(9)nofree,noalll
      ENDIF
      
      IF (noall .GT. nbig) THEN
         WRITE(6,*)'Error stopping because noall ', noall,
     $        ' gt nbig ', nbig
         STOP
      END IF
      IF(iconq.NE.icon)THEN
         IF (noalll .GT. nbig) THEN
            WRITE(6,*)'Error stopping because noalll ', noalll,
     $           ' gt nbig ', nbig
            STOP
         END IF
      ENDIF
      
      WRITE(6,*)'number of open channels=',nopen
      WRITE(6,*)'number of bound + free functions=',notot
      WRITE(6,*)'number of bound + free+q-space functions=',noall
      IF(nbtot.EQ.npvec)THEN
         READ(7)((hdenom(i,j),i=1,notot),j=iw,notot)
         DO 114 i=1,notot
            DO 114 j=iw,notot
 114        hdenom(j,i)=hdenom(i,j)   
      ELSE
         READ(7)((hdenom(i,j),i=1,ntot),j=1,ntot)
      ENDIF
      CALL rdbinsqr(hpp,nbig,npvec,8)
      IF(istat.EQ.istatic)THEN
         CALL rdbinsqr(hpp,nbig,nbtot,99)
      ENDIF
      IF(nbtot.EQ.npvec.OR.istat.EQ.istatic)THEN
         DO 2 i=1,nbtot
            DO 2 j=1,nbtot
  2      hdenom(i,j)=hpp(i,j)
      ENDIF
      IF(npdim.GT.0)THEN
         CALL rdbinsqr(hqq,nbig,npdim,8)
c!!!!!!!!!!!!!!

         IF(istat.EQ.idump)THEN
            WRITE(6,*)'first energy'
            WRITE(6,667)(hqq(i,i),i=1,npdim)
 667        FORMAT(10e12.4)
            STOP
         ENDIF
         DO 3 i=1,npdim
            DO 3 j=1,npdim
               ii=i+notot
               jj=j+notot
  3         hdenom(ii,jj)=hqq(i,j)
            IF(iopen.EQ.1)THEN
c      hfix(1)=hqq(1,1)+etarget+kchan(1)**2/2.
               hfix(1)=hqq(1,1)+dele
               IF(npdim.GT.1)THEN
                  DO 7 i=2,npdim
 7                hfix(i)=hqq(1,i)
               ENDIF
            ENDIF
       DO 4 i=1,npvec
          READ(8)(hpq(j,i),j=1,npdim)
          DO 5 j=1,npdim
             hdenom(i,j+notot)=hpq(j,i)
 5        hdenom(j+notot,i)=hpq(j,i)
  4    CONTINUE
      ENDIF
      IF(iprint.NE.0) THEN
         WRITE(6,107)
 107     FORMAT(//' denominator matrix of (h-e)')
         DO 60 i=1,noall
 60      WRITE(6,108) i,(hdenom(j,i),j=1,noall)
 108     FORMAT(1x,i3,6e12.5,/,(4x,6e12.5))
      ENDIF
      READ(7)((htopp(i,j),i=1,notot),j=1,nofree)
      READ(7)nread
      WRITE(6,*)' no. occupied plus scattering terms', nread
      READ(7)((cov(i,j),i=1,nread),j=1,nfree)
      DO 6 k=1,nbtot
         DO 6 i=1,npdim
            DO 6 j=1,nofree
  6   htopp(i+notot,j)=htopp(i+notot,j)+hpq(i,k)*cov(k,j)
c
c add extra terms to RHS for correlated TARGET CASE
c
      IF(iopen.EQ.1)THEN
       DO 8 i=1,npdim
       DO 8 j=1,nofree
  8       htopp(i+notot,j)=htopp(i+notot,j)+coef*hfix(i)
     # *cov(nread,j)
c**********NEW**********
       DO 88 i=1,nbtot
       DO 88 j=1,nofree
  88      htopp(i,j)=htopp(i,j)+coef*hpq(1,i)*cov(nread,j)
c**********************
      ENDIF
      IF(iprint.NE.0) THEN
         WRITE(6,109)
 109     FORMAT(//' numerator matrix of (h-e)')
         DO 61 i=1,nofree
 61      WRITE(6,108) i,(htopp(j,i),j=1,noall)
      ENDIF
c      CALL second(tnow)
      elapse=tnow-time
      WRITE(6,777)elapse
 777  FORMAT(" elapsed time for setup is ",e12.4)
      time=tnow      
      IF(ie.GE.iestart)THEN
c      CALL matinv(hdenom,nall,htopp,nfree,det,nbig)
c      CALL csifa(hdenom,nbig,nall,kpvt,info)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         nosolve=noall
c strike nkill columns from hdenom
      IF(istat.EQ.ikill)THEN
c         WRITE(6,*)'before, hdenom(n,j)',(hdenom(noall,i),i=1,noall)
         nosolve=noall-nkill
          DO i=1,noall
             jcnt=1
             jkeep=1
             DO j=1,noall
                IF(j.EQ.jkill(jcnt))THEN
                   IF(jcnt.LT.nkill)THEN
                      jcnt=jcnt+1
                   ENDIF
                ELSE
                   hdenom(i,jkeep)=hdenom(i,j)
                   jkeep=jkeep+1
                ENDIF
             ENDDO
          ENDDO
c strike nkill rows from hdenom
          DO j=1,noall-nkill
             icnt=1
             ikeep=1
             DO i=1,noall
                IF(i.EQ.jkill(icnt))THEN
                   IF(icnt.LT.nkill)THEN
                      icnt=icnt+1
                   ENDIF
                ELSE
                   hdenom(ikeep,j)=hdenom(i,j)
                   ikeep=ikeep+1
                ENDIF
             ENDDO
          ENDDO
c         WRITE(6,*)'after, hdenom(n,j)',(hdenom(nosolve,i),i=1,nosolve)
c         WRITE(6,*)'after, hdenom(j,n)',(hdenom(i,nosolve),i=1,nosolve)
c strike nkill rows from htopp
c          WRITE(6,*)'before, htopp(j,n)',(htopp(i,nofree),i=1,noall)
          DO j=1,nofree
             icnt=1
             ikeep=1
             DO i=1,noall
                IF(i.EQ.jkill(icnt))THEN
                   IF(icnt.LT.nkill)THEN
                      icnt=icnt+1
                   ENDIF
                ELSE
                   htopp(ikeep,j)=htopp(i,j)
                   ikeep=ikeep+1
                ENDIF
             ENDDO
          ENDDO
       ENDIF
c         WRITE(6,*)'after, htopp(j,n)',(htopp(i,nofree),i=1,nosolve)
         CALL cgefs(hdenom,nbig,nosolve,htopp(1,1),1,ind,work,kpvt)
         DO 666 i=2,nofree
 666        CALL cgefs(hdenom,nbig,nosolve,htopp(1,i),2,ind,work,kpvt)
         ENDIF
         IF(iprint.NE.0) THEN
            WRITE(6,187)
 187        FORMAT(//' solution matrix ')
            DO 79 i=1,nosolve
c            DO 79 i=1,npdim
 79            WRITE(6,108) i,(htopp(i,j),j=1,nofree)
c 79            WRITE(6,108) i,(htopp(notot+i,j),j=1,nofree)
         ENDIF
c!!!!!!!!!!!!!!!!
c pad htopp WITH zeros IF necessary
      IF(istat.EQ.ikill)THEN
c         WRITE(6,*)'sol, htopp(j,n)',(htopp(i,nofree),i=1,nosolve)
         DO i=1,nosolve
            DO j=1,nofree
               hdenom(i,j)=htopp(i,j)
            ENDDO
         ENDDO
         DO j=1,nofree
            icnt=1
            ifill=1
            DO i=1,noall
               IF(i.EQ.jkill(icnt))THEN
                  htopp(i,j)=0.0
                  IF(icnt.LT.nkill)THEN
                     icnt=icnt+1
                  ENDIF
               ELSE
                  htopp(i,j)=hdenom(ifill,j)
                  ifill=ifill+1
               ENDIF
            ENDDO
         ENDDO
c         WRITE(6,*)'new sol, htopp(j,n)',(htopp(i,nofree),i=1,noall)
      ENDIF
      WRITE(6,*)'iconq,icon: ',iconq,icon
      IF(iconq.NE.icon)THEN 
         WRITE(9)((htopp(i,j),i=1,noall),j=1,nofree)
      ELSE
         ij=0
         DO i=1,nofree
            DO j=1,nqvecs
               ij=ij+1
               store(ij)=htopp(notot+j,i)
            ENDDO
         ENDDO
c         WRITE(6,*)'q-space coefficients'
c         CALL cmatout(store,nqvecs,nofree,nqvecs,nofree,6)
         CALL expand(htopp,store,tq,nbig,nsmall,nofree,
     $    notot,nqvecs,nq)
c            WRITE(6,*)' solution matrix in original basis'
c            DO 799 i=1,noall
c 799           WRITE(6,108) i,(htopp(i,j),j=1,nofree)
         WRITE(9)((htopp(i,j),i=1,noalll),j=1,nofree)
      ENDIF
c      CALL second(tnow)
      elapse=tnow-time
      WRITE(6,778)elapse
      time=tnow
 778  FORMAT(" time for matinv is ",e12.4)
      DO 103 i=1,nbig
      DO 103 j=1,nsmall
 103  htop(i,j)=0.
      READ(7)((htop(i,j),i=1,notot),j=1,nofree)
      IF(iprint.NE.0) THEN
      WRITE(6,110)
110   FORMAT(//' conjugate numerator matrix of (h-e)')
      DO 62 i=1,nofree
62    WRITE(6,108) i,(htop(j,i),j=1,noall)
      ENDIF
      DO 500 ilm=1,nofree
      DO 500 jlm=1,nofree
      tmat(ilm,jlm)=zdotu(noall,htopp(1,jlm),1,htop(1,ilm),1)
500   CONTINUE
      READ(7)((smat(i,j),i=1,nofree),j=1,nofree)
      IF(ie.LT.iestart)go to 1
      IF(iprint.NE.0) THEN
      WRITE(6,113)
 113  FORMAT(//' first Born matrix')
      DO 63 i=1,nofree
63    WRITE(6,108) i,(smat(j,i),j=1,nofree)
      ENDIF
      DO 501 ilm=1,nofree
      DO 501 jlm=1,nofree
      tmat(ilm,jlm)=-2.*(-tmat(ilm,jlm)+smat(ilm,jlm))
 501  CONTINUE
      WRITE(6,115)
 115  FORMAT(//' T-Matrix')
      DO 64 i=1,nofree
64    WRITE(6,108) i,(tmat(j,i),j=1,nofree)
      istart = 0
      DO 879 ic=1,nopen
      jstart = 0
      ni = nlm(ic)
      DO 878 jc=1,nopen
      nj=nlm(jc)
      ist=istart+1
      ifin=istart+ni
      jst=jstart+1
      jfin=jstart+nj
      WRITE(77)ic,jc,ni,nj,kchan(ic),kchan(jc)
c      WRITE(88,177)ic,jc,ni,nj,kchan(ic),kchan(jc)
 177  FORMAT(4i5,2f20.10)
      WRITE(77)((tmat(ii,jj),ii=ist,ifin),jj=jst,jfin)
c      WRITE(88,277)((tmat(ii,jj),ii=ist,ifin),jj=jst,jfin)
 277  FORMAT(4e20.10)
      summod = 0.
      DO 877 ilm=1,ni
      DO 877 jlm=1,nj
      isub = istart + ilm
      jsub = jstart + jlm
877   summod = summod + ABS(tmat(isub,jsub))**2
      xsecmat(ic,jc) = 4.0*3.141592654*summod/kchan(ic)**2
      jstart = jstart + nlm(jc)
878   CONTINUE
      istart = istart + nlm(ic)
879   CONTINUE
       WRITE(6,874)
874   FORMAT(//,' total cross sections: row index = initial chnl,',
     # ' column index = final chnl')
      WRITE(6,871) (i,i=1,nopen)
871   FORMAT(4x,6(6x,i2,4x))
      DO 873 i=1,nopen
      WRITE(6,872) i, (xsecmat(i,j),j=1,nopen)
873   CONTINUE
872   FORMAT(1x,i2,1x,6e12.5/(4x,6e12.5))
c      CALL second(tnow)
      elapse=tnow-time
      WRITE(6,779)elapse
      time=tnow
 779  FORMAT("  time at end is ",e12.4)
 1    CONTINUE
      IF(iconq.EQ.icon)THEN
         CALL iosys('close kohn',0,0,0,' ')
         CALL iosys('close all',0,0,0,' ')
      ENDIF

c      CALL EXIT
      END
