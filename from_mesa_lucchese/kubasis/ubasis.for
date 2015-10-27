      PROGRAM basis
c
c computes contracted gaussian basis functions over a grid of points
c
c **** changed to READ input file called geobas made by mesa in m950 ***
c  all basis set and geometry input resides on that (binary) file
c             CWM  7/6/89
c
c
      PARAMETER (maxpts=500,maxpts4=4*maxpts,maxnbf=#maxnbfkohn)
      PARAMETER (maxprim=#maxprimkohn,maxprim1=20)
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (mxcen=maxprim1)
      COMMON/parms/ngrid,marg,nbf
      COMMON // buff(maxpts4)
      DIMENSION grid(maxpts,3),eta(maxprim,5),eta1(maxprim1,5)
      DIMENSION ll(maxprim),mm(maxprim),nn(maxprim),cent(4,mxcen)
     1,nfirst(maxnbf)
      DIMENSION l1(maxprim1),m1(maxprim1),n1(maxprim1)
      DIMENSION  nlast(maxnbf)
      DIMENSION vbuf(maxnbf*maxpts),vpp(maxpts,maxnbf)
     1,vpbuf(maxnbf*maxpts)
      DIMENSION v(maxpts,maxnbf),it(3),rit(3)
c..unicos
c      CALL link ("unit6=(outbas,create,hc),unit5=inbas
c     x ,unit9=(grid,open,abs),
c     x unit8=geobas,
c     x unit59=terminal//")
c..unicos
      OPEN(6,file='outbas',status='unknown')
      OPEN(5,file='inbas',status='unknown')
      OPEN(8,file='geobas',form='unformatted',status='old')
      OPEN(10,file='vbas',form='unformatted',status='unknown')
c      istime=time_()
      READ(5,*)npt
      WRITE(6,*)" npt=",npt
      irecl=4*npt*8
      CALL openabs(9,'grid',irecl)
      WRITE(6,*)" npt=",npt
      CALL rdabs(9,rit,3,0)
      it(1)=rit(1)
      it(2)=rit(2)
      it(3)=rit(3)
      WRITE(6,*)" npt=",npt
      WRITE(6,*)it
      ngrid=it(3)
      IF(it(1).NE.npt)THEN
         WRITE(6,*)" stopping because grid buffer length is wrong"
         WRITE(6,*)it(1),it(2),it(3)
         WRITE(6,*)npt,irecl
         STOP
      ENDIF
c
      maxtyp = 8
c
c  READ basis set DATA
c
      CALL inputs(nbf,ll,mm,nn,eta,nfirst,nlast,
     x nnuc,cent,maxnbf,mxcen,maxprim)
      WRITE(10)ngrid,npt,nbf
c
c READ in a BLOCK of grid points and transfer to a temporary location
c
      iset=1
      iread=0
      iquit=0
      marg=min0(ngrid,npt)
      nread=4*marg
      CALL rdabs(9,buff(1),nread,iset)
      iset=iset+1
      iread=iread+marg
 32    CONTINUE
      narg=marg
      CALL scopy(narg,buff(1),4,grid(1,1),1)
      CALL scopy(narg,buff(2),4,grid(1,2),1)
      CALL scopy(narg,buff(3),4,grid(1,3),1)
c      CALL scopy(narg,buff(4),4,wt(1),1)
      iremn=ngrid-iread
      IF(iremn.EQ.0)THEN
        iquit=1
        go to 34
      ENDIF
      marg=min0(iremn,npt)
      nread=4*marg
      CALL rdabs(9,buff(1),nread,iset)
      iset=iset+1
      iread=iread+marg
34    CONTINUE
c
c  OPEN loops on basis functions to evaluate their
c   contributions at the current BLOCK
c  of grid points
c
c      CALL second(time0)
      DO 5 j=1,nbf
      DO 5 i=1,narg
      vpp(i,j) = 0.
 5    v(i,j) = 0.e0
      ij=0
      DO 500 i=1,nbf
      is=nfirst(i)
      IF=nlast(i)
      icn= IF-is +1
      IF(icn .GT. maxprim1) THEN
         WRITE(6,*)'Error too many primatives in one contraction',
     $        icn, maxprim1
         STOP 'Too many primitives'
      END IF
      DO 50 ii=is,IF
      DO 50 iii=1,5
      indx  = ii-is+1
      l1(indx) = ll(ii)
      m1(indx)=mm(ii)
      n1(indx)=nn(ii)
50    eta1(indx,iii)=eta(ii,iii)
c      IF(i.EQ.19)WRITE(*,"(5f8.4)")((eta1(iz,jz),jz=1,5),iz=1,is-IF+1)
      CALL gaussval(i,maxpts,narg,grid,icn,eta1,l1,m1,n1,v,vpp,
     $     maxnbf,maxprim1)
  500 CONTINUE
c      CALL second(time1)
c      time = time1-time0
c      WRITE(6,143) time
143   FORMAT(' gaussval time for this block of grid points = ',f10.5)
      DO 601 ii=1,nbf
      idum=(ii-1)*narg + 1
      CALL scopy(narg,v(1,ii),1,vbuf(idum),1)
      CALL scopy(narg,vpp(1,ii),1,vpbuf(idum),1)
601   CONTINUE
      WRITE(10)(vbuf(i),i=1,narg*nbf)
      WRITE(10)(vpbuf(i),i=1,narg*nbf)
c
c RETURN to start a fetch another BLOCK of points
c
      IF(iquit.EQ.0)go to 32
c      iftime=time_()
c      itime=iftime-istime
c      WRITE(6,884)itime
884   FORMAT(' times: cpu ',i10)
      CALL closeabs(9)
      STOP
      END
