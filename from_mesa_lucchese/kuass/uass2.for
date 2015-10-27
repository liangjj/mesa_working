      PROGRAM uass2
c
c  computes ylms on grid
c
c  REAL valued ylms normed properly to 1 - 10/6/88
c ****  tnr 11/20/2013 replaced plm routine WITH more stable lpmn routine for
c  computing associated Legendre functions
c ****
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (ngrid=500,ng7=7*ngrid)
      PARAMETER (llmax=#maxltop,jmax=ngrid*2*llmax,imax=4*ngrid)
      PARAMETER (fourpi=4.d0*3.14159265358979d0)
      DIMENSION it(3),rit(3)
      COMPLEX*16 cmp
      COMMON /parm/nbig,npt,lmax,mumax
      COMMON/stuff/scr(ng7)
      DIMENSION buff(jmax),biff(imax)
      DIMENSION p(ngrid,0:llmax),x(ngrid),pml(0:llmax,0:llmax,ngrid)
      DIMENSION xx(ngrid),yy(ngrid),zz(ngrid),r(ngrid),sphi(ngrid)
     1 ,cmp(ngrid),cmphi(ngrid),smphi(ngrid),cphi(ngrid)
      DIMENSION pp(ngrid,0:llmax),pm(ngrid,0:llmax),fac(0:60)
      DATA tol/1.e-8/
      OPEN(6,file='outass',status='unknown')
      OPEN(5,file='inass',status='unknown')
      OPEN(10,file='ylms',form='unformatted',status='unknown')
      READ(5,*,END=999)npt,mumax,lmax
      irecl=4*npt*8
      CALL openabs(9,'grid',irecl)
      CALL rdabs(9,rit,3,0)
      it(1)=rit(1)
      it(2)=rit(2)
      it(3)=rit(3)
      WRITE(6,*)it
      nbig=it(3)
      IF(it(1).NE.npt)THEN
         WRITE(6,*)" Error stopping because grid buffer length is wrong"
         WRITE(6,*)it(1),it(2),it(3)
         WRITE(6,*)npt,irecl
         STOP
      ENDIF
c
      lmmax=lmax+mumax
      fac(0)=1.
      DO 2 i=1,lmmax
2     fac(i)=i*fac(i-1)
c
c READ in a BLOCK of grid points and transfer to a temporary location
c
      iset=1
      iread=0
      iquit=0
c
      WRITE(10)nbig,npt,lmax,mumax
      marg=min0(nbig,npt)
      nread=4*marg
      CALL rdabs(9,biff(1),nread,iset)
      iset=iset+1
      iread=iread+marg
 32    CONTINUE
      narg=marg
      CALL scopy(narg,biff(1),4,xx(1),1)
      CALL scopy(narg,biff(2),4,yy(1),1)
      CALL scopy(narg,biff(3),4,zz(1),1)
      iremn=nbig-iread
      IF(iremn.EQ.0)THEN
        iquit=1
        go to 34
      ENDIF
      marg=min0(iremn,npt)
      nread=4*marg
      CALL rdabs(9,biff(1),nread,iset)
      iset=iset+1
      iread=iread+marg
34    CONTINUE
      DO 3 i=1,narg
      r(i)=SQRT(xx(i)**2+yy(i)**2+zz(i)**2)
      x(i)=zz(i)/r(i)
      xsq=x(i)**2
      xqq=ABS(1.-xsq)
      IF(xqq.LT.tol)THEN
      cphi(i)=1.   
      sphi(i)=0.
      IF(x(i))66,67,68
 66   x(i)=-1.
      go to 3
 67   WRITE(6,*)" Error screwup in ass1"
      STOP
 68   x(i)=1.
      ELSE
      cphi(i)=xx(i)/r(i)/SQRT(1.-xsq)
      sphi(i)=yy(i)/r(i)/SQRT(1.-xsq)
      ENDIF
3     CONTINUE
      CALL lpmn(llmax,mumax,lmax,pml,narg,x)
      DO 1 mu=0,mumax
         DO j=1,narg
         DO i=0,lmax
            p(j,i)=pml(mu,i,j)
         ENDDO
         ENDDO
c      CALL plm(x,narg,mu,lmax,p)
      IF(mu.EQ.0)THEN
       DO 10 i=0,lmax
       const=SQRT((2*i+1)/fourpi)
       DO 10 j=1,narg
       p(j,i)=p(j,i)*const
10     CONTINUE
       DO 11 i=0,lmax
       CALL scopy(narg,p(1,i),1,buff(1+i*narg),1)
c       WRITE(*,*)"ylm,l=,m=",i,mu
c       WRITE(*,"(2i5,4f10.4)")(i,mu,xx(iz),yy(iz),zz(iz),
c     &     buff(i*narg+iz),iz=1,narg)
11     CONTINUE
       jbuf=narg*(lmax+1)
       WRITE(10)(buff(i),i=1,jbuf)
      ELSE
       DO 4 i=1,narg
       cmp(i)=(CMPLX(cphi(i),sphi(i)))**mu
       cmphi(i)=DBLE(cmp(i))
       smphi(i)=dimag(cmp(i))
4      CONTINUE
       DO 12 i=mu,lmax
       const=SQRT((2*i+1)/fourpi*fac(i-mu)/fac(i+mu))
c  SQRT(2) factor added to normalize "real valued" ylms
       const=const*SQRT(2.0)
       DO 12 j=1,narg
       pp(j,i)=p(j,i)*const*cmphi(j)
       pm(j,i)=p(j,i)*const*smphi(j)
12     CONTINUE
       DO 13 i=mu,lmax
       CALL scopy(narg,pp(1,i),1,buff(1+(i-mu)*2*narg),1)
       CALL scopy(narg,pm(1,i),1,buff(1+narg+(i-mu)*2*narg),1)
c       WRITE(*,*)"ylm,l=,m=",i,2*mu-1 
c       WRITE(*,"(2i5,4f10.4)")(i,2*mu-1,xx(iz),yy(iz),zz(iz),
c     &     buff((i-mu)*2*narg+iz),iz=1,narg)
c
c       WRITE(*,*)"ylm,l=,m=",i,2*mu
c       WRITE(*,"(2i5,4f10.4)")(i,2*mu,xx(iz),yy(iz),zz(iz),
c     &     buff((i-mu)*2*narg+narg+iz),iz=1,narg)
13     CONTINUE
       jbuf=2*narg*(1+lmax-mu)
       WRITE(10)(buff(i),i=1,jbuf)
      ENDIF
1     CONTINUE
c
c RETURN to start a fetch another BLOCK of points
c
      IF(iquit.EQ.0)go to 32
c      DO 2 i=1,narg
c2     WRITE(6,100)x(i),mu,(p(i,j),j=mu,lmax)
100   FORMAT(" x=",f10.5,"  mu = ",i3/(10e12.4))
999   CONTINUE
      CALL closeabs(9)
      STOP
      END
