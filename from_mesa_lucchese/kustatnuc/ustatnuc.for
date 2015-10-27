      PROGRAM ustatnuc
      PARAMETER (maxchan=#maxchan)
      PARAMETER (maxp=maxchan*(maxchan+1)/2)
      PARAMETER (mxbuf=1000)
      PARAMETER (maxnbf=#maxnbfkohn)
      PARAMETER (maxprim=#maxprimkohn)
      PARAMETER (maxprim1=20)
      PARAMETER (mxcen=20)
      IMPLICIT REAL*8 (a-h,o-z)
c     
c     changed 6/11/92 by BHL to speed up transformation in vints
c
c vectorized computation of static potential*quadrature weights
c on a set of grid points
c
c  this version also works for coulomb CASE, for which it substracts
c  zresidual/r from static potential
c
c  changed output file to INCLUDE parameters at start CWMc 7-22-88
c
c   ***** changed to READ density matrices and basis FUNCTION information
c         from files denmat and geobas which are made by the mesa code m950
c              CWM 7/7/89
c
c  changed to eliminate basis functions for which the sum of
c   ABS(rho(i,j)) summed over j and channel pairs is less than a
c   given tolerance.   CWMc 7/10/89
c
c changed 7/28/91 by TNR to speed up vectorized computation of vints
c
c  variable DIMENSION for number of channels added via PARAMETER statement
c
c  1-25-90 ERROR in elim corrected, computation of ndel taken out of loop
c
      DIMENSION dmat(maxnbf, maxnbf)
      DIMENSION rit(3)
      INTEGER*4 it(3)
      DIMENSION grid(mxbuf,3),eta(maxprim,5)
      DIMENSION eta1(maxprim1,5),eta2(maxprim1,5)
      DIMENSION ll(maxprim),mm(maxprim),nn(maxprim),cent(4,mxcen)
      DIMENSION nfirst(maxnbf)
      DIMENSION l1(maxprim1),l2(maxprim1),m1(maxprim1),m2(maxprim1)
      DIMENSION n1(maxprim1),n2(maxprim1)
      DIMENSION   nlast(maxnbf)
      DIMENSION vbuf(maxp*mxbuf),wt(mxbuf)
      DIMENSION v(mxbuf,maxp),den(maxp)
      DIMENSION allrho((maxnbf*(maxnbf+1))/2,maxp),labels(maxp)
      DIMENSION soo(maxp)
      DIMENSION sumpot(maxp)
      DIMENSION buff(40000)
      INTEGER ivparms(2)
      COMMON /nmbrs /  pi, piterm, pitern, acrcy, scale, icanon
c..unicos      COMMON /nmbrs /  pi, piterm, pitern
c
c.     CALL link ("unit6=(outstat,create,hc),unit5=instat
c.    x ,unit9=(grid,open,abs),
c.    x unit8=geobas,
c.    x unit12=denmat,
c.    x unit59=terminal//")
c
      OPEN(5,file='instat',status='unknown')
      OPEN(66,file='outstat',status='unknown')
      OPEN(8,file='geobas',form='unformatted',status='old')
      OPEN(12,file='denmat',form='unformatted',status='old')
      OPEN(10,file='vstat',form='unformatted',status='unknown')
c.bhl routine to mimic ltss
      READ(5,*) zres
      WRITE(6,*)' ZRes  ',zres
      READ(5,*)npt
      WRITE(6,*)' NPT ',npt
      irecl=4*npt*8
      CALL openabs(9,'grid',irecl)
      CALL rdabs(9,rit,3,0)
      it(1)=rit(1)
      it(2)=rit(2)
      it(3)=rit(3)
      ngrid=it(3)
      IF(it(1).NE.npt)THEN
         WRITE(6,*)" Error stopping because grid buffer length is wrong"
         WRITE(6,*)it(1),it(2),it(3)
         WRITE(6,*)npt,irecl
         STOP
      ENDIF
c
      pi = 3.14159265358979e0
      piterm=2.e0/pi**0.5e0
      pitern=pi**1.5e0
      maxtyp = 8
      CALL generf (maxtyp,maxrng)
      DO 2 j=1,maxp
   2   sumpot(j) = 0.0
c
c  READ basis set DATA
c
      CALL inputs(nbf,ll,mm,nn,eta,nfirst,nlast,
     x nnuc,cent,maxnbf,mxcen,maxprim)
c
c  READ in residual charge
c

      IF(ABS(zres).GT.1.e-8) THEN
      nnuc=nnuc+1
      IF (nnuc .GT. mxcen) THEN
         WRITE(6,*) 'Error nnuc > mxcen'
         STOP 'Bad nnuc'
      END IF
      cent(1,nnuc)=0.
      cent(2,nnuc)=0.
      cent(3,nnuc)=0.
      cent(4,nnuc)=-zres
      ENDIF
c
c  READ density matrix -- in this version it refers
c  to contracted gaussians and comes from mesa m950
c
      CALL rho(istate,labels,allrho,nbf,maxp,maxnbf,dmat)
c
c  eliminate basis functions which DO not contribute to static
c  potential for any channel pair
c
      CALL elim(nbf,ll,mm,nn,eta,nfirst,nlast,
     x allrho,istate,maxp,maxnbf,maxprim)
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
      ivparms(1)=istate
      ivparms(2)=ngrid
      WRITE(6,*)' NGrid ',ngrid
      WRITE(10)(ivparms(i),i=1,2)
 32    CONTINUE 
      narg=marg
      CALL scopy(narg,buff(1),4,grid(1,1),1)
      CALL scopy(narg,buff(2),4,grid(1,2),1)
      CALL scopy(narg,buff(3),4,grid(1,3),1)
      CALL scopy(narg,buff(4),4,wt(1),1)
c      DO 859 i=1,narg
c 859      WRITE(66,858)i,(grid(i,j),j=1,3),wt(i)
 858      FORMAT(i5,4d12.4)
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
c  contributions to the static potential at the current BLOCK
c  of grid points
c
      DO 5 j=1,istate
      DO 5 i=1,narg
 5    v(i,j) = 0.e0
      ij=0
      DO 500 i=1,nbf
      is=nfirst(i)
      IF=nlast(i)
      icn= IF-is +1
      IF (icn .GT. maxprim1) THEN
         WRITE(6,*) 'Error too man primitives', icn, maxprim1
         STOP 'Too many primitives'
      END IF
      DO 50 ii=is,IF
      DO 50 iii=1,5
      indx  = ii-is+1
      l1(indx) = ll(ii)
      m1(indx)=mm(ii)
      n1(indx)=nn(ii)
      eta1(indx,iii)=eta(ii,iii)
 50   CONTINUE
      DO 500 j=1,i
      ij=ij+1
      DO 53 ists=1,istate
53    den(ists) = allrho(ij,ists)
      js=nfirst(j)
      jf=nlast(j)
      jcn= jf-js +1
      DO 60 jj=js,jf
      DO 60 jjj=1,5
      jndx  = jj-js+1
      l2(jndx) = ll(jj)
      m2(jndx)=mm(jj)
      n2(jndx)=nn(jj)
      eta2(jndx,jjj)=eta(jj,jjj)
 60   CONTINUE
c
      CALL vints(l1,l2,m1,m2,n1,n2,icn,jcn,eta1,eta2,grid,narg
     x ,v,den,soo,istate,maxp,maxprim1)
  500 CONTINUE
c
c  add in the contribution of the nuclear attraction potentials
c
      DO 602 ii=1,istate
      i1=labels(ii)/1000
      j1=labels(ii)-1000*i1
      IF(i1.NE.j1)go to 602
      DO 600 i=1,nnuc
      DO 599 j=1,narg
      dist=(cent(1,i)-grid(j,1))**2+(cent(2,i)-grid(j,2))**2
     x + (cent(3,i)-grid(j,3))**2
      potn=  cent(4,i)/SQRT(dist)
599   v(j,ii)=v(j,ii) - potn
600   CONTINUE
  602 CONTINUE

      DO 39 j=1,istate
      DO 39 i=1,narg
39    v(i,j)=v(i,j)*wt(i)
c      DO i=1,istate
c      WRITE(66,*)"grid positions and grid vstat for state ",i
c      WRITE(66,"(4f10.5)")((grid(izhang,jzhang),jzhang=1,3),
c     &                        v(izhang,i), izhang=1,narg)
c      ENDDO
      DO 9333 j=1,istate
      DO 9333 isum=1,narg
9333  sumpot(j)=sumpot(j)+v(isum,j)
c
c transfer this BLOCK of potential points to a scratch array and buffer
c out to disk.
c
      DO 601 ii=1,istate
      idum=(ii-1)*narg + 1
      CALL scopy(narg,v(1,ii),1,vbuf(idum),1)
601   CONTINUE
      WRITE(10)(vbuf(i),i=1,narg*istate)
c
c RETURN to start a fetch another BLOCK of points
c
      IF(iquit.EQ.0)go to 32
      WRITE(66,991) (iisum,sumpot(iisum),iisum=1,istate)
      WRITE(6,991) (iisum,sumpot(iisum),iisum=1,istate)
991   FORMAT(' sumpot ',i3,'   = ',e16.8)
      CALL closeabs(9)
      STOP
      END
