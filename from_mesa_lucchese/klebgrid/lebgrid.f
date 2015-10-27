      program lebgrid
c     
c     new grid program to generate 3-D grid using Becke's scheme.
c     This version uses lebedev's quadrature of order 17
c     for the angular quadrature
c
      parameter (mxcen=20,npmx=1000000,mxshell=30)
      parameter (mxlebpts=6000)
      implicit real*8 (a-h,o-z)
      character*8 itheta, iphi
      dimension cth(mxlebpts),sth(mxlebpts),cph(mxlebpts)
      dimension sph(mxlebpts),ww(mxlebpts)
      dimension dummy(2),scr(64),x(64),w(64),rad(mxcen)
      dimension a(3,mxcen),rr(mxshell),nr(mxshell)
      dimension grid(4,npmx,mxcen),iatom(mxcen),arad(mxcen,mxcen)
      dimension gg(4,npmx*mxcen),it(3),rit(3)
      dimension r(mxcen,mxcen),d(mxcen),p(npmx),wtt(mxcen)
      data pi/3.14159265358/,miter/3/
      open(5,file='ingrid',status='unknown')
      open(6,file='outgrid',status='unknown')
      read(5,*)nbuff
      read(5,*)ncent,miter,irad
      read(5,*)lmaxk
c     
c if irad .ne. zero, atomic radii will be read in and used to weight polyhedra
c
      IF (ncent .gt. mxcen) THEN
         WRITE (6,*) 'Error ncent > mxcen'
         stop 'Error Bad dimension'
      end if

      read(5,*)((a(iii,i),iii=1,3),i=1,ncent)
 100  format(3f10.5)
      write(6,*)' atomic centers'
      write(6,100)((a(iii,i),iii=1,3),i=1,ncent)
      if(irad.ne.0)then
         write(6,*)' atomic radii'
         read(5,*)(rad(i),i=1,ncent)
         write(6,103)(rad(i),i=1,ncent)
 103     format(10f10.5)
      endif
c
c open a loop on centers
c
      itotal=0
      do 1 i=1,ncent
      read(5,*)nshell,wtt(i)
      ngr=1+nshell
      if(wtt(i).gt.1.d-8)then
      read(5,*)(rr(j),j=1,ngr)
      read(5,*)(nr(j),j=1,nshell)
      read(5,110)itheta,iphi
 110  format(2a8)
      wt=wtt(i)
      write(6,*)
      write(6,101)i,(rr(j),j=1,ngr)
 101  format(//' atom',i3/' radial shells:',(10f8.3))
      write(6,102)(nr(j),j=1,nshell)
 102  format(' radial quad. orders:',(10i5))
      write(6,110)itheta,iphi
      write(6,104)wt
 104  format(' the weight for this atom is ',f10.5)
      call ldpw(cth,sth,cph,sph,ww,max(2*LMaxK+3, 17),nleb)
      if(nleb.gt.mxlebpts)then
         write(6,*)'Error mxlebpts not large enough', nleb, mxlebpts
         stop 'Error stopping bad dimension'
      end if
      write(6,*)'Using ', nleb, ' Lebedev points'
c      write(*,*)"sth"
c      write(*,*)sth
c      write(*,*)"cth"
c      write(*,*)cth
c      write(*,*)"sph"
c      write(*,*)sph
c      write(*,*)"cph"
c      write(*,*)cph
c      write(*,*)"za_leb"
c      do jj = 1,nleb
c         write(*,*)jj,sth(jj)*cph(jj),sth(jj)*sph(jj),cth(jj)
c      enddo
c
c loop on radial points
c
      irun=0
      do 2 j=1,nshell
       if(j.gt.1.and.nr(j).eq.nr(j-1))go to 11
       nrr=nr(j)
       call gaussq(1,nrr,0.,0.,0,dummy,scr,x,w)
c       write(*,*)"unweighted"
c       do index = 1,nrr
c          write(*,"(5x, 2f16.7)")x(index),w(index)
c       enddo
 11    continue
       ampr=(rr(j+1)-rr(j))/2.
       abpr=(rr(j+1)+rr(j))/2.
c       write(*,*)"weighted"
       do 52 l=1,nrr
        zz=ampr*x(l)+abpr
        wz=zz*zz*ampr*w(l)
c        write(*,"(5x, 2f16.7)")zz,wz
c
c loop on angular points
c
c        write(*,*)"before transform"
        do 51 jj=1,nleb
          isub=irun+jj+nleb*(l-1)
          grid(1,isub,i)=zz*sth(jj)*cph(jj)+a(1,i)
          grid(2,isub,i)=zz*sth(jj)*sph(jj)+a(2,i)
          grid(3,isub,i)=zz*cth(jj)+a(3,i)
          grid(4,isub,i)=ww(jj)*wz*wt
c          write(*,"(i5,4f12.8)")isub, (grid(igrid, isub,i), igrid=1,4)
 51      continue
   52   continue
        irun=irun+nleb*nrr
        if(irun.gt.npmx)then
         write(6,*)' Error npmx exceeded',irun, npmx
         stop 'Error stopping npmx exceeded'
        endif
 2     continue
       itotal=itotal+irun
       iatom(i)=irun
       write(6,666)irun
 666   format(' number of points on this atom is',i7)
c
c compute sum of quadrature weights
c
        sum=0.
        do 53 ii=1,irun
         sum=sum+grid(4,ii,i)
 53     continue
        exact=wtt(i)*4./3.*pi*rr(ngr)**3
        write(6,200)sum,exact
 200    format(' sum of quadrature weights for this atom is '
     $  ,d16.8/' exact answer is ',d16.8)
      else
      write(6,223)i
223   format(//' atom',i3,' has no weight'/)      
      endif
 1    continue
c
c perform voronoi transformation
c
      call voronoi(miter,ncent,a,iatom,grid,npmx,d,r,p,wtt,irad
     $  ,rad,arad)
c
c write grid to disk
c
      ii=0
      do 443 i=1,ncent
         if(wtt(i).lt.1.d-8)go to 443
         npt=iatom(i)
c         write(*,*)"center, points", i, npts
         do 4421 j=1,npt
            ij=j+ii
            do 442 k=1,4
               gg(k,ij)=grid(k,j,i)  
 442        continue
c            write(*,"(4f10.5)")(gg(k,ij), k=1,4)
 4421    continue
         ii=ii+iatom(i)
 443  continue
      irecl=4*nbuff*8
      iwrit=nbuff*4
      call openabs(7,'grid',irecl)
      it(1)=nbuff
      it(2)=irecl
      it(3)=itotal
      rit(1)=it(1)
      rit(2)=it(1)
      rit(3)=it(3)
      write(6,*)it(1),it(2),it(3)
c
      call wrabs(7,rit,3,0)
      irectot=4*itotal*8
      ipass=irectot/irecl
      if(irecl*ipass.ne.irectot)ipass=ipass+1
      ij=1
      do 555 i=1,ipass
        call wrabs(7,gg(1,ij),iwrit,i)
        ij=ij+nbuff
555     continue
      write(6,444)itotal
444   format(' total number of grid points is ',i7)    
      call closeabs(7)
      stop
      end
