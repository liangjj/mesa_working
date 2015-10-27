*deck %W% %G% 
      subroutine twoint(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                 i1,j1,ex,alpha,c,ainv,xyza,expon,xyz,
     #                 nv,nprim,nat,nbtype,mmax,kappa,t,f,pi252)
c***begin prologue     %M%
c***date written       920427 (yymmdd) 
c***revision date      %G%
c
c***keywords           
c***author             martin, richard (lanl) 
c***source             %W% %G%
c***purpose            
c      module to form the two-dimensional integrals over primitives
c      for the overlap integrals
c***description
c     
c    
c
c***references
c      obara and saika, j.chem,phys., 84,3963(1986).
c
c***routines called
c
c***end prologue       %M%
      implicit integer (a-z)
c
c     input arrays(unmodified)
      real*8 ex(nprim),c(3,nat)
c
c     output arrays
c     the actual two-electron integrals are accumulated in xyz(,,,,0)
      real*8 xyz(nv,0:imax,0:jmax,3,0:mmax)
c
c     scratch arrays
      real*8 temp(nv,0:imax,0:jmax+1)
      real*8 alpha(nv,2),ainv(nv,5),xyza(nv,3),expon(nv)
      real*8 kappa(nv),t(nv),f(nv,0:mmax)
      real*8 ix,jx,scalar,scalri,scalrj
      real*8 tmp(0:16)
c
      real*8 pi252,half,one,two
      logical debug
      parameter (half=0.5d+00,one=1.0d+00,two=2.0d+00,debug=.true.)
c
      common/io/inp,iout
c
cccccc right now this uses too much core.  and use the alpha arrays for kappa,etc.
      i2=i1+nprimi-1
      j2=j1+nprimj-1
c
      n=0
      do 20 jprim=j1,j2
         jx=ex(jprim)
         do 10 iprim=i1,i2
            ix=ex(iprim)
            n=n+1
            alpha(n,1)=ix
            alpha(n,2)=jx
   10    continue
   20 continue
      if(debug) then
         write(iout,*) 'alpha'
         write(iout,*) (alpha(i,1),i=1,n)
         write(iout,*) (alpha(i,2),i=1,n)
      endif
c
      if (n.ne.nv) then
         call lnkerr('twoint: error in vector lengths')
      end if
c
      call vinv(ainv(1,1),alpha(1,1),n)
      call vinv(ainv(1,2),alpha(1,2),n)
      call vadd(ainv(1,3),alpha(1,1),alpha(1,2),n)
      call vinv(ainv(1,3),ainv(1,3),n)
      call vmul(ainv(1,4),alpha(1,1),ainv(1,3),n)
      call vmul(ainv(1,5),alpha(1,2),ainv(1,3),n)
      write(iout,*) 'ainv'
      write(iout,*) (ainv(ij,1),ij=1,n)
      write(iout,*) (ainv(ij,2),ij=1,n)
      write(iout,*) (ainv(ij,3),ij=1,n)
      write(iout,*) (ainv(ij,4),ij=1,n)
      write(iout,*) (ainv(ij,5),ij=1,n)
c
c     ----- form the natural center W, (ix*xi + jx*xj)/(ix+jx) -----
c
      do 30 coord=1,3
         call smul(xyza(1,coord),alpha(1,1),c(coord,iatom),n)
         call saxpy(n,c(coord,jatom),alpha(1,2),1,xyza(1,coord),1)
         call vmul(xyza(1,coord),xyza(1,coord),ainv(1,3),n)
   30 continue
      write(iout,*) 'natural center'
      write(iout,*) (xyza(ij,1),ij=1,n)
      write(iout,*) (xyza(ij,2),ij=1,n)
      write(iout,*) (xyza(ij,3),ij=1,n)
c
c     ----- form exponential parameter T -----
c
      call vmul(expon,alpha(1,1),alpha(1,2),n)
      call vmul(t,expon,ainv(1,3),n)
      scalar=((c(1,iatom)-c(1,jatom))**2+(c(2,iatom)-c(2,jatom))**2+
     $         (c(3,iatom)-c(3,jatom))**2)
      call smul(t,t,scalar,n)
      if(debug) then
         write(iout,*) 't:'
         call matout(t,nprimi,nprimj,nprimi,nprimj,iout)
      endif
      write(iout,*) 'ainv3',(ainv(ij,3),ij=1,n)
c
c     ----- form exponential prefactor -----
      call vinv(kappa,expon,n)
      call vsqrt(expon,ainv(1,3),n)
      call vmul(kappa,kappa,expon,n)
      call smul(kappa,kappa,pi252,n)
      if(debug) then
         write(iout,*) 'kappa:'
         write(iout,*) (kappa(ij),ij=1,n)
      endif
c
c     ----- form two-dimensional integrals -----
c
      call rzero(xyz,n*(imax+1)*(jmax+1)*3*(mmax+1))
c
c     ----- get the f(m) table for this primitive block. -----
c
c     the dimensioning here is awkward.
      do 50 i=1,n
         call fmoft(tmp,mmax,t(i))
         do 40 m=0,mmax
            f(i,m)=tmp(m)
   40    continue
   50 continue
      if(debug) then
         write(iout,*) 'f(m,t):'
         do 60 m=0,mmax
            write(iout,*) 'm:',m
            call matout(f(1,m),nprimi,nprimj,nprimi,nprimj,iout)
   60    continue
      endif
c
c     ----- form necessary (s|s)(m) integrals to begin recursion
c           begin the x and y 2-d integrals at 1.0
c           the z component has the prefactor. 
c
      do 70 m=0,mmax
         call vmul(t,kappa,f(1,m),n)
         call vfill(xyz(1,0,0,1,m),one,n)
         call vfill(xyz(1,0,0,2,m),one,n)
         call vmove(xyz(1,0,0,3,m),t,n)
   70 continue
         if(debug) then
            write(iout,*) 'sss 2e integrals'
            write(iout,*) (xyz(ij,0,0,3,0),ij=1,n)
         endif
c
c     ----- begin recursion. -----     
      call loops(imax,jmax)
c
      do 210 coord=1,3
c
c        ----- bump the angular momentum on center a -----
c              this generates all (a,0) integrals
         call ssub(expon,xyza(1,coord),c(coord,iatom),n)
         li=0
         do 170 lj=0,jmax
            do 160 m=0,mmax-li-lj-1
               call vwxy(temp(1,li,lj+1),temp(1,li,lj+1),
     $                   expon,temp(1,li,lj),+1,n)
               if(lj.gt.0) then
                  scalrj=float(j)*half
                  call smul(t,ainv(1,2),scalrj,n)
                  call vwxy(temp(1,li,lj+1),
     $                      temp(1,li,lj+1),
     $                      t,temp(1,li,lj-1),+1,n)
                  call vmul(t,t,ainv(nv,4),n)
                  call vwxy(temp(1,li,lj+1),
     $                      temp(1,li,lj+1),
     $                      t,temp(1,li-1,lj),-1,n)
               endif
  160       continue
  170    continue
c
c        ----- bump the angular momentum on center b -----
c              this generates all (a,b) integrals
         call ssub(expon,xyza(1,coord),c(coord,jatom),n)
         do 200 li=0,imax-1
            do 190 lj=0,jmax+1
c
c     ****************         i got to about her and quit
                  call vwxy(temp(1,li+1,lj),
     $                      terp(1,li+1,lj),
     $                      expon,temp(1,li,lj+1),+1,n)
                  if(lj.gt.0) then
                     scalrj=float(j)*half
                     call smul(t,ainv(1,2),scalrj,n)
                     call vwxy(temp(1,li,lj+1,coord,m),
     $                         temp(1,li,lj+1,coord,m),
     $                         t,temp(1,li,lj-1,coord,m),+1,n)
                     call vmul(t,t,ainv(nv,5),n)
                     call vwxy(xyz(1,li,lj+1,coord,m),
     $                         xyz(1,li,lj+1,coord,m),
     $                         t,xyz(1,li,lj-1,coord,m+1),-1,n)
                  endif
                  if(li.gt.0) then
                     scalri=float(i)*half
                     call smul(t,ainv(1,3),scalri,n)
                     call vwxy(xyz(1,li,lj+1,coord,m),
     $                         xyz(1,li,lj+1,coord,m),
     $                         t,xyz(1,li-1,lj,coord,m+1),+1,n)
                  endif
  180          continue
  190       continue
  200    continue
  210 continue
      if(debug) then
         write(iout,*) 'twoint:'
         do 230 li=0,imax
            do 220 lj=0,jmax
               write(iout,*) 'li,lj:',li,lj
               write(iout,*) 'x:'
               write(iout,*) (xyz(i,li,lj,1,0),i=1,n)
               write(iout,*) 'y:'
               write(iout,*) (xyz(i,li,lj,2,0),i=1,n)
               write(iout,*) 'z:'
               write(iout,*) (xyz(i,li,lj,3,0),i=1,n)
  220       continue
  230    continue
      endif
c
      return
      end
