*deck @(#)sint3.f	1.3  1/19/93
      subroutine sint3(iatom,jatom,katom,imax,jmax,kmax,
     #                 nprimi,nprimj,nprimk,
     #                 i1,j1,k1,ex,exaux,c,caux,
     #                 alpha,ainv,xyza,expon,kappa,xyz,
     #                 nv,nprim,npaux,nat,pi)
c
c***begin prologue     sint3
c***date written       920402    (yymmdd) 
c***revision date      yymmdd
c
c***keywords           
c***author             martin, richard (lanl) 
c***source             @(#)sint3.f	1.3   1/19/93
c***purpose            
c   module to form the two-dimensional integrals over primitives
c   for three-center overlap integrals
c***description
c   this routine forms two-dimensional 3-center overlap integrals
c   by the recursion technique described by obara and saika. 
c
c***references
c   obara and saika, j.chem. phys., 84, 3963(1986).
c***routines called
c
c***end prologue       sint3
c
      implicit integer (a-z)
c
c     ----- input arrays (unmodified) -----
      real*8 ex(nprim),exaux(npaux)
      real*8 c(3,nat),caux(3,nat)
c
c     ----- output arrays; two-dimensional 3-center overlap integrals -----
c     the first index spans the primitives on the three centers.
c     the "auxiliary" basis set is denoted by center k.
      real*8 xyz(nv,0:imax,0:jmax,0:kmax,3)
c
c     ----- scratch arrays -----
      real*8 alpha(nv,3),xyza(nv,3),ainv(nv),expon(nv),kappa(nv)
c
c     ----- local variables -----
      real*8 ix,jx,kx,scalri,scalrj,scalrk,one,half,pi
      logical debug
c
      parameter (one=1.0d+00,half=0.5d+00)
      parameter (debug=.false.)
c
      common/io/inp,iout
c
c     ----- determine lengths for primitive exponents -----
      i2=i1+nprimi-1
      j2=j1+nprimj-1
      k2=k1+nprimk-1
c
c     ----- load the primitive exponents into temporary work space -----
      n=0
      do 30 kprim=k1,k2
         kx=exaux(kprim)
         do 20 jprim=j1,j2
            jx=ex(jprim)
            do 10 iprim=i1,i2
               ix=ex(iprim)
               n=n+1
               alpha(n,1)=ix
               alpha(n,2)=jx
               alpha(n,3)=kx
   10       continue
   20    continue
   30 continue
c
      if (n.ne.nv) then
         call lnkerr('m352: sint3 error in vector lengths')
      end if
c
c     ----- form the inverse sum:  1/(ix+jx+kx) -----
      call vadd(ainv,alpha(1,1),alpha(1,2),n)
      call vadd(ainv,ainv,alpha(1,3),n)
      call vinv(ainv,ainv,n)
      if(debug) then
         write(iout,*) ' inverse sum:'
         write(iout,*) (ainv(i),i=1,n)
      endif
c
c     ----- form the natural centers G -----
c           Gx = (ix*xi +jx*xj +kx*xk) / (ix+jx+kx)
      do 40 coord=1,3
         call smul(xyza(1,coord),alpha(1,1),c(coord,iatom),n)
         call saxpy(n,c(coord,jatom),alpha(1,2),1,xyza(1,coord),1)
         call saxpy(n,caux(coord,katom),alpha(1,3),1,xyza(1,coord),1)
         call vmul(xyza(1,coord),xyza(1,coord),ainv,n)
   40 continue
      if(debug) then
         write(iout,*) 'coords:'
         write(iout,*) (c(coord,iatom),c(coord,jatom),
     $                  caux(coord,katom),coord=1,3)
      endif
c
c     ----- form exponential prefactor -----
      scalri=c(1,iatom)*c(1,iatom) +c(2,iatom)*c(2,iatom)
     $      +c(3,iatom)*c(3,iatom) 
      scalrj=c(1,jatom)*c(1,jatom) +c(2,jatom)*c(2,jatom)
     $      +c(3,jatom)*c(3,jatom) 
      scalrk=caux(1,katom)*caux(1,katom) +caux(2,katom)*caux(2,katom)
     $      +caux(3,katom)*caux(3,katom) 
      call smul(expon,alpha(1,1),scalri,n)
      call saxpy(n,scalrj,alpha(1,2),1,expon,1)
      call saxpy(n,scalrk,alpha(1,3),1,expon,1)
      call vmul(expon,expon,ainv,n)
c
c     ----- get norm of G -----
      call rzero(alpha(1,1),n)
      do 50 coord=1,3
         call vwxy(alpha(1,1),alpha(1,1),xyza(1,coord),xyza(1,coord),
     $             1,n)
   50 continue
      if(debug) then
         write(iout,*) 'G2:',(alpha(i,1),i=1,n)
      endif
c
c     ----- form the sss 3-center overlap integrals -----
      call vsub(expon,alpha(1,1),expon,n) 
      call vdiv(expon,expon,ainv,n)
      call vexp(kappa,expon,n)
c
      call smul(alpha(1,1),ainv,pi,n)
      call vsqrt(alpha(1,2),alpha(1,1),n)
      call vmul(alpha(1,1),alpha(1,2),alpha(1,1),n)
      call vmul(kappa,kappa,alpha(1,1),n)
      if(debug) then
         write(iout,*) 'kappa:'
         write(iout,*) (kappa(i),i=1,n)
      endif
c
c     ----- form two-dimensional integrals for higher "l" by recursion --- 
c           begin the x and y 2-d integrals at 1.0
c           the z component contains sss overlap.
c
      call rzero(xyz,n*(imax+1)*(jmax+1)*(kmax+1)*3)
      call vfill(xyz(1,0,0,0,1),one,n)
      call vfill(xyz(1,0,0,0,2),one,n)
      call vmove(xyz(1,0,0,0,3),kappa,n)
      if(debug) then
         write(iout,*) 'sss overlap:'
         write(iout,*) ('x',xyz(i,0,0,0,1),i=1,n)
         write(iout,*) ('y',xyz(i,0,0,0,2),i=1,n)
         write(iout,*) ('z',xyz(i,0,0,0,3),i=1,n)
      endif
c
      do 200 coord=1,3
         call ssub(alpha(1,1),xyza(1,coord),c(coord,iatom),n)
         call ssub(alpha(1,2),xyza(1,coord),c(coord,jatom),n)
         call ssub(alpha(1,3),xyza(1,coord),caux(coord,katom),n)
c
c        ----- bump the angular momentum on center c -----
c              this generates all (0,0,c) integrals.
         li=0
         lj=0
         do 140 lk=0,kmax-1
            call vmul(xyz(1,li,lj,lk+1,coord),alpha(1,3),
     $                xyz(1,li,lj,lk,coord),n)
            if(lk.gt.0) then
               scalrk=float(lk)*half
               call smul(expon,ainv,scalrk,n)
               call vwxy(xyz(1,li,lj,lk+1,coord),
     $                   xyz(1,li,lj,lk+1,coord),
     $                   expon,xyz(1,li,lj,lk-1,coord),+1,n)
            endif
  140    continue
c
c        ----- bump the angular momentum on center b -----
c              this generates all (0,b,c) integrals
         li=0
         do 160 lj=0,jmax-1
            do 150 lk=0,kmax
               call vmul(xyz(1,li,lj+1,lk,coord),alpha(1,2),
     $                   xyz(1,li,lj,lk,coord),n)
               if(lj.gt.0) then
                  scalrj=float(lj)*half
                  call smul(expon,ainv,scalrj,n)
                  call vwxy(xyz(1,li,lj+1,lk,coord),
     $                      xyz(1,li,lj+1,lk,coord),
     $                      expon,xyz(1,li,lj-1,lk,coord),+1,n)
               endif
               if(lk.gt.0) then
                  scalrk=float(lk)*half
                  call smul(expon,ainv,scalrk,n)
                  call vwxy(xyz(1,li,lj+1,lk,coord),
     $                      xyz(1,li,lj+1,lk,coord),
     $                      expon,xyz(1,li,lj,lk-1,coord),+1,n)
               endif
  150       continue
  160    continue
c
c        ----- bump the angular momentum on center a -----
c              this completes the generation of all (a,b,c)
         do 190 li=0,imax-1
            do 180 lj=0,jmax
               do 170 lk=0,kmax
                  call vmul(xyz(1,li+1,lj,lk,coord),alpha(1,1),
     $                                 xyz(1,li,lj,lk,coord),n)
                  if(li.gt.0) then
                     scalri=float(li)*half
                     call smul(expon,ainv,scalri,n)
                     call vwxy(xyz(1,li+1,lj,lk,coord),
     $                         xyz(1,li+1,lj,lk,coord),
     $                         expon,xyz(1,li-1,lj,lk,coord),+1,n)
                  endif
                  if(lj.gt.0) then
                     scalrj=float(lj)*half
                     call smul(expon,ainv,scalrj,n)
                     call vwxy(xyz(1,li+1,lj,lk,coord),
     $                         xyz(1,li+1,lj,lk,coord),
     $                         expon,xyz(1,li,lj-1,lk,coord),+1,n)
                  endif
                  if(lk.gt.0) then
                     scalrk=float(lk)*half
                     call smul(expon,ainv,scalrk,n)
                     call vwxy(xyz(1,li+1,lj,lk,coord),
     $                         xyz(1,li+1,lj,lk,coord),
     $                         expon,xyz(1,li,lj,lk-1,coord),+1,n)
                  endif
  170          continue
  180       continue
  190    continue
  200 continue
c
c
      return
      end
