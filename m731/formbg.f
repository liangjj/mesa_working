*deck @(#)formbg.f	5.1  11/6/94
      subroutine formbg(nz,ianz,iz,bl,alpha,beta,nparm,b,ib,
     $                 g,xm,cz,iscr,scr,dump,toang,det)
c***begin prologue     formbg.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)formbg.f	5.1   11/6/94
c***purpose            
c***description
c     given a z-matrix, this routine will form the wilson b and g
c     matrices.  these may be subsequently used to transform
c     cartesian first derivatives to internal coordinates.
c
c     arguments
c
c     nz     ... number of entries in the z-matrix.
c     ianz   ... vector of length nz containing integer atomic
c                numbers.
c     iz     ... integer connectivity matrix of dimension
c                (nz*4).
c     bl     ... vector of length nz containing bond-lengths.
c     alpha  ... vector of length nz containing first angles.
c     beta   ... vector of length nz containing second angles.
c     nparm  ... number of degrees of freedom (3*nz-6).
c     b      ... output b-matrix (3*4*nparm).
c     ib     ... integer portion of b-matrix (4*nparm).
c     g      ... output g-matrix (nparm*nparm).
c     xm     ... scratch array of length (nz*5).
c     cz     ... scratch array of length (3*nz).
c     scr    ... scratch array of length (4*nparm)
c     iscr   ... integer scratch array of length (nz)
c     dump   ... dump flag.
c     toang  ... bohr to angstrom conversion factor
c     det    ... determinant of the inverse matrix.
c
c
c     "molecular vibrations. the theory of infrared and raman vibrational
c     spectra", e.b.wilson,jr., j.c.decius, and p.c.cross, mcgraw-hill, 1955.
c
c***references
c
c***routines called
c
c***end prologue       formbg.f
      implicit none
c     --- input variables -----
      integer nz,nparm
      logical dump
c     --- input arrays (unmodified) ---
      integer ianz(nz),iz(4,nz)
      real*8 bl(nz),alpha(nz),beta(nz)
c     --- input arrays (scratch) ---
      integer iscr(nparm)
      real*8 xm(nz,5),cz(3*nz),scr(4*nparm)
c     --- output arrays ---
      integer ib(4,nparm)
      real*8 b(3,4,nparm),g(nparm,nparm)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,j,k,l,i1,j1,ij,ibi,iparm
      real*8 toang,det,r,zero,one
c
      parameter (zero=0.0d+00,one=1.0d+00)
c
      common/io/inp,iout
c
 1000 format(' formbg:  ztoc coordinates.')
 1010 format(' ib matrix:')
 1020 format(' b matrix:')
 1030 format(' g matrix:')
 1040 format(' g-inverse matrix:')
c
c     --- print the coordinates?
      if(dump) then
         write(iout,1000)
         call corpr1(iout,nz,ianz,cz,toang)
      endif
c
c     --- prepare xm.
      do 10 i=1,nz
         do 5 j=1,3
            xm(i,j)=one
    5    continue
   10 continue
      xm(1,1)=zero
      xm(1,2)=zero
      xm(1,3)=zero
      xm(2,1)=zero
      xm(2,2)=zero
      xm(3,2)=zero
c
c     --- form the b-matrix.
      call izero(ib,4*nparm)
      call rzero(b,3*4*nparm)
c
c     --- loop over all rows of the z-matrix.
      do 100 i=2,nz
c
c        --- first element is a bond stretch.
         call str(i-1,i,iz(1,i),b,ib,cz)
         if(bl(i).lt.zero) then
            do 20 j=1,3
               do 15 k=1,2
                  b(j,k,i-1)=-b(j,k,i-1)
   15          continue
   20       continue
         endif
c
c        --- angle bend (alpha).
         if(i.gt.2) then
            iparm=nz-3+i
            call bend(iparm,i,iz(1,i),iz(2,i),b,ib,cz)
            if(alpha(i).lt.zero) then
               do 40 j=1,3
                  do 30 k=1,3
                        b(j,k,iparm)=-b(j,k,iparm)
   30             continue
   40          continue
            endif
            if(i.gt.3) then
               iparm=nz+nz-6+i
               if(iz(4,i).eq.0) then
                  call tors(iparm,i,iz(1,i),iz(2,i),iz(3,i),b,ib,cz)
               else
                  call bend(iparm,i,iz(1,i),iz(3,i),b,ib,cz)
                  if(beta(i).lt.zero) then
                     do 60 j=1,3
                        do 50 k=1,3
                           b(j,k,iparm)=-b(j,k,iparm)
   50                   continue
   60                continue
                  endif
               endif
            endif
         endif
  100 continue
c
c     --- apply mask to b.
      do 130 i=1,nparm
         do 120 i1=1,4
            ibi=ib(i1,i)
            if(ibi.ne.0) then
               do 110 l=1,3
                 b(l,i1,i)=b(l,i1,i)*xm(ibi,l)
  110          continue
            endif
  120    continue
  130 continue
c
c     --- possibly print b and ib.
      if(dump) then
         write(iout,1010)
c        --- float ib into scr for printing.
         ij=0
         do 140 j=1,nparm
            do 140 i=1,4
               ij=ij+1
               scr(ij)=float(ib(i,j))
  140    continue
         call matout(scr,4,nparm,4,nparm,iout)
         write(iout,1020)
         call matout(b,12,nparm,12,nparm,iout)
      endif
c
c     --- form g-matrix.
      do 200 i=1,nparm
         do 190 j=1,i
            r=zero
            do 180 i1=1,4
               ibi=ib(i1,i)
               if(ibi.ne.0) then
                  do 170 j1=1,4
                     if(ibi.eq.ib(j1,j)) then
                        do 160 l=1,3
                           r=r+b(l,i1,i)*b(l,j1,j)
  160                   continue
                     endif
  170             continue
               endif
  180       continue
            g(i,j)=r
            g(j,i)=r
  190    continue
  200 continue
c
c     --- print the g-matrix?
      if(dump) then
         write(iout,1030)
         call matout(g,nparm,nparm,nparm,nparm,iout)
      endif
c
c     --- form fi=g(-1)*b*fx .
      call minvrt(g,nparm,nparm,r,iscr,scr)
c
c     --- print g**(-1)?
      if(dump) then
         write(iout,1040)
         call matout(g,nparm,nparm,nparm,nparm,iout)
      endif
c
      det=abs(r)
c
c
      return
      end
