*deck @(#)putf.f	5.1  11/6/94
      subroutine putf(nz,lbl,lalpha,lbeta,nparm,nvar,f,force,dump)
c***begin prologue     putf.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)putf.f	5.1   11/6/94
c***purpose            
c***description
c   routine to add force contributions over internal coordinates
c   into the vector over actual variables.
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       putf.f
      implicit none
c     --- input variables -----
      integer nz,nparm,nvar
      logical dump
c     --- input arrays (unmodified) ---
      integer lbl(nz),lalpha(nz),lbeta(nz)
      real*8 f(nparm)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 force(nvar)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,j,ibl,ialpha,ibeta
      real*8 dx
c
      common/io/inp,iout
c
 1000 format(' putf: contents of force.')
 1010 format(1x,i3,e20.10)
c
c     --- initialize the force array.
      call rzero(force,nvar)
      if(nz.ge.2) then
         do 10 i=2,nz
            ibl=abs(lbl(i))
            if(ibl.ne.0) then
                  dx=f(i-1)
                  if(lbl(i).lt.0) dx=-dx
                  force(ibl)=force(ibl)+dx
            endif
   10    continue
c
         if(nz.ge.3) then
            j=nz-3
            do 20 i=3,nz
               ialpha=abs(lalpha(i))
               if(ialpha.ne.0) then
                     dx=f(i+j)
                     if(lalpha(i).lt.0) dx=-dx
                     force(ialpha)=force(ialpha)+dx
               endif
   20       continue
c
            if(nz.ge.4) then
               j=nz+nz-6
               do 30 i=4,nz
                  ibeta=abs(lbeta(i))
                  if(ibeta.ne.0) then
                        dx=f(i+j)
                        if(lbeta(i).lt.0) dx=-dx
                        force(ibeta)=force(ibeta)+dx
                  endif
   30          continue
            endif
         endif
      endif
c
c     --- possibly print the internal forces.
      if(dump) then
         write(iout,1000)
         write(iout,1010) (i,force(i),i=1,nvar)
      endif
c
c
      return
      end
