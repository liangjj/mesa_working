*deck @(#)maket.f	5.1  11/6/94
      subroutine maket(t,npf,nbf,ptcont,natoms,nbtype,noprim,nocont,
     $                 cf,numcf,nocart,minmom,maxmom,start,pstart)
c
c***begin prologue     maket.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe, paul(lanl)
c***source             @(#)maket.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       maket.f
      implicit none
c     --- input variables -----
      integer natoms,nbtype,npf,nbf,numcf
c     --- input arrays (unmodified) ---
      integer ptcont(natoms,nbtype),noprim(natoms,nbtype)
      integer nocont(natoms,nbtype),nocart(nbtype),minmom(nbtype)
      integer maxmom(nbtype),start(natoms,nbtype),pstart(natoms,nbtype)
      real*8 t(npf,nbf),cf(numcf)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer atom,type
c
c     --- zero the t matrix ---
      call rzero(t,npf*nbf)
c
c     --- loop through atoms and types, filling in the transformation
      do 100 atom=1,natoms
         do 90 type=1,nbtype
            if (noprim(atom,type).gt.0) then
               call maket1(t,npf,nbf,noprim(atom,type),
     $                     nocont(atom,type),cf(ptcont(atom,type)),
     $                     minmom(type),maxmom(type),nocart,
     $                     start(atom,type),pstart(atom,type))
            end if
   90    continue
  100 continue
c
c
      return
      end
