*deck @(#)get1dm.f	1.1  4/25/95
      subroutine get1dm(d,nnp,ndmat,dij,ni,nj,nfi,nfj,is,js)
c***begin prologue     get1dm.f
c***date written       870702   (yymmdd)  
c***revision date      11/6/94      
c
c   december 8, 1993   rlm at lanl
c     modifying for use in direct scf code.
c***keywords           density matrices, direct
c***author             saxe, paul (lanl)
c***source             @(#)get1dm.f	1.1   4/25/95
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       get1dm.f
      implicit none
c     --- input variables -----
      integer nnp,ndmat
      integer ni,nj
      integer nfi,nfj
      integer is,js
c     --- input arrays (unmodified) ---
      real*8 d(nnp,ndmat)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 dij(ni,nj,nfi,nfj,ndmat)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i,j,if,jf,jpos,ipos
      integer ii,ij,jj,dmat
c
c
c     -- extract just the local ij portion of the full density ---
      do 50 dmat=1,ndmat
         do 40 jf=1,nfj
            jpos=js+jf-nfj
            do 30 if=1,nfi
               ipos=is+if-nfi
               do 20 j=1,nj
                  jj=jpos+j*nfj
                  do 10 i=1,ni
                     ii=ipos+i*nfi
                     ij=max(ii,jj)*(max(ii,jj)-1)/2+min(ii,jj)
                     dij(i,j,if,jf,dmat)=d(ij,dmat)
 10               continue
 20            continue
 30         continue
 40      continue
 50   continue
c
c
      return
      end
