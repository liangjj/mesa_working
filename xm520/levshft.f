*deck @(#)levshft.f	5.1 11/6/94
      subroutine levshft(f,nnp,nbf,maxocc,shift)
c***begin prologue     levshft.f
c***date written       930515   (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           level-shift 
c***author             martin, richard (lanl) 
c***source             @(#)levshft.f	5.1   11/6/94
c***purpose            adds a constant to the diagonal elements of
c                      virtual block of the fock matrix
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       levshft.f
      implicit none
c     --- input variables -----
      integer nnp,maxocc,nbf
      real*8 shift
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 f(nnp)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer ii,i
c
      common/io/inp,iout
c
c
c     --- add the level shift to the diagonal elements of the
c         virtual block of the fock matrix
      ii=maxocc*(maxocc+1)/2
      do 10 i=maxocc+1,nbf
         ii=ii+i
         f(ii)=f(ii)+shift
   10 continue
c
c
      return
      end
