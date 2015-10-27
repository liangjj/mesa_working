*deck @(#)gdmat.f	1.1  4/25/95
      subroutine gdmat(d,c,nbf,nnp,mink,maxk)
c***begin prologue     gdmat.f
c***date written       840906  (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe,paul (lanl) 
c***source             @(#)gdmat.f	1.1   4/25/95
c***purpose            forms the density matrices.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       gdmat.f
      implicit none
c     --- input variables -----
      integer nbf,nnp,mink,maxk
c     --- input arrays (unmodified) ---
      real*8 c(nbf,nbf)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 d(nnp)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i,j,k,ij
c
c
      ij=0
      call rzero(d,nnp)
      do 5 i=1,nbf
         do 4 j=1,i
            ij=ij+1
               do 1 k=mink,maxk
                  d(ij)=d(ij)+c(i,k)*c(j,k)
    1          continue
    4    continue
    5 continue
c
c
      return
      end
