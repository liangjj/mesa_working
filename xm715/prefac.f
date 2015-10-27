*deck %W%  %G%
      subroutine prefac(ar,a,nij,br,b,nkl,expon,index,ijindx,klindx,
     $                  kl,len,lenv,thresh,pi252,t1)
c***begin prologue     prefac.f
c***date written       931215  
c***revision date      11/6/94
c
c   18 december, 1993  rlm at lanl
c      modifying for use in direct j-matrix
c***keywords           
c***author             lengsfield, byron 
c***source             %W%   %G%
c***purpose            eliminates integrals from the list of those to be
c                      computed based on a test of the prefactor.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       prefac.f
      implicit none
c     --- input variables -----
      integer nij,nkl,len
      real*8 thresh,pi252
c     --- input arrays (unmodified) ---
      integer ijindx(nij,2),klindx(nkl,2)
      real*8 ar(nij),a(nij),br(nkl),b(nkl),expon(len)
c     --- input arrays (scratch) ---
      real*8 t1(nij)
c     --- output arrays ---
      integer index(len,6)
c     --- output variables ---
      integer kl,lenv
c     --- scratch arrays ---
c     --- local variables ---
      integer ij,k,l,lenvsv
      integer inp,iout
      real*8 brkl,bkl
c
      common/io/inp,iout
c
c     --- loop over the primitives in the kl block ---
      lenv=0
    1 continue
         kl=kl+1
         brkl=br(kl)
         bkl=b(kl)
c        --- compute prefactor for all ij primitives with this
c            kl pair and test against the threshhold ---
         do 2 ij=1,nij
            t1(ij)=exp(-ar(ij)-brkl)*pi252/(a(ij)*bkl*sqrt(a(ij)+bkl))
    2    continue
c
         lenvsv=lenv
         do 3 ij=1,nij
            if (abs(t1(ij)).gt.thresh) then
               lenv=lenv+1
               expon(lenv)=t1(ij)
               index(lenv,5)=ij
            end if
    3    continue
c
c        --- save the indices for those we must compute ---
         k=klindx(kl,1)
         l=klindx(kl,2)
         do 4 ij=lenvsv+1,lenv
            index(ij,3)=k
            index(ij,4)=l
            index(ij,6)=kl
    4    continue
c
c        --- if there is room for another full ij block and we have not
c            screened all the kl primitives, then do it again.
      if (len-lenv.ge.nij.and.kl.lt.nkl) go to 1
c
c     --- fill in the i and j indices in index ---
      do 5 ij=1,lenv
         index(ij,1)=ijindx(index(ij,5),1)
         index(ij,2)=ijindx(index(ij,5),2)
    5 continue
c
c
      return
      end
