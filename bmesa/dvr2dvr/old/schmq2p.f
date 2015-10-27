*deck schmq2p.f
      subroutine schmq2p(vqq,vpp,thresh,n,nqin,nqout,npin,npout,
     1                   schmdt,sets)
c***begin prologue     schma2b
c***date written       960801  (yymmdd)
c***revision date              (yymmdd)
c
c***keywords           gram-schmidt, orthogonalization
c***author             barry schneider(nsf)
c***source
c***purpose            gram-schmidt orthogonalization.
c***description        a set of non-orthonormal vectors are input and 
c                      orthonormalized using a gram-schmidt process that 
c                      checks for linear dependencies. the routine outputs 
c                      the orthonormal functions in the same locations as the
c                      input functions.  the number of output functions is
c                      less than or equal to the number of input functions. 
c
c                      vpp(n,*) = input vectors.  assumed to orthonormal
c                                                 in p-space
c                      vqq(n,*) = input vectors.  assumed to orthonormal
c                                                 in q-space    
c                      thresh   = acceptance tolerance on overlap
c                      n        = length of vector space
c                      np       = number of p-space vectors
c                      nq       = number q-space vectors input
c                      nout     = number q-space vectors outputted
c                      schmdt   = perform either one or two
c                                 orthonormalizations on set
c***references
c***routines called    saxpy(clams), sdot(clams), sscal(clams)
c
c***end prologue       schma2b
      implicit integer (a-z)
c
      real*8 vpp(n,*), vqq(n,*)
      real*8 sdot, norm, thresh, ovrlap
      logical schmdt
      character*(*) sets
c
      common /io/inp, iout
c
      if(sets.eq.'same') then
         call gschmt(vpp,thresh,n,1,npin,npout,schmdt)
      else         

c     schmidt orthogonalize the bigger set to the smaller set
         ntrial=0
         do 10 i=1,nqout
            norm=sqrt(sdot(n,vqq(1,i),1,vqq(1,i),1))
            if(norm.lt.thresh) go to 10
            call sscal(n,1.d0/norm,vqq(1,i),1)
            do 20 j=1,npout
               ovrlap=sdot(n,vpp(1,j),1,vqq(1,i),1)
               call saxpy(n,-ovrlap,vpp(1,j),1,vqq(1,i),1)
               norm=sqrt(sdot(n,vqq(1,i),1,vqq(1,i),1))
               if (norm.lt.thresh) go to 10
   20       continue
            call sscal(n,1.0d+00/norm,vqq(1,i),1)
            ntrial=ntrial+1                
            call copy(vqq(1,i),vqq(1,ntrial),n)
   10    continue
c        now orthonormalize the remaining vectors   
         call gschmt(vqq,thresh,n,1,ntrial,nqout,schmdt)
      endif               
c
      return
      end
