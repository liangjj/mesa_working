*deck mgschm.f
      subroutine mgschm(v,thresh,idrop,n,nstart,nfin,nout,schmdt)
c***begin prologue     mgschm
c***date written       960801  (yymmdd)
c***revision date              (yymmdd)
c
c***keywords           modified gram-schmidt, orthogonalization
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
c                          v(n,*) = input vectors
c                          n      = dimension of vector space
c                          nstart = beginning vector
c                          nfin   = ending vector
c                          nout   = number vectors outputted
c                          thresh = threshold for dropping vectors
c                          schmdt = perform either one or two
c                                   orthonormalizations on set
c                          idrop = drop array based on thresh
c***references
c***routines called    saxpy(clams), sdot(clams), sscal(clams)
c
c***end prologue       gschmt
      implicit integer (a-z)
c
      real*8 v(n,*)
      dimension idrop(*)
      real*8 sdot, norm, ovrlap, thresh
      logical schmdt
c
      common /io/inp, iout
c
      ntimes=1
      if(schmdt) then
         ntimes=2
      endif
      do 1000 trips=1,ntimes
         count=0
         do 10 i=nstart,nfin
            count=count+1
            norm=sqrt(sdot(n,v(1,i),1,v(1,i),1))
            if(norm.ge.thresh) then
               idrop(count)=1
               call sscal(n,1.d0/norm,v(1,i),1)
               do 20 j=i+1,nfin
                  ovrlap=sdot(n,v(1,j),1,v(1,i),1)
                  call saxpy(n,-ovrlap,v(1,j),1,v(1,i),1)
   20          continue
            else
               idrop(count)=0
            endif
   10    continue
 1000 continue   
      count=0
      nout=0
      do 30 i=nstart,nfin
         count=count+1
         if(idrop(count).ne.0) then
            nout=nout+1
            call copy(v(1,i),v(1,nout),n)
         endif
 30   continue   
c
c
      return
      end
