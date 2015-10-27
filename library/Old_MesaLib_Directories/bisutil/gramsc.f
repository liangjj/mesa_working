*deck gramsc.f
      subroutine gramsc(v,thresh,n,nin,nout,prnt)
c***begin prologue     gramsc
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
c                          v(n,*) = input vectors
c                          thresh = acceptance tolerance on overlap
c                          n      = dimension of vector space
c                          nin    = number of input vectors
c                          nout   = number output vectors
c***references
c***routines called    saxpy(clams), sdot(clams), sscal(clams)
c
c***end prologue       gramsc
      implicit integer (a-z)
c
      real*8 v(n,*)
      real*8 sdot, norm, thresh, ovrlap
      logical prnt
      character*80 title
c
      common /io/inp, iout
c     
      nout=0
      do 10 i=1,nin
         norm=sqrt(sdot(n,v(1,i),1,v(1,i),1))
         if(norm.lt.thresh) go to 10
         call sscal(n,1.d0/norm,v(1,i),1)
         if(nout.ne.0) then
            do 20 j=1,nout
               ovrlap=sdot(n,v(1,j),1,v(1,i),1)
               call saxpy(n,-ovrlap,v(1,j),1,v(1,i),1)
               norm=sqrt(sdot(n,v(1,i),1,v(1,i),1))
               if (norm.lt.thresh) go to 10
   20       continue
            call sscal(n,1.0d+00/norm,v(1,i),1)
         endif
         nout=nout+1
         call copy(v(1,i),v(1,nout),n)
   10 continue
c
      if(prnt) then
         write(iout,1) nin, nout
      endif
      if(nout.eq.0) then
         call lnkerr('no vectors from gramsc')
      endif
      return
 1    format(/,1x,'number of vectors in to schmidt '
     1            'orthogonalization = ',i4,/,1x,
     2            'number of vectors out after schmidt '
     3            'orthogonalization = ',i4)

      end
