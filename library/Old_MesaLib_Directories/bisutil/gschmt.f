*deck gschmt.f
      subroutine gschmt(v,thresh,n,nstart,nfin,nout,schmdt,prnt)
c***begin prologue     gschmt
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
c                          nstart = beginning vector
c                          nfin   = ending vector
c                          nout   = number vectors outputted
c                          schmdt = perform either one or two
c                                   orthonormalizations on set
c***references
c***routines called    saxpy(clams), sdot(clams), sscal(clams)
c
c***end prologue       gschmt
      implicit integer (a-z)
c
      real*8 v(n,*)
      real*8 sdot, norm, thresh, ovrlap
      logical schmdt, prnt
      character*80 title
c
      common /io/inp, iout
c     
      if(prnt) then
         title='trial vectors before orthonormalization'
         call prntrm(title,v,n,nfin,n,nfin,iout)
      endif
      ntrial=nstart-1
      nout=0      
      do 10 i=nstart,nfin
         norm=sqrt(sdot(n,v(1,i),1,v(1,i),1))
         if(norm.lt.thresh) go to 10
         call sscal(n,1.d0/norm,v(1,i),1)
         if(ntrial.ne.0) then
            do 20 j=1,ntrial
               ovrlap=sdot(n,v(1,j),1,v(1,i),1)
               call saxpy(n,-ovrlap,v(1,j),1,v(1,i),1)
               norm=sqrt(sdot(n,v(1,i),1,v(1,i),1))
               if (norm.lt.thresh) go to 10
   20       continue
            call sscal(n,1.0d+00/norm,v(1,i),1)
         endif
         ntrial=ntrial+1                
         nout=nout+1
         call copy(v(1,i),v(1,ntrial),n)
   10 continue
      if(schmdt) then  
         do 30 i=nstart,ntrial
            norm=sqrt(sdot(n,v(1,i),1,v(1,i),1))
            call sscal(n,1.d0/norm,v(1,i),1)
            do 40 j=1,i-1
               ovrlap=sdot(n,v(1,j),1,v(1,i),1)
               call saxpy(n,-ovrlap,v(1,j),1,v(1,i),1)
               norm=sqrt(sdot(n,v(1,i),1,v(1,i),1))
   40       continue
            call sscal(n,1.0d+00/norm,v(1,i),1)
   30     continue
      endif
      if(prnt) then
         write(iout,1) nstart, nfin                     
         write(iout,2) nout
         title='trial vectors after orthonormalization'
         call prntrm(title,v,n,nout,n,nout,iout)
      endif 	 
 1    format(/,5x,'orthogonalizing vectors = ',i5,' to vectors = ',i5)
 2    format(/,5x,'number of output vectors = ',i5)          
c
c
      return
      end
