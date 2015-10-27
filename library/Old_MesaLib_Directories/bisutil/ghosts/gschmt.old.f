*deck gschmt.f
      subroutine gschmt(va,vb,thresh,idrop,n,na,nb,nout,schmdt,type)
c***begin prologue     gschmt
c***date written       960801  (yymmdd)
c***revision date              (yymmdd)
c
c***keywords           gram-schmidt, orthogonalization
c***author             barry schneider(nsf)
c***source
c***purpose            gram-schmidt orthogonalization of two sets.
c***description        a set of orthonormal vectors, va, and a set of
c                      non-orthonormal vectors, vb, are input.  the set vb
c                      is schmidt orthogonalized to va and then optionally
c                      internally orthonormalized using a gram-schmidt 
c                      process that checks for linear dependencies. 
c                      the routine outputs the new orthonormal functions 
c                      vb in the same locations as the input functions.
c                      the number of output functions is less than or 
c                      equal to nb the number of input functions. 
c
c                          va(n,na) = input vector set a
c                          vb(n,nb) = input vector set b
c                          thresh = acceptance tolerance on overlap
c                          idrop  = scratch integer array
c                          n      = dimension of vector space
c                          na     = size of vector set a
c                          nb     = size of vector set b
c                          nout   = number vectors outputted
c                          type   = gram-schmidt, modified-gram-schmidt or
c                                   simple projection of the space spanned
c                                   by vector set a out of vector set b.
c                          schmdt = perform either one or two
c                                   orthonormalizations on set
c***references
c***routines called    saxpy(clams), sdot(clams), sscal(clams)
c
c***end prologue       gschmt
      implicit integer (a-z)
c
      real*8 va(n,na), vb(n,nb)
      dimension idrop(nb)
      real*8 sdot, norm, thresh, ovrlap
      character*(*) type
      logical schmdt
c
      common /io/inp, iout
c
c     orthogonalize the set vb to the set va
c
      write(iout,11) type
      ntimes=1
      if(schmdt) then
         ntimes=2
      endif   
      write(iout,1) na, nb
      nout0=nb  
      if(na.ne.0) then
         do 1000 trips=1,ntimes           
            count=0
            do 10 i=1,nout0
               do 20 j=1,na
                  ovrlap=sdot(n,va(1,j),1,vb(1,i),1)
                  call saxpy(n,-ovrlap,va(1,j),1,vb(1,i),1)
 20            continue
               norm=sqrt(sdot(n,vb(1,i),1,vb(1,i),1))
               if(norm.ge.thresh) then
                  count=count+1
                  call sscal(n,1.0d+00/norm,vb(1,i),1)
                  call copy(vb(1,i),vb(1,count),n)
               endif
 10         continue
            nout0=count 
 1000    continue
      endif
      write(iout,2) nout0
      if(type.eq.'orthogonalize set a to set b') then
         return
      endif            
c
c     orthonormalize the set vb using either gram-schmidt
c                           or
c                 modified gram-schmidt
c
      write(iout,11) type
 11   format(a80)
      nout1=nout0 
      if(type.eq.'gram-schmidt') then
         do 2000 trips=1,ntimes
            nout=0
            write(iout,*) nout, nout1 
            do 30 i=1,nout1
               norm=sqrt(sdot(n,vb(1,i),1,vb(1,i),1))
               if(norm.lt.thresh) go to 30
               call sscal(n,1.d0/norm,vb(1,i),1)
               if(nout.ne.0) then
                  do 40 j=1,nout
                     ovrlap=sdot(n,vb(1,j),1,vb(1,i),1)
                     call saxpy(n,-ovrlap,vb(1,j),1,vb(1,i),1)
                     norm=sqrt(sdot(n,vb(1,i),1,vb(1,i),1))
                     if (norm.lt.thresh) go to 30
   40             continue
                  call sscal(n,1.0d+00/norm,vb(1,i),1)
               endif
               nout=nout+1
               call copy(vb(1,i),vb(1,nout),n)
   30       continue
            nout1=nout
 2000    continue                
      elseif(type.eq.'modified-gram-schmidt') then 
         do 3000 trips=1,ntimes
            do 50 i=1,nout1
               norm=sqrt(sdot(n,vb(1,i),1,vb(1,i),1))
               if(norm.ge.thresh) then
                  idrop(i)=1
                  call sscal(n,1.d0/norm,vb(1,i),1)
                  do 60 j=i+1,nout1
                     ovrlap=sdot(n,vb(1,j),1,vb(1,i),1)
                     call saxpy(n,-ovrlap,vb(1,j),1,vb(1,i),1)
   60             continue
               else
                  idrop(i)=0
               endif
   50       continue
 3000    continue   
         nout=0
         do 70 i=1,nout1
            if(idrop(i).ne.0) then
               nout=nout+1
               call copy(vb(1,i),vb(1,nout),n)
            endif
 70      continue
      endif            
      write(iout,3) nout                          
c
c
      return
 1    format(/,1x,'schmidt orthogonalization of two sets of vectors',
     1      /,1x,
     2            'set a already assumed to be orthonormal',/,1x, 
     3            'number of vectors in set a = ',i4,/,1x,
     4            'number of vectors in set b = ',i4)
 2    format(/,1x,'number of set b vectors after orthogonalization '
     1            'to set a vectors = ',i4)
 3    format(/,1x,'final number of set b vectors = ',i4)      
      end
