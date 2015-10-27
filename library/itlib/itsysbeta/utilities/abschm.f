*deck abschm.f
      subroutine abschm(va,vb,thresh,n,na,nb,nout,schmdt,prnt)
c***begin prologue     abschm
c***date written       960801  (yymmdd)
c***revision date              (yymmdd)
c
c***keywords           gram-schmidt, orthogonalization
c***author             barry schneider(nsf)
c***source
c***purpose            gram-schmidt orthogonalization.
c***description        a set of non-orthonormal vectors, vb, are orthogonalized 
c                      to another set of vectors, va, using a gram-schmidt process 
c                      that checks for linear dependencies. the set of vectors,
c                      va, are already assumed to be internally orthonormal. 
c
c                          va(n,*) = input vectors
c                          vb(n,*) = input as non-orthogonal vectors and output
c                                    as orthonormal set.
c                          thresh  = acceptance tolerance on overlap
c                          n       = dimension of vector space
c                          nstart  = beginning vector
c                          nfin    = ending vector
c                          nout    = number vectors outputted
c                          schmdt  = perform either one or two
c                                    orthonormalizations on set
c***references
c***routines called    saxpy(clams), sdot(clams), sscal(clams)
c
c***end prologue       gschmt
      implicit integer (a-z)
c
      real*8 va(n,na), vb(n,nb)
      real*8 ddot, norm, thresh, ovrlap
      logical schmdt, prnt
c
      common /io/inp, iout
c

      ntimes=1
      if(schmdt) then
         ntimes=2
      endif
      upper=nb
      do 1000 trips=1,ntimes
         nout=0      
         do 10 i=1,upper
            do 20 j=1,na
               ovrlap=ddot(n,va(1,j),1,vb(1,i),1)
               call daxpy(n,-ovrlap,va(1,j),1,vb(1,i),1)
   20       continue
            norm=sqrt(ddot(n,vb(1,i),1,vb(1,i),1))
            if(norm.gt.thresh) then
               call dscal(n,1.0d+00/norm,vb(1,i),1)
               nout=nout+1
               call copy(vb(1,i),vb(1,nout),n)  
            endif
   10    continue
         upper=nout
 1000 continue
      if(prnt) then        
         write(iout,1) na, nb   
         write(iout,2) nout
      endif	 	 
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
      end
