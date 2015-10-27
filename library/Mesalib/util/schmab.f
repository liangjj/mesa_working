*deck schmab.f
      subroutine schmab(va,vb,thresh,n,na,nb,nout,prnt)
c***begin prologue     schmab
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
c                          na      = number of va vectors
c                          nb      = number of initial vb vectors
c                          nout    = number of final vb vectors
c***references
c***routines called    saxpy(clams), sdot(clams), sscal(clams)
c
c***end prologue       gschmt
      implicit integer (a-z)
c
      real*8 va(n,na), vb(n,nb)
      real*8 sdot, norm, thresh, ovrlap
      logical prnt
c
      common /io/inp, iout
c
      nout=0      
      do 10 i=1,nb
         do 20 j=1,na
            ovrlap=sdot(n,va(1,j),1,vb(1,i),1)
            call saxpy(n,-ovrlap,va(1,j),1,vb(1,i),1)
   20    continue
         norm=sqrt(sdot(n,vb(1,i),1,vb(1,i),1))
         if(norm.gt.thresh) then
            call sscal(n,1.0d+00/norm,vb(1,i),1)
            nout=nout+1
            call copy(vb(1,i),vb(1,nout),n)  
         endif
   10 continue
      if(prnt) then        
         write(iout,1) na, nb   
         write(iout,2) nout
      endif	 	 
      if(nout.eq.0) then
         call lnkerr('no vectors from schmab')
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
