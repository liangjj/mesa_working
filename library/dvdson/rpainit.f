*deck rpainit.f
c***begin prologue     rpainit
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           restart davidson for rpa matrices
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       rpainit
      subroutine rpainit(vold,hvold,vnew,hvnew,eigvec,tmp1,tmp2,
     1                   btmp,b,thresh,n,nend,maxvec,prnt)
      implicit integer (a-z)
      real*8 vold, hvold, vnew, hvnew, eigvec, tmp1, tmp2, b 
      real*8 btmp
      real*8 thresh
      complex*16 cdum
      logical prnt
      character*80 title
      dimension vold(n,maxvec), hvold(n,maxvec)
      dimension vnew(n,maxvec), hvnew(n,maxvec)
      dimension eigvec(maxvec,nend), tmp1(nend,nend), tmp2(nend,nend)
      dimension b(maxvec,nend), btmp(maxvec,nend)
      common/io/inp, iout
      write(iout,1)
c
c     first step is to orthogonalize the current set of nend vectors.
c
c     calculate the transformation coefficients relating
c     the non-orthogonal and orthogonal sets.
c
      do 10 i=1,nend
         do 20 j=1,nend
            tmp1(i,j) = eigvec(i,j)
 20      continue
 10   continue
      if (prnt) then
          title='old transformation matrix before initialization'
          call prntrm(title,tmp1,nend,nend,nend,nend,iout)
      endif       
      call gschmt(tmp1,thresh,nend,1,nend,nout,.true.,.false.)
      if (prnt) then
          title='new transformation coefficients for initialization'
          call prntrm(title,tmp1,nend,nend,nend,nend,iout)
      endif          
      if(nout.ne.nend) then
         call lnkerr('quit in rpainit')
      endif
c
c     use that to transform vectors and hamiltonian on vectors
c     to the new set.
c
      call ebc(vnew,vold,tmp1,n,nend,nend)
      call ebc(hvnew,hvold,tmp1,n,nend,nend)
      call copy(vnew,vold,n*nend)
      call copy(hvnew,hvold,n*nend)
      if (prnt) then
          title='new davidson vectors for initialization'
          call prntrm(title,vnew,n,nend,n,nend,iout)
          title='new hamiltonian on davidson vectors for initialization'
          call prntrm(title,hvnew,n,nend,n,nend,iout)
      endif    
c
c
c     transform the small matrices to that basis.
c
      if (prnt) then
          title='new transformation coefficients for initialization'
          call prntrm(title,tmp1,nend,nend,nend,nend,iout)
      endif          
      do 30 i=1,nend
         do 40 j=1,nend
            tmp2(i,j)=b(i,j)
 40      continue
 30   continue       
      call udagmu(tmp2,cdum,tmp1,cdum,eigvec,cdum,nend,nend,'real')
c
c     copy tmp2 to the working array b.
c
      do 50 i=1,nend
         do 60 j=1,nend
            b(i,j) = tmp2(i,j)
 60      continue
 50   continue
      if(prnt) then  
         title='transformed initial small matrix'
         call prntrm(title,b,nend,nend,maxvec,maxvec,iout)
      endif    
      return
 1    format(/,10x,'***** maximum number of vectors exceeded *****',
     1        /10x,'***** transform to best current basis *****')
      end
