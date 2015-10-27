*deck creinit.f
c***begin prologue     creinit
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           restart davidson for non-symmetric matrices
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       creinit
      subroutine creinit(vold,hvold,vnew,hvnew,eigvec,tmp1,tmp2,
     1                   btmp,b,thresh,n,nend,m,maxvec,prnt)
      implicit integer (a-z)
      complex*16 vold, hvold, vnew, hvnew, eigvec, tmp1, tmp2, b 
      complex*16 btmp
      real*8 thresh, dum
      logical prnt
      character*80 title
      dimension vold(n,maxvec), hvold(n,maxvec)
      dimension vnew(n,maxvec), hvnew(n,maxvec)
      dimension eigvec(maxvec,nend), tmp1(nend,nend), tmp2(nend,nend)
      dimension b(maxvec,nend), btmp(maxvec,nend)
      common/io/inp, iout
      write(iout,1) m
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
          call prntcm(title,tmp1,nend,nend,nend,nend,iout)
      endif       
      call cschmt(tmp1,thresh,nend,1,nend,nout,.true.)
      if (prnt) then
          title='new transformation coefficients for initialization'
          call prntcm(title,tmp1,nend,nend,nend,nend,iout)
      endif          
      if(nout.ne.nend) then
         call lnkerr('quit in creinit')
      endif
c
c     use that to transform vectors and hamiltonian on vectors
c     to the new set.
c
      call cebc(vnew,vold,tmp1,n,nend,nend)
      call cebc(hvnew,hvold,tmp1,n,nend,nend)
      call cc2opy(vnew,vold,n*nend)
      call cc2opy(hvnew,hvold,n*nend)
      if (prnt) then
          title='new davidson vectors for initialization'
          call prntcm(title,vnew,n,nend,n,nend,iout)
          title='new hamiltonian on davidson vectors for initialization'
          call prntcm(title,hvnew,n,nend,n,nend,iout)
      endif    
c
c
c     transform the small matrices to that basis.
c
      if (prnt) then
          title='new transformation coefficients for initialization'
          call prntcm(title,tmp1,nend,nend,nend,nend,iout)
      endif          
      do 30 i=1,nend
         do 40 j=1,nend
            tmp2(i,j)=b(i,j)
 40      continue
 30   continue       
      call udagmu(dum,tmp2,dum,tmp1,dum,eigvec,nend,nend,'complex')
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
         call prntcm(title,b,nend,nend,maxvec,maxvec,iout)
      endif    
      return
 1    format(/,10x,'***** maximum number of vectors exceeded *****',
     1        /10x,'      contract back to ',1x,i4,' vectors')
      end
