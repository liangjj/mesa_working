*deck reinit.f
c***begin prologue     reinit
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           restart davidson
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       reinit
      subroutine reinit(vold,hvold,vnew,hvnew,trmat,b,btmp,n,
     1                  nend,m,maxvec,prnt)
      implicit integer (a-z)
      real*8 vold, hvold, vnew, hvnew, b, btmp, trmat
      complex*16 cdum
      logical prnt
      character*80 title
      dimension vold(n,maxvec), hvold(n,maxvec)
      dimension vnew(n,maxvec), hvnew(n,maxvec)
      dimension b(maxvec,nend), btmp(maxvec,nend)
      dimension trmat(maxvec,nend)
      common/io/inp, iout
      if (prnt) then
          title='old transformation matrix before initialization'
          call prntrm(title,trmat,nend,nend,maxvec,nend,iout)
      endif    
c
c     copy the current best vectors and hamiltonian on vectors from the
c     scratch to the working arrays.
c
      call copy(vold,vnew,n*maxvec)
      call copy(hvold,hvnew,n*maxvec)
      if (prnt) then
          title='new davidson vectors for initialization'
          call prntrm(title,vnew,n,nend,n,nend,iout)
          title='new hamiltonian on davidson vectors for initialization'
          call prntrm(title,hvnew,n,nend,n,nend,iout)
      endif    
c
c     transform the small matrices to the best basis.
c
c      call ebcxx(btmp,b,trmat,nend,nend,nend,maxvec,maxvec,maxvec)
c      call ebtcxx(b,trmat,btmp,nend,nend,nend,maxvec,maxvec,maxvec)
      call udagmu(b,cdum,trmat,cdum,btmp,cdum,nend,maxvec,'real')
      if(prnt) then  
         title='transformed initial small matrix'
          call prntrm(title,b,nend,nend,maxvec,maxvec,iout)
      endif    
c
c     copy b to the working array btmp.
c
      do 10 i=1,nend
         do 20 j=1,i
            btmp(i,j) = b(i,j)
            btmp(j,i) = b(j,i)
 20      continue
 10   continue   
      return
      end
