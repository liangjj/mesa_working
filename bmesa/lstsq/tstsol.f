*deck tstsol.f
c***begin prologue     tstsol
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            test of iterative solution.
c***                   
c***references         
c
c***routines called    
c***end prologue       tstsol
      subroutine tstsol(ham,rhs,vec,tvec,n,nrhs)
      implicit integer (a-z)
      real*8 ham, rhs, vec, tvec, error, sdot
      character*80 title
      dimension ham(n,n), rhs(n,nrhs), vec(n,nrhs), tvec(n,nrhs) 
      common/io/inp, iout
      do 10 i=1,n
         ham(i,i) = 1.d0 + ham(i,i)
 10   continue   
      call ebc(tvec,ham,vec,n,n,nrhs)
      do 20 i=1,nrhs
         call vsub(tvec(1,i),tvec(1,i),rhs(1,i),n)
         error=sqrt( sdot(n,tvec(1,i),1,tvec(1,i),1) )
         write(iout,1) i, error
         call copy(tvec(1,i),rhs(1,i),n)
 20   continue   
      title='error vectors'
      call prntrm(title,rhs,n,nrhs,n,nrhs,iout)
      return
 1    format(/,1x,'solution = ',i4,1x,'final rms error = ',e15.8)
      end       



