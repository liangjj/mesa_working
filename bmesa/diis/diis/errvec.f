*deck errvec.f
c***begin prologue     errvec
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       errvec
      subroutine errvec(matrix,rhs,xold,err,errmax,iter,n,prnt)
      implicit integer (a-z)
      real*8 matrix, rhs, xold, err
      real*8 rms, errmax, sdot
      logical prnt
      character*80 title
      character*3 itoc
      dimension matrix(n,n), rhs(n), err(n), xold(n)
      common/io/inp, iout
      call copy(rhs,err,n)
      call ambc(err,matrix,xold,n,n,1)
      rms = sqrt (sdot(n,err,1,err,1) )
      errmax=0.d0
      do 10 i=1,n
         errmax=max(errmax,abs(err(i)))
 10   continue   
      write(iout,1) iter, rms, errmax   
      if(prnt) then
         title='negative of error matrix for iteration = '//itoc(iter)
         call prntrm(title,err,n,1,n,1,iout)
      endif
      return
 1    format(/,5x,'iteration         = ',i3,/,5x,
     1            'rms error         = ',e15.8,/,5x,
     2            'max. abs. element = ',e15.8)
      end       
