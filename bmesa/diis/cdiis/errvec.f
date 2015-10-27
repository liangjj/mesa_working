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
      complex*16 matrix, rhs, xold, err
      real*8 rms, errmax
      logical prnt
      character*80 title
      character*3 itoc
      dimension matrix(n,n), rhs(n), err(n), xold(n)
      common/io/inp, iout
      call cc2opy(rhs,err,n)
      call cambc(err,matrix,xold,n,n,1)
      rms = 0.d0
      do 10 i=1,n
         rms = rms + conjg(err(i))*err(i)
 10   continue   
      rms = sqrt (rms )
      errmax=0.d0   
      do 20 i=1,n
         errmax=max(errmax,abs(err(i)))
 20   continue   
      write(iout,1) iter, rms, errmax   
      if(prnt) then
         title='negative of error matrix for iteration = '//itoc(iter)
         call prntcm(title,err,n,1,n,1,iout)
      endif
      return
 1    format(/,5x,'iteration         = ',i3,/,5x,
     1            'rms error         = ',e15.8,/,5x,
     2            'max. abs. element = ',e15.8)
      end       
