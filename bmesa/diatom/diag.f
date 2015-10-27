*deck ker.f
c***begin prologue     ker
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            kinetic energy matrix elements in r coordinate.
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       ker
      subroutine ker(p,dp,ddp,x,wt,tmat,eig,work,dum,n,npt,typ,prn)
      implicit integer (a-z)
      real*8 p, dp, ddp, x, wt, tmat, eig, work, dum
      logical prn, typ
      character*80 title
      dimension p(npt,0:n-1), dp(npt,0:n-1), ddp(npt,0:n-1)
      dimension x(npt), wt(npt), tmat(n,n), eig(n), work(n), dum(n)
      common/io/inp, iout 
      call rzero(tmat,n*n)
      do 10 i=1,n
         do 20 j=1,i
            do 30 k=1,npt
               tmat(i,j)=tmat(i,j) + p(k,i-1) * wt(k) * ddp(k,j-1)
 30         continue
            tmat(i,j)=tmat(i,j) - p(npt,i-1)*dp(npt,j-1)
            tmat(i,j)=-.5d0*tmat(i,j)
            tmat(j,i)=tmat(i,j)
 20      continue
 10   continue   
      if (prn) then
          title='kinetic energy matrix elements in r coordinate'
          call prntrm(title,tmat,n,n,n,n,iout)
      endif
      if(typ) then
         call tred2(n,n,tmat,eig,work,ham)
         call tql2(n,n,eig,work,tmat,ierr)
      else
         call tred1(n,n,tmat,eig,work,dum)
         call tql1(n,eig,work,ierr)
      endif                  
      title='eigenvalues of radial kinetic energy'
      call prntrm(title,eig,n,1,n,1,iout)
      return
      end       
