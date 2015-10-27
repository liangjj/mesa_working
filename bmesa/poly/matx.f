*deck matx.f
c***begin prologue     matx
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            matrix elements of local functions of coordinate
c***                   in diagonal representation.
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       matx
      subroutine matx(pold,pnew,x,wts,f,mat,c,dum,n,npts)
      implicit integer (a-z)
      real*8 pold, pnew, x, wts, f, mat, c, dum, sdot
      character*80 title
      dimension pold(npts,0:n-1), pnew(npts,0:n-1), x(npts)
      dimension wts(npts), f(npts), c(n)
      dimension mat(*), dum(npts,n)
      common/io/inp, iout 
      do 10 i=1,npts
         f(i)=wts(i)*x(i)*exp(-x(i))
c         f(i)=wts(i)*x(i)*x(i)
   10 continue
      do 20 i=1,n
         mat(i)=sdot(npts,f,1,pold(1,i-1),1)
   20 continue
      title='expansion coefficients in original basis'
      call prntrm(title,mat,n,1,n,1,iout) 
       do 30 i=1,n
         mat(i)=sdot(npts,f,1,pnew(1,i-1),1)
   30 continue
      title='expansion coefficients in new basis'
      call prntrm(title,mat,n,1,n,1,iout)                               
c      call vmmul(f,pnew,dum,npts,n)
c      do 40 i=1,n
c         mat(i)=sdot(npts,pnew(1,i-1),1,dum(1,i),1)
c   40 continue              
c      title='matrix of local co-ordinate function'
c      call prntrm(title,mat,n,1,n,1,iout)
c      do 50 i=1,n
c         c(i)=c(i)*mat(i)
c   50 continue
c      title='expansion coefficients in nonstandard fashion'
c      call prntrm(title,c,n,1,n,1,iout)         
      return
      end       
