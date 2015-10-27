*deck prdiff
      subroutine prdiff(x,y,type,neq,nstp)
      real*8 x, y
      logical type
      character*80 title
      character*15 fptoc
      dimension x(nstp), y(neq,nstp)
      common/io/ inp, iout
      if(type) then
         n=neq/2
         do 10 i=1,nstp
            title='real and imaginary values of y at x = '//fptoc(x(i))
            call prntrm(title,y(1,i),n,2,n,2,iout)
            write(iout,*) cos(x(i)), -sin(x(i))
 10      continue
      else
         do 20 i=1,nstp
            title='value of y at x = '//fptoc(x(i))
            call prntrm(title,y(1,i),neq,1,neq,q,iout)
 20      continue             
      endif
      return
      end


