*deck prdiff
      subroutine prdiff(p1,p2,p3,q1,q2,q3,t,psi,dim,n1,n2,n3,n,nstp)
      implicit integer(a-z)
      real*8 p1, p2, p3, q1, q2, q3, t, psi
      character*80 title
      character*15 fptoc
      dimension q1(n1), q2(n2), q3(n3)
      dimension p1(n1,n1), p2(n2,n2), p3(n3,n3), t(nstp), psi(n,2,nstp)
      common/io/ inp, iout
      if(dim.eq.1) then
         do 10 i=1,nstp
            write(iout,1) t(i)
            write(iout,*) cos(.5d0*t(i)*t(i)), -sin(.5d0*t(i)*t(i))
            do 20 j=1,n1
               write(iout,2) q1(j), psi(j,1,i), psi(j,2,i)
 20         continue   
 10      continue
      elseif(dim.eq.2) then
         do 30 i=1,nstp
            write(iout,1) t(i)
            count=0
         do 40 j=1,n1
            do 50 k=1,n2
               count=count+1
               write(iout,3) q1(j), q2(k), psi(count,1,i), 
     1                                     psi(count,2,i) 
 50         continue
 40      continue
 30   continue   
      elseif(dim.eq.3) then
         do 60 i=1,nstp
            write(iout,1) t(i)
            count=0
            do 70 j=1,n1
               do 80 k=1,n2
                  do 90 l=1,n3
                     count=count+1
                     write(iout,4) q1(j), q2(k), q3(l), psi(count,1,i), 
     1                                        psi(count,2,i) 
 90               continue
 80            continue
 70         continue   
 60      continue   
      else   
         call lnkerr('error in dimension')
      endif
      return
 1    format(/,5x,'time = ',e15.8)
 2    format(/,5x,'x        = ',e15.8,/,5x,
     1            'real psi = ',e15.8,1x,'imag psi = ',e15.8)
 3    format(/,5x,'x        = ',e15.8,/,5x,
     1            'y        = ',e15.8,/,5x,
     2            'real psi = ',e15.8,1x,'imag psi = ',e15.8)
 4    format(/,5x,'x        = ',e15.8,/,5x,
     1            'y        = ',e15.8,/,5x,
     2            'z        = ',e15.8,/,5x,
     2            'real psi = ',e15.8,1x,'imag psi = ',e15.8)
      end


