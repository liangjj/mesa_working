*deck h1e.f
c***begin prologue     h1e
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            energy scale a 1-D vector.
c***                                      
c***references         
c
c***routines called    
c***end prologue       h1e
      subroutine h1e(vecin,eig1,eig,n1,nvc)
      implicit integer (a-z)
      real*8 vecin, eig1, eig
      real*8 test, zero, nrzero, one
      dimension vecin(n1,nvc), eig1(n1), eig(nvc)
      common/io/inp, iout
      data zero, nrzero, one / 0.d0, 1.0d-06, 1.d0 /
      do 10 i=1,nvc
         do 20 j=1,n1
            test = eig(i) - eig1(j)         
            if(abs(test).ge.nrzero) then
               vecin(j,i) = vecin(j,i)/test
            else
               vecin(j,i) = one           
            endif
 20      continue
 10   continue   
      return
      end       
