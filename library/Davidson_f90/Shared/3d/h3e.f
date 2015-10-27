*deck h3e.f
c***begin prologue     h3e
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            energy scale a 2-D vector.
c***                                      
c***references         
c
c***routines called    
c***end prologue       h3e
      subroutine h3e(vecin,eig1,eig2,eig3,eig,n1,n2,n3,nvc)
      implicit integer (a-z)
      real*8 vecin, eig1, eig2, eig3, eig
      real*8 test, zero, nrzero, one
      dimension vecin(n3,n2,n1,nvc), eig1(n1), eig2(n2), eig3(n3)
      dimension eig(nvc)
      common/io/inp, iout
      data zero, nrzero, one / 0.d0, 1.0d-06, 1.d0 /
      do 10 i=1,nvc
         do 20 j=1,n1
            do 30 k=1,n2
               do 40 l=1,n3
                  test = eig(i) - eig1(j) - eig2(k) - eig3(l)
                  if(abs(test).ge.nrzero) then
                     vecin(l,k,j,i) = vecin(l,k,j,i)/test
                  else
                     vecin(l,k,j,i) = one
                  endif  
 40            continue
 30         continue
 20      continue
 10   continue   
      return
      end       
