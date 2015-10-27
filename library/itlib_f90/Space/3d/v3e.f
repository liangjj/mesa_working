*deck v3e.f
c***begin prologue     v3e
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
c***end prologue       v3e
      subroutine v3e(vecin,eig1,eig2,eig3,n1,n2,n3,nc,nvc)
      implicit integer (a-z)
      real*8 vecin, eig1, eig2, eig3
      real*8 test, zero, nrzero, one
      dimension vecin(n3,n2,n1,nc,nvc), eig1(n1), eig2(n2), eig3(n3)
      common/io/inp, iout
      data zero, nrzero, one / 0.d0, 1.0d-06, 1.d0 /
      do 10 i=1,n1
         do 20 j=1,n2
            do 30 k=1,n3
               test = -1.d0/( eig1(i) + eig2(j) + eig3(k) )
               if(abs(test).ge.nrzero) then
                  do 40 ic=1,nc
                     do 50 l=1,nvc
                       vecin(k,j,i,ic,l) = vecin(k,j,i,ic,l)*test
 50                  continue   
 40               continue
               endif
 30         continue
 20      continue
 10   continue   
      return
      end       
