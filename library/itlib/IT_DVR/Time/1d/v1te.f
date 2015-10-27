*deck v1te.f
c***begin prologue     v1te
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
c***end prologue       v1te
      subroutine v1te(vecin,eig1,n1,nt,nc,nvc)
      implicit integer (a-z)
      real*8 vecin, eig1
      real*8 test, zero, nrzero, one
      dimension vecin(n1,nt,nc,2,nvc), eig1(n1)
      common/io/inp, iout
      data zero, nrzero, one / 0.d0, 1.0d-06, 1.d0 /
      do 10 i=1,n1
         test = -1.d0/eig1(i)         
         if(abs(test).ge.nrzero) then
            do 20 j=1,nt
               do 30 ic=1,nc
                  do 40 k=1,nvc
                     vecin(i,j,ic,1,k) = vecin(i,j,ic,1,k)*test
                     vecin(i,j,ic,2,k) = vecin(i,j,ic,2,k)*test
 40               continue   
 30            continue
 20         continue
         endif
 10   continue   
      return
      end       
