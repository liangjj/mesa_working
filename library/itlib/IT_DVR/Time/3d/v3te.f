*deck v3te.f
c***begin prologue     v3te
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
c***end prologue       v3te
      subroutine v3te(vecin,eig1,eig2,eig3,n1,n2,n3,nt,nc,nvc)
      implicit integer (a-z)
      real*8 vecin, eig1, eig2
      real*8 test, zero, nrzero, one
      dimension vecin(n3,n2,n1,nt,nc,2,nvc)
      dimension eig1(n1), eig2(n2), eig3(n3)
      common/io/inp, iout
      data zero, nrzero, one / 0.d0, 1.0d-06, 1.d0 /
      do 10 i=1,n1
         do 20 j=1,n2
            do 30 k=1,n3
               test = -1.d0/( eig1(i) + eig2(j) + eig3(k) )
               if(abs(test).ge.nrzero) then
                  do 40 l=1,nt
                     do 50 ic=1,nc
                        do 60 m=1,nvc
                           vecin(k,j,i,l,ic,1,m) = 
     1                     vecin(k,j,i,l,ic,1,m)*test
                           vecin(k,j,i,l,ic,2,m) = 
     1                     vecin(k,j,i,l,ic,2,m)*test

 60                     continue
 50                  continue   
 40               continue
               endif
 30         continue
 20      continue
 10   continue   
      return
      end       
