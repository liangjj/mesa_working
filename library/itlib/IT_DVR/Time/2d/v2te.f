*deck v2te.f
c***begin prologue     v2te
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
c***end prologue       v2te
      subroutine v2te(vecin,eig1,eig2,n1,n2,nt,nc,nvc)
      implicit integer (a-z)
      real*8 vecin, eig1, eig2, fac
      real*8 test, zero, nrzero, one
      dimension vecin(n2,n1,nt,nc,2,nvc), eig1(n1), eig2(n2)
      common/io/inp, iout
      data zero, nrzero, one / 0.d0, 1.0d-06, 1.d0 /
      do 10 i=1,n1
         do 20 j=1,n2
	    test = -1.d0/( eig1(i) + eig2(j) )
            if(abs(test).ge.nrzero) then
               do 30 k=1,nt
                  do 40 ic=1,nc
	             do 50 l=1,nvc
                        vecin(j,i,k,ic,1,l) = vecin(j,i,k,ic,1,l)*test
                        vecin(j,i,k,ic,2,l) = vecin(j,i,k,ic,2,l)*test
 50                  continue   
 40               continue
 30            continue
            endif
 20      continue
 10   continue   
      return
      end       
