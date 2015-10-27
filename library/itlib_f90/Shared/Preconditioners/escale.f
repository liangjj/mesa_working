*deck escale.f
c***begin prologue     escale
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            energy scale vectors.
c***                                      
c***references         
c
c***routines called    
c***end prologue       escale
      subroutine escale(vecin,vecout,eig,eig0,n,nvc,m)
      implicit integer (a-z)
      real*8 vecin, vecout, eig0, eig
      real*8 test, zero, nrzero, one
      dimension vecin(m,nvc), vecout(m,nvc), eig0(n), eig(nvc)
      common/io/inp, iout
      data zero, nrzero, one / 0.d0, 1.0d-06, 1.d0 /
      do 10 i=1,nvc
         do 20 j=1,n
            test = eig(i) - eig0(j)
            if(abs(test).ge.nrzero) then
               vecout(j,i) = vecin(j,i)/test
            else
               vecout(j,i) = one           
            endif            
 20      continue
 10   continue   
      return
      end       


