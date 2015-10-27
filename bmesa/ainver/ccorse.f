*deck ccorse.f
c***begin prologue     ccorse
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       ccorse
      subroutine ccorse(fine,corse,pf,pc,wt,nf,nc)
      implicit integer (a-z)
      real*8 fine, corse, pf, pc, wt
      dimension fine(nf), corse(nc), pf(nf,nf), pc(nf,nc), wt(nf)
      common/io/inp, iout
      do 10 i=1,nc
         corse(i) = 0.d0
         do 20 j=1,nf
            corse(i) = corse(i) + fine(j)*wt(j)*pc(j,i)*pf(j,j)
 20      continue
 10   continue   
      return
      end
