*deck cvexp.f
c***begin prologue     cvexp
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            exponentiation
c***                   
c***references         
c
c***routines called    
c***end prologue       cvexp
      subroutine cvexp(vexp,v,s,n)
      implicit integer (a-z)
      complex*16 vexp, eye
      real*8 v, s
      dimension vexp(n), v(n)
      data eye / (0.d0,1.d0) /
      common/io/inp, iout
      do 10 i=1,n
         vexp(i) = exp(-eye*s*v(i))
 10   continue
      return   
      end       






