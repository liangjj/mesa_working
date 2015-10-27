*deck rdiag.f
c***begin prologue     rdiag
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           diagonalization
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for real diagonalization.
c***                   
c***references         
c
c***routines called    
c***end prologue       rdiag
      subroutine rdiag(v,eig,work,n,m)
      implicit integer (a-z)
      real*8 v, eig, work
      dimension v(n,*), eig(*), work(*)
      common/io/inp, iout 
      call tred2(n,m,v,eig,work,v)
      call tql2(n,m,eig,work,v,ierr)
      return
      end       
