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
      subroutine rdiag(v,eig,work,n)
      implicit integer (a-z)
      real*8 v, eig, work
      dimension v(n,n), eig(n), work(*)
      common/io/inp, iout 
      call tred2(n,n,v,eig,work,v)
      call tql2(n,n,eig,work,v,ierr)
      return
      end       
