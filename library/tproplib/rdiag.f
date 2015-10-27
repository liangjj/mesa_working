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
      subroutine rdiag(v,eig,vec,work,n,m)
      implicit integer (a-z)
      character*80 title
      real*8 v, eig, vec, work
      dimension v(n,*), eig(*), work(*), vec(n,*)
      common/io/inp, iout
      call tred2(n,m,v,eig,work,vec)
      call tql2(n,m,eig,work,vec,ierr)
      if(ierr.ne.0) then
         call lnkerr('error from direct diagonalization routine')
      endif   
      return
      end       
