*deck subham.f
c***begin prologue     subham
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            time-dependent hamiltonian not including
c***                   non-linear term
c***                   
c***description        total hamiltonian in space contructed explicitly.  
c***references         
c
c***routines called    
c***end prologue       subham
      subroutine subham(ham,hsub,n)
      implicit integer (a-z)
      real*8 ham, hsub
      logical prn
      character*80 title
      dimension ham(n,n), hsub(*)
      common/io/inp, iout      

      if (prn) then
          title='full space and time-dependent hamiltonian'
          call prntrm(title,ham,n,n,n,n,iout)
      endif     
      return
      end       
