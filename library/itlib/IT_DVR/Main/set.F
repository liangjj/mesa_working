*deck set.f
c***begin prologue     set
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           set locations
c***author             schneider, barry (nsf)
c***source             
c***routines called    
c***end prologue      set
      subroutine set(ham,v,eigv,s,eig,n)
      implicit integer (a-z)
      common/io/inp, iout
      ham=1
      v=ham+n*n
      eigv=v+n
      s=eigv+n*n
      eig=s+2*n
      return
      end       
