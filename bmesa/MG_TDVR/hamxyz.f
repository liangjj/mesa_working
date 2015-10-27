*deck hamxyz.f
c***begin prologue     hamxyz
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           1-d hamiltonian
c***author             schneider, barry (nsf)
c***source             
c***purpose            calculate one body hamiltonian matrix.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       hamxyz
      subroutine hamxyz(zham,zv,ham,work,n,ug,mgr,maxd,test,parray,
     1                  prnt,spac)
      implicit integer (a-z)
      real*8 array
      real*8 ham, work 
      character*80 title
      logical test, prnt, spac
      dimension n(mgr), zham(mgr), zv(mgr), ham(*) 
      dimension work(1)
      common/io/inp, iout
      pointer(parray,array(1))
      if(prnt) then
         title='hamiltonian'
         call prntrm(title,array(zham(ug)),n(ug),n(ug),n(ug),n(ug),iout)
      endif
      e=1
      w=e+maxd      
      if(test) then
         call copy(array(zham(ug)),ham,n(ug)*n(ug))
         call addd(ham,ham,array(zv(ug)),n(ug))
         call tred2(n(ug),n(ug),ham,work(e),work(w),ham)
         call tql2(n(ug),n(ug),work(e),work(w),ham,ierr)
         title='eigenvalues'
         call prntrm(title,work(e),n(ug),1,n(ug),1,iout) 
      endif
      if(spac) then
         call rzero(array(zham(ug)),n(ug)*n(ug))
         call rzero(array(zv(ug)),n(ug))
      endif
      return         
      end       


