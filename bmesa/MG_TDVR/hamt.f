*deck hamt.f
c***begin prologue     hamt
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           1-d hamiltonian
c***author             schneider, barry (nsf)
c***source             
c***purpose            calculate one body time matrix.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       hamt
      subroutine hamt(zham,zv,n,ug,mgr,test,parray,prnt)
      implicit integer (a-z)
      real*8 array
      real*8 work 
      character*80 title
      logical test, prnt
      dimension n(mgr), zham(mgr), zv(mgr)
      common/io/inp, iout
      pointer(parray,array(1))
      pointer(pwork,work(1))
      matrix=1
      eigr=matrix+n(ug)*n(ug)
      if(test) then
         eigi=eigr+n(ug)
         vecl=eigi+n(ug)
         vecr=vecl+n(ug)*n(ug)
         need=wpadti(vecr+n(ug)*n(ug))
      else
         need=wpadti(eigr)
      endif
      call memory(need,pwork,ngot,'work',0)
      call copy(array(zham(ug)),work(matrix),n(ug)*n(ug))
      if(prnt) then
         title='hamiltonian'
         call prntrm(title,work(matrix),n(ug),n(ug),n(ug),n(ug),iout)
      endif
      if(test) then
         call eigslv(work(matrix),work(eigr),work(eigi),
     1               work(vecl),work(vecr),n(ug))
      endif
      call memory(-ngot,pwork,idum,'work',idum)
      return         
      end       


