*deck vadpr.f
c***begin prologue     vadpr
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           pair potential
c***author             schneider, barry (nsf)
c***source             
c***purpose            potential matrix
c***                   
c***description        calculate the potential matrix
c***                   in coordinate space.
c***                   
c***references         
c
c***routines called    
c***end prologue       vadpr
      subroutine vadpr(pgrid,pv,type,a,b,vword,n,zeroit,prnt)
      implicit integer (a-z)
      integer*8 pgrid, pv
      real*8 q1, q2, v, a, b
      character*(*) type
      logical prnt, zeroit
      dimension pgrid(2), n(2), a(2,2), b(2,2)      
      common/io/inp, iout
      pointer (pq1,q1(1))
      pointer (pq2,q2(1))
      pointer (pv,v(1))
      need=wptoin(n(1)*n(2))
      call memory(need,pv,vword,'vint',0)
      if(zeroit) then
         call rzero(v,n(1)*n(2))
      endif
      if(type.eq.'none') then
         return
      endif 
c
c     assign variables and locations
c
      pq1=pgrid(1)
      pq2=pgrid(2)
      if(type.eq.'well') then
         call vwelhy(v,q1,q2,a,b,n)
      elseif(type.eq.'exponential') then
         call vexphy(v,q1,q2,a,b,n)        
      elseif(type.eq.'coulomb') then
         call vcouhy(v,q1,q2,a,n)        
      elseif(type.eq.'harmonic-oscillator') then
         call vhmohy(v,q1,q2,a,b,n)
      else
         call lnkerr('error in potential')
      endif
      return
 1    format(/,5x,'interaction potential = ',a32)
      end       














