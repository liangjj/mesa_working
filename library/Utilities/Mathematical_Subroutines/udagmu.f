*deck udagmu.f
c***begin prologue     udagmu
c***date written       970430   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            perform a unitary transformation on a matrix
c                           +
c***description        m = u m u
c***
c                                m = m(n,n) = input and output
c                                u = u(n,n) = unitary transformation matrix
c                                s = s(n,n) = scratch matrix
c***references         
c
c***routines called    
c***end prologue       udagmu
      subroutine udagmu(m,cm,u,cu,s,cs,n,dim,type)
      implicit integer (a-z)
      real*8 m, u, s
      complex*16 cm, cu, cs
      character*(*) type
      dimension m(dim,n), u(dim,n), s(dim,n)
      dimension cm(dim,n), cu(dim,n), cs(dim,n)
      common/io/inp, iout 
      if(type.eq.'real') then
         call ebcxx(s,m,u,n,n,n,dim,dim,dim)
         call ebtcxx(m,u,s,n,n,n,dim,dim,dim)
      elseif(type.eq.'complex') then
         call cebcx(cs,dim,cm,dim,cu,dim,n,n,n)
         call cehbtcx(cm,dim,cu,dim,cs,dim,n,n,n)
      endif
      return
      end       
