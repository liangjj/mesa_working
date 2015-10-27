*deck mkh3d.f
c***begin prologue     mkh3d
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            right hand side of inhomogeneous time-dependent
c***                   hamiltonian for 3-dimensional spatial hamiltonian.
c***                   
c***references         
c
c***routines called    
c***end prologue       mkh3d
      subroutine mkh3d(h1,h2,h3,vecout,vecin,tmp,nd,nt,m)
      implicit integer (a-z)
      real*8 h1, h2, h3, vecout, vecin, tmp
      dimension nd(*)
      dimension h1(nd(1),nd(1)), h2(nd(2),nd(2)), h3(nd(3),nd(3))
      dimension vecout(nd(3),nd(2),nd(1),nt,m)
      dimension vecin(nd(3),nd(2),nd(1),m)
      dimension tmp(nd(3),nd(2),nd(1),m)
      common/io/inp, iout
      call ebc(tmp,h3,vecin,nd(3),nd(3),nd(2)*nd(1)*m)
      call filvt(vecout,tmp,nd(1)*nd(2)*nd(3),nt,m)
      do 10 i=1,nd(1)
         do 20 j=1,m
            call ebc(tmp(1,1,i,j),vecin(1,1,i,j),h2,nd(3),nd(2),nd(2))
 20      continue   
 10   continue   
      call filvt(vecout,tmp,nd(1)*nd(2)*nd(3),nt,m)
      do 30 i=1,m
         call ebc(tmp(1,1,1,i),vecin(1,1,1,i),h1,nd(3)*nd(2),
     1                                           nd(1),nd(1))
 30   continue   
      call filvt(vecout,tmp,nd(1)*nd(2)*nd(3),nt,m)
      return
      end       


