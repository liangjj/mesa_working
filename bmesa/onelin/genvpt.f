c $Header: genvpt.f,v 1.2 92/12/12 09:34:14 bis Exp $
*deck genvpt
      subroutine genvpt (diag,sudiag,spdiag,v,f,x,s4,rdel,energy,
     1                   refe,temp,veco,vecn,value,bcond,points,
     2                   nfd,urefe)
      implicit integer(a-z)
      real*8 diag, sudiag, spdiag, v, f, x, s4, rdel, energy, temp
      real*8 veco, vecn, value, refe
      character*(*) bcond
      logical urefe
      dimension diag(0:*), sudiag(0:*), spdiag(0:*)
      dimension v(0:points), f(*), x(*), s4(4), temp(0:points,6)
      dimension veco(points), vecn(points)
      common /io/ inp, iout
*
*          subroutine to generate a new iteration in the vector
*          sequence used to solve the differential equations after they
*          have been reduced to a set of matrix equations via a
*          finite difference technique.
*
*          this subroutine is specifically designed for a 
*          tridiagonal system and based on perturbation theory.
*
*
*
      do 10 i=1,points
         temp(i,4) = ( 2.d0*v(i) - ( energy-refe ))*veco(i)
   10 continue
      call copy(temp(0,1),sudiag,points+1)
      call copy(temp(0,2),diag,points+1)
      call copy(temp(0,3),spdiag,points+1)
      if (urefe) then
           call bndry(diag,sudiag,spdiag,f,temp(0,4),temp(0,5),
     1                x,s4(4),rdel,refe,points,points,nfd,
     2                bcond,value)      
      else
           call bndry(diag,sudiag,spdiag,f,temp(0,4),temp(0,5),
     1                x,s4(4),rdel,energy,points,points,nfd,
     2                bcond,value)     
      endif 
      call sgtsl(nfd,sudiag(1),diag(1),spdiag(1),temp(1,5),info)
      temp(points,5)=(s4(4)-s4(1)*temp(points-2,5)
     1                              -s4(2)*temp(points-1,5))/s4(3)
      call copy(temp(1,5),vecn,points)
      if (info.ne.0) then
          call lnkerr('error in call to sgtsl')
      endif
      return
      end


