c $Header: genvec.f,v 1.2 92/12/12 09:34:11 bis Exp $
*deck genvec
      subroutine genvec (diag,sudiag,spdiag,veco,vecn,points,nfd,type)
      implicit integer(a-z)
      real*8 diag, sudiag, spdiag, veco, vecn
      character *(*) type
      dimension diag(0:points), sudiag(0:points), spdiag(0:points)
      dimension veco(points), vecn(points)
      common /io/ inp, iout
      data first /0/
      save first
*
*          subroutine to generate a new iteration in the vector
*          sequence used to solve the differential equations after they
*          have been reduced to a set of matrix equations via a
*          finite difference technique.
*
*          this subroutine is specifically designed for a 
*          tridiagonal system.
*
*
*
      if (first.eq.0) then
          write(iout,*) '     preconditioner = ',type
          first=1
      endif
      call rzero(vecn,points)
      do 10 i=2,nfd-1
         vecn(i) = vecn(i) + sudiag(i)*veco(i-1) + diag(i)*veco(i) + 
     1                       spdiag(i)*veco(i+1)
   10 continue
      vecn(1) = vecn(1) + diag(1)*veco(1) + spdiag(1)*veco(2)
      vecn(nfd) = vecn(nfd) + diag(nfd)*veco(nfd) + 
     1                            sudiag(nfd)*veco(nfd-1)
c      write(iout,*) vecn
      if (type.eq.'unit-matrix') then
          write(iout,*) 'unit matrix'
          return
      elseif (type.eq.'diagonal') then
          do 20 i=1,nfd
             vecn(i)=vecn(i)/diag(i)
   20     continue
      elseif (type.eq.'lower-triangular') then
          vecn(1)=vecn(1)/diag(1) 
          do 30 i=2,nfd
             vecn(i)=(vecn(i)-sudiag(1)*vecn(i-1))/diag(i)
   30     continue
      else
          call lnkerr('you screwed up in genvec')
      endif
      return
      end
