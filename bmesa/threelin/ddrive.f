*deck ddrive
      subroutine ddrive(v,x,g,f,diag,sudiag,spdiag,rhs,rmin,rmax,
     1                  rdel,energy,value,m,last,pottyp,drive,bcond)
      implicit integer(a-z)
      real*8 v, x, g, f, diag, sudiag, spdiag, rhs, energy, value
      real*8 s4, rmin, rmax, rdel
      character*(*) pottyp, drive, bcond
      dimension v(0:*), x(0:*), g(0:*), f(0:*), diag(0:*), sudiag(0:*)
      dimension spdiag(0:*), rhs(0:*), s4(4)
      common /io/ inp, iout
*
      call mkgrd(x,rmin,rmax,rdel,m)
      call potntl(v,x,m,last,pottyp)
      call mkrhs(v,x,g,energy,m,last,drive)
      call mkfx(v,f,energy,0.d0,m,'with v',.false.)
      rhs(0)=0.d0
      call numerv(diag,sudiag,spdiag,f,s4,rdel,m,last,nfd)
      call bndry(diag,sudiag,spdiag,f,g,rhs,x,s4(4),rdel,energy,m,
     1           last,nfd,bcond,value)
      call sgtsl(nfd,sudiag(1),diag(1),spdiag(1),rhs(1),info)
      rhs(last)=(s4(4)-s4(1)*rhs(last-2)-s4(2)*rhs(last-1))/s4(3)
      if (info.ne.0) then
          call lnkerr('solution to tridiagonal equations bad')
      endif
      call phase(rhs,x,energy,last)
      return
      end











