c $Header: vdrive.f,v 1.2 92/12/12 09:35:12 bis Exp $
*deck vdrive
      subroutine vdrive(v,x,g,f,diag,sudiag,spdiag,rhs,guess,exvc,
     1                  exiter,aold,anew,bold,bnew,temp,ipvt,
     2                  rmin,rmax,rdel,energy,refe,convg,ovtol,value,
     3                  points,m,iter,pottyp,drive,bcond,gtype,mgrid,
     4                  slndir,urefe,ops)
      implicit integer(a-z)
      real*8 v, x, g, f, diag, sudiag, spdiag, rhs, guess, exvc
      real*8 exiter, aold, anew, bold, bnew, temp, rmin, rmax, rdel
      real*8 energy, convg, ovtol, value, s4, twodel, refe
      character*(*) pottyp, drive, bcond, gtype, mgrid, slndir, ops
      logical urefe
      dimension v(0:m), x(0:m), g(0:m), f(0:m), diag(0:m)
      dimension sudiag(0:m), spdiag(0:m), rhs(0:m), guess(0:*)
      dimension exvc(*), exiter(*), aold(*), anew(*), bold(*), bnew(*)
      dimension temp(0:points,6), ipvt(*), s4(4)
      common /io/ inp, iout
c**********************************************************************c
c**********************************************************************c
c          as stated in the main code ignore any multi-grid            c
c                              or                                      c
c          reference energy type calculations.  they do not yet work   c
c**********************************************************************c
c**********************************************************************c
      write(iout,*) '     approach         = ',mgrid
      if (urefe) then
          write(iout,*) '     method           = reference energy'//
     1                                           ' method'
          write(iout,*) '     reference energy = ',refe
      endif
      pnts=points
      mhold=m
c**********************************************************************c
c                 for now just jump to the else part of the next       c
c                 statement                                            c
c**********************************************************************c
      if (mgrid.ne.'single-grid') then
*
*     calculate x, v, g and f on the doubled grid
*
          twodel=2.d0*rdel
          points=(rmax-rmin)/twodel
          write(iout,*) '     initial number of points = ',points
          m=points+10
          call mkgrd(x,rmin,rmax,twodel,m)
          call potntl(v,x,m,points,pottyp) 
          call mkrhs(v,x,g,energy,m,points,drive)
          call mkfx(v,f,energy,0.d0,m,'with v',.false.)
          call numerv(diag,sudiag,spdiag,f,s4,twodel,m,points,nfd)
c**********************************************************************c
c     on exit from numerv s4 contains the subdiagonal, diagonal
c     and superdiagonal values of the numerov factors at the
c     next to last point.
c**********************************************************************c
          call bndry(diag,sudiag,spdiag,f,g,rhs,x,s4(4),twodel,
     1               energy,m,points,nfd,bcond,value)
c**********************************************************************c
c     on exit from bndry s4(4) is filled in with proper numerov
c     right hand side at the next to last point.
c**********************************************************************c
*     solve the wavefunction on the doubled grid
          call sgtsl(nfd,sudiag(1),diag(1),spdiag(1),rhs(1),info)
          rhs(0)=0.d0
          rhs(points)=(s4(4)-s4(1)*rhs(points-2)-s4(2)*rhs(points-1))
     1                           /s4(3)
          if (info.ne.0) then
              call lnkerr('solution to tridiagonal equations bad')
          endif
*     interpolate this solution on the fine grid
          points=pnts
          m=mhold
          count=0
          do 10 i=0,points,2
             guess(i)=rhs(count)
             guess(i+1)=(rhs(count+1)+rhs(count))*.5d0
             count=count+1
   10     continue
          call copy(guess(0),rhs(0),points)
          call mkgrd(x,rmin,rmax,rdel,m)
          call potntl(v,x,m,points,pottyp) 
          call mkrhs(v,x,g,energy,m,points,drive)
          call mkfx(v,f,energy,0.d0,m,'with v',.false.)
          call numerv(diag,sudiag,spdiag,f,s4,rdel,m,points,nfd)
          call bndry(diag,sudiag,spdiag,f,g,rhs,x,s4(4),rdel,energy,
     1                m,points,nfd,bcond,value)
      elseif (mgrid.eq.'single-grid') then
          call mkgrd(x,rmin,rmax,rdel,m)
          call potntl(v,x,m,points,pottyp)
          if ( slndir.eq.'matrix') then
                  call mkrhs(v,x,g,energy,m,points,drive)
                  call mkfx(v,f,energy,0.d0,m,'with v',.false.)
                  call numerv(diag,sudiag,spdiag,f,s4,rdel,m,points,nfd)
                  call bndry(diag,sudiag,spdiag,f,g,rhs,x,s4(4),rdel,
     1                       energy,m,points,nfd,bcond,value)
          elseif (slndir.eq.'one-plus-matrix'.or.
     1            slndir.eq.'one-minus-matrix') then
                  gtype='unit-matrix'
*      calculate the difference equation matrix elements with v and g zero
*      and no boundary condition specified. we will fix this later.
*      this is tantamount to the solution of the first order wavefunction.
                  call mkfx(v,f,energy,refe,m,'without v',urefe)
                  call numerv(diag,sudiag,spdiag,f,s4,rdel,m,points,nfd)
                  call copy(sudiag,temp(0,1),points+1)
                  call copy(diag,temp(0,2),points+1)
                  call copy(spdiag,temp(0,3),points+1)
*
*      put in the desired inhomogeneity then calculate correct
*      rhs and boundary condition.
                  if (urefe) then
                      call mkrhs(v,x,g,refe,m,points,drive)
                      call newtrm(x,g,energy,refe,m,points)
                      call bndry(diag,sudiag,spdiag,f,g,rhs,x,
     1                           s4(4),rdel,refe,m,points,nfd,
     2                           bcond,value)
                  else
                      call mkrhs(v,x,g,energy,m,points,drive)
                      call bndry(diag,sudiag,spdiag,f,g,rhs,x,
     1                           s4(4),rdel,energy,m,points,nfd,
     2                           bcond,value)
                  endif
*      solve the zeroth order equation.
                  call sgtsl(nfd,sudiag(1),diag(1),spdiag(1),
     1                       rhs(1),info)
*      calculate the last point.
                  rhs(0)=0.d0
                  rhs(points)=(s4(4)-s4(1)*rhs(points-2)
     1                              -s4(2)*rhs(points-1))/s4(3)
                  call phase(rhs,x,energy,points)
          endif
      endif
      call gen0(diag,sudiag,spdiag,rhs,guess,m,gtype)
      call slvit(diag,sudiag,spdiag,v,f,rhs,x,s4,rdel,energy,
     1           refe,temp,guess,exvc,exiter,aold,anew,bold,bnew,
     2           ipvt,points,nfd,iter,convg,ovtol,gtype,slndir,
     3           bcond,urefe,ops)
      if ( slndir.eq.'matrix') then
           rhs(0)=0.d0
           rhs(points)=(s4(4)-s4(1)*rhs(points-2)
     1                       -s4(2)*rhs(points-1))/s4(3)
      endif
      call phase(rhs,x,energy,points)
      return
      end









