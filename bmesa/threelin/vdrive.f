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

      end









