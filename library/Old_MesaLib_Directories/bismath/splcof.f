*deck splcof
      subroutine splcof(n,ydata,k,nbreak,break,c,sc)
      implicit real*8(a-h,o-z)
      dimension break(nbreak+1), ydata(n), c(k,nbreak), sc((4*n+1)*k)   
      common /io/ inp, iout       
c
c  solve for coefficients in b-spline representation.
c
      ia=n+k
      call copy(ydata,sc(ia+1),n)
      iq=ia+n
      idummy=iq+n*(3*k-2)
      call banmat(n,k-1,k-1,1,1,sc(iq+1),n,sc(ia+1),n,det,sc(idummy+1))
c
c  convert b-spline representation to ppiecewise polynomial
c  representation.
c
      call bsplpp(sc,sc(ia+1),n,k,sc(idummy+1),break,c,nbreak)
      return
      end
