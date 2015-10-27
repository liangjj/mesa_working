*deck mkspln
c***begin prologue     mkspln
c***date written       xxxxxx   (yymmdd)
c***revision date      890422   (yymmdd)
c***keywords           m6004, link 6004, bessel, spline
c***author             rescigno, t. n.(llnl)
c***source             m6004
c***purpose            spline fit free functions
c***description        regular and irregular functions spline fit
c***references       
c
c***routines called
c***end prologue       mkspln
      subroutine mkspln(x,xinv,j,jp,y,yp,hs,hsder,cj,cy,d,toim,toim2,
     1                  rdel,alpha,gamma,ncut,ipow,nr,lmax,ltop)
      implicit integer (a-z)
      real *8 x, alpha, gamma, toim, toim2, term, term2, cr, crp
      real *8 crpp, gr, grp, grpp, gc, gcp, gcpp, j, jp, y, yp
      real *8 rdel, xinv, rdeli
      complex *16 hs, hsder, yp1, yp2, ai, precon, cj, d, cy
      dimension x(nr), xinv(nr), toim(nr), toim2(nr), j(nr,0:ltop)
      dimension jp(nr,0:ltop), y(nr,0:ltop), yp(nr,0:ltop)
      dimension hs(nr,0:lmax), hsder(nr,0:lmax), cj(nr,0:lmax)
      dimension cy(nr,0:lmax), d(nr), ipow(0:lmax)
      ai=dcmplx(0.d0,1.d0)
      precon = -0.5d0*ai
      rdeli=1.d+00/rdel
      do 10 i=1,nr
         toim(i) = exp(-alpha*x(i))
         toim2(i)=exp(-gamma*x(i))
   10 continue
      do 20 l=0,lmax
         do 30 i=1,nr
            term = 1.0d0-toim(i)
            term2=1.0d0-toim2(i)
            cr = term**(2*l+3)
            crp = (2*l+3)*alpha*toim(i)*term**(2*l+2)
            crpp = crp*((2*l+2)*alpha*toim(i)/term - alpha)
            gr = term2**ncut
            grp=ncut*toim2(i)*gamma*term2**(ncut-1)
            grpp=grp*((ncut-1)*gamma*toim2(i)/term2 - gamma)
            gc=cr*gr
            gcp=cr*grp+crp*gr
            gcpp= cr*grpp+2.d0*grp*crp + gr*crpp
            hs(i,l) = gr*(ai*j(i,l) -cr*y(i,l))*xinv(i)
            hsder(i,l)=precon*(j(i,l)*grpp + 2.d0*grp*jp(i,l) +
     1                ai*(y(i,l)*gcpp + 2.d0*gcp*yp(i,l)))*
     2                                  xinv(i)**ipow(l)
   30    continue
         yp1=(hs(2,l)-hs(1,l))*rdeli
         yp2 = (hs(nr,l)-hs(nr-1,l))*rdeli
         call spline(x,hs(1,l),nr,yp1,yp2,d,cj(1,l))
         yp1=(hsder(2,l)-hsder(1,l))*rdeli
         yp2 = (hsder(nr,l)-hsder(nr-1,l))*rdeli
         call spline(x,hsder(1,l),nr,yp1,yp2,d,cy(1,l))
   20 continue
      return
      end
