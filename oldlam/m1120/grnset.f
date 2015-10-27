      subroutine grnset(gr,cj2,en,rk,charge,pt,ptend,lval,nptmx,type)
      implicit integer(a-z)
      logical logkey
      real *8 gr, en, rk, charge, pt, ptend, cj1, cj2
      real *8 f, fp, g, gp
      character *10 type
      common /io/ inp, iout
      dimension gr(nptmx), pt(nptmx)
*
*
*               this routine returns either the regular function
*               f if type is function or the function g defined
*               as the linear combination of regular and irregular
*               with zero derivative at the point ptend.
*               to within a sign the product of these will be
*               the greens function.
      if (type.eq.'function') then
          cj1=0.e+00
          cj2=1.e+00
      elseif (type.eq.'derivative') then
          cj1=1.e+00
          call grncal(ptend,lval,en,charge,f,fp,g,gp)
          cj2=-gp/fp
      elseif (type.eq.'standard') then
          call grncal(ptend,lval,en,charge,f,fp,g,gp)
          gr(1)=f
          gr(2)=fp
          gr(3)=g
          gr(4)=gp
          cj1=(gp*f-fp*g)
          cj2=f*g/cj1
          g=g-(gp/fp)*f
          cj2=-f*g/cj1
          return
      else
          call lnkerr ('error in grncal call')
      endif
*
      do 10 ir=1,nptmx
         call grncal (pt(ir),lval,en,charge,f,fp,g,gp)
         gr(ir)=(cj1*g+cj2*f)/sqrt(rk)
   10 continue
*
*
      return
*
c
      end
