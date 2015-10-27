      subroutine phase(rmat,ek,rad,charge,typ)
      implicit integer (a-z)
      dimension gr(4)
      common /io/ inp, iout
      character *8 typ
      character *10 type
      real *8 rmat, ek, rad, pshift, energy, sn, cn, coef, charge
      real *8 ratio, ptend, pt, gr
      lval=0
      pt=rad
      ptend=rad
      nptmx=4
      energy=ek*ek
      type='standard'
      call grnset(gr,ratio,energy,ek,charge,pt,ptend,lval,
     1            nptmx,type)
      coef=(rmat*gr(2)-gr(1))/(gr(3)-rmat*gr(4))
      pshift=atan (coef)
      write (iout,10) energy,pshift
      write (iout,20) rmat, ratio
   10 format (/,10x,'phase shift at energy',1x,e15.8,/,10x,
     1        e15.8)
   20 format (/,10x,'r matrix',1x,e15.8,2x,'ratio',1x,e15.8)
      return
      end
