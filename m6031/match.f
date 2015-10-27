*deck match
      subroutine match(hfbox,eig,rbox,energy,l,v0,ncon,nen)
      implicit integer (a-z)
      real*8 hfbox, eig, energy, rmta, drmta, rmte, drmte, v0
      real*8 rbox, k, reg, ireg , dreg, direg, sn, cn
      real*8 tandel, phasea, phasee
      dimension eig(ncon), hfbox(ncon), energy(nen)
      common /io/ inp, iout
      write(iout,1)
      do 10 ien=1,nen
         call rmat(rmta,drmta,rmte,drmte,hfbox,eig,rbox,energy(ien),
     1             v0,l,ncon,1,.false.)          
   20    continue
         k=sqrt(2.d0*energy(ien))
         sn=sin(k*rbox)
         cn=cos(k*rbox)
         reg=sn
         dreg=k*cn
         ireg=-cn
         direg=k*sn
         if (l.ne.0) then
             reg=sn/(k*rbox) - cn
             ireg=-sn - cn/(k*rbox)
             dreg=cn/rbox - sn/(k*rbox*rbox) +k*sn
             direg=-k*cn +sn/rbox +cn/(k*rbox*rbox)
         endif
         tandel=-(reg-rmte*dreg)/(direg*rmte-ireg)   
         phasee=atan(tandel)
         tandel=-(reg-rmta*dreg)/(direg*rmta-ireg)
         phasea=atan(tandel) 
         write(iout,2) energy(ien), phasee, phasea, tandel
   10 continue            
    1 format(/,3x,'    energy    ','   ex. phase  ','     ap. phase   ',
     1       '   tan del   ')
    2 format(1x,4e15.8) 
      return
      end
