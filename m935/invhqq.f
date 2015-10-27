*deck @(#)invhqq.f	5.1  11/6/94
      subroutine invhqq(hqq,t,grad,npvec,npdim,mdim,energy,filtyp)
      implicit integer (a-z)
      real*8 hqq(npdim,*),t(npdim,npvec),grad(*),energy
      real*8 det(2)
      character*(*) filtyp
      common /io/ inp,iout
c
c     ----- direct inversion of h(qq) in core -----  
c
      twalks=mdim-npdim
c
      ix=1
      jx=1
      do 1 i=1,npvec
         call scopy(npdim,grad(jx),1,grad(ix),1)
         ix=ix+npdim
         jx=jx+mdim
 1    continue
c
      call iosys('read real hqq from '//filtyp,npdim*npdim,hqq,
     $            0,' ')
      do 2 i=1,npdim
         hqq(i,i)=hqq(i,i)-energy
 2    continue
c
      call scopy(npdim*npvec,grad,1,hqq(1,npdim+1),1)
c
      call rminv(npdim,npdim,hqq,npvec,hqq(1,npdim+1),det,'all')
      call scopy(npvec*npdim,hqq(1,npdim+1),1,t,1)
c
c
      return
      end
