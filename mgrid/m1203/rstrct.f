*deck rstrct
c***begin prologue     rstrct
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            restrict fine grid to coarse grid solution
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       rstrct
      subroutine rstrct(uc,uf,nc,nf)
c
      implicit integer (a-z)
      real*8 uc, uf, half, twelth
      character*80 title
      dimension uc(nc,nc,nc), uf(nf,nf,nf)
      common/io/inp, iout
      data half, twelth/.5d0,.083333333333333333333333d0/
c     do all interior points-that is all points not on face of cube      
      do 10 kc=2,nc-1
         kf=2*kc-1
         do 20 jc=2,nc-1
            jf=2*jc-1
            do 30 ic=2,nc-1
               if=2*ic-1
               uc(ic,jc,kc)=half*uf(if,jf,kf) +
     1                twelth*( uf(if+1,jf,kf) + uf(if-1,jf,kf) + 
     2                         uf(if,jf+1,kf) + uf(if,jf-1,kf)    +
     3                         uf(if,jf,kf+1) + uf(if,jf,kf-1) )
 30         continue   
 20      continue
 10   continue
      nfx=2*nc-1
c     do lower and upper (x,y) planes      
      do 40 ic=1,nc
         if=2*ic-1
         do 50 jc=1,nc
            jf=2*jc-1
            uc(ic,jc,1)=uf(if,jf,1)
            uc(ic,jc,nc)=uf(if,jf,nfx)
 50      continue
 40   continue
c     do left and right (x,z) plane
      do 60 ic=1,nc
         if=2*ic-1
         do 70 kc=1,nc
            kf=2*kc-1
            uc(ic,1,kc)=uf(if,1,kf)
            uc(ic,nc,kc)=uf(if,nfx,kf)
 70      continue
 60   continue
c     do back and front (y,z) plane
      do 80 jc=1,nc
         jf=2*jc-1
         do 90 kc=1,nc
            kf=2*kc-1
            uc(1,jc,kc)=uf(1,jf,kf)
            uc(nc,jc,kc)=uf(nfx,jf,kf)
 90      continue
 80   continue                                
c      title='coarse solution after injection'
c      call prntrm(title,uc,nc,nc,nc,nc,iout)
      return
      end





