*deck @(#)bfprnt.f	1.1 9/8/91
c***begin prologue     bfprnt
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            print bound-free matrix elements
c***references
c
c***routines called
      subroutine bfprnt(ovbf,hpvb,nlm,nchan,dimc,maxlm,ncon,exdim)
      implicit integer (a-z)
      common /io/ inp, iout
      real *8 rowv
      character *80 title
      character *3 itoc
      character *8 rowt, colt
      complex *16 ovbf, hpvb
      dimension ovbf(maxlm,ncon,nchan), hpvb(maxlm,ncon,nchan,exdim)
      dimension nlm(dimc)
      colv=-99
      rowv=-99.d0
      do 10 ch1=1,nchan
         do 20 ch2=1,exdim
            title='hpvb matrix ch1-'//itoc(ch1)//' ch2-'//itoc(ch2)
            call cmprir(hpvb(1,1,ch1,ch2),rowv,colv,nlm(ch1),
     1                  ncon,maxlm,ncon,title,rowt,colt,iout)
   20    continue
   10 continue
      do 50 ch1=1,nchan
         title='bound free overlap matrix ch1-'//itoc(ch1)
         call cmprir(ovbf(1,1,ch1),rowv,colv,nlm(ch1),ncon,
     1               maxlm,ncon,title,rowt,colt,iout)
   50 continue
      return
      end
