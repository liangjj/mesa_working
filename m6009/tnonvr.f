*deck @(#)tnonvr.f	1.1 9/8/91
c***begin prologue     tnonv
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           tnonv, link 1106, kohn variational
c***author             schneider, barry (lanl), rescigno, tom(llnl)  
c***source             m1106
c***purpose            find non-variational t-matrix
c***description        
c***                   
c***                   
c***                  
c***references         schneider and rescigno, physical review
c
c***routines called    iosys, util and mdutil
c***end prologue       tnonv
      subroutine tnonvr (vff,tmat,seig,svec,dum,dir,ntchn)
      complex*16 vff, tmat, seig, svec, catan, cfac
      character *80 title
      logical dir
      real*8 dum, phase, epsum, rtest, atan2, impart, repart 
      dimension vff(ntchn,ntchn), tmat(ntchn,ntchn), seig(ntchn)
      dimension svec(ntchn,ntchn), dum(3*ntchn)
      common /io/ inp, iout
      cfac=dcmplx(0.d0,2.d0)
      write (iout,220)
      if (.not.dir) then
          title='t-matrix'
         else
          title='k-matrix'
      endif
      call prntcm(title,vff,ntchn,ntchn,ntchn,ntchn,iout)
      call cc2opy(vff,tmat,ntchn*ntchn)
      if (.not.dir) then
         do 150 i=1,ntchn
            do 140 j=1,ntchn
               tmat(i,j)=cfac*tmat(i,j)
  140       continue
  150    continue
         do 160 j=1,ntchn
            tmat(j,j)=tmat(j,j)+1.d+00
  160    continue   
      endif
      job=1
      call cgeev (tmat,ntchn,ntchn,seig,svec,ntchn,dum,job,info)
      write (iout,200)
      if (.not.dir) then
         epsum=0.d+00
         do 170 i=1,ntchn
            impart=imag(seig(i))
            repart=real(seig(i))
            phase=atan2(impart,repart)
            phase=phase*.5d+00
            rtest=seig(i)*conjg(seig(i))
            write (iout,210) i,phase,rtest
            epsum=epsum+phase
  170    continue
      else
         epsum=0.d+00
         do 180 i=1,ntchn
            phase=catan(seig(i))
            write (iout,210) i,phase,rtest
            epsum=epsum+phase
  180    continue
      endif
      write (iout,300) epsum
      return
  220 format (/,10x,'variationally uncorrected results')
  200 format (/,10x,' eigenphases of s matrix')
  210 format (/,2x,'phase no.',1x,i3,2x,'phase =',2x,e15.8,2x,'modulus =
     1 ',2x,f10.5)
  300 format(/,5x,'eigenphase sum:',1x,e15.8)
      end
