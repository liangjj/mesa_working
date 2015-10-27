*deck @(#)wrtpot.f	1.1 9/7/91
c***begin prologue     wrtpot
c***date written       xxxxxx   (yymmdd)
c***revision date      890409   (yymmdd)
c***keywords           m6002, link 6003, v integral
c***authors            schneider, barry (lanl)
c***source             m6002
c***purpose            write and print out static potential
c***references       
c***routines called    iosys
c***end prologue       wrtpot
      subroutine wrtpot(valint,grid,npnts,nstri,nptmx,nwrite,
     1                 nwds,ipass,prnt)
      implicit integer (a-z)
      real *8 valint, grid
      logical prnt
      dimension valint(nptmx,nstri), grid(4,nptmx)
      common /io/ inp, iout
      do 10 ist=1,nstri
         call iosys ('write real "static potential" to vstat without '//
     1               'rewinding',npnts,valint(1,ist),0,' ')
         write (iout,100) ipass, ist, npnts
         nwds=nwds+npnts
         nwrite=nwrite+1
   10 continue
      if (prnt) then
          call prntv(valint,grid,nstri,npnts,ipass)
      endif
  100 format(5x,i3,11x,i2,10x,i8,8x,i8)
      return
      end
