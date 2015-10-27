*deck gngr.f
      subroutine gngr(rmin,rcol,rmtch,stp,n,nptcol,nptmax)
      implicit integer (a-z)
      real*8 rmin, rcol, rmtch, stp, rtst
      common /io/ inp, iout
      nptcol=0
      rtst=rmin
      do 10 i=1,n
         if(rtst.ge.rcol) then
            go to 20
         else
            nptcol=nptcol+1
            rtst=rtst+stp  
         endif
 10   continue
 20   nptmax=0
      rtst=rmin
      do 30 i=1,n
         if (rtst.gt.rmtch) then
            go to 40
         else
            nptmax=nptmax+1
            rtst=rtst+stp  
         endif
 30   continue
 40   nptmax=nptmax-1
      write(iout,*) nptmax
      return
      end
