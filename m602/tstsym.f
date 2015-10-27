*deck @(#)tstsym.f	1.2  7/30/91
      subroutine tstsym(hold,cvec,s,temp,tt,num,ndeg,test)
      implicit integer(a-z)
      real*8 hold(num,ndeg),cvec(num,ndeg),tt(ndeg,ndeg),
     1 s(num,num),temp(num,ndeg),test,xx
c
      common /io/ inp,iout
c
      call ebc(temp,s,hold,num,num,ndeg)
c
      call ebtc(tt,cvec,temp,ndeg,num,ndeg)
c
      do 1 i=1,ndeg
      xx=0.d0
       do 2 j=1,ndeg
       xx=xx+tt(i,j)**2
  2    continue
       if(abs(xx-1.d0).gt.test) then
        write(iout,*)' projection error '
        call lnkerr(' m604: tstsym ')
       end if
  1   continue
c
      return
      end
