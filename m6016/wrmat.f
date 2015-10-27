      subroutine wrmat(a,n,tit,title)
      implicit integer (a-z)
      real *8  a
      character *24 tit
      character *3 ans
      character *80 title
      dimension a(*)
      common /io/ inp,iout
      call iosys ('does '//tit//' exist on tmat',0,0,0,ans)
      if (ans.eq.'yes') then
          call iosys ('read real '//tit//' from tmat',2*n*n,a,0,' ')
          call prntcm(title,a,n,n,n,n,iout)
      else
          write (iout,10) tit
      endif
      return
   10 format(/,5x,'file:',a16,' not on tmat')
      end
