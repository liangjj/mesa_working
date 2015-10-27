*deck prntcmn
      subroutine prntcmn(title,a,n,m,ia,ja,iout,frmt)
      implicit integer (a-z)
      complex *16 a
      character *80 title
      character *1 frmt
      dimension a(ia,ja)
      write(iout,1) title
    1 format(/,a80,/)
      ibeg=0
      do 10 i=1,m,2
      iend=min(ibeg+2,m)
      ibeg=ibeg+1
      write (iout,20) (ii,ii=ibeg,iend)
   20 format(1x,'col',2(15x,i6,19x))
      if (frmt.eq.'f') then
          do 30 j=1,n
             write (iout,40) (a(j,k),k=ibeg,iend)
   30     continue
      else
	  do 50 j=1,n
	     write (iout,60) (a(j,k),k=ibeg,iend)
   50     continue
      endif
      ibeg=iend
   10 continue
   40 format(4x,2(2f20.12))
   60 format(4x,2(2e20.12))
      return
      end
