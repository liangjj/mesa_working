*deck prntfmn
      subroutine prntfmn(title,a,n,m,ia,ja,iout,frmt)
      implicit integer (a-z)
      real*8 a
      character *80 title
      character *1 frmt
      dimension a(ia,ja)
      write(iout,1) title
    1 format(/,a80,/)
      ibeg=0
      do 10 i=1,m,5
      iend=min(ibeg+5,m)
      ibeg=ibeg+1
      write (iout,20) (ii,ii=ibeg,iend)
   20 format(1x,'col',5(5x,i6,9x))
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
   40 format(4x,5f20.12)
   60 format(4x,5e20.12)
      return
      end
