*deck wrvec
      subroutine wrvec(title,vec,n,iout)
      implicit integer (a-z)
      character *80 title
      real *8 vec
      dimension vec(n)
      write(iout,1) title
      write(iout,2) (vec(i),i=1,n)
    1 format(/,a80)
    2 format( (/,5x,4(e15.8,1x) ) )
      return
      end
