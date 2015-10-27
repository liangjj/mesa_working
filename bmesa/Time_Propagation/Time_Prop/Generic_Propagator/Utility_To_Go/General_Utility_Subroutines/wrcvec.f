*deck wrcvec
      subroutine wrcvec(title,vec,n,iout)
      implicit integer (a-z)
      character *80 title
      complex *16 vec
      dimension vec(n)
      write(iout,1) title
      write(iout,2) (vec(i),i=1,n)
    1 format(/,a80,/)
    2 format(2(5x,'(',e15.8,1x,',',e15.8,')' ) )
      return
      end
