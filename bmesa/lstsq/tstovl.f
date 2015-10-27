*deck tstovl.f
      subroutine tstovl(vec,n,nvc)
      implicit integer (a-z)
      character*80 title
      real*8 vec, ovl, sdot
      dimension vec(n,nvc)
      common/io/inp, iout
      do 10 i=1,nvc
         do 20 j=1,i
            ovl = sdot(n,vec(1,i),1,vec(1,j),1)
            write(iout,1) i, j, ovl
   20    continue
   10 continue 
      return
    1 format(1x,'i = ',i3,1x,'j = ',i3,1x,'overlap = ',e15.8)
      end       
