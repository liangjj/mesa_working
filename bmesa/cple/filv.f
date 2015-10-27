*deck filv
      subroutine filv(v,nic,nio,njc,njo,ni,nj,n)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 v
      dimension v(n,n)
      character*80 title
      character*1 itoc
      nitot=nic+nio
      njtot=njc+njo
      do 10 i=1,nitot
         do 20 j=1,njtot
            v(i,j)=.5d0
   20    continue
   10 continue
      title='v matrix for channels = ('//itoc(ni)//','//itoc(nj)//')'
      call prntrm(title,v,nitot,njtot,n,n,iout)
      return
      end
