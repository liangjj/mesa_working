*deck @(#)vecout.f	5.1  11/6/94
      subroutine vecout(a,nad,nbd)
      implicit real*8 (a-h,o-z)
      common /io/ inp,iout
c
      dimension a(nad,nbd)
      character*1 line(126)
      data line /126*'-'/
      save line
    1 format(2x,10(7x,i5))
    2 format(2x,21a6)
    3 format(2x,i2,2x,10f12.7)
    4 format(/,' eval  ',10f12.7)
   44 format(/,' occ   ',10f12.7)
    5 format(/)
c
      n=nbd
      ii=0
      jj=0
  200 ii=ii+1
      jj=jj+1
      kk= 5*jj
      nn=n
      if(n.gt.kk) nn=kk
      ll=2*(nn-ii+1)+1
      write(iout,1) (i,i=ii,nn)
cc    write(6,4) (b(j),j=ii,nn)
cc    write(6,44) (c(j),j=ii,nn)
      write(iout,2) (line(i),i=1,ll)
      do 101 i=1,nad
  101 write(iout,3) i,(a(i,j),j=ii,nn)
      if(n.le.kk) go to 201
      write(iout,5)
      ii=kk
      go to 200
  201 return
      end
