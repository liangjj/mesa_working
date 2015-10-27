      subroutine mcrt1a(nft,nsym,nob,ltrb1,buf)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
c     implicit real*8(a-h,o-p,r-z),integer*2(q)
cc
cmp   extended dummy buf
cc
      dimension nob(2),buf(2),ihd(3)
c
      call lnkerr(' mcrt1a entered ')
c
c.bl
c.bl      do 30 i=1,nsym
c.bl  nt=(nob(i)*(nob(i)+1))/2
c.bl  if (nt.eq.0) go to 30
c.bl  len=ltrb1
c.bl  do 20 j=1,nt,ltrb1
c.bl  n=j+ltrb1-1
c.bl  if(n.gt.nt) len=nt-j+1
cc.ibm    read (nft) ia,ib,ic,(buf(ix+l),l=1,len)
c.bl  call sread(nft,ihd,3)
c.bl  call sread(nft,buf(ix+1),intowp(len))
c.bl   20 ix=ix+len
c.bl   30 continue
c
      return
      end
