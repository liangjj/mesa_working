*deck %W%  %G%
      subroutine mcrt1a(nft,nsym,nob,ltrb1,buf)
C
C***Begin prologue
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
C
C***References
C
C***Routines called    (none)
C
C***End prologue
C
      implicit real*8(a-h,o-z)
cc
cmp   extended dummy buf
cc
      dimension nob(2),buf(2),ihd(3)
c
      call lnkerr(' MCRT1A entered ')
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
