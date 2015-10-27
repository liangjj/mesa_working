*deck @(#)lonep.f	5.1  11/6/94
      subroutine lonep(ndip,nsite,site,zl,j)
      implicit none
c
      integer j,nsite
      integer ndip(5,nsite)
      real*8 site(3,nsite),zl(3,nsite)
c
      integer ntype,i,i1,i2,i3,i4
      real*8 xwork(3,10),half
c
      parameter (half=0.5d+00)
c
      ntype=ndip(1,j)
      i1=ndip(2,j)+j
      i2=ndip(3,j)+j
      i3=ndip(4,j)+j
      i4=ndip(5,j)+j
      call vs(site(1,i2),site(1,i1),xwork(1,1))
      call norm(xwork(1,1),xwork(1,1))
      call vs(site(1,i3),site(1,i1),xwork(1,2))
      call norm(xwork(1,2),xwork(1,2))
c     --- ntype=1 is trigonal, ntype=2 is tetrahedral
      if(ntype.eq.1) then
c        --- trigonal type; three atoms are given
         do 20 i=1,3
            zl(i,j)=-half*(xwork(i,1)+xwork(i,2))
   20    continue
         call norm(zl(1,j),zl(1,j))
      else if(ntype.eq.2) then
c        --- tetrahedral type; four atoms are given
         call vs(site(1,i4),site(1,i1),xwork(1,3))
         call norm(xwork(1,3),xwork(1,3))
c        --- find plane 1,2,3
         call vs(xwork(1,2),xwork(1,1),xwork(1,4))
         call vs(xwork(1,3),xwork(1,1),xwork(1,5))
         call vxv(xwork(1,4),xwork(1,5),xwork(1,6))
         call norm(xwork(1,6),xwork(1,6))
c        --- reverse the plane normal vector
         do 40 i=1,3
            zl(i,j)=-xwork(i,6)
   40    continue
      endif
c
c
      return
      end
