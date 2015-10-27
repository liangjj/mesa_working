*deck @(#)bond.f	5.1  11/6/94
      subroutine bond(jz,site,zl,nsite,ident,nbegin,tsite)
      implicit none
c
      integer nsite,nbegin
      integer jz(nsite)
      real*8 site(3,nsite),zl(3,nsite),tsite(3,nsite)
      integer ident(2,nsite)
c
      integer i,j,nend
      real*8 bmax,dist
      parameter (bmax=1.6d+00)
c
      nend=nbegin+nsite-1
c     ---  nsite will be incremented
      nsite=nbegin-1
c     --- look for bonds; examine site(i,j) from j=nbegin to j=nend
      do 20 i=nbegin,nend-1
         do 10 j=i+1,nend
            call fdist(site(1,i),site(1,j),dist)
            if (dist.le.bmax) then
c              --- no h-h centers
               if ((jz(i).ne.1).or.(jz(j).ne.1)) then
                  nsite=nsite+1
                  call center(site(1,i),site(1,j),tsite(1,nsite),
     $                        zl(1,nsite))
                  ident(1,nsite)=i
                  ident(2,nsite)=j
               endif
            endif
   10    continue
   20 continue
c
c     --- copy back to site
      do 40 j=nbegin,nsite
         do 30 i=1,3
            site(i,j)=tsite(i,j)
   30    continue
   40 continue
c
c
      return
      end
