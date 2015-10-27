*deck @(#)sitout.f	5.1  11/6/94
      subroutine sitout(site,zl,nsite,ident,nbegin,label,nbond)
      implicit none
c
      integer nsite,nbegin,nbond
      integer ident(2,nsite)
      character*5 label(nsite)
      real*8 site(3,nsite),zl(3,nsite)
c
      integer j,ii,jj,nend
      integer inp,iout
      common/io/inp,iout
c
 1000 format(x/' listing of atomic and bond sites')
 1010 format(' atomic sites')
 1020 format(x,a5,3f10.5)
 1040 format(' bond sites and bond direction cosines'/,
     $                 12x,'     x         y         z'
     $               //'         lx        ly        lz')
 1060 format(x,a5,x,a5,6f10.5)
c
      write(iout,1000)
c     --- go to bond listing if nbond=1
      if(nbond.eq.0) then
c        --- list atomic sites
         write(iout,1010)
         do 10 j=1,nsite
            write(iout,1020) label(j),site(1,j),site(2,j),site(3,j)
   10    continue
      else if(nbond.eq.1) then
         write(iout,1040)
         do 20 j=nbegin,nsite
            ii=ident(1,j)
            jj=ident(2,j)
            write(iout,1060) label(ii),label(jj),site(1,j),site(2,j),
     $                       site(3,j),zl(1,j),zl(2,j),zl(3,j)
   20    continue
      else if(nbond.eq.2) then
c        --- if nbond=2 set nend for end of saved atom list
         nend=nbegin-1
c        --- list atomic sites
         write(iout,1010)
         do 30 j=1,nend
            write(iout,1020) label(j),site(1,j),site(2,j),site(3,j)
   30    continue
         write(iout,1040)
         do 40 j=nbegin,nsite
            ii=ident(1,j)
            jj=ident(2,j)
            write(iout,1060) label(ii),label(jj),site(1,j),site(2,j),
     $                       site(3,j),zl(1,j),zl(2,j),zl(3,j)
   40    continue
      endif
c
c
      return
      end
