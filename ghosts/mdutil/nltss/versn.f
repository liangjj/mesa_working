*deck %W%  %G%
      subroutine versn(iout)
c***begin prologue     versn
c***date written       850601  (yymmdd)
c***revision date      861201  (yymmdd)
c
c   1 december 1986   pws at lanl
c        adding machine 4 and setting version to mesa 1.0
c
c***keywords           mesa, version
c***author             martin, richard (lanl)
c***source             %W%   %G%
c***purpose            prints the date and current version of the mesa system.
c***description
c                      call versn(iout)
c                         iout   output file.
c
c***references
c***routines called    (none)
c***end prologue       versn
      implicit integer(a-z)
      parameter (nmx=8)
      character today*24
      character*80 curver,copyr,authrz,line
      character*1 mx,mch(nmx)
      character*16 machin,mchtyp(nmx)
      character site*32
      character*8 user
      data curver/' mesa(%I%);%G%;cray-nltss.'/
      data copyr/'     (c) 1990, the university of california.'/
      data authrz
     $/'     p.w. saxe, b.h. lengsfield iii, r.l. martin, and m. page.'/
      data site/'livermore national laboratory;'/
      data mch/'v','w','x','y','1','2','3','4'/
      data mchtyp/'cray-1','cray-1','cray-1','cray-1','cray-xmp/24',
     $            'cray-xmp/48','cray-xmp/48','cray-xmp/416'/
c
 1000 format(a80)
 1010 format(5x,80a1)
c
c     write the version information.
      write(iout,1000) curver
      write(iout,1000) copyr
      write(iout,1000) authrz
c
c     accumulate and write operating conditions.
      call dattim(today)
      write(iout,1010) (today(i:i),i=1,len(today))
      line=site
c     call mach(mx)
c     call locase(mx,mx)
c     do 10 i=1,nmx
c        if(mx.eq.mch(i)) machin=mchtyp(i)
c  10 continue
c..bhl
      mx='e'
      machin=mchtyp(6)
c..bhl
      pos=cskipb(line,' ')+1
      line(pos:)=machin
      pos=cskipb(line,' ')+1
      line(pos:)='('//mx//');'
      call usrnam(user)
      pos=cskipb(line,' ')+1
      line(pos:)=user
      write(iout,1010) (line(i:i),i=1,len(line))
c
c
      return
      end
