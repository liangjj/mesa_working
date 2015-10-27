*deck %W%  %G%
      subroutine versn(iout)
c
c***begin prologue     versn
c***date written       850601  (yymmdd)
c***revision date      870207  (yymmdd)
c
c   1 march    1990   rlm at lanl
c        making bsd 4.3 unix version for sun 4/260 workstations.
c
c   7 february 1987   pws at lanl
c        making bsd 4.2 unix version on sun 3/50 and 3/160 workstations
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
c
      implicit integer(a-z)
c
      parameter (nmx=6)
      character today*24
      character*80 curver,copyr,authrz,line
      character*32 machin,mch(nmx)
      character*32 mchtyp(nmx),mx
      character site*32
      character*8 user
      data curver/' mesa(%I%);%G%;cos'/
      data copyr/'     (c) 1990, the university of california.'/
      data authrz
     $/'     p.w. saxe, b.h. lengsfield iii, r.l. martin, and m. page.'/
      data site/'los alamos national laboratory;'/
c
      data mch /'t12-server','t12-c1','t12-c2','t12-c3','t12-c4',
     $          'vivaldi'/
      data mchtyp /'sun 3/160',4*'sun 3/50','sun 4/280s'/
c
 1000 format(a80)
 1020 format(a80,/)
 1010 format(5x,80a1)
c
c     write the version information.
c
      write(iout,1000) curver
      write(iout,1000) copyr
      write(iout,1020) authrz
c
c     accumulate and write operating conditions.
c
      call dattim(today)
      write(iout,1010) (today(i:i),i=1,len(today))
      line=site
      if (hostnm(mx).ne.0) mx='unknown'
c
      machin='unknown'
      do 10 i=1,nmx
         if(mx.eq.mch(i)) machin=mchtyp(i)
   10 continue
      pos=cskipb(line,' ')+1
      line(pos:)=machin
      pos=cskipb(line,' ')+1
      line(pos:)='('//mx
      pos=cskipb(line,' ')+1
      line(pos:)=');'
      call usrnam(user)
      pos=cskipb(line,' ')+1
      line(pos:)=user
      write(iout,1010) (line(i:i),i=1,len(line))
c
c
      return
      end
