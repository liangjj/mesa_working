*deck @(#)versn.f	5.1 11/6/94
      subroutine versn(iout)
c
c***begin prologue     versn
c***date written       850601  (yymmdd)
c***revision date      920317  (yymmdd)
c
c   17 march   1992   rlm at lanl
c        making hp-ux 8.07 unix version for hp 700 series workstations.
c
c***keywords           mesa, version
c***author             martin, richard (lanl)
c***source             @(#)versn.f	5.1 11/6/94
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
      parameter (nmx=7)
      character today*24
      character*80 curver,copyr,line
      character*72 authrz(2)
      character*32 machin,mch(nmx)
      character*32 mchtyp(nmx),mx
      character site*32
      character*8 user
      data curver/' mesa(5.1);11/6/94;hp-ux 9.01 unix.'/
      data copyr/'     (c) 1990, the university of california.'/
      data authrz(1)/
     $  'b.h. lengsfield iii,r.l. martin,p.w. saxe,t.v. russo,m. page,'/
      data authrz(2)/
     $  'g.j. tawa,b. schneider,m.o. braunstein,p.j. hay,a.k. rappe'/
      data site/'los alamos national laboratory;'/
c
      data mch /'lithium','sodium','rubidium','cesium','kalium',
     $          'francium','platinum'/
      data mchtyp /'hp 735','hp710','hp 705','hp735','hp705','hp735',
     $             'hp735'/
      save curver,copyr,authrz,site,mch,mchtyp
c
 1000 format(a80)
 1010 format(5x,80a1)
 1020 format(5x,a72,/,5x,a72,/)
c
c     write the version information.
      write(iout,1000) curver
      write(iout,1000) copyr
      write(iout,1020) authrz(1),authrz(2)
c
c     accumulate and write operating conditions.
      call dattim(today)
      write(iout,1010) (today(i:i),i=1,len(today))
      line=site
      if (hostnm(mx).ne.0) mx='unknown'
c
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
