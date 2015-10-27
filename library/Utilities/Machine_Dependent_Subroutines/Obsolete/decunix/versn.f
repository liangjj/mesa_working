*deck @(#)versn.f	1.2  5/30/91
      subroutine versn(iout)
c
c***begin prologue     versn
c***date written       850601  (yymmdd)
c***revision date      870207  (yymmdd)
c
c   2 february 1992   bis at nsf
c        finalizing hpux unix version for hp 9000/700 workstations
c        use being made of optimizing pre-processor and vector library
c        where possible. some links will not run properly and the makefiles
c        reflect that.
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
c***source             @(#)versn.f	1.2   5/30/91
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
      parameter (nmx=3)
      character today*24
      character*80 curver,copyr,authrz,line
      character*32 machin,mch(nmx)
      character*32 mchtyp(nmx),mx
      character*32 site
      character*16 user
      data curver/' mesa(1.2);5/30/91;unix.'/
      data copyr/'     (c) 1990, the university of california.'/
      data site /'National Science Foundation;'/
c
c
 1000 format(a80)
 1020 format(a80,/)
 1010 format(5x,80a1)
      authrz ='     P.W. Saxe, B.H. Lengsfield iii, R.L. Martin, '//
     $              'M. Page and B. Schneider.'
      mch(1) = 'bohr'
      mch(2) = 'noether'
      mch(3) = 'feynman'
      mchtyp(1) = 'Dec Alpha 600MHz'
      mchtyp(2) = 'Dec Alpha 500MHz'
      mchtyp(3) = 'HP 735 100MHz'
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
c     if (hostnm(mx).ne.0) mx='unknown'
c
      mx='bohr'
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
