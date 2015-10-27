*deck @(#)pm812.f	5.1  11/6/94
      subroutine pm812(z,a)
c***begin prologue     m812
c***date written       871117   (yymmdd)
c***revision date      880110   (yymmdd)
c   10 january 1988    bhl at brl
c   mcscf defaulted to .false. for multi-reference ci run
c
c
c***keywords           m812, link 812, density matrix, transformation,
c                      one-electron, two-electron
c***author             saxe, paul (lanl)
c***source             @(#)pm812.f	5.1   11/6/94
c***purpose            transforms the one- and two-particle density
c  matrices from the molecular-orbital to the atomic-orbital basis.
c
c***description
c
c***references
c
c***routines called
c***end prologue       m812
c
      implicit integer (a-z)
c
      character*4096 ops
      character chrkey*32,xform*32
      character*32 occ
      character*128 namint,moden,aoden
      logical logkey
      logical ci
      logical mcscf
      real*8 z(*)
      integer a(*)
c
      common /io/     inp,ioutpt
c
      data maxcor /20000/
      save maxcor
c
c
 1000 format(1x,'m812: density matrix transformation')
 1030 format(5x,'memory use        ',16x,i9)
cdir$ fastmd
c
c     ----- recover the options string -----
c
      ops=' '
      call iosys('read character options from rwf',-1,0,0,ops)
c
c
      ci=logkey(ops,'ci',.false.,' ')
      mcscf=.false.
      if(.not.ci) then
         mcscf=logkey(ops,'mcscf',.false.,' ')
      end if

c
c     process the other options.
c
      xform=chrkey(ops,'ci=tvector',' ',' ')
      if(xform.eq.' ') then
         call iosys('read character "transformation vector" from rwf',
     $        -1,0,0,xform)
      endif
      call iosys('write character "transformation vector" to rwf',
     $            0,0,0,xform)
c
c     ----- work out the name of the corresponding orbital eigenvalues
c           or occupations
c
      call locase(xform,xform)
      if (xform(1:10).eq.'"no vector') then
         occ='"no occ'//xform(10:)
      else if (xform.eq.'"mcscf vector"') then
         occ='"mcscf orbital energies"'
      else
         occ='"orbital energies"'
      end if
c
c     ----- read dimensions etc from wherever they are -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbf,0,' ')
      nnp=(nbf+1)*nbf/2
c
      if (mcscf) then
         call iosys('read integer mc_ncore from rwf',1,ncore,0,' ')
         call iosys('read integer mc_nactive from rwf',1,nactiv,0,' ')
      else
         ncore=0
         nactiv=nbf
      end if
      nnpact=nactiv*(nactiv+1)/2
c
      write(ioutpt,1000)
c
      nsym=1
      write (ioutpt,1050) xform,nbf
 1050 format (5x,'vector used:                   ',a16,/,
     $        5x,'number of basis functions:     ',9x,i4)
c
c     ----- now get the transformation vector -----
c
      c=1
      eigval=c+nbf**2
      e=eigval+nbf
      need=wpadti(e+nbf)
c
      call getscm(need,z,ngot,'m812 xform vector',0)
c
      call iosys('read real '//xform//' from rwf',-1,z(c),0,' ')
      call iosys('read real '//occ//' from rwf',nbf,z(e),0,' ')
c
c     ----- open the integral and density files and check the orthonormality 
c                   of the orbitals
c
      call iosys('read character "mo density filename" from rwf',
     $            0,0,0,moden)
      call iosys('open moden as old',0,0,0,moden)
c
c     create the output density matrix file.
c
      call iosys('read character "ao density filename" from rwf',
     $            0,0,0,aoden)
      call iosys('open aoden as new',0,0,0,aoden)
c
      s=iadtwp(need)
      need=wpadti(s+nnp)
c
      call getscm(need,z,ngot,'m812 one-electron',0)
c
      call iosys('read real "overlap integrals" from rwf',
     $            nnp,z(s),0,' ')
c
      t1=iadtwp(need)
      t2=t1+nbf**2
      t3=t2+nbf**2
      ao1dm=t3+nbf**2
      mo1dm=ao1dm+nnp
      need=wpadti(mo1dm+nnpact)
c
      call getscm(need,z,ngot,'m812 chknrm',0)
c
      call chknrm(z(c),z(s),z(t1),z(t2),z(t3),nbf,nnp)
c
c     ----- transform the one-particle density matrix-----
c
      if (mcscf) then
         call iosys('read real "mcscf mo 1pdm" from moden',
     $        nnpact,z(mo1dm),0,' ')
      else
         call iosys('read real "mo 1pdm" from moden',
     $        nnpact,z(mo1dm),0,' ')
      end if
c
      call trn1dm(z(ao1dm),z(c),nbf,nnp,z(t1),z(t2),z(mo1dm),nactiv,
     $     nnpact,ncore)
c
      if (mcscf) then
         call iosys('write real "mcscf active ao 1pdm" to rwf',
     $        nnp,z(ao1dm),0,' ')
      else
         call iosys('write real "ci ao 1pdm" to rwf',
     $        nnp,z(ao1dm),0,' ')
      end if
c
c     ----- reallocate core for the two-particle transformation -----
c
      t1=c+nbf**2
      t2=t1+nbf**2
      val=t2+nbf**2
      lenbin=max(1920,nnp)
      lab=wpadti(val+lenbin)
      bin=lab+lenbin
      out=iadtwp(bin+lenbin)
      ntriang=min(nnp,30000/nnp+1)
      in=out+nnp*ntriang
      asort=in+nnp*ntriang
      call getscm(0,z,maxcor,'possible',0)
c
      lnsort=min(wptoin(nnp**2),maxcor-wpadti(asort)-100)
      need=wpadti(asort)+lnsort
c
      call getscm(need+100,z,ngot,'m812 two-particle',0)
c
      write(ioutpt,1030) need
      if(ntriang.le.0) call lnkerr('not enough memory')
c
      continue
c
      call trn2dm(z(out),nnp,ntriang,z(c),nbf,z(t1),z(t2),
     #     z(asort),lnsort,nsym,z(val),a(lab),a(bin),lenbin,
     $     z(in),nactiv,nnpact,ncore,ci,mcscf)
c
c     ----- check the energy, if desired -----
c
      if (logkey(ops,'m812=check-energy',.false.,' ')) then
c
         call iosys('read character "integral filename" from rwf',
     $               0,0,0,namint)
         call iosys('open ints as old',0,0,0,namint)
c
         ntriang=min(nnp,30000/nnp+1)
         t1=1
         t2=t1+nbf**2
         dm=t2+nbf**2
         ints=dm+ntriang*nnp
         need=wpadti(ints+ntriang*nnp)
c
         call check(z(t1),z(t2),z(dm),z(ints),nbf,nnp,ntriang,mcscf)
c
         call iosys('close ints',0,0,0,' ')
c
      end if
c
c     ----- and exit with grace -----
c
      call iosys('close moden',0,0,0,' ')
      call iosys('close aoden',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end
