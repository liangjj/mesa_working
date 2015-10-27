*deck @(#)pm602.f	1.5 9/3/91 
      program m602
c***begin prologue     pm602
c***date written       901230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m602, link 601, rbitals
c***author             schneider barry (lanl)
c***source             @(#)pm602.f	1.3 7/30/91
c***purpose            project symmetry orbitals.
c***
c***routines called
c***end prologue       pm602
      implicit integer(a-z)
c
      real*8  z, z1
      integer a, a1
      logical scat
      character*4096 ops
      character*32 xform
      character*128 nmfile
      character*128 chk
      real*8 tstmo, smllmo, fpkey, smleig, xx, maxerr, tol
      character*8 reduced
      character*8 chrkey, scatyp, filtyp
      character*3 ans 
      logical logkey, prn, symtre
      dimension prn(10)
      common /io/ inp, iout
      pointer(p,z(1)), (p,a(1))
      pointer(p1,z1(1)), (p1,a1(1))
c
    2 format(/,1x,'m602: symmetry projection')
    3 format(5x,'memory use',18x,i9)
c
c     ----- recover the options string -----
c
c      call manmem(0,idum,idum,'m602',idum)
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
c
      prn(1)=logkey(ops,'print=sym=old-salcs',.false.,' ')
      prn(2)=logkey(ops,'print=sym=new-salcs',.false.,' ')
      prn(3)=logkey(ops,'print=sym=vectors',.false.,' ')
      prn(4)=logkey(ops,'print=sym=overlap-eigenvalues',.false.,' ')
      prn(5)=logkey(ops,'print=sym=overlap-matrix',.false.,' ')
      prn(6)=logkey(ops,'print=sym=all',.false.,' ')
      if(prn(6)) then
         call setprn(prn,5)
      endif
      smllmo=fpkey(ops,'sym=smallmo',1.d-02,' ')
      tstmo=fpkey(ops,'sym=testmo',1.d-06,' ')
      smleig=fpkey(ops,'sym=smleig',1.d-05,' ')
      symtre=logkey(ops,'sym=reorder-orbitals',.false.,' ')
      tol=fpkey(ops,'sym=overlap-tolerance',1.d-08,' ')
      maxerr=fpkey(ops,'sym=tol',1.d-06,' ')
      call iosys ('read integer "number of basis functions" from rwf',1,
     1             num,0,' ')
      write(iout,2)
      if(symtre) then
         nsmall=intkey(ops,'internal-orbitals',0,' ')
         if(nsmall.eq.0) then
            call iosys('read integer "number of alpha electrons" from'//
     1                 ' rwf',1,nsmall,0,' ')
         endif
         nvirt=num-nsmall 
         write(iout,909) nsmall, nvirt
      endif
      scat=.false.
      scatyp=chrkey(ops,'scattering','none',' ') 
c
c     we need to ignore trailing characters in determining which
c     type of calculation we are performing.
c
      if(scatyp(1:4).eq.'kohn') then
         call iosys('read character "kohn filename" from rwf',
     #               0,0,0,nmfile)
         filtyp='kohn'
         filtyp=filtyp(1:4)
         scatyp='kohn'
         scat=.true.
      else if(scatyp(1:8).eq.'r-matrix') then
         call iosys('read character "r-matrix filename" from rwf',
     #               0,0,0,nmfile)
         filtyp='rmtrx'
         filtyp=filtyp(1:5)
         scatyp='r-matrix'
         scat=.true.
      endif    
c
c     ----- get the dimensions, etc -----
c

      call iosys ('does "old nbf" exist on rwf',0,0,0,ans)
      if (ans.ne.'no') then
          call iosys('read integer "old nbf" from rwf',1,oldnum,0,' ')
      else
          oldnum=num
      endif
      call iosys('read integer "number of irreducible'//
     1           ' representations" from rwf',1,nsym,0,' ')
c
      if(nsym.gt.20) then
         call lnkerr('m602: nsym must be lt 20')
      endif
c
      nnp=(num+1)*num/2
c
c
c     ----- allocate core -----
c
c
c
c
c     ----- get the core needed -----
      write(iout,800) oldnum, num
c
c     allocate core for main part of calculation
c
      salcn=1
      cvec=salcn+num*num
      s=cvec+num*num
      slcsm=wpadti(s+num*num)
      nummo=slcsm+nsym+nsym
      molst=nummo+nsym
      aosym=molst+nsym*num+nsym*num
      aolst=aosym+num
      eval=iadtwp(aolst+num)
      t1=eval+num+num 
      t2=t1+num*num
      t3=t2+num*num
      t4=t3+num*num
      t5=t4+num*num
      need=wpadti(t5+num*num)
      call getmem(need,p,ngot,'m602',0)
c
c     ----- calculate the salcs in the reduced basis -----
c     -----       if needed -----
c
      ans='no' 
      call iosys ('does reduced exist on rwf',-1,0,0,ans)
      if (logkey(ops,'drop',.false.,' ')) then
          if (ans.eq.'no') then
              salco=1
              pkindx=wpadti(salco+oldnum*oldnum) 
              wds=pkindx+oldnum
              call getmem(wds,p1,ngot1,'tm602',0)
              call nsalc(z1(salco),z(salcn),a1(pkindx),nsym,a(slcsm),
     1                   a(slcsm+nsym),oldnum,num,prn)
              call getmem(-ngot1,p1,idum,'tm602',idum)
          else
              call iosys('read integer "number of symmetry orbitals" '//
     1                   'from rwf',nsym,a(slcsm),0,' ')
              call iosys('read real "salc transformation matrix" '//
     1                   'from rwf',num*num,z(salcn),0,' ')
          endif
      else 
          call iosys('read integer "number of symmetry orbitals" from'//
     1               ' rwf',nsym,a(slcsm),0,' ')
          call iosys('read real "salc transformation matrix" from rwf',
     1                num*num,z(salcn),0,' ')
      endif
      if (scat) then
          call iosys('read character "checkpoint filename" from rwf',
     1               -1,0,0,chk)
          call iosys('open chk as unknown',0,0,0,chk)
          call detsym(z(salcn),nsym,a(aosym),a(slcsm),a(aolst),num)
          call iosys ('close chk',0,0,0,' ')
      endif
c
c     ----- pick up the transformation vectors from the desired -----  
c     -----                    file             -----
c
      if(logkey(ops,'sym='//scatyp,.false.,' ')) then
         kw=1
         write(iout,910)
         call iosys('open '//filtyp//' as old',0,0,0,nmfile)
         call iosys('read character "transformation vector" from '
     #               //filtyp,-1,0,0,xform)
         call iosys('read real '//xform//' from '//filtyp,
     #               num*num,z(cvec),0,' ')
         call iosys('close '//filtyp,0,0,0,' ')
      elseif(logkey(ops,'sym=chk',.false.,' ')) then
         kw=-1
         call iosys('read character "checkpoint filename" from rwf',
     1              -1,0,0,chk)
         call iosys('open chk as old',0,0,0,chk)
         call iosys ('read real "scf vector" from chk',num*num,
     1                z(cvec),0,' ')
         call iosys('read real "orbital energies" from chk',
     1              num,z(eval),0,' ')
         call iosys('close chk',0,0,0,' ')
      else
         kw=0
         write(iout,911)
         call iosys('read character "checkpoint filename" from rwf',
     1              -1,0,0,chk)
         call iosys('open chk as unknown',0,0,0,chk)
         call iosys ('read real "scf vector" from rwf',num*num,
     1                z(cvec),0,' ')
         call iosys('read real "orbital energies" from rwf',
     1               num,z(eval),0,' ')
      endif
c
c----------------------------------------------------------------------c
c           calculate the symmetry orbitals                            c
c----------------------------------------------------------------------c
      call symorb(z(cvec),z(salcn),z(t1),z(s),z(eval),z(t2),z(t3),z(t4),
     1            z(t5),tol,maxerr,a(slcsm),a(nummo),a(molst),nsym,num,
     2            symtre,nsmall,prn(3))
c      call ebc(z(t1),z(s),z(cvec),num,num,num)
c      call ebtc(z(t2),z(cvec),z(t1),num,num,num)
c      call matout(z(t2),num,num,num,num,iout)
c
      if(.not.logkey(ops,'sym=nowrite',.false.,' ')) then
         if(kw.le.0) then
            write(iout,801)
            call iosys ('write real "scf vector" to rwf',num*num,
     #                   z(cvec),0,' ')
            call iosys('write real "orbital energies" to rwf',num,
     1                  z(eval),0,' ')
            call iosys('write real "guess vector" to rwf',num*num,
     1                  z(cvec),0,' ')
            call iosys('write real "guess orbital energies" to rwf',
     1                  num,z(eval),0,' ')
         else
            write(iout,802)
            call iosys('open '//filtyp//' as old',0,0,0,nmfile)
            call iosys('write real '//xform//' to '//filtyp,
     #                  num*num,z(cvec),0,' ')
            call iosys('close '//filtyp,0,0,0,' ')
         end if
      end if
c
      call getmem(-ngot,p,idum,'m602',idum)
c
c     ----- and exit gracefully -----
c
      call chainx(0)
c
c      stop
 800  format(/,1x,'number of old atomic orbitals = ',i4,/1x,
     1            'number of new atomic orbitals = ',i4)
 801  format(/,' vectors have been written to the rwf file')
 802  format(/,' vectors have been written to the scattering file')
 909  format(/,5x,'number of internal orbitals = ',i4,/,5x,
     1         'number of virtual orbitals  = ',i4,/,2x,
     2         'they will be symmetry ordered ')
 910  format(/,5x,' reading the scattering file ')
 911  format (/,5x,' reading the rwf file ')
      end














