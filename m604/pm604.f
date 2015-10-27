*deck @(#)pm604.f	5.1 11/6/94 
      subroutine pm604(z,a,maxcor)
c***begin prologue     m604
c***date written       901230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m604, link 601, rbitals
c***author             schneider barry (lanl)
c***source             @(#)pm604.f	1.3 7/30/91
c***purpose            modification of bhl code to project symmetry
c***                   orbitals but in a basis containing dropped 
c***                   functions.
c***
c***routines called
c***end prologue       m604
      implicit integer(a-z)
c
      real*8  z(1)
      integer a(1)
      character*4096 ops
      character*32 xform
      character*128 kohn
      character*128 chk
      real*8 tstmo, smllmo, fpkey, smleig, xx, maxerr
      character*3 ans
      logical logkey, prvec, prold, prnew
      common /io/ inp, iout
c
    2 format(1x,'m604:')
    3 format(5x,'memory use',18x,i9)
c
c     ----- recover the options string -----
c
      call iosys ('read character options from rwf',-1,0,0,ops)
c
      prold=logkey(ops,'ksym=print=old-salcs',.false.,' ')
      prnew=logkey(ops,'ksym=print=new-salcs',.false.,' ')
      prvec=logkey(ops,'ksym=print=vectors',.false.,' ')
      smllmo=fpkey(ops,'ksym=smallmo',1.d-02,' ')
      tstmo=fpkey(ops,'ksym=testmo',1.d-06,' ')
      smleig=fpkey(ops,'ksym=smleig',1.d-05,' ')
      symtre=intkey(ops,'ksym=symmetry',0,' ')
      maxerr=fpkey(ops,'ksym=tol',1.d-06,' ')
      if(symtre.ne.0) then
         nsmall=intkey(ops,'kohn=nsmall',0,' ')
         if(nsmall.eq.0) then
            call iosys('read integer "number of alpha electrons" from'//
     1                 ' rwf',1,nsmall,0,' ')
         endif
         write(iout,909) nsmall+1
      endif
 
c
c     ----- get the dimensions, etc -----
c
      call iosys ('read integer "number of basis functions" from rwf',1,
     1             num,0,' ')
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
         call lnkerr('m604: nsym must be lt 20')
      endif
c
      nnp=(num+1)*num/2
c
c
c     ----- allocate core -----
c
c
      call iosys ('read integer maxsiz from rwf',1,maxcor,0,' ')
c
c
c     ----- get the core needed -----
c     ----- some space is wasted when no basis functions are dropped -----
c
      need=oldnum*oldnum+(7+nsym)*num*num+num+2+nsym*nsym
      need=wpadti(need)+3*nsym+2*num*num+num+oldnum
      need=iadtwp(need)
      call iosys ('write integer maxsiz to rwf',1,need,0,' ')
      call getscm(need,z,ngot,'m604',0)
c
c     allocate core
      salco=1
      xm=salco+oldnum*oldnum
      salcn=xm+num*num
      cvec=salcn+num*num
      proj=cvec+num*num
      temp=proj+num*num*nsym
      hold=temp+num*num
      tt=hold+num*num
      eval=tt+num*num
      tsym=eval+num+2
      s=tsym+nsym*nsym
      mosym=wpadti(s+num*num)
      mosmo=mosym+nsym
      numsym=mosmo+nsym
      mosave=numsym+nsym
      aosave=mosave+num*num
      mov=aosave+num*num
      aov=mov+num
      pkindx=aov+oldnum
      aolst=pkindx
c
c     ----- calculate the salcs in the reduced basis -----
c     -----       if needed -----
c
      ans='no' 
      call iosys ('does character ran402 exist on rwf',0,0,0,ans)
      if (logkey(ops,'drop',.false.,' ')) then
          if (ans.eq.'no') then
              call nsalc(z(salco),z(salcn),a(pkindx),nsym,a(mosmo),
     1                   a(mosym),oldnum,num,prold,prnew)
          else
              call iosys('read integer "number of symmetry orbitals" '//
     1                   'from rwf',nsym,a(mosym),0,' ')
              call iosys('read real "salc transformation matrix" '//
     1                   'from rwf',num*num,z(salcn),0,' ')
          endif
      else 
          call iosys('read integer "number of symmetry orbitals" from'//
     1               ' rwf',nsym,a(mosym),0,' ')
          call iosys('read real "salc transformation matrix" from rwf',
     1                oldnum*oldnum,z(salcn),0,' ')
      endif
      if (logkey(ops,'kohn',.false.,' ')) then
          call iosys('read character "checkpoint filename" from rwf',
     1               -1,0,0,chk)
          call iosys('open chk as unknown',0,0,0,chk)
          call detsym(z(salcn),nsym,a(numsym),a(mosym),a(aolst),num)
          call iosys ('close chk',0,0,0,' ')
      endif
c
c     ----- pick up the transformation vectors from the desired -----  
c     -----                    file             -----
c
      if(logkey(ops,'ksym=kohn',.false.,' ')) then
         kw=1
         write(iout,910)
         call iosys('read character "kohn filename" from rwf',
     1             -1,0,0,kohn)
         call iosys('open kohn as old',0,0,0,kohn)
         call iosys('read character "transformation vector" from kohn',
     1              -1,0,0,xform)
         call iosys('read real '//xform//' from kohn',num*num,z(cvec),
     1               0,' ')
         call iosys('close kohn',0,0,0,' ')
      elseif(logkey(ops,'ksym=chk',.false.,' ')) then
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
         call iosys('open chk as old',0,0,0,chk)
         call iosys ('read real "scf vector" from rwf',num*num,
     1                z(cvec),0,' ')
         call iosys('read real "orbital energies" from rwf',
     1               num,z(eval),0,' ')
      endif
c
      notest=intkey(ops,'ksym=notest',0,' ')
      if(logkey(ops,'ksym=noeval',.false.,' ')) then
         nnn=intkey(ops,'ksym=noeval',num,' ')
         xx=1.d+00
         do 99 i=1,nnn
            z(eval-1+i)=xx
            xx=xx+1.d+00
  99     continue
      endif
c
c----------------------------------------------------------------------c
c           calculate the symmetry orbitals                            c
c----------------------------------------------------------------------c
      call symorb(z(cvec),z(salcn),z(temp),a(mosave),a(aosave),
     1            z(hold),z(tt),z(proj),z(s),z(xm),a(mov),a(aov),
     2            z(eval),z(tsym),a(mosym),a(numsym),nsym,num,smllmo,
     3            tstmo,smleig,symtre,maxerr,nsmall,kw,notest,prvec)
c
      if(.not.logkey(ops,'ksym=nowrite',.false.,' ')) then
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
            call iosys('open kohn as old',0,0,0,kohn)
            call iosys('write real '//xform//' to kohn',num*num,z(cvec),
     1                  0,' ')
            call iosys('close kohn',0,0,0,' ')
         end if
      end if
c
      call iosys ('write integer maxsiz to rwf',1,maxcor,0,' ')
c
c     ----- and exit gracefully -----
c
      call chainx(0)

c
      return
 801  format(/,' vectors have been written to the rwf file')
 802  format(/,' vectors have been written to the kohn file')
 909  format(/,' the virtual orbitals start at orbital ',i4,
     1' and will be symmetry ordered ')
 910  format(/,5x,' reading the kohn file ')
 911  format (/,5x,' reading the rwf file ')
      end
