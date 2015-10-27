*deck @(#)merge.f	5.1  11/6/94
      subroutine merge(c,t1,t2,eig,ev1,ev2,nbf,ops)
      implicit integer(a-z)
      real*8 c(nbf,nbf),t1(nbf,*),t2(nbf,*),eig(*),ev1(*),ev2(*)
      integer ind(3,50)
      character*(*) ops
      logical logkey
      character*128 namchk,namsav
      character *4 itoc
      common /io/ inp,iout
      data ind/150*-1/
      save ind
c
c     prepare the checkpoint file for action.
      call iosys('read character "checkpoint filename" from rwf',
     $        0,0,0,namchk)
      call iosys('open chk as old',0,0,0,namchk)
c
c     retrieve the old basis parameters from the chk.
c
      call iosys('read integer "number of basis functions" from chk',
     $        1,nbasj,0,' ')
c
      if(nbasj.ne.nbf) then
         write(iout,*)' nbf nbasj ',nbf,nbasj
         call lnkerr(' merge: basis mismatch')
      end if
c
      n2=nbf*nbf
c
      isave=0
c
      if(logkey(ops,'set1=chk',.false.,' ')) then
         if(logkey(ops,'set1=chk=mcscf',.false.,' ')) then
            call iosys('read real "mcscf vector" from chk',n2,t1,0,' ')
            call iosys('read real "mcscf orbital energies" from chk',
     $                  nbf,ev1,0,' ')
         else if(logkey(ops,'set1=chk=no',.false.,' ')) then
            root=intkey(ops,'set1=chk=no',1,' ')
            call iosys('read real "no vector '//itoc(root)//'"'
     $                  //' from chk',n2,t1,0,' ')
            call iosys('read real "no occ '//itoc(root)//'" from chk',
     $                  nbf,ev1,0,' ')
         else
            call iosys('read real "scf vector" from chk',n2,t1,0,' ')
            call iosys('read real "orbital energies" from chk',
     $                 nbf,ev1,0,' ')
         endif
      else if(logkey(ops,'set1=sav',.false.,' ')) then
         isave=1
c..unicos
         namsav='sav'
c..unicos
         write(iout,*)'    opening unit sav as ',namsav
         call iosys('open sav as old',0,0,0,namsav)
         call iosys('read integer "number of basis functions"'
     $              //' from sav',1,nbasj,0,' ')
c
         if(nbasj.ne.nbf) then
            write(iout,*)' set1-sav nbf nbasj ',nbf,nbasj
            call lnkerr(' merge: basis mismatch')
         end if
c
         if(logkey(ops,'set1=sav=mcscf',.false.,' ')) then
            call iosys('read real "mcscf vector" from sav',n2,t1,0,' ')
            call iosys('read real "mcscf orbital energies" from sav',
     $                nbf,ev1,0,' ')
         else if(logkey(ops,'set1=sav=no',.false.,' ')) then
            root=intkey(ops,'set1=sav=no',1,' ')
            call iosys('read real "no vector '//itoc(root)//'"'
     $                  //' from sav',n2,t1,0,' ')
            call iosys('read real "no occ '//itoc(root)//'" from sav',
     $                  nbf,ev1,0,' ')
         else
            call iosys('read real "scf vector" from sav',n2,t1,0,' ')
            call iosys('read real "orbital energies" from sav',
     $                 nbf,ev1,0,' ')
         endif
      else
         call lnkerr(' merge: set1 not found in option string')
      end if
c
c
      if(logkey(ops,'set2=chk',.false.,' ')) then
         if(logkey(ops,'set2=chk=mcscf',.false.,' ')) then
            call iosys('read real "mcscf vector" from chk',n2,t2,0,' ')
            call iosys('read real "mcscf orbital energies" from chk',
     $                nbf,ev2,0,' ')
         else if(logkey(ops,'set2=chk=no',.false.,' ')) then
            root=intkey(ops,'set2=chk=no',1,' ')
            call iosys('read real "no vector '//itoc(root)//'"'
     $                  //' from chk',n2,t2,0,' ')
            call iosys('read real "no occ '//itoc(root)//'" from chk',
     $                  nbf,ev2,0,' ')
         else
            call iosys('read real "scf vector" from chk',n2,t2,0,' ')
            call iosys('read real "orbital energies" from chk',
     $                nbf,ev2,0,' ')
         endif
      else if(logkey(ops,'set2=sav',.false.,' ')) then
c
         if(isave.eq.0) then
c..unicos
            namsav='sav'
c..unicos
            write(iout,*)'    opening unit sav as ',namsav
            call iosys('open sav as old',0,0,0,namsav)
            call iosys('read integer "number of basis functions"'
     $                 //' from sav',1,nbasj,0,' ')
c
            if(nbasj.ne.nbf) then
               write(iout,*)' set1-sav nbf nbasj ',nbf,nbasj
               call lnkerr(' merge: basis mismatch')
            end if
c
         end if
         if(logkey(ops,'set2=sav=mcscf',.false.,' ')) then
            call iosys('read real "mcscf vector" from sav',n2,t2,0,' ')
            call iosys('read real "mcscf orbital energies" from sav',
     $                  nbf,ev2,0,' ')
         else if(logkey(ops,'set2=sav=no',.false.,' ')) then
            root=intkey(ops,'set2=sav=no',1,' ')
            call iosys('read real "no vector '//itoc(root)//'"'
     $                  //' from sav',n2,t2,0,' ')
            call iosys('read real "no occ '//itoc(root)//'" from sav',
     $                  nbf,ev2,0,' ')
         else
            call iosys('read real "scf vector" from sav',n2,t2,0,' ')
            call iosys('read real "orbital energies" from sav',
     $                nbf,ev2,0,' ')
         endif
      else
         call lnkerr(' merge: set2 not found in option string')
      end if
c
      if(logkey(ops,'guess=ind',.false.,' ')) then
         call intarr(ops,'guess=ind',ind,150,' ')
      else
         call lnkerr(' merge: ind array not found ')
      end if
c
      call rzero(c,n2)
c
      ix=1
      iv=1
      ic=0
c
      write(iout,101)
 101  format(/,'  merge input parameters',/,
     #         '  vector set   starting vector   ending vector')
 102  format(2x,i6,7x,i6,14x,i6)
      do 1 i=1,150
         if(ind(1,i).eq.-1) go to 2
         write(iout,102)(ind(k,i),k=1,3)
         ix=ind(2,i)
         nvec=ind(3,i)-ind(2,i)+1
         lcop=nvec*nbf
         ic=ic+nvec
         if(ic.gt.nbf) then
            write(iout,*)' merge: too many vectors'
            call lnkerr(' m401: input errors in merge')
         end if
         if(ind(1,i).eq.1) then
            call scopy(nvec,ev1(ix),1,eig(iv),1)
            call scopy(lcop,t1(1,ix),1,c(1,iv),1)
         else if (ind(1,i).eq.2) then
            call scopy(nvec,ev2(ix),1,eig(iv),1)
            call scopy(lcop,t2(1,ix),1,c(1,iv),1)
         else
           call lnkerr(' error in set no. in ind')
         end if
         iv=iv+nvec
   1  continue
c
  2   continue
c
      if(ic.ne.nbf) then
         write(iout,*)' nvec ne nbf ',ic,nbf
         call lnkerr(' input error in merge')
      end if
c
c

      return
      end
