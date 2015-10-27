*deck @(#)vecin.f	1.1 9/7/91
c***begin prologue     vecin
c***date written       880814   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           vecin, link 2703
c***author             schneider, barry (lanl)
c***source             m2703
c***purpose            input scf vectors
c***                   
c***references         none
c
c***routines called
c***end prologue       vecin
      subroutine vecin (trans,nmo,ncon,prnt)
      implicit integer (a-z)
      logical tstvec, logkey, prnt, posinp
      character *8 cpass
      character *3 itoc 
      character *32 xform
      character *800 card 
      character *80 title
      real *8 trans
      common/io/ inp, iout
      dimension trans(ncon,nmo)
      if ( posinp('$vectors',cpass) ) then
           call cardin(card)
           tstvec=logkey(card,'unit-matrix',.false.,' ')
           if (.not.tstvec) then
                do 10 i=1,nmo
                    call fparr(card,'vector-'//itoc(i),trans(1,i),
     1                         ncon,' ')
   10           continue
           else
                call rzero(trans,ncon*nmo)
                do 20 i=1,nmo
                   trans(i,i)=1.d+00
   20           continue
          endif
          if (prnt) then
              title='transformation matrix'
              call prntrm(title,trans,ncon,nmo,ncon,nmo,iout)
          endif
          return
      endif
      call iosys('read character "transformation vector" from '//
     1           'kohndt',0,0,0,xform)
      call iosys('read real '//xform//' from kohndt',nmo*ncon,
     1            trans,0,' ')      
      if (prnt) then
          title='transformation matrix'
          call prntrm(title,trans,ncon,nmo,ncon,nmo,iout)
      endif
      return
      end
