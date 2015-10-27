*deck @(#)reordr.f	1.1  11/30/90
      subroutine reordr(tempc,c,cdens,iout,nbf,norbs,nnp,ncore,ops,
     $     bflabl,tempe,e)
c
c***begin prologue
c***date written       yymmdd   (yymmdd)
c***revision date      871105   (yymmdd)
c
c     5 november 1987  pws at lanl
c         adding option to not reorder the vector.
c
c    13 october 1987   pws at lanl
c         adding options to print the vector before and after
c         reordering.
c
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)reordr.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
c
c***purpose: to reorder the scf vector from the scf to the drt
c            ordering and also to form the frozen core density
c            matrix (cdens).
c
c paul saxe                 20 august 1984                  lanl
c
      implicit integer (a-z)
c
      real*8 tempc(nbf,nbf),c(nbf,norbs),cdens(nnp)
      real*8 tempe(nbf)
      real*8 e(norbs)
      character*(*) ops
      character*16 bflabl(nbf)
      logical logkey
      integer iout(nbf)
c
      common /io/     inp,ioutpt
c
c
 1000 format(1x,'reordr: drtorb.ne.norbs',2i5)
c
      if (logkey(ops,'transformation=scf',.false.,' ').or.
     #    logkey(ops,'kohn',.false.,' ')) then
         do 450 bf=1,nbf
            iout(bf)=bf
 450     continue
      end if
      ncore=0
      if(logkey(ops,'ci=full',.false.,' ')) then
         do 460 bf=1,nbf
            if(iout(bf).eq.-1) ncore=ncore+1
            iout(bf)=bf
  460    continue
         do 461 bf=1,ncore
            iout(bf)=-1
  461    continue
         do 462 bf=ncore+1,nbf
             iout(bf)=bf-ncore
  462    continue
      end if
      write(ioutpt,*) 'in reorder'
      write(ioutpt,*) (iout(bf),bf=1,nbf)
c
      if (logkey(ops,'print=transformation=vector=input',.false.,
     $     ' ')) then
         write (ioutpt,100)
 100     format (t6,'transformation vector in input order:')
         call iosys('read character "basis function labels" from rwf',
     $               -1,0,0,bflabl)
c
         call wvec(tempc,tempe,nbf,nbf,bflabl,' ')
      end if
c
      call rzero(cdens,nnp)
      ncore=0
      drtorb=0
      do 10 bf=1,nbf
         if (iout(bf).lt.0) then
c
c     ----- frozen core orbitals -----
c
            ncore=ncore+1
            ij=0
            do 2 i=1,nbf
               do 1 j=1,i
                  ij=ij+1
                  cdens(ij)=cdens(ij)+tempc(i,bf)*tempc(j,bf)*2.0d+00
    1          continue
    2       continue
c
         else if (iout(bf).gt.0) then
c
c     ----- orbitals to transform -----
c
            drtorb=drtorb+1
            e(iout(bf))=tempe(bf)
            do 3 i=1,nbf
               c(i,iout(bf))=tempc(i,bf)
    3       continue
         end if
   10 continue
c
      if (logkey(ops,'print=transformation=vector=reordered',.false.,
     $     ' ')) then
         write (ioutpt,110)
 110     format (t6,'transformation vector in drt order:')
         call iosys('read character "basis function labels" from rwf',
     $               -1,0,0,bflabl)
c
         call wvec(c,e,nbf,norbs,bflabl,' ')
      end if
c
      if (drtorb.ne.norbs) then
         write(ioutpt,1000) drtorb,norbs
         call lnkerr(' ')
      end if
c
c
      return
      end
