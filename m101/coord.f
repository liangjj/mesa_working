*deck  @(#)coord.f	5.2  2/5/95
      subroutine coord(natoms,ian,z,c,radius,
     $                 bastyp,grdtyp,atname,el,inau,toang)
c***begin prologue     coord.f
c***date written       850601  yymmdd
c***revision date      2/5/95
c
c    5 october 1987   pws at lanl
c      allowing for return from ffnext of 'replacement' for cases
c      such as 'pt1(basis=dz-ecp1) ...'
c
c***keywords           cartesian, coordinates
c***author             martin, richard (lanl)
c***source             @(#)coord.f	5.2  2/5/95
c***purpose            reads atoms, basis sets, and geometry in cartesian
c                      coordinate form.
c***description
c     call coord(natoms,ian,z,c,bastyp,grdtyp,atname,inau,toang)
c***references
c***routines called    ffnext(chr), ctoi(chr), isubst(m101),
c                      lnkerr(mdutil), ctofp(chr)
c***end prologue       coord.f
      implicit none
c     --- input variables -----
      integer natoms
      real*8 toang
      logical inau
c     --- input arrays (unmodified) ---
      integer ian(*)
      real*8 z(*)
      character*(*) bastyp(*),grdtyp(*),atname(*)
      character*(*) el(0:*)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 c(*),radius(*)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer iatom,cind,icur,bcur,bbcur,ecur,eecur
      integer posb,endb,iicur,i,j
      integer ctoi,isubst
      real*8 ctofp, fpkey, chg
      real*8 zero
      character*16 ffnext, found, symbol
      character*80 card
      character work*80,chrkey*16
      logical logkey
c
      parameter(zero=0.0d0)
c
      common/io/inp,iout
c
 1000 format(a)
 1010 format(' looking for an atomic number or symbol and found: ',80a1)
 1020 format(' looking for a cartesian coordinate and found: ',80a1)
c
c     --- module to read the geometry in cartesian coordinates.
c         loop over atoms until an end of section is found.
      iatom=0
      cind=0
   20 read(inp,1000,end=9000) card
         if(index(card,'$').ne.0) goto 40
         call locase(card,card)
         icur=0
         found=ffnext(card,icur,bcur,ecur)
c        --- skip blank cards.
         if(ecur.eq.0) goto 20
         iatom=iatom+1
         if(found.eq.'integer') then
            ian(iatom)=ctoi(card(bcur:ecur))
            z(iatom)=ian(iatom)
            radius(iatom)=zero
            bastyp(iatom)=' '
            grdtyp(iatom)=' '
            atname(iatom)=el(ian(iatom))
         else if(found.eq.'string'.or.found.eq.'replacement') then
c           --- look for a possible basis,point charge,or radius specification.
            posb=index(card,'(')
            endb=index(card,')')
            if(posb.eq.0) then
               ian(iatom)=isubst(card(bcur:ecur))
               z(iatom)=ian(iatom)
               radius(iatom)=zero
               bastyp(iatom)=' '
               grdtyp(iatom)=' '
               atname(iatom)=card(bcur:ecur)
            else
c              --- center can be either atomic number or atomic symbol.
               iicur=0
               ecur=endb
               icur=endb
               symbol=ffnext(card(bcur:posb-1),iicur,bbcur,eecur)
               if(symbol.eq.'integer') then
                  ian(iatom)=ctoi(card(bcur:posb-1))
                  z(iatom)=ian(iatom)
                  radius(iatom)=zero
                  atname(iatom)=el(ian(iatom))
               else if(symbol.eq.'string') then
                  ian(iatom)=isubst(card(bcur:posb-1))
                  z(iatom)=ian(iatom)
                  radius(iatom)=zero
                  atname(iatom)=card(bcur:posb-1)
               else
                  write(iout,1010) (card(j:j),j=bcur,posb-1)
                  call lnkerr(' ')
               endif
c              --- possibly override the basis and nuclear charge
c                  or radius specifications.
               work=card(posb+1:endb-1)
               bastyp(iatom)=chrkey(work,'basis',' ',' ')
               if(logkey(work,'nobasis',.false.,' '))
     $            bastyp(iatom)='nobasis'
               grdtyp(iatom)=chrkey(work,'grid',' ',' ')
               if(logkey(work,'nogrid',.false.,' '))
     $            grdtyp(iatom)='nogrid'
               chg=ian(iatom)
               z(iatom)=fpkey(work,'z',chg,' ')
               radius(iatom)=fpkey(work,'r',zero,' ')
            endif
         else
            write(iout,1010) (card(j:j),j=bcur,ecur)
            call lnkerr(' ')
         endif
c
         do 30 i=1,3
            found=ffnext(card,icur,bcur,ecur)
            if(found.eq.'floating point') then
               cind=cind+1
               c(cind)=ctofp(card(bcur:ecur))
            else
               write(iout,1020) (card(j:j),j=bcur,ecur)
               call lnkerr(' ')
            endif
   30   continue
      goto 20
c
c     --- have found all coordinate input.
c         set natoms, possibly convert to atomic units, and return.
c
   40 natoms=iatom
      if(.not.inau) then
         do 50 i=1,cind
            c(i)=c(i)/toang
            radius(i)=radius(i)/toang
   50    continue
      endif
c
c
      return
 9000 call lnkerr(' no terminus found for the $geom section.')
      end
