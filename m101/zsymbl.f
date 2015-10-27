*deck @(#)zsymbl.f	5.2  2/5/95
      subroutine zsymbl(ianz,z,iz,bl,alpha,beta,lbl,lalpha,lbeta,
     $                 nsymbl,nz,bastyp,grdtyp,radius,symbls,namcnt)
c***begin prologue     zsymbl.f
c***date written       850601  yymmdd
c***revision date      2/5/95
c
c     12 april 1988   rlm at lanl
c     updating ccur cursor from ffnext so it can parse cards with
c     basis or point charge specifications.
c
c     14 august 1987   pws at lanl
c     allowing the atomic name to be a string or replacement pattern
c     for handling 'cr(basis=dz-ecp2)'
c***keywords           z-matrix
c***author             binkley, et.al., guassian 82.
c***                   martin, richard (lanl)
c***source             @(#)zsymbl.f	5.2   2/5/95
c***purpose            module to parse the z-matrix.
c***description
c     call zsymbl(ianz,z,iz,bl,alpha,beta,lbl,lalpha,lbeta,
c                nsymbl,nz,bastyp,grdtyp,radius,symbls,namcnt)
c
c     module to parse the z-matrix.
c     the input arrays are explained in m101.
c
c***references
c***routines called    ffnext(chr), isubst(m101), lnkerr(mdutil),
c                      zcentr(m101), zparm(m101), ctoi(chr)
c***end prologue       zsymbl.f
      implicit none
c     --- input variables -----
      character*80 card
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer ianz(*),iz(4,*),lbl(*),lalpha(*),lbeta(*)
      character*(*) bastyp(*),grdtyp(*),namcnt(*)
      character*80 symbls(*)
      real*8 z(*),bl(*),alpha(*),beta(*),radius(*)
c     --- output variables ---
      integer nz,nsymbl
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer ccur,start,end,posb,endb,ocur,scur
      integer isubst,ctoi
      character ffnext*16,found*16,center*16,outstr*80
      character itoc*16,work*80,chrkey*16,fptoc*16,temp*16
      real*8 fpkey
      real*8 chg
      real*8 zero
      logical logkey
c
      parameter (zero=0.0d0)
c
      common/io/inp,iout
c
 1000 format(a)
 1010 format(80a1)
c
c     --- let's go.
      nz=0
      nsymbl=0
c
c     --- top of the loop for parsing the z-matrix.
   40 read(inp,1000,end=9000) card
         if(card.eq.' '.or.index(card,'$').ne.0) return
         call locase(card,card)
         nz=nz+1
c        --- look for the center/basis set specification.
         ccur=0
         found=ffnext(card,ccur,start,end)
         if(found.eq.'string'.or.found.eq.'replacement') then
c           look for a basis set specification.
            posb=index(card,'(')
            endb=index(card,')')
            if(posb.ne.0) then
c              parse basis set/atomic charge/solvent radius.
               work=card(posb+1:endb-1)
               bastyp(nz)=chrkey(work,'basis',' ',' ')
               if (logkey(work,'nobasis',.false.,' '))
     $              bastyp(nz)='nobasis'
               grdtyp(nz)=chrkey(work,'grid',' ',' ')
               if (logkey(work,'nogrid',.false.,' '))
     $              grdtyp(nz)='nogrid'
               center=card(start:start+posb-2)
               ianz(nz)=isubst(center)
               chg=ianz(nz)
               z(nz)=fpkey(work,'z',chg,' ')
               radius(nz)=fpkey(work,'r',zero,' ')
               ccur=endb
            else
               bastyp(nz)=' '
               grdtyp(nz)=' '
               center=card(start:end)
               ianz(nz)=isubst(center)
               z(nz)=ianz(nz)
               radius(nz)=zero
            endif
c
            namcnt(nz)=center
            temp=fptoc(z(nz))
            outstr=' '//center(1:4)//' '//bastyp(nz)(1:16)//' '
     $             //temp(1:8)
            ocur=31
         else
            call lnkerr(' original center specification must be a '
     $                //'string. i found: '//card(start:end))
         endif
c
         if(nz.gt.1) then
c           --- get the name of the center to which this is attached.
            call zcentr(card,ccur,outstr,ocur,iz(1,nz),namcnt,nz)
c           --- get the bond length.
            scur=1
            symbls(nz)=' '
            call zparm(card,ccur,outstr,ocur,bl(nz),lbl(nz),
     $                 symbls(nz),scur,nsymbl)
         endif
c
         if(nz.gt.2) then
c           --- get the third center and the internuclear angle.
            call zcentr(card,ccur,outstr,ocur,iz(2,nz),namcnt,nz)
            call zparm(card,ccur,outstr,ocur,alpha(nz),lalpha(nz),
     $                 symbls(nz),scur,nsymbl)
         endif
c
         if(nz.gt.3) then
c           --- get the fourth center,the dihedral angle,
c               and the last integer.
            call zcentr(card,ccur,outstr,ocur,iz(3,nz),namcnt,nz)
            call zparm(card,ccur,outstr,ocur,beta(nz),lbeta(nz),
     $                 symbls(nz),scur,nsymbl)
c
            found=ffnext(card,ccur,start,end)
            if(found.eq.'eos') then
               iz(4,nz)=0
            else if(found.eq.'integer') then
               iz(4,nz)=ctoi(card(start:end))
            else
               call lnkerr(' unrecognized symbol in the z-matrix: '
     $                     //card(start:end))
            endif
            outstr(ocur+1:)=itoc(iz(4,nz))
            ocur=ocur+4
         endif
c        write(iout,1010) (outstr(k:k),k=1,ocur)
      goto 40
c
c
 9000 call lnkerr(' end of file while parsing the z-matrix. ')
      return
      end
