*deck @(#)setcom.f	5.1  11/6/94
      subroutine setcom(rttyp,route,lenrt,ops,lenops)
c***begin prologue     setcom.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c   october 16, 1992  rlm at lanl
c     adding new list of physical constants
c        (the command oldconstants in the route will load the old ones). 
c***keywords           rwf, common blocks
c***author             martin, richard (lanl)
c***source             @(#)setcom.f	5.1   11/6/94
c***purpose            initializes certain system common blocks.
c***description
c     call setcom(rttyp,route,lenrt,ops,lenops)
c       rttyp   the type of route.
c       route   the route in it's nonstandard form.
c       lenrt   the length of the route string.
c       ops     the options string.
c       lenops  the length of the options string.
c
c     this module initializes some mesa common blocks.
c     it writes archive information, physical constants,
c     the route in nonstandard(internal) form, and route status
c     information to the rwf.
c***references
c***routines called    iosys(io), phyfil(m1),
c                      izero(math), rzero(math)
c***end prologue       setcom.f
      implicit none
c     --- input variables -----
      integer lenrt,lenops
      character*(*) rttyp,route,ops
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer mxlnk,mxlops
      parameter (mxlnk=200,mxlops=4096)
c
      integer ncards,curcrd,curlnk,ov,links,nlnks
      integer nenter(mxlnk),lnknos(mxlnk),coruse(mxlnk)
      integer beg,end,lenst
      real*8 times(3,mxlnk)
      real*8 phycon(30)
      logical logkey
      character*4096 info
      character jump*8, temp*80
c
      common/arcinf/info
      common/rtstat/ncards,curcrd,curlnk,ov,links(mxlnk),nlnks
c
c
      if(rttyp.eq.'restart') then
      else
c
c        --- write the archive information to the rwf.  currently unused.
         call iosys('write character "archive information" to rwf',
     $               0,0,0,info)
c
c        --- load and write the physical constants to the rwf.
         if(logkey(ops,'oldconstants',.false.,' ')) then
            call phyold(30,phycon)
         else
            call phyfil(30,phycon)
         endif
         call iosys('write real angstrom/bohr to rwf',1,phycon(1),0,
     $              ' ')
         call iosys('write real kg/amu to rwf',1,phycon(2),0,' ')
         call iosys('write real esu/e- to rwf',1,phycon(3),0,' ')
         call iosys('write real planck to rwf',1,phycon(4),0,' ')
         call iosys('write real avogadro to rwf',1,phycon(5),0,' ')
         call iosys('write real j/cal to rwf',1,phycon(6),0,' ')
         call iosys('write real m/bohr to rwf',1,phycon(7),0,' ')
         call iosys('write real j/hartree to rwf',1,phycon(8),0,' ')
         call iosys('write real lightspeed to rwf',1,phycon(9),0,' ')
         call iosys('write real boltzmann to rwf',1,phycon(10),0,' ')
         call iosys('write real fine-structure to rwf',
     $               1,phycon(11),0,' ')
         call iosys('write real pi to rwf',1,phycon(12),0,' ')
         call iosys('write real "electron mass" to rwf',
     $               1,phycon(13),0,' ')
         call iosys('write real "ideal molar volume" to rwf',
     $               1,phycon(14),0,' ')
c
c        --- prepare the route information for chainx.
         call rmvnb(route,route)
         lenrt=index(route,' ')
c        --- blow this up to card size.
c            change chainx someday so it can read compressed route.
         ncards=0
         beg=1
         call iosys('create character route on rwf',-1,0,0,' ')
   10    end=index(route(beg:),';')
         if(end.ne.0) then
            ncards=ncards+1
            temp=route(beg:beg+end-1)
            call iosys('write character route to rwf without'
     $                 //' rewinding',0,0,0,temp)
            beg=beg+end
            goto 10
         endif
         call iosys('endfile route on rwf',0,0,0,' ')
c
c        --- pack the options string and write it to the rwf.
         call pakstr(ops,lenops)
         if(lenops.gt.mxlops) call lnkerr('option string too long.')
         call iosys('write character options to rwf',mxlops,0,0,ops)
c
c        --- prepare data needed for the chaining process.
         curcrd=0
         curlnk=1
         ov=0
         call izero(links,mxlnk)
         nlnks=1
         links(1)=1
         jump=' '
         lenst=5+mxlnk
         call iosys('write integer "route status" to rwf',
     $               lenst,ncards,0,' ')
         call iosys('write character "route jump" to rwf',8,0,0,jump)
c
c        --- initialize the job statistics blocks and write to rwf.
         call izero(nenter,mxlnk)
         call izero(lnknos,mxlnk)
         call izero(coruse,mxlnk)
         call rzero(times,3*mxlnk)
         call iosys('write integer "links:nenter" to rwf',
     $               mxlnk,nenter,0,' ')
         call iosys('write integer "links:lnknos" to rwf',
     $               mxlnk,lnknos,0,' ')
         call iosys('write integer "links:coruse" to rwf',
     $               mxlnk,coruse,0,' ')
         call iosys('write real "links:times" to rwf',
     $               3*mxlnk,times,0,' ')
c
      endif
c
c
      return
      end
