*deck @(#)zparm.f	5.2  4/17/95
      subroutine zparm(card,icur,outstr,ocur,x,lx,symbls,scur,nsymbl)
c***begin prologue     zparm.f
c***date written       850601  yymmdd
c***revision date      4/17/95
c***keywords           z-matrix
c***author             binkley, et.al., gaussian 82.
c***                   martin, richard (lanl)
c***source             @(#)zparm.f	5.2   4/17/95
c***purpose            reads a z-matrix parameter (variable or constant).
c***description
c     call zparm(card,icur,outstr,ocur,x,lx,symbls,scur,nsymbl)
c
c     module to read a z-matrix parameter (variable or constant).
c
c     if a constant is found, this is returned in x, and lx is set to 2.
c     in this case a name of '0' is stored in symbls.
c
c     if a variable (name) is found, then this name is appended to symbls,
c     nsymbl is incremented, x is set to 0.0, and lx is set to 3
c     (or -3 if -name is found).
c
c     card is the current z-matrix card being parsed.
c     it is parsed beginning at card(icur+1:).  outstr is an output
c     string which will eventually be sent to the printer.
c***references
c***routines called    ffnext(chr), ctofp(chr), lnkerr(mdutil)
c***end prologue       zparm.f
      implicit none
c     --- input variables -----
      integer icur,scur,ocur
      integer nsymbl
c     --- input arrays (unmodified) ---
      character*(*) card
c     --- input arrays (scratch) ---
c     --- output arrays ---
      character*(*) outstr
      character*80 symbls
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer lx,start,end
      integer inp,iout
      real*8 x,zero,ctofp
      character*16 found,ffnext
      character*24 parm,fptoc
c
      data zero/0.0d+00/
      save zero
c
      common/io/inp,iout
c
c     --- see what's next on this card. 
      found=ffnext(card,icur,start,end)
      parm=card(start:end)
      if(found.eq.'floating point') then
         x=ctofp(card(start:end))
         lx=2
         symbls(scur+1:)='0'
         outstr(ocur+1:)=fptoc(x)
      else if(found.eq.'string') then
         x=zero
         nsymbl=nsymbl+1
         if(parm(1:1).eq.'-') then
            lx=-3
            symbls(scur+1:)=parm(2:)
         else
            lx=3
            symbls(scur+1:)=parm
         endif
         outstr(ocur+1:)=parm
      else
         call lnkerr(' z-matrix parameter must be either a string '
     $               //'or a floating point constant:'//parm)
      endif
c
      scur=scur+16
      ocur=ocur+10
      outstr(ocur-1:)=' '
c
c
      return
      end
