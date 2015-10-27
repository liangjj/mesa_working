*deck @(#)zvar.f	5.1  11/6/94
      subroutine zvar(bl,alpha,beta,lbl,lalpha,lbeta,values,
     $                fpvec,intvec,nsymbl,nz,nvar,anames,symbls)
c***begin prologue     zvar.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c***keywords           z-matrix
c***author             binkley, et.al., gaussian 82.
c***                   martin, richard (lanl)
c***source             @(#)zvar.f	5.1   11/6/94
c***purpose            reads variables section of the z-matrix.
c***description
c     call zvar(bl,alpha,beta,lbl,lalpha,lbeta,values,
c               fpvec,intvec,nsymbl,nz,nvar,anames,symbls)
c
c     module to read variables section of the z-matrix.
c     the input arrays are explained in m101.
c***references
c***routines called    ffnext(chr), ctofp(chr), lnkerr(mdutil), streqc(chr),
c                      locase(chr), azmatch(m101), crjust(chr)
c***end prologue       zvar.f
      implicit none
c     --- input variables -----
      integer nsymbl,nz
c     --- input arrays (unmodified) ---
      integer lbl(*),lalpha(*),lbeta(*),intvec(*)
      character*80 symbls(*)
      character*(*) anames(*)
      real*8 bl(*),alpha(*),beta(*),values(*),fpvec(*)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
      integer nvar
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer ccur,start,end,poseq,i
      logical streqc,ifvar
      character ffnext*16,found*16,sname*16,svalue*16,card*80
      character fcname*16,fcval*16
      real*8 ctofp
c
      common/io/inp,iout
c
 1000 format(a)
 1010 format(a80)
 1070 format(' undefined z-matrix symbol, bond length: ',i4)
 1072 format(' undefined z-matrix symbol, angle alpha: ',i4)
 1074 format(' undefined z-matrix symbol, angle beta: ',i4)
c
c     --- top of the loop for parsing the variables.
      nvar=0
   40 read(inp,1000) card
         if(card.eq.' '.or.index(card,'$').ne.0) goto 100
         call locase(card,card)
         nvar=nvar+1
         ccur=0
c        --- look for a variable replacement.
         found=ffnext(card,ccur,start,end)
         if(found.eq.'replacement') then
            poseq=index(card(start:end),'=')
            anames(nvar)=card(start:start+poseq-2)
            values(nvar)=ctofp(card(start+poseq:end))
            sname=anames(nvar)
            svalue=card(start+poseq:end)
         else
            call lnkerr(' incorrect variable replacement syntax.'
     $             //' i found: '//card(start:end))
         endif
c
c        --- look for force constant flag.
         fcname=' '
         fcval=' '
         found=ffnext(card,ccur,start,end)
         if(found.eq.'replacement') then
            poseq=index(card(start:end),'=')
            fcname=card(start:start+poseq-2)
            fcval=card(start+poseq:end)
            if(.not.streqc(fcname(1:3),'d2e'))
     $        call lnkerr(' unrecognized force constant flag: '//fcname)
            if(streqc(fcval(1:5),'numer')) then
               intvec(nvar)=3
            else
               intvec(nvar)=1
               fpvec(nvar)=ctofp(fcval)
            endif
         endif
c
c        --- look for matches in the z-matrix.
c            this card was just processed as if it were truly a variable.
c            maybe it's not.  check this out and adjust appropriately.
         call locase(card,card)
         if(index(card,'constant').ne.0) then
            ifvar=.false.
            call zmatch(symbls,anames(nvar),values(nvar),ifvar,nz,nvar,
     $                  nsymbl,bl,lbl,alpha,lalpha,beta,lbeta)
            nvar=nvar-1
         else
            ifvar=.true.
            call zmatch(symbls,anames(nvar),values(nvar),ifvar,nz,nvar,
     $                  nsymbl,bl,lbl,alpha,lalpha,beta,lbeta)
         endif
c
      goto 40
c
  100 continue
c
c     --- make sure that all symbols in the z-matrix have been defined.
      if(nsymbl.eq.0) then
         return
      else if(nsymbl.lt.0) then
         call lnkerr(' more variable replacements than variables.')
      else
         do 120 i=1,nz
            if(iabs(lbl(i)).eq.3000) write(iout,1070)
            if(iabs(lalpha(i)).eq.3000) write(iout,1072)
            if(iabs(lbeta(i)).eq.3000) write(iout,1074)
  120    continue
         call lnkerr(' ')
      endif
c
c
      return
      end
