*deck @(#)zmatch.f	5.1  11/6/94
      subroutine zmatch(symbls,name,value,ifvar,nz,nvar,nsymbl,bl,lbl,
     $                 alpha,lalpha,beta,lbeta)
c***begin prologue     zmatch.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c***keywords           z-matrix
c***author             binkley, et.al., gaussian 82
c***                   martin, richard (lanl)
c***source             @(#)zmatch.f	5.1   11/6/94
c***purpose            looks for a specific symbol in the symbolic z-matrix.
c***description
c     call zmatch(symbls,name,value,ifvar,nz,nvar,nsymbl,bl,lbl,
c                 alpha,lalpha,beta,lbeta)
c
c     module to search for a match in the symbolic z-matrix for
c     the symbol 'name'.
c        symbls ... the list of symbols found on each z-matrix card.
c        name   ... the name of the symbol to be matched.
c        value  ... the value of the symbol.
c        ifvar  ... .true./.false.  the symbol is a variable/constant.
c        nsymbl  ... the number of symbols remaining unmatched in the
c                   z-matrix. this is decremented as symbols are found
c                   and replaced with value.
c        bl     ... the bond lengths.
c        lbl    ... bond lengths flags.
c        alpha  ... first bond angle.
c        lalpha ... first angle flag.
c        beta   ... dihedral(second) angle.
c        lbeta  ... second angle flag.
c***references
c***routines called    cskipb(chr), cskipf(chr), streqc(chr), lnkerr(mdutil)
c***end prologue       zmatch.f
c***end prologue       zmatch.f
      implicit none
c     --- input variables -----
      integer nz,nvar,nsymbl
      character*16 name
      logical ifvar
c     --- input arrays (unmodified) ---
      character*80 symbls(*)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer lbl(*),lalpha(*),lbeta(*)
      real*8 bl(*),alpha(*),beta(*)
c     --- output variables ---
      real*8 value
c     --- scratch arrays ---
c     --- local variables ---
      integer i,start,end,cskipf
      character zsyms*16
      logical ok,streqc
c
c     --- loop over the z-matrix searching for a match to 'name'.
c         on each card, check to see if any of the parameters match.
      ok=.false.
      do 100 i=2,nz
c        --- check the bond length.
         start=cskipf(symbls(i),' ')
         end=index(symbls(i)(start:),' ')+start-1
         zsyms=symbls(i)(start:end)
         if(streqc(zsyms,name)) then
            ok=.true.
            nsymbl=nsymbl-1
            if(ifvar) then
               lbl(i)=sign(nvar,lbl(i))
            else
               bl(i)=value
               if(lbl(i).lt.0) bl(i)=-bl(i)
               lbl(i)=0
            endif
         endif
c
c        --- check the first angle.
         if(i.gt.2) then
            start=cskipf(symbls(i)(end:),' ')+end-1
            end=index(symbls(i)(start:),' ')+start-1
            zsyms=symbls(i)(start:end)
            if(streqc(zsyms,name)) then
               ok=.true.
               nsymbl=nsymbl-1
               if(ifvar) then
                  lalpha(i)=sign(nvar,lalpha(i))
               else
                  alpha(i)=value
                  if(lalpha(i).lt.0) alpha(i)=-alpha(i)
                  lalpha(i)=0
               endif
            endif
         endif
c
c        --- check the second(dihedral) angle.
         if(i.gt.3) then
            start=cskipf(symbls(i)(end:),' ')+end-1
            end=index(symbls(i)(start:),' ')+start-1
            zsyms=symbls(i)(start:end)
            if(streqc(zsyms,name)) then
               ok=.true.
               nsymbl=nsymbl-1
               if(ifvar) then
                  lbeta(i)=sign(nvar,lbeta(i))
               else
                  beta(i)=value
                  if(lbeta(i).lt.0) beta(i)=-beta(i)
                  lbeta(i)=0
               endif
            endif
         endif
c
  100 continue
c
      if(.not.ok)
     $   call lnkerr('symbol not found in z-matrix: '
     $               //name)
c
c
      return
      end
