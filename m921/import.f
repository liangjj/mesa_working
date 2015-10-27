*deck @(#)import.f	5.1  11/6/94
      subroutine import(c,nwks,coef,walk,nimprt,cutoff,nlarge,
     #                  lgcoef,lgwk)
c
c***begin prologue  import
c***date written   851003   (yymmdd)
c***revision date  880422   (yymmdd)
c
c   22 april 1988  bhl at   llnl
c     read mcscf orbital energies for mrci case
c
c***keywords most important configurations, breakdown
c
c***author  saxe, paul,    (lanl)
c***purpose  to find the most important configurations in a ci vector
c
c***description  import finds the 'nimprt' most important
c       configurations in a ci vector, and stores the coefficients and
c       configuration numbers ordered in coef and walk, respectively.
c
c       on input:
c
c            c          real (nwks)
c                       the ci vector.
c
c            nwks       integer
c                       the number of configurations.
c
c            nimprt     integer
c                       the number of important configurations to find.
c
c
c        on output:
c
c            coef       real (nimprt)
c                       the ordered largest coefficients.
c
c            walk       integer (nimprt)
c                       the configuration number associated with the
c                       coefficients in 'coef'.
c
c
c***references
c
c***routines called  rzero (util), izero (util)
c***end prologue import
c
c
      implicit integer (a-z)
c
      real*8 cutoff
      real*8 c(nwks),coef(nimprt),lgcoef(nlarge)
      integer walk(nimprt),lgwk(nlarge)
c
c     ----- zero the arrays for string the most important coefficients
c           and walk numbers
c
      call rzero(lgcoef,nlarge)
      call izero(lgwk,nlarge)
      call rzero(coef,nimprt)
      call izero(walk,nimprt)
c
c     ----- extract the coefficients larger than cutoff from the vector
c
      nlg=0
      do 20 wk=nwks,1,-1
         if (abs(c(wk)).gt.cutoff) then
            nlg=nlg+1
            if (nlg.gt.nlarge) call lnkerr('nlarge too small')
            lgwk(nlg)=wk
            lgcoef(nlg)=c(wk)
         end if
   20 continue
c
c     ----- pass through the vector, seeing if the current coefficient
c           is larger than the smallest of our current most important
c           list. we start at the end of the list simply for efficiency
c           as many times the most important configurations are near
c           end of the list.
c
      nimprt=min(nimprt,nlg)
      do 4 wk=1,nlg
         if (abs(lgcoef(wk)).lt.abs(coef(nimprt))) go to 4
c
c        ----- new element for the list, so go up the list, moving
c              elements down one bin, till we can insert the new one
c
         do 2 i=nimprt,2,-1
            if (abs(lgcoef(wk)).lt.abs(coef(i-1))) go to 3
            coef(i)=coef(i-1)
            walk(i)=walk(i-1)
    2    continue
c
         i=1
    3    continue
         coef(i)=lgcoef(wk)
         walk(i)=lgwk(wk)
    4 continue
c
c
      return
      end
