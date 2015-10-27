*deck @(#)onepdm.f	5.1  11/6/94
      subroutine onepdm(arc,wt,nlwks,orbsym,b,srow,
     #                c,s,d1,
     #                nrows,norbs,nlevs,orbfrm,nsym,nwks,nnp,
     #                irowsv,segsv,pagesv,iwtsv,traksv,jrowsv,
     #                jwtsv,hrowsv,hwt,harcsv,
     #                acoef,bcoef,trans)
c
c***begin prologue  onepdm
c***date written   851003   (yymmdd)
c***revision date  870225   (yymmdd)
c
c   25 february 1987     pws at lanl
c     adding transition density matrices.
c
c***keywords one-particle density matrix, density matrices
c
c***author  saxe, paul,    (lanl)
c***purpose  to compute a one-particle density matrix from a ci vector.
c
c***description  onepdm computes the one-particle density matrix given
c       a ci vector and some of the drt arrays. the density matrix in
c       this context is the entity which, when dotted with the full
c       list of integrals, gives the energy.
c
c       on input:
c
c          arc      integer (4,nrows)
c                   the arc array of the drt.
c
c          wt       integer (4,nrows)
c                   the weights of the arcs in the drt.
c
c          nlwks    integer (nrows)
c                   the number of lower walks from the rows in the drt.
c
c          orbsym   integer (norbs)
c                   the symmetry of the ci orbitals in guga ordering.
c                   symmetries start at 0.
c
c          b        integer (nrows)
c                   the 'b' values of the rows in the drt.
c
c          srow     integer
c                   the row to consider as the top of the graph. this
c                   can be used to get smaleer ci's out of larger ones.
c
c          c        real (nwks)
c                   the ci vector.
c
c          s        real (nwks)
c                   currently not used, but will be the second vector
c                   for transition density matrices.
c
c          nrows    integer
c                   the number of rows in the drt.
c
c          norbs    integer
c                   the number of orbitals in the ci.
c
c          nlevs    integer
c                   the number of levels in the drt. (norbs+1)
c
c          orbfrm   integer
c                   the fermi level in terms of orbitals, ie. the
c                   number of virtual orbitals.
c
c          nsym     integer
c                   the number of symmetries in the point group used
c                   for the ci.
c
c          nwks     integer
c                   the number of configurations.
c
c          nnp      integer
c                   the size of a triangular matrix 'norbs' on a side.
c                   (norbs+1)*norbs/2
c
c
c       on output:
c
c          d1       real (nnp)
c                   the one-particle density matrix.
c
c
c       scratch:
c
c          irowsv   integer (nlevs)
c          segsv    integer (nlevs)
c          pagev    integer (nlevs)
c          iwtsv    integer (nlevs)
c          traksv   integer (nlevs)
c          jrowsv   integer (nlevs)
c          jwtsv    integer (nlevs)
c          hrowsv   integer (nlevs)
c          hwt      integer (nlevs)
c          harcsv   integer (nlevs)
c          acoef    real (nlevs)
c          bcoef    real (nlevs)
c
c
c***references
c
c***routines called  sdot (calmath), rzero (math), lnkerr (mdutil)
c***end prologue onepdm
c
      implicit integer (a-z)
c
      real*8 c(nwks),acoef(nlevs),bcoef(nlevs),s(nwks)
      real*8 d1(norbs,norbs)
      integer arc(4,nrows),wt(4,nrows),nlwks(nrows)
      integer orbsym(norbs),b(nrows)
      integer irowsv(nlevs),segsv(nlevs),pagesv(nlevs),iwtsv(nlevs)
      integer traksv(nlevs),jrowsv(nlevs),jwtsv(nlevs)
      integer hrowsv(nlevs),hwt(nlevs),harcsv(nlevs)
      logical trans
c
c     ----- local arrays -----
c
      real*8 cfs(420),coeffs(20,21)
      integer segmax(22),trak(35),jcond(35),kcond(35),nxtpag(35)
      integer iarc(35),jarc(35),segmin(22)
c
c     ----- local variables -----
c
      real*8 crite,root2,rootn2,toor2,toorn2
      real*8 d,loopij,t,loopji
c
c     ----- external functions -----
c
      real*8 sdot
c
      common /io/ inp,iout
c
      equivalence (segmin(2),segmax(1)),(coeffs,cfs)
c
      data segmin(1) /0/
      data segmax/ 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,14,
     a            21,28,35,35,35,35,35/
      data   trak
     #  /  3,  3,  2,  7,  2,  7, 10,  0,  0,  0,
     #     0,  0,  8,  0,  0,  0,  0,  0,  0,  8,
     #     0,  0,  0,  0,  0,  0,  6,  0,  0,  0,
     #     0,  0,  0,  6,  0/
      data  jcond
     #  /  1,  1,  1,  1,  1,  1,  1,  0,  0,  0,
     #     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     #     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     #     0,  0,  0,  0,  0/
      data  kcond
     #  /  1,  1,  1,  1,  1,  1,  1,  0,  0,  0,
     #     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     #     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     #     0,  0,  0,  0,  0/
      data nxtpag
     #  / 18, 17,  0, 15,  0, 16,  0, 15, 15,  0,
     #    16, 15,  0, 15, 16,  0, 16, 15, 16,  0,
     #    16, 17, 17,  0, 18, 17,  0, 17, 18,  0,
     #    18, 17, 18,  0, 18/
      data   jarc
     #  /  2,  3,  2,  4,  3,  4,  4,  1,  2,  1,
     #     2,  3,  2,  4,  1,  1,  2,  3,  3,  3,
     #     4,  1,  2,  1,  2,  3,  2,  4,  1,  1,
     #     2,  3,  3,  3,  4/
      data   iarc
     #  /  1,  1,  2,  2,  3,  3,  4,  1,  2,  3,
     #     3,  3,  4,  4,  1,  2,  2,  2,  3,  4,
     #     4,  1,  2,  3,  3,  3,  4,  4,  1,  2,
     #     2,  2,  3,  4,  4/
c
      save segmin,segmax,trak,jcond,kcond,nxtpag,jarc,iarc
c
c     ----- build the coefficient table -----
c
      call rzero(coeffs,20*21)
c
      do 700 i=3,20
         t = float(i-2)
         coeffs(i,1) = sqrt(t/(t+1.0d+00))
         coeffs(i,2) = -coeffs(i,1)
         coeffs(i,3) = coeffs(i,1)/sqrt(2.0d+00)
         coeffs(i,4) = -coeffs(i,3)
         coeffs(i,5) = sqrt((t+1.0d+00)/t)
         coeffs(i,6) = -coeffs(i,5)
         coeffs(i,7) = coeffs(i,5)/sqrt(2.0d+00)
         coeffs(i,8) = -coeffs(i,7)
         coeffs(i,9) = sqrt((t+2.0d+00)/(t*2.0d+00))
         coeffs(i,10) = -coeffs(i,9)
         coeffs(i,11) = sqrt(t/(2.0d+00*(t+2.0d+00)))
         coeffs(i,12) = -coeffs(i,11)
         coeffs(i,13) = sqrt(2.0d+00/(t*(t+1.0d+00)))
         coeffs(i,14) = sqrt(t*(t+2.0d+00))/(t+1.0d+00)
         coeffs(i,15) = -sqrt(t*(t+2.0d+00))/(t+1.0d+00)
         coeffs(i,16) = sqrt((t-1.0d+00)*(t+2.0d+00)/(t*(t+1.0d+00)))
         coeffs(i,17) = -coeffs(i,16)
         coeffs(i,18)=-sqrt(2.0d+00/(t*(t+2.0d+00)))
         coeffs(i,19) = 1.0d+00/t
         coeffs(i,20) = -1.0d+00/t
         coeffs(i,21) = -sqrt(2.0d+00)/t
  700 continue
c
c     ----- set up some constants -----
c
      crite=0.00001d+00
      root2=sqrt(2.0d+00)
      rootn2=-root2
      toor2=1.0d+00/root2
      toorn2=-toor2
c
c     ----- zero storage for density matrix -----
c
      call rzero(d1,norbs**2)
c
c     ----- loop over opening levels of loops, starting at top of graph
c           irow is the distinct row at this level which the single j
c           walk passes through.
c
      do 1000 level=nlevs,2,-1
c
c        ----- search for head segments from top of graph to level
c              of opening of loops
c
         hrow=srow
         hlev=nlevs
         hwt(hlev)=0
         harc=0
         if (hlev.eq.level) go to 403
c
  401    continue
         harc=harc+1
         if (harc.gt.4) then
            hlev=hlev+1
            if (hlev.gt.nlevs) go to 1000
            hrow=hrowsv(hlev-1)
            harc=harcsv(hlev-1)
            go to 401
         end if
c
         hnxt=arc(harc,hrow)
         if (hnxt.eq.0) go to 401
c
         if (hlev-1.eq.level) go to 402
c
         hlev=hlev-1
         hrowsv(hlev)=hrow
         harcsv(hlev)=harc
         hwt(hlev)=hwt(hlev+1)+wt(harc,hrow)
         harc=0
         hrow=hnxt
         go to 401
c
  402    continue
         irow=hnxt
         iwtsv(level)=hwt(hlev)+wt(harc,hrow)
         go to 404
c
  403    continue
         irow=hrow
         iwtsv(level)=0
c
  404    continue
         jrow=irow
         jwtsv(level)=iwtsv(level)
c
         page=1
         lev=level
         levm1=lev-1
         i=levm1
         ia=i*(i-1)/2
         is=orbsym(i)
         seg=segmin(page)
         segmx=segmax(page)
c
c        ----- initialize the stacks -----
c
         traksv(lev)=0
         acoef(lev)=1.0d+00
c
c        ----- test next segment in this page to see if has correct j arc
c
  201    continue
         seg=seg+1
         if (seg.gt.segmx) then
c
c           ----- have tried all segments at this level, so back up one
c
            lev=lev+1
            if (lev.gt.level) then
               go to 401
            end if
c
            levm1=lev-1
            seg=segsv(levm1)
            page=pagesv(levm1)
            segmx=segmax(page)
            irow=irowsv(levm1)
            jrow=jrowsv(levm1)
            go to 201
         end if
c
         jnxt=arc(jarc(seg),jrow)
         if (jnxt.eq.0) go to 201
c
c        ----- check if i arc leads to a valid row -----
c
         inxt=arc(iarc(seg),irow)
         if (inxt.eq.0) go to 201
c
c        ----- check for a change in tracks -----
c
         if (trak(seg).ne.0) then
            traksv(levm1)=trak(seg)
         else
            traksv(levm1)=traksv(lev)
         end if
c
      go to
     #   ( 1,  1,  1,  3,  1,  4, 50,  1,  6,  1,
     #     9,  2,  5,  2,  1,  1,  2, 36, 11, 10,
     #     2,  1,  6,  1,  9,  2,  5,  2,  1,  1,
     #     2, 36, 11, 10,  2), seg
c
  1      acoef(levm1)=acoef(lev)
         goto 120
  2      acoef(levm1)=-acoef(lev)
         goto 120
   3     junk = b(jrow) + 2
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
   4     junk = b(jrow) + 83
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
   5     junk = b(jrow) + 82
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
   6     junk = b(jrow) + 261
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
   9     junk = b(jrow) + 362
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  10     junk = b(jrow) + 3
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  11     junk = b(jrow) + 263
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  36     junk = b(jrow) + 384
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  50     acoef(levm1) = acoef(lev) + acoef(lev)
         d=0.5d+00
         go to 120
  120    continue
c
c        ----- update stacks -----
c
         newpag=nxtpag(seg)
         if (newpag.eq.0) then
            j=levm1
c
            if (orbsym(j).ne.is) go to 201
c
c           ----- have finished a loop, so check if the rows are
c                 the same. (because of 't' values from interacting
c                 spaces, the rows may not be the same, though their
c                 a,b and s values must be.
c
            if (inxt.ne.jnxt) then
c
c              ----- push row, etc on stack and search down for closing -----
c
               tlev=lev-1
               irowsv(tlev)=irow
               jrowsv(tlev)=jrow
               iwtsv(tlev)=iwtsv(tlev+1)+wt(iarc(seg),irow)
               jwtsv(tlev)=jwtsv(tlev+1)+wt(jarc(seg),jrow)
               tarc=0
               itrow=inxt
               jtrow=jnxt
c
  400          continue
               tarc=tarc+1
               if (tarc.gt.4) then
                  tlev=tlev+1
                  if (tlev.ge.lev) go to 201
                  itrow=irowsv(tlev-1)
                  jtrow=jrowsv(tlev-1)
                  tarc=segsv(tlev-1)
                  go to 400
               end if
c
               inxt=arc(tarc,itrow)
               if (inxt.eq.0) go to 400
c
               jnxt=arc(tarc,jtrow)
               if (jnxt.eq.0) go to 400
c
               if (inxt.eq.jnxt) then
                  iwt=iwtsv(tlev)+wt(tarc,itrow)+1
                  jwt=jwtsv(tlev)+wt(tarc,jtrow)+1
c
                  loopij=sdot(nlwks(inxt),c(iwt),1,s(jwt),1)
c..bhl
                  loopji=sdot(nlwks(inxt),s(iwt),1,c(jwt),1)
c..bhl
c..bug                  loopji=sdot(nlwks(inxt),s(iwt),1,c(iwt),1)
c..bhl
                  if(iwt.eq.jwt) loopji=0.0d+00
c
                  go to (411,412,413,414,415,416,417,418,419,420),
     #                                                 traksv(levm1)
                  call lnkerr('bad track')
c
  411             continue
                     call lnkerr('411')
                     go to 430
  412             continue
                     d1(i,j)=d1(i,j)+loopij*acoef(levm1)
                     d1(j,i)=d1(j,i)+loopji*acoef(levm1)
                     go to 430
  413             continue
                     d1(i,j)=d1(i,j)+loopij*acoef(levm1)
                     d1(j,i)=d1(j,i)+loopji*acoef(levm1)
                     go to 430
  414             continue
                     call lnkerr('414')
                     go to 430
  415             continue
                     call lnkerr('415')
                     go to 430
  416             continue
                     d1(i,j)=d1(i,j)+loopij*acoef(levm1)
                     d1(j,i)=d1(j,i)+loopji*acoef(levm1)
                     go to 430
  417             continue
                     d1(i,j)=d1(i,j)+loopij*acoef(levm1)
                     d1(j,i)=d1(j,i)+loopji*acoef(levm1)
                     go to 430
  418             continue
                     d1(i,j)=d1(i,j)+loopij*acoef(levm1)
                     d1(j,i)=d1(j,i)+loopji*acoef(levm1)
                     go to 430
  419             continue
                     call lnkerr('419')
                     go to 430
  420             continue
                     d1(i,j)=d1(i,j)+loopij*acoef(levm1)
                     d1(j,i)=d1(j,i)+loopji*acoef(levm1)
  430             continue
                  go to 400
c
               else
                  tlev=tlev-1
                  if (tlev.le.1) call lnkerr('error in tail segments')
                  irowsv(tlev)=itrow
                  jrowsv(tlev)=jtrow
                  iwtsv(tlev)=iwtsv(tlev+1)+wt(tarc,itrow)
                  jwtsv(tlev)=jwtsv(tlev+1)+wt(tarc,jtrow)
                  segsv(tlev)=tarc
                  tarc=0
                  itrow=inxt
                  jtrow=jnxt
                  go to 400
               end if
c
            else
c------------------------------------------------------------------
c
               iwt=iwtsv(lev)+wt(iarc(seg),irow)+1
               jwt=jwtsv(lev)+wt(jarc(seg),jrow)+1
c
               loopij=sdot(nlwks(inxt),c(iwt),1,s(jwt),1)
               loopji=sdot(nlwks(inxt),s(iwt),1,c(jwt),1)
               if(iwt.eq.jwt) loopji=0.0d+00
c
               go to (211,212,213,214,215,216,217,218,219,220),
     #                                                 traksv(levm1)
               call lnkerr('bad track')
c
  211          continue
                  call lnkerr('211')
                  go to 230
  212          continue
                  d1(i,j)=d1(i,j)+loopij*acoef(levm1)
                  d1(j,i)=d1(j,i)+loopji*acoef(levm1)
                  go to 230
  213          continue
                  d1(i,j)=d1(i,j)+loopij*acoef(levm1)
                  d1(j,i)=d1(j,i)+loopji*acoef(levm1)
                  go to 230
  214          continue
                  call lnkerr('214')
                  go to 230
  215          continue
                  call lnkerr('215')
                  go to 230
  216          continue
                  d1(i,j)=d1(i,j)+loopij*acoef(levm1)
                  d1(j,i)=d1(j,i)+loopji*acoef(levm1)
                  go to 230
  217          continue
                  d1(i,j)=d1(i,j)+loopij*acoef(levm1)
                  d1(j,i)=d1(j,i)+loopji*acoef(levm1)
                  go to 230
  218          continue
                  d1(i,j)=d1(i,j)+loopij*acoef(levm1)
                  d1(j,i)=d1(j,i)+loopji*acoef(levm1)
                  go to 230
  219          continue
                  call lnkerr('219')
                  go to 230
  220          continue
                  d1(i,j)=d1(i,j)+loopij*acoef(levm1)
                  d1(j,i)=d1(j,i)+loopji*acoef(levm1)
  230          continue
            end if
c
c           ----- go back up and try next segment -----
c
            go to 201
         end if
c
c        ----- update stacks and then go down one level -----
c
         if (lev.le.2) go to 201
c
         irowsv(levm1)=irow
         jrowsv(levm1)=jrow
         pagesv(levm1)=page
         iwtsv(levm1)=iwtsv(lev)+wt(iarc(seg),irow)
         jwtsv(levm1)=jwtsv(lev)+wt(jarc(seg),jrow)
         segsv(levm1)=seg
         irow=inxt
         jrow=jnxt
         page=newpag
c
         lev=lev-1
         levm1=lev-1
         seg=segmin(page)
         segmx=segmax(page)
         go to 201
 1000 continue
c
c     ----- halve the off-diagonals to account for ij ji symmetry -----
c
c      ij=0
c      do 1020 i=1,norbs
c         do 1010 j=1,i-1
c            ij=ij+1
c            d1(ij)=d1(ij)*0.5d+00
c 1010    continue
c         ij=ij+1
c 1020 continue
c
c
      return
      end
