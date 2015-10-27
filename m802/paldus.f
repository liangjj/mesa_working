*deck @(#)paldus.f	5.1  11/6/94
      subroutine paldus(nlevs,na,nb,ns,nrowmx,norbs,ncode,cor,vir,
     #                  levfrm,excita,nrows,nwks,out,
     #                  levpt,levnr,a,b,s,x,orbcod,delele,orbsym,
     #                  arc,nlwks,numel,symtp,spin,ntypes)
c
c***begin prologue   paldus
c***date written     841212   (yymmdd)
c***revision date    850820   (yymmdd)
c
c    20 august 1985 modified by pws at lanl to handle special groups of
c                   orbitals (type>10). this added numel,symtp, spin
c                   and ntypes to the argument list.
c
c***keywords         drt,distinct row table,unitary group ci,
c                    configuration interaction
c***author           saxe,paul, (lanl)
c***purpose          to form the rough distinct row table (drt)
c***description
c
c   on entry:
c
c      nlevs    integer
c               the number of levels in the graph (norbs+1)
c
c      na       integer
c               the 'a' value of the top of the graph
c
c      nb       integer
c               the 'b' value of the top of the graph
c
c      ns       integer
c               the symmetry (0, 1, ...) of the top of the graph
c
c      nrowmx   integer
c               the maximum number of rows the graph can have
c
c      norbs    integer
c               the number of orbitals spanned by the graph
c
c      ncode    integer
c               the number of possible different orbital types
c
c      cor      integer
c               the code for a restricted core orbital
c
c      vir      integer
c               the code for a restricted virtual orbital
c
c      levfrm   integer
c               the fermi-level, ie. the level above the virtual space
c
c      excita   integer
c               the excitation level allowed from the reference
c
c      out      integer
c               unit for writing output messages, warnings, etc.
c
c      ntypes   number of possible orbital types. those >10 are special
c               groups of orbitals for which occupancies are specified
c
c   arrays on input:
c
c      orbcod   integer (norbs)
c               the reference walk in terms of orbital codes
c
c      delele   integer (ncode)
c               the chaonge in the number of electrons for each
c               possible orbital type
c
c      orbsym   integer (norbs)
c               the symmetry (0, 1, ....) of the orbitals
c
c      numel    integer (ntypes)
c               number (cumulative) of electrons for type>10
c
c      symtp    integer (ntypes)
c               cumulative spin for type>10, <0 if no specification
c
c      spin     integer (ntypes)
c               cumulative spin (0,1,2...) for type>10, <0 if no
c               specification
c
c   scalars returned:
c
c      nrows    integer
c               the actual number of rows in the graph
c
c      nwks     integer
c               the number of walks (configurations) in this drt
c
c   arrays returned:
c
c      levpt    integer (nlevs)
c               the offset for the rows of each level
c
c      levnr    integer (nlevs)
c               the number of rows in each level
c
c      a        integer (nrows)
c               the 'a' values of the rows
c
c      b        integer (nrows)
c               the 'b' values of the rows
c
c      s        integer (nrows)
c               the symmetry associated with the rows
c
c      x        integer (nrows)
c               the excitation value associated with the rows
c
c      arc      integer (4,nrows)
c               the connections between rows on adjacent levels
c
c      nlwks    integer (nrows)
c               the number of lower walks from each row
c
c***references
c***routines called  pakdrt, lwrwks
c***end prologue  paldus
c
      implicit integer (a-z)
c
      integer levpt(nlevs),levnr(nlevs),a(nrowmx),b(nrowmx),s(nrowmx)
      integer x(nrowmx),orbcod(norbs),delele(ncode),orbsym(norbs)
      integer arc(4,nrowmx),nlwks(nrowmx),numel(ntypes),symtp(ntypes)
      integer spin(ntypes)
      logical spini
c
      data spini /.false./
      save spini
c
c
c     ----- set up the top row of the graph -----
c
      levpt(nlevs)=0
      levnr(nlevs)=1
      a(1)=na
      b(1)=nb
      s(1)=ns
      x(1)=0
      arc(1,1)=0
      arc(2,1)=0
      arc(3,1)=0
      arc(4,1)=0
      rs1=ns
      nelerf=2*na+nb
c
c     ----- check for group orbitals if their restrictions meet the
c           top row
c
      code=orbcod(nlevs-1)
      if (code.gt.10) then
         if (spin(code-10).lt.0) spin(code-10)=nb
         if (symtp(code-10).lt.0) symtp(code-10)=ns
         if (numel(code-10).ne.nelerf.or.nb.ne.spin(code-10).or.
     #       ns.ne.symtp(code-10)) then
            write(out,*) ' problem with group orbital specifications:'
            write(out,*) '   numel ',numel(code-10),nelerf
            write(out,*) '   spin ',spin(code-10),nb
            write(out,*) '   symtyp ',symtp(code-10),ns
            call lnkerr(' ')
         end if
      end if
c
c     ----- now descend down through graph, adding rows and arcs -----
c
      do 1000 level=nlevs,2,-1
         nrowm1=0
         code=orbcod(level-1)
         nrowlv=levnr(level)
         ptlv=levpt(level)
         ptlvm1=ptlv+nrowlv
         levpt(level-1)=ptlvm1
c
c        ----- if singly occupied, update symmetry -----
c
         if (code.le.10) then
            if (delele(code).eq.1) rs1=xor(rs1,orbsym(level-1)-1)
c
c           ----- determine the number of electrons left in the ref
c
            nelerf=nelerf-delele(code)
         end if
c
c        ----- decide the range of a and b values possible for the
c               next level down -----
c
         amax=0
         amin=999999
         bmax=0
         bmin=999999
         do 1 row=ptlv+1,ptlv+nrowlv
            amax=max(amax,a(row))
            amin=min(amin,a(row))
            bmax=max(bmax,b(row))
            bmin=min(bmin,b(row))
    1    continue
         amin=max(0,amin-1)
         bmin=max(0,bmin-1)
         bmax=bmax+1
c
c        ----- loop through the possible rows in this next level -----
c
         do 23 atest=amax,amin,-1
            do 22 btest=bmax,bmin,-1
c
c              ----- loop through the rows on the level above, and
c                    see if any connect to the row under consideration -----
c
               do 21 row=ptlv+1,ptlv+nrowlv
                  do 20 case=1,4
                     if (code.eq.cor.and.case.ne.4) go to 19
                     if (code.eq.vir.and.case.ne.1) go to 19
                     ia=a(row)
                     ib=b(row)
                     is=s(row)
                     ic=level-1-ia-ib
c
c                    ----- determine which row this case takes us to
c
                     go to (6,7,8,9), case
                     call lnkerr('bad case value')
c
c                    ----- unoccupied arc -----
c
    6                continue
                        ic=ic-1
                        dele=0
                        go to 10
c
c                    ----- alpha occupied arc -----
c
    7                continue
                        ib=ib-1
                        is=xor(is,orbsym(level-1)-1)
                        dele=1
                        go to 10
c
c                    ----- beta occupied arc -----
c
    8                continue
                        ia=ia-1
                        ib=ib+1
                        ic=ic-1
                        is=xor(is,orbsym(level-1)-1)
                        dele=1
c                       if (code.eq.alp.and.spini) dele=2
                        go to 10
c
c                    ----- doubly occupied arc -----
c
    9                continue
                        ia=ia-1
                        dele=2
   10                continue
c
c                    ----- check if this is the row of interest -----
c
                     if (ia.ne.atest.or.ib.ne.btest.or.ic.lt.0) go to 19
c
c                    ----- if we are down to the virtual space, check
c                          that the number of electrons left is not
c                          greater than the excitation level -----
c
                     nele=2*ia+ib
                     if (code.gt.10) then
c
c                       ----- check for groups of orbitals constraints
c
                        if (level.gt.2.and.code.ne.orbcod(level-2))
     #                                                             then
                           lcode=orbcod(level-2)-10
                           if (nele.ne.numel(lcode)) go to 19
                           if (symtp(lcode).ge.0.and.
     #                                   symtp(lcode).ne.is) go to 19
                           if (spin(lcode).ge.0.and.spin(lcode).ne.ib)
     #                                                         go to 19
                        end if
                     else
                        if (level-1.le.levfrm.and.nele.gt.excita)
     #                                              go to 19
c
c                       ----- check that we are within the excitation
c                          level number of electrons from the reference
c
                        if (nele.gt.nelerf+excita.or.
     #                      nele.lt.nelerf-excita) go to 19
                     end if
c
c                    ----- determine if this is a new row to add to the
c                          graph, or one already in the graph
c
                     do 13 rowm1=ptlvm1+1,ptlvm1+nrowm1
                        if (ia.ne.a(rowm1).or.ib.ne.b(rowm1).or.
     #                      is.ne.s(rowm1)) go to 13
c
c                       ----- if we are in the virtual space, or using
c                             groups of orbitals, a b and s
c                             are the only attributes of the row, so
c                             we have a match, and must add the arc ---
c
                        if (level-1.lt.levfrm.or.code.gt.10) go to 17
c
c                       ----- determine the excitation level of this
c                                                      row
c
                        ix=x(row)
                        if (dele.gt.delele(code)) ix=ix+dele-
     #                                                  delele(code)
c
c                       ----- check again if too far from the
c                                                         reference -----
c
                        if (level-1.eq.levfrm) then
                           frmx=ix
                        else
                           frmx=0
                        end if
                        if (nele.gt.nelerf+excita-ix.or.
     #                      nele.lt.nelerf-excita+frmx) go to 19
c
c                       ----- check if ix value equals that of this
c                                                            row -----
c
                        if (level-1.eq.levfrm) ix=0
                        if (ix.eq.x(rowm1)) go to 17
   13                continue
c
c                    ----- at this point, the new row does not match
c                          any in the graph, but before adding this row,
c                          we must again check the excitation levels to
c                          see if it is indeed a valid addition.
c
                     if (code.le.10) then
                        ix=x(row)
                        if (dele.gt.delele(code)) ix=ix+dele-
     #                                               delele(code)
                        if (level-1.eq.levfrm) then
                           frmx=ix
                        else
                           frmx=0
                        end if
                        if (nele.gt.nelerf+excita-ix.or.
     #                      nele.lt.nelerf-excita+frmx) go to 19
                        if (level-1.le.levfrm) ix=0
                    else
                        ix=0
                    end if
c
c                    ----- finally, this new row is not already in the
c                          graph, and is a valid row from all points of
c                          view. therefore, add it to the graph.
c
                     nrowm1=nrowm1+1
                     rowm1=ptlvm1+nrowm1
c
c                    ----- check that have enough space... -----
c
                     if (rowm1.gt.nrowmx) then
                        write (out,15) nrowmx,level-1
   15                   format (//,' ##### drt: paldus -- not enough ',
     #                          'space to finish drt...have maximum ',
     #                          'of',i6,' rows at level',i4,//)
                        call lnkerr('not enough space for drt')
                     end if
c
c                    ----- add row -----
c
                     a(rowm1)=ia
                     b(rowm1)=ib
                     s(rowm1)=is
                     x(rowm1)=ix
                     arc(1,rowm1)=0
                     arc(2,rowm1)=0
                     arc(3,rowm1)=0
                     arc(4,rowm1)=0
c
c                    ----- add the arc from the upper level to this row
c                          n.b. this section is shared with cases where
c                          the row was already in the graph
c
   17                continue
                     arc(case,row)=rowm1-ptlvm1
c
   19                continue
   20             continue
   21          continue
   22       continue
   23    continue
c
c        ----- save the number of rows at this new level -----
c
         levnr(level-1)=nrowm1
 1000 continue
c
c     ----- set the total number of rows -----
c
      nrows=levpt(1)+levnr(1)
c
c     ----- combine redundant rows -----
c
      call pakdrt(nlevs,nrows,a,b,s,x,arc,levpt,levnr,nlwks,.false.)
c
c     ----- find the number of lower walks from each row and
c           eliminate useless rows
c
      call lwrwks(a,b,s,x,nlwks,arc,levpt,levnr,nlevs,nrows,nwks,out)
c
c
      return
      end
