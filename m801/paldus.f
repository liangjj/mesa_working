*deck @(#)paldus.f	5.1  11/6/94
      subroutine paldus(ops,bfsym,bfkey,bfcode,orbtbf,a,b,s
     #,                 levpt,levnr,arc,nlwks,x,nelecs,spc,
     $                  exv,exref,spini,symi,nopen)
c
c***revision date      870907   (yymmdd)
c
c      7 september 1987  pws at lanl
c                    adding the ability to restrict configurations
c                    with a minumim and maximum number of open shells
c                    these options are drt=nopen=(min=n,max=m)
c
c***********************************************************************
c     this subroutine generates the distinct row table by searching    *
c     from the top of the graph down possible paths. note that there   *
c     is a program loop over orbitals running from statement labeled   *
c     '4'. the idea is for each level, search across all possible      *
c     a and b combinations (atest and btest) then across all points    *
c     on the row above (row) testing whether a particular arc (case)   *
c     gets down to the points given by atest,btest. the complication   *
c     lies in computing the excitation level when electrons are        *
c     excited into occupied or multi-reference orbitals. this is       *
c     accomplished using the x array which contains the number of      *
c     electrons excited into orbitals above the present point for      *
c     any walk leading to the particular point in question. therefore, *
c     points may be the same in all respects except for x value and    *
c     finally, at the fermi level the x value determines if the walk   *
c     is possible or not. for simple cases such as high spin open      *
c     shells the interacting space is implimented by counting spin-    *
c     flips as excitations into orbitals.                              *
c                                                                      *
c     the variables amax, amin, bmax, and bmin limit the portion       *
c     of the search to only those points that can be reached from      *
c     the present row. special orbitals are handled by input           *
c     explicitly the x values for all possible partial walks in the    *
c     special orbitals. originally the excitation level is set         *
c     to the sum of the general value and the value in the valence     *
c     (%) space. upon leaving the valence space it is reset to the     *
c     general value, thus giving say all singles and doubles           *
c     from references of all singles and doubles in the valence        *
c     space.                                                           *
c                                                                      *
c     maxb keeps track of the largest b value encountered. this        *
c     can then be used in the ci program to compute sufficient         *
c     coefficients for segment values. if this is done, there is       *
c     no limit to the value of the spin possible. finally the last     *
c     portion of the routine eliminates all points and arcs from       *
c     walks that dont make it from head to tail.                       *
c***********************************************************************
c
      implicit integer (a-z)
      integer numint
      logical spini,symi
      character*1 multrf,valenc,bfkey,key
      character*3 codes,words*18
      logical count
      logical logkey
      integer intkey
      character*(*) ops
c
      common /tapes/  out,errout,input,drttap
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     #,               levocc,spec,sspesh,val
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     #,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     #,               valvir,opensh,multi,speshl
      common /drtchr/ codes(9),words(9),multrf,valenc
      common /code/   fzc, fzv, cor, vir, doc, uoc, alp, bet, spe
      common /dimens/ nbf,nsym,norbs,nrowsp,nrows4p,nrows,nrows4
     #,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     #,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc,nvref
     #,               nijvir
c
      dimension bfsym(nbf),bfkey(nbf),bfcode(nrefs,nbf),orbtbf(norbs)
      dimension a(nrowsp),b(nrowsp),s(nrowsp),nlwks(nrowsp)
      dimension levpt(nlevs),levnr(nlevs),arc(nrows4p)
      dimension x(nrefs,nrowsp),nelecs(nrefs),spc(*)
      integer nopen(nrowsp)
c
c     array spc is of dimension nspc. the asterick is
c     used because the vax rejects a dimension of zero
c
      data spcnt /0/
      save spcnt
c
   92 format(5x,'number of interacting walks',i11)
   94 format(5x,'total number of walks      ',i11)
c
c
      count=logkey(ops,'drt=nopen',.false.,' ')
      minopn=intkey(ops,'drt=nopen=min',0,' ')
      maxopn=intkey(ops,'drt=nopen=max',nlevs,' ')
c
      if (count) then
         write (out,1093) minopn,maxopn
 1093    format (5x,'restricting number of open shells to ',
     $        i3,' to ',i3)
      end if
c
      maxb=-999999
c     default to doubles. move this into calling routine someday.
c     is this stuff working right?
      if(exv.eq.0) exv=2
      exref=exv+exref
      excita=exref
      do 1 ref=1,nrefs
         nelecs(ref)=2*na+nb
         x(ref,1)=0
    1 continue
      levm1=nlevs
      levpt(levm1)=0
      levnr(levm1)=1
      a(1)=na
      b(1)=nb
      s(1)=ns
      x(1,1)=0
      nopen(1)=0
      rs1=ns
c
c     ----- here begins the loop over orbitals (levels) -----
c
    4 if (levm1.le.1) go to 25
      lev=levm1
      levm1=lev-1
      if (levm1.le.levocc) excita=exref
      if (levm1.lt.levval) excita=exv
      nrowm1=0
      bf=orbtbf(levm1)
      key=bfkey(bf)
      code=bfcode(1,bf)
      nrowlv=levnr(lev)
      pontlv=levpt(lev)
      pontm1=pontlv+nrowlv
      levpt(levm1)=pontm1
      if (code.ne.spe) go to 5
      if (spcnt.eq.0) nelecs(1)=nelecs(1)-nespec
      spcnt=spcnt+1
      go to 11
    5 continue
      if (delele(bfcode(1,bf)).eq.1) rs1=xor(rs1,bfsym(bf))
      do 2 ref=1,nrefs
         nelecs(ref)=nelecs(ref)-delele(bfcode(ref,bf))
    2 continue
   11 continue
c
      amax=0
      amin=999999
      bmax=0
      bmin=999999
      do 3 row=pontlv+1,pontlv+nrowlv
         if (a(row).gt.amax) amax=a(row)
         if (a(row).lt.amin) amin=a(row)
         if (b(row).gt.bmax) bmax=b(row)
         if (b(row).lt.bmin) bmin=b(row)
    3 continue
      if (amin.gt.0) amin=amin-1
      if (bmin.gt.0) bmin=bmin-1
      if (bmax.gt.maxb) maxb=bmax
      bmax=bmax+1
c
      do 23 ajunk=1,amax-amin+1
         atest=amax-ajunk+1
         do 22 bjunk=1,bmax-bmin+1
            btest=bmax-bjunk+1
            do 21 row=pontlv+1,pontlv+nrowlv
               do 20 case=1,4
                  if (code.eq.cor.and.case.ne.4) go to 19
                  if (code.eq.vir.and.case.ne.1) go to 19
                  ia=a(row)
                  ib=b(row)
                  is=s(row)
                  iopen=nopen(row)
                  ic=levm1-ia-ib
                  go to (6,7,8,9),case
c
                  write (errout,18) case
   18             format (//,' impossible case value in paldus:',i5)
                  call lnkerr('impossible case value in paldus')
c
    6             ic=ic-1
                  dele=0
                  go to 10
c
    7             ib=ib-1
                  is=xor(is,bfsym(bf))
                  dele=1
                  if (count) iopen=iopen+1
                  go to 10
c
    8             ia=ia-1
                  ib=ib+1
                  ic=ic-1
                  is=xor(is,bfsym(bf))
                  dele=1
                  if (count) iopen=iopen+1
                  if (code.eq.alp.and.spini) dele=2
                  go to 10
c
    9             continue
                  ia=ia-1
                  dele=2
c
   10             continue
                  if (ia.lt.0.or.ib.lt.0.or.ic.lt.0) go to 19
                  if (ia.ne.atest.or.ib.ne.btest) go to 19
                  if (2*ia+ib.gt.excita.and.levm1.le.levfrm) go to 19
                  if (count.and.iopen.gt.maxopn) go to 19
                  if (count.and.levm1.eq.1.and.iopen.lt.minopn) go to 19
                  if (levm1.eq.1) iopen=0
                  possbl=0
                  do 60 ref=1,nrefs
                     if ((2*ia+ib.le.nelecs(ref)+excita).and
     #               .   (2*ia+ib.ge.nelecs(ref)-excita)) possbl=1
   60             continue
                  if (possbl.eq.0.) go to 19
c
                  do 13 rowm1=pontm1+1,pontm1+nrowm1
                     if (ia.ne.a(rowm1).or.ib.ne.b(rowm1).or.is.ne.s(
     *               rowm1).or.(count.and.iopen.ne.nopen(rowm1)))
     $                    go to 13
                     if (levm1.lt.levfrm) go to 17
                     diff=0
                     possbl=0
                     if (levval-levfrm.eq.1) possbl=1
                     if (code.ne.spe) go to 62
                     if (spcnt.eq.spec) go to 61
                     possbl=1
                     ix=4*(x(1,row)+case-1)
                     if (ix.ne.x(1,rowm1)) diff=1
                     go to 63
   61                continue
                     ref=1
                     ix=spc(x(1,row)+case)
                     frmx=0
                     if (levm1.eq.levfrm) frmx=ix
                     if ((2*ia+ib.le.nelecs(ref)+excita-ix).and
     #               .   (2*ia+ib.ge.nelecs(ref)-excita+frmx)) possbl=1
                     if (levm1.eq.levfrm) ix=0
                     if (ix.ne.x(1,rowm1)) diff=1
                     go to 63
   62                continue
                     do 40 ref=1,nrefs
                        ix=x(ref,row)
                        if (dele.gt.delele(bfcode(ref,bf)).and.key.ne.
     *                  valenc)ix=ix+dele-delele(bfcode(ref,bf))
                        if (lev.ne.levval) go to 70
                        ixsv=nelecs(ref)-(2*ia+ib)
                        ix=ix+ixsv
                        if (ix.lt.0) then
                           if (spini) ix=(ix+ib)/2
                        else
                           if (spini.and.ib.gt.ix) ix=ix+(ib-ix)/2
                        end if
                        if (ix.lt.0) ix=0
                        if (symi.and.ix.eq.0.and.is.ne.rs1.and.ixsv.eq.0
     *                  ) ix=1
                        if (2*ia+ib.le.excita-ix) ix=0
   70                   continue
                        frmx=0
                        if (levm1.eq.levfrm) frmx=ix
                        if ((2*ia+ib.le.nelecs(ref)+excita-ix).and
     #                  .   (2*ia+ib.ge.nelecs(ref)-excita+frmx))
     *                   possbl=1
                        if (levm1.eq.levfrm) ix=0
   40                if (ix.ne.x(ref,rowm1)) diff=1
   63                continue
                     if (possbl.eq.0) go to 19
                     if (diff.eq.0) go to 17
   13             continue
c
c     ----- check that this is indeed a possible new point,     -----
c           especially, if got here because a, b or s different
c
                  possbl=0
                  if (levval-levfrm.eq.1) possbl=1
                  if (code.ne.spe) go to 162
                  if (spcnt.eq.spec) go to 161
                  possbl=1
                  go to 163
  161             continue
                  ref=1
                  ix=spc(x(1,row)+case)
                  frmx=0
                  if (levm1.eq.levfrm) frmx=ix
                  if ((2*ia+ib.le.nelecs(ref)+excita-ix).and
     #            .   (2*ia+ib.ge.nelecs(ref)-excita+frmx)) possbl=1
                  go to 163
  162             continue
                  do 140 ref=1,nrefs
                  ix=x(ref,row)
                  if (dele.gt.delele(bfcode(ref,bf)).and.key.ne.valenc)
     #            ix=ix+dele-delele(bfcode(ref,bf))
                  if (lev.ne.levval) go to 170
                  ixsv=nelecs(ref)-(2*ia+ib)
                  ix=ix+ixsv
                  if (ix.lt.0) then
                     if (spini) ix=(ix+ib)/2
                  else
                     if (spini.and.ib.gt.ix) ix=ix+(ib-ix)/2
                  end if
                  if (ix.lt.0) ix=0
                  if (symi.and.ix.eq.0.and.is.ne.rs1.and.ixsv.eq.0) ix=1
                  if (2*ia+ib.le.excita-ix) ix=0
  170             continue
                  frmx=0
                  if (levm1.eq.levfrm) frmx=ix
                  if ((2*ia+ib.le.nelecs(ref)+excita-ix).and
     #            .   (2*ia+ib.ge.nelecs(ref)-excita+frmx)) possbl=1
  140             continue
  163             continue
                  if (possbl.eq.0) go to 19
c
                  nrowm1=nrowm1+1
                  rowm1=pontm1+nrowm1
                  if (rowm1.lt.nrowsp) go to 16
                  write (errout,15) nrowsp,levm1
   15             format (//,' not enough space to make drt, have only r
     *oom for',           i6,' rows and are exceeding that at level',i4)
                  call lnkerr('not enough core while making drt')
   16             a(rowm1)=ia
                  b(rowm1)=ib
                  s(rowm1)=is
                  if (count) then
                     nopen(rowm1)=iopen
                  else
                     nopen(rowm1)=0
                  end if
c
                  if (levm1.le.levfrm) go to 50
                  if (lev.ne.levval) go to 42
                  do 77 ref=1,nrefs
                  ixsv=nelecs(ref)-(2*ia+ib)
                  ix=ixsv
                  if (ix.lt.0) then
                     if (spini) ix=(ix+ib)/2
                  else
                     if (spini.and.ib.gt.ix) ix=ix+(ib-ix)/2
                  end if
                  if (ix.lt.0) ix=0
                  if (symi.and.ix.eq.0.and.is.ne.rs1.and.ixsv.eq.0) ix=1
                  if (2*ia+ib.le.excita-ix) ix=0
                  x(ref,rowm1)=x(ref,row)+ix
   77             continue
                  go to 17
c
   42             continue
                  if (code.ne.spe) go to 44
                  if (spcnt.eq.spec) go to 43
                  x(1,rowm1)=4*(x(1,row)+case-1)
                  go to 17
   43             continue
                  x(1,rowm1)=spc(x(1,row)+case)
                  go to 17
   44             continue
                  do 41 ref=1,nrefs
                     ix=x(ref,row)
                     if (dele.gt.delele(bfcode(ref,bf)).and.key.ne.
     *               valenc)ix=ix+dele-delele(bfcode(ref,bf))
                     if (2*ia+ib.le.excita-ix) ix=0
                     x(ref,rowm1)=ix
   41             continue
                  go to 17
   50             continue
                  do 51 ref=1,nrefs
                     x(ref,rowm1)=0
   51             continue
   17             continue
                  arc((row-1)*4+case)=rowm1-pontm1
   19             continue
   20          continue
   21       continue
   22    continue
   23 continue
      levnr(levm1)=nrowm1
      go to 4
   25 continue
c
c     ----- set weight of bottom of shavitt graph to one, -----
c                  eliminate all other bottoms
c
      nroot=0
      do 27 root=levpt(1)+1,levpt(1)+levnr(1)
      if (a(root).ne.0.or.b(root).ne.0.or.s(root).ne.0) go to 26
      nroot=nroot+1
      nlwks(root)=1
   26 continue
   27 continue
c
      if (nroot.eq.1) go to 29
      write (errout,28) nroot
   28 format (//,' invalid number of bottoms to graph:',i4,//)
      call lnkerr('invalid number of bottoms to graph')
   29 continue
c
c     ----- generate the weights (nlwks) of all rows -----
c
      if (lev.lt.2) go to 91
      do 32 lev=2,nlevs
         levm1=lev-1
         pontm1=levpt(levm1)
         do 31 row=levpt(lev)+1,levpt(lev)+levnr(lev)
            nlwk=0
            do 30 case=1,4
               arcpt=(row-1)*4+case
               if (arc(arcpt).gt.0) nlwk=nlwk+nlwks(arc(arcpt)+pontm1)
   30       continue
            nlwks(row)=nlwk
   31    continue
   32 continue
c
      nwks=nlwks(1)
      if (spini) then
         write (out,92) nwks
      else
         write (out,94) nwks
      end if
c
c     ----- remove all rows with zero weights -----
c
      pont=1
      do 37 junk=1,nlevs-1
         lev=nlevs-junk
         nrowlv=levnr(lev)
         pontlv=levpt(lev)
         levpt(lev)=pont
         do 36 row=pontlv+1,pontlv+nrowlv
         if (nlwks(row).eq.0) go to 35
         pont=pont+1
         a(pont)=a(row)
         b(pont)=b(row)
         s(pont)=s(row)
         nlwks(pont)=nlwks(row)
         do 34 case=1,4
            arc((pont-1)*4+case)=arc((row-1)*4+case)
   34    continue
         do 55 rowp1=levpt(lev+1)+1,levpt(lev+1)+levnr(lev+1)
         do 54 case=(rowp1-1)*4+1,(rowp1-1)*4+4
            if (arc(case).eq.row-pontlv) arc(case)=pont-levpt(lev)
   54    continue
   55    continue
         go to 36
   35    continue
         do 53 rowp1=levpt(lev+1)+1,levpt(lev+1)+levnr(lev+1)
            do 52 case=(rowp1-1)*4+1,(rowp1-1)*4+4
               if (arc(case).eq.row-pontlv) arc(case)=0
   52       continue
   53    continue
   36    continue
         levnr(lev)=pont-levpt(lev)
   37 continue
c
   91 continue
      nwks=nlwks(1)
      nrows=pont
      nrows4=nrows*4
      nrowoc=0
      do 38 lev=levfrm,nlevs
         nrowoc=nrowoc+levnr(lev)
   38 continue
      nrow4o=nrowoc*4
      return
      end
