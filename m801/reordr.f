*deck @(#)reordr.f	5.1  11/6/94
      subroutine reordr(ntype,bfnum,numsym,bfsym,bfkey,bfcode,orbtbf
     #,                 iout,orbsym)
c
c***********************************************************************
c     this subroutine determines the order of the orbitals in the ci   *
c     portion of the calculation. the order is set in the block data   *
c     drtcod -- 1 appears at the bottom of te graph, highest number    *
c     at the top. iout(n) gives ci number of nth scf orbital, with     *
c     a -1 for frozen cores and 0 for frozen virtual. orbtbf(n) gives  *
c     the scf number for the nth ci orbital.                           *
c***********************************************************************
c
      implicit integer (a-z)
      integer numint
      character*1 multrf,valenc,bfkey
      character*3 codes,words*18
c
      common /dimens/ nbf,nsym,norbs,nrowsp,nrows4p,nrows,nrows4
     #,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     #,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc,nvref
     #,               nijvir
      common /tapes/  out,errout,input,drttap
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     #,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     #,               valvir,opensh,multi,speshl
      common /drtchr/ codes(9),words(9),multrf,valenc
      common /code/   fzc, fzv, cor, vir, doc, uoc, alp, bet, spe
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     #,               levocc,spec,sspesh,val
c
      dimension ntype(ntypes),bfnum(ntypes,nbf),numsym(ntypes,nsym)
      dimension bfsym(nbf),bfcode(nrefs,nbf),orbtbf(norbs)
      dimension bfkey(nbf),iout(nbf),orbsym(norbs)
c
      orb=0
      na=0
      nb=0
      ns=0
      levfrm=0
      levocc=-999999
      levopn=999999
      levmul=999999
      levval=999999
c
      do 7 type=1,ntypes
         ntp=ntype(type)
         if (ntp.le.0) go to 6
         if (type.eq.speshl.and.levopn.eq.999999) levopn=orb+1
         if (type.eq.opensh) levopn=orb+1
         if (type.eq.multi ) levmul=orb+1
         if (type.eq.valocc) levocc=orb+1
         if (type.eq.valvir) levval=orb+2
         do 5 oldsym=1,nsym
            if ((type.eq.multi).and.oldsym.ne.1) go to 5
            sym=nsym-oldsym
            if (type.eq.virtul) sym=oldsym-1
            do 4 junk=1,ntp
               num=ntp-junk+1
               bf=bfnum(type,num)
               if (type.eq.multi) sym=bfsym(bf)
               if (sym.ne.bfsym(bf)) go to 3
               code=bfcode(1,bf)
               if (type.ne.frozen) go to 1
               if (code.eq.fzc) iout(bf)=-1
               if (code.eq.fzv) iout(bf)=0
               go to 2
    1          orb=orb+1
               iout(bf)=orb
               orbsym(orb)=sym+1
               na=na+dela(code)
               nb=nb+delb(code)
               if (delele(code).eq.1) ns=xor(ns,sym)
               orbtbf(orb)=bf
    2          continue
    3          continue
    4       continue
    5    continue
    6    continue
         if (type.eq.virtul) levfrm=orb+1
    7 continue
c
      if (levocc.eq.-999999) levocc=levfrm
      orbfrm=levfrm-1
c
      return
      end
