*deck @(#)reordr.f	5.1  11/6/94
      subroutine reordr(ntypes,nsym,norbs,nbf,speshl,opensh,multi,
     #                  valocc,valvir,virtul,frozen,fzc,fzv,ncodes,
     #                  ntype,bfnum,bfsym,bfcode,dela,delb,delele,
     #                  na,nb,ns,levfrm,levocc,levopn,levmul,levval,
     #                  orbfrm,
     #                  iout,orbsym,orbtbf,vtype)
c
c***begin prologue    reordr
c***date written      850102   (yymmdd)
c***revision date     910101   (yymmdd)
c
c     1 january 1991  bhl at llnl
c        fixing so that m901 can run off an m802 drt list
c    20  august 1985  modified by pws at lanl to handle groups of
c                     orbitals with code>10
c
c***keywords          drt, distinct row table
c
c***author            saxe, paul,    (lanl)
c***purpose           to reorder the orbitals from the scf to the ci (drt)
c                     ordering.
c***description
c
c
c***references
c
c***routines called   (none)
c***end prologue      reordr
c
      implicit integer (a-z)
c
      integer ntype(ntypes),bfnum(ntypes,nbf),bfsym(nbf),bfcode(nbf)
      integer iout(nbf),orbsym(norbs),dela(ncodes),delb(ncodes)
      integer delele(ncodes),orbtbf(norbs)
c
      orb=0
      levfrm=0
      levocc=-999999
      levopn=999999
      levmul=999999
      levval=999999
c
c     ----- loop through the various types of orbitals, noting key
c           levels such as the fermi-level as we reach them
c
      do 7 type=1,ntypes
         ntp=ntype(type)
         if (ntp.le.0) go to 6
c
c        ----- check for boundary levels -----
c
         if (type.eq.speshl.and.levopn.eq.999999) levopn=orb+1
         if (type.eq.opensh) levopn=orb+1
         if (type.eq.multi)  levmul=orb+1
         if (type.eq.valocc) levocc=orb+1
         if (type.eq.valvir) levval=orb+2
c
c        ----- sort orbitals within a type into ascending, or
c              for virtuals descending, symmetries. note that
c              multi reference orbitals are not sorted -----
c
         do 5 junk=1,nsym
            if (type.eq.virtul.or.type.eq.vtype) then
               sym=junk-1
            else
               sym=nsym-junk
            end if
            if (type.eq.multi.and.junk.ne.1) go to 5
            do 4 num=ntp,1,-1
               bf=bfnum(type,num)
               if (type.eq.multi) sym=bfsym(bf)
               if (sym.ne.bfsym(bf)) go to 3
               code=bfcode(bf)
               if (type.eq.frozen) then
                  if (code.eq.fzc) iout(bf)=-1
                  if (code.eq.fzv) iout(bf)=0
               else
                  orb=orb+1
                  iout(bf)=orb
                  orbsym(orb)=sym+1
                  if (code.le.10) then
                     na=na+dela(code)
                     nb=nb+delb(code)
                     if (delele(code).eq.1) ns=xor(ns,sym)
                  end if
                  orbtbf(orb)=bf
               end if
    3          continue
    4       continue
    5    continue
    6    continue
         if (type.eq.virtul.or.type.eq.vtype) levfrm=orb+1
    7 continue
c
      if (levocc.eq.-999999) levocc=levfrm
      orbfrm=max(0,levfrm-1)
c
      return
      end
