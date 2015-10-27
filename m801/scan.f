*deck @(#)scan.f	5.1  11/6/94
      subroutine scan(ntype,bfnum,numsym,bfsym,bfkey,bfcode,mosym,
     $                symrwf,noinp,spini)
c
c***********************************************************************
c     reads in the orbital codes using four function subroutines--     *
c     getcnt, getkey, getcod and getsym and fills the arrays ntype     *
c     with no. orbitals of a type,numsym with no. of symmetry and type,*
c     bfnum with scf no. of orbital of type and no. within type for use*
c     in reordr to reorder the orbitals. since the input subroutines   *
c     ignore blanks, the input is free format and codes may be split   *
c     over lines, etc. the last portion of this routine reads in multi-*
c     reference orbitals.                                              *
c***********************************************************************
c
      implicit integer (a-z)
      integer numint
      character*1 multrf,valenc,getkey,clef,junkc
      character*(*) bfkey(nbf)
      character*3 codes
      character*18 words
      logical noinp,spini,symrwf
c
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     #,               levocc,spec,sspesh,val
      common /dimens/ nbf,nsym,norbs,nrowsp,nrows4p,nrows,nrows4
     #,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     #,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc,nvref
     #,               nijvir
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     #,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     #,               valvir,opensh,multi,speshl
      common /drtchr/ codes(9),words(9),multrf,valenc
      common /code/   fzc, fzv, cor, vir, doc, uoc, alp, bet, spe
      common /tapes/  out,errout,input,drttap
c
      dimension ntype(ntypes),bfnum(ntypes,nbf),numsym(ntypes,nsym)
      dimension bfsym(nbf),bfcode(nrefs,nbf),mosym(nbf)
      dimension numtyp(10),tpcode(10)
c
    1 format (t10,i3,'-',i3,2x,i2,4x,a1,1x,a18,8x,i1)
   11 format (5x,'orbital information:',
     $       /7x,'functions',3x,'#',2x,'key',1x,'orbital type',
     $        7x,'symmetry')
 1000 format(5x,'orbital information:',
     $      /7x,'functions',3x,'#',2x,'key',1x,'orbital type')
 1010 format(t10,i3,'-',i3,2x,i2,4x,a1,1x,a18)
      if (symrwf) then
         write(out,1000)
      else
         write (out,11)
      end if
      bf=0
      norbs=0
      if (noinp) then
         call iosys('read integer "number of shells" from rwf',
     $        1,nshell,0,' ')
         if (nshell.lt.2) call lnkerr('scan problems with default '//
     #                         'orbital types on rwf')
         call iosys('read integer "number in shell" from rwf',-1,
     #              numtyp,0,' ')
         tpcode(1)=doc
         tpcode(nshell)=uoc
         if (nshell.eq.3) tpcode(2)=alp
         if (nshell.eq.4) then
            tpcode(2)=spe
            tpcode(3)=spe
         end if
      else
         call getlin
      end if
      tp=0
c
    2 if (bf.ge.nbf) go to 4
      tp=tp+1
      if (noinp) then
         repcnt=numtyp(tp)
         clef=' '
         code=tpcode(tp)
         sym=0
      else
         repcnt=getcnt()
         clef=getkey()
         code=getcod()
         if(symrwf) then
            sym=0
         else
            sym=getsym()
         end if
      end if
c
      if (symrwf) then
         write(out,1010) bf+1,bf+repcnt,repcnt,clef,words(code)
      else
         write(out,1)bf+1,bf+repcnt,repcnt,clef,words(code),sym+1
      end if
c
      do 3 junk=1,repcnt
         bf=bf+1
         if (symrwf) then
c           note that the symmetry indices used here are biased by -1.
            sym=mosym(bf)-1
         end if
         if (code.ne.fzc.and.code.ne.fzv) norbs=norbs+1
         type=occupd
         if (code.eq.spe) type=speshl
         if (code.eq.spe) sspesh=xor(sspesh,sym)
         if (code.eq.uoc) type=virtul
         if (code.eq.fzc.or.code.eq.fzv) type=frozen
         if (code.eq.cor) type=rescor
         if (code.eq.vir) type=resvir
         if (code.eq.alp.or.code.eq.bet) type=opensh
         if (clef.eq.valenc.and.code.ne.uoc) type=valocc
         if (clef.eq.valenc.and.code.eq.uoc) type=valvir
ctemp
         if (clef.eq.valenc) type=valvir
cend
         if (clef.eq.multrf) type=multi
         ntype(type)=ntype(type)+1
         bfnum(type,ntype(type))=bf
         numsym(type,sym+1)=numsym(type,sym+1)+1
         bfsym(bf)=sym
         bfkey(bf)=clef
         bfcode(1,bf)=code
    3 continue
      go to 2
c
    4 continue
c
c     check for nonsense in the input.
      if(norbs.gt.nbf) call lnkerr('drt: too many orbitals defined.')
c     check that we are in d2h or lower
      maxsym=0
      do 1200 i=1,nbf
         maxsym=max(maxsym,bfsym(i))
 1200 continue
      if(maxsym.gt.7) call lnkerr('ci can only handle d2h and lower')
c
c     ----- check for multi reference of valence-space cis -----
c                  to disable the interacting space
      if (ntype(multi).ne.0) spini=.false.
c
      spec=ntype(speshl)
c
      if (nrefs.lt.2) return
c
      do 10 ref=2,nrefs
         write (out,12) ref
   12    format (//,t30,'reference #',i3,//,t25,'basis fct',t35,'type of
     * ',       'orbital',t54,'sym',/,t25,'---------',t35,'--------'
     #   ,       '-------',t54,'---')
         do 9 bf=1,nbf
            if (bfkey(bf).eq.multrf) go to 5
            bfcode(ref,bf)=bfcode(1,bf)
            go to 8
    5       continue
            junkc=getkey()
            bfcode(ref,bf)=getcod()
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c original  if (bfsym(bf).eq.getsym()) go to 7
c           qetsym=getsym()
c           if (bfsym(bf).eq.qetsym) go to 7
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            if(symrwf) then
            else
               qetsym=getsym()
               if(bfsym(bf).ne.qetsym) then
                  write (errout,6) ref,bf
    6             format (//,' symmetry error in reference',i4,
     *                    ', basis fnct',i4,//)
                  call lnkerr('symmetry error in multi-reference')
               end if
            end if
            write (out,13) bf,words(bfcode(ref,bf)),bfsym(bf)+1
   13       format (t25,i5,5x,a18,2x,i1)
    8       continue
    9    continue
   10 continue
      return
      end
