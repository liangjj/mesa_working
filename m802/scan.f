*deck @(#)scan.f	5.1  11/6/94
      subroutine scan(nbf,fzc,fzv,spe,uoc,cor,vir,alp,bet,occupd,
     #                speshl,virtul,frozen,rescor,resvir,opensh,
     #                valocc,valvir,multi,valenc,multrf,words,
     #                norbs,ntypes,nsym,nrefs,inp,out,
     #                ntype,bfnum,numsym,bfsym,bfkey,bfcode,symrwf,
     #                mosym)
c
c***begin prologue    scan
c***date written      841219   (yymmdd)
c***revision date     900415   (yymmdd)
c
c   15 april,  1990   rlm at lanl
c      modified to pick up symmetry information from rwf if requested.
c   19 august, 1985   pws at lanl
c      modified to handle orbital codes typnn;s
c
c***keywords          distinct row table, drt, guga-ci
c***author            saxe, paul,    (lanl)
c***purpose           to read and assimilate the drt codes form the input file
c***description
c
c
c***references
c
c***routines called   getlin, getcnt, getkey, getcod, getsym, lnkerr
c
c***end prologue  scan
c
      implicit integer (a-z)
c
      character*1 multrf,valenc,getkey,clef,bfkey(nbf)
      character*18 words(*)
      integer ntype(ntypes),bfnum(ntypes,nbf),numsym(ntypes,nsym)
      integer bfsym(nbf),bfcode(nbf,nrefs),mosym(nbf)
      logical symrwf
c
c     ----- prepare for pretty output .... -----
c
      if (symrwf) then
         write(out,101)
      else
         write (out,1)
      end if
    1 format (t6,'orbital information:',/,
     #        t10,'basis fcts',t22,'#',t24,'key',t31,'type of orbital',
     #        t48,'sym',/,t10,'----------',t22,'-',t24,'---',t31,
     #        '---------------',t48,'---')
  101 format (t6,'orbital information:',/,
     #        t10,'basis fcts',t22,'#',t24,'key',t31,'type of orbital',
     #        /,t10,'----------',t22,'-',t24,'---',t31,
     #        '---------------')
c
c     ----- do some initialising -----
c
      bf=0
      norbs=0
      sspesh=0
      call getlin
      do 30 i=1,ntypes
         ntype(i)=0
         do 29 j=1,nsym
            numsym(i,j)=0
   29    continue
   30 continue
c
c     ----- read in each code, including a repetition count, and
c                                                          process
c
    2 continue
         repcnt=getcnt()
         clef=getkey()
         code=getcod()
         if(symrwf) then
            sym=0
         else
            sym=getsym()
         end if
c
         if (code.le.10) then
            if(symrwf) then
               write(out,103) bf+1,bf+repcnt,repcnt,clef,words(code)
            else
               write (out,3) bf+1,bf+repcnt,repcnt,clef,words(code),
     $                       sym+1
            end if
         else
            if (symrwf) then
               write(out,1033) bf+1,bf+repcnt,repcnt,clef,code-10
            else
               write (out,33) bf+1,bf+repcnt,repcnt,clef,code-10,sym+1
            end if
         end if
    3    format (t12,i3,'-',i3,i4,2x,a1,t31,a18,t48,i3)
  103    format (t12,i3,'-',i3,i4,2x,a1,t31,a18)
   33    format (t12,i3,'-',i3,i4,2x,a1,t31,'type',i3,t48,i3)
 1033    format (t12,i3,'-',i3,i4,2x,a1,t31,'type',i3)
c
c        ----- loop through repetition count -----
c
         do 5 junk=1,repcnt
            bf=bf+1
c
            if (bf.gt.nbf) then
               write (out,4) junk,repcnt,nbf
    4          format (//,' ##### drt:scan -- error in orbital codes.',
     #                 /,10x,i5,'th repetition of',i5,' exceeds nbf:',
     #                  i5,//)
               call lnkerr('scan: too many orbital codes')
            end if
c
            if(symrwf) then
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
            if (clef.eq.multrf) type=multi
            if (code.gt.10) type=code
            if (type.gt.ntypes) call lnkerr('too many types')
c
            ntype(type)=ntype(type)+1
            bfnum(type,ntype(type))=bf
            numsym(type,sym)=numsym(type,sym)+1
            bfsym(bf)=sym
            bfkey(bf)=clef
            bfcode(bf,1)=code
    5    continue
      if (bf.lt.nbf) go to 2
c
c
      if (ntype(11).gt.0) return
c
c
c
c     ----- if there is more than one reference, read in the
c           additional codes.
c
      do 20 ref=2,nrefs
         write (out,6) ref
    6    format (//,t30,'reference #',i3,//,t22,'number',t35,
     #           'type of orbital',t54,'sym',/,t22,'----------',t35,
     #           '---------------',t54,'---')
c
         repcnt=0
         do 19 bf=1,nbf
            if (bfkey(bf).ne.multrf) then
               bfcode(bf,ref)=bfcode(bf,1)
            else
               if (repcnt.le.0) then
                  repcnt=getcnt()
                  clef=getkey()
                  code=getcod()
                  if (symrwf) then
                     sym=mosym(bf)-1
                  else
                     sym=getsym()
                  end if
                  write (out,8) repcnt,words(code),sym
    8             format (t22,i8,5x,a18,2x,i1)
               end if
c
               if (sym.ne.bfsym(bf)) then
                  write (out,7) ref,bf,bfsym(bf),sym
    7             format (//,' drt: scan -- error in symmetry for ',
     #                    'reference',i4,/,t10,'basis function',i4,
     #                    ' has symmetry ',i1,', code has symmetry ',i1,
     #                    //)
                  call lnkerr('drt: scan symmetry error in multi-ref')
               end if
c
               repcnt=repcnt-1
               bfcode(bf,ref)=code
            end if
   19    continue
   20 continue
c
c
      return
      end
