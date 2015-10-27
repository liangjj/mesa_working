*deck @(#)scan.f	1.2  7/30/91
      subroutine scan(nbf,fzc,fzv,spe,uoc,cor,vir,alp,bet,occupd,
     #                speshl,virtul,frozen,rescor,resvir,opensh,
     #                valocc,valvir,multi,valenc,multrf,words,
     #                norbs,ntypes,nsym,nrefs,inp,out,
     #                ntype,bfnum,numsym,bfsym,bfkey,bfcode,symchk)
c
c***begin prologue  scan
c***date written   841219   (yymmdd)
c***revision date  850819   (yymmdd)
c***keywords  distinct row table, drt, guga-ci
c***author  saxe, paul,    (lanl)
c***purpose  to read and assimilate the drt codes form the input file
c***description
c
c
c***references
c
c***routines called  getlin, getcnt, getkey, getcod, getsym, lnkerr
c
c    modified 19 august 1985 by pws at lanl to handle orbital codes
c              typnn;s
c
c***end prologue  scan
c
      implicit integer (a-z)
c
      character*1 multrf,valenc,getkey,clef,bfkey(nbf)
      character*18 words(*)
      integer ntype(ntypes),bfnum(ntypes,nbf),numsym(ntypes,nsym)
      integer bfsym(nbf),bfcode(nbf,nrefs)
      logical symchk
c
c     ----- prepare for pretty output .... -----
c
      write (out,1)
    1 format (//,t28,'**** orbital information ****',//,
     #        t21,'basis fcts',t33,'#',t35,'key',t42,'type of orbital',
     #        /,t21,'----------',t33,'-',t35,'---',
     #          t42,'---------------')
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
         sym=getsym()
         if(.not.symchk) then 
            call ifill(bfsym(bf+1),sym,repcnt)
         endif
         write(out,*)
         if (code.le.10) then
             write (out,3) bf+1,bf+repcnt,repcnt,clef,words(code)
         else
             write (out,33) bf+1,bf+repcnt,repcnt,clef,code-10
         end if
    3    format (t22,i3,'-',i3,i4,2x,a1,t43,a18)
   33    format (t22,i3,'-',i3,i4,2x,a1,t43,'orbital type',i3)
c
c
c        ----- loop through repetition count -----
c
         begin=bf+1
         do 5 junk=1,repcnt
            bf=bf+1
            sym=bfsym(bf)
c
            if (bf.gt.nbf) then
               write (out,4) junk,repcnt,nbf
    4          format (//,' ##### drt:scan -- error in orbital codes.',
     #                 /,10x,i5,'th repetition of',i5,' exceeds nbf:',
     #                  i5,//)
               call lnkerr('scan: too many orbital codes')
            end if
c
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
            bfkey(bf)=clef
            bfcode(bf,1)=code
    5    continue
         end=bf
         write(out,55) (bfsym(i),i=begin,end)
 55      format(/,t43,'symmetry information = ',(t66,10(i1,1x)))
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
                  sym=getsym()
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
