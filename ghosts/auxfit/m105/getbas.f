*deck @(#)getbas.f	1.1  11/20/92
      subroutine getbas(atsymb,atomno,atomz,basnam,typeno,nprim,ncont,
     $                  ptprim,ptcont,ex,cf,shname,nshell,mxsh,mxpr,
     $                  mxcf,typnam,nctype,ntypes,namdat,inp,iout)
c***begin prologue     getbas
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)  .
c***keywords           basis set, integrals
c***author             saxe, paul (lanl).
c***source
c***purpose            recursive search for basis set information.
c***description
c                      call getbas(atsymb,atomno,atomz,basnam,typeno,nprim,
c                                  ncont,ptprim,ptcont,ex,cf,shname,nshell,
c                                  mxsh,mxpr,mxcf,typnam,nctype,ntypes,
c                                  namdat,inp,iout)
c
c
c***references
c***routines called    lnkerr(mdutil),ffnext(chr),ctofp(chr),ctoi(cchr),
c                      skipln(chr)
c***end prologue       getbas
      implicit integer (a-z)
c
      parameter (stakmx=100)
c
      logical streqc
      real*8 atomz
      character*(*) atsymb,basnam,typnam(ntypes),shname(mxsh)
      character*80 line,token
      character*16 ffnext
      character*3 nxtlin,skipln
      character*8 type
      real*8 ctofp,scalef
      real*8 ex(mxpr),cf(mxcf)
      integer typeno(mxsh),nprim(mxsh),ncont(mxsh),ptprim(mxsh)
      integer ptcont(mxsh),nctype(ntypes),stakln(stakmx),stakfl(stakmx)
c
c     ----- initialise parameters -----
c
      nshell=0
      totpr=0
      totco=0
      stakpt=0
      token=basnam
      do 501 i=1,mxsh
         typeno(i)=0
         nprim(i)=0
         ncont(i)=0
         ptprim(i)=0
         ptcont(i)=0
  501 continue
c
c
c     ----- search for the key, which must be on a line starting
c           with '$'
c
c     first search through the input file.
  101 continue
      namfil=inp
   99 continue
      rewind namfil
      lineno=0
    1 continue
         if (nxtlin(line,namfil,lineno).eq.'eof') then
            if(namfil.eq.inp) then
c              token not found in input.  look on the data file.
               namfil=namdat
               goto 99
            else
               call lnkerr(' unexpected eof while reading data'
     $                    //' in m102')
            endif
         endif
         if (line(1:1).ne.'$') go to 1
         pos=1
         if (ffnext(line,pos,start,end).ne.'string') go to 1
         if (line(start:end).ne.token) go to 1
c
c        ----- found the right basis, so check if the atom is correct
c
         if (ffnext(line,pos,start,end).ne.'string') go to 1
         if (.not.streqc(line(start:end),atsymb)) go to 1
c
c     ----- found the right location in the file, so get the basis
c           set type etc.
c      type=xx (name=xx) (nprimitives=n) (ncontracted=n) (scalefactor=f)
c
      if (nxtlin(line,namfil,lineno).eq.'eof')
     $   call lnkerr(' unexpected eof while reading data'
     $              //' in m102')
  100 continue
         if (line(1:1).eq.'$') then
c
c            ----- either include or end of basis information
c                  at this level
c
            pos=1
             token=ffnext(line,pos,start,end)
             if (line(start:end).ne.'include') then
c
c                ----- end of section, so 'pop' stack if needed -----
c
                 if (stakpt.le.0) go to 30
c
                 namfil=stakfl(stakpt)
                 lineno=stakln(stakpt)
                 if (skipln(namfil,lineno).eq.'eof')
     $              call lnkerr(' unexpected eof while reading data'
     $                         //' in m102')
                 stakpt=stakpt-1
                 if (nxtlin(line,namfil,lineno).eq.'eof')
     $              call lnkerr(' unexpected eof while reading data'
     $                         //' in m102')
                 go to 100
              else
c
c                ----- new include, so 'push' stack -----
c
                 stakpt=stakpt+1
                 if (stakpt.gt.stakmx) then
                    write (iout,102) stakmx
  102               format (//,' ##### m102: too deep a recursion ',
     #                      'on definition of basis:',i5,//)
                    call lnkerr('m102 recursion on basis too deep')
                 end if
                 stakfl(stakpt)=namfil
                 stakln(stakpt)=lineno
                 if (ffnext(line,pos,start,end).ne.'string') then
                    write (iout,103) line
  103               format (//,' ##### m102: faulty $include ',
     #                      'statement:',/,1x,a80,//)
                    call lnkerr('m102 faulty $include')
                 end if
                 token=line(start:end)
                 go to 101
              end if
            end if
c
c     ----- increment shell counter and check we fit still -----
c
      nshell=nshell+1
      if (nshell.gt.mxsh) then
         write (iout,2) mxsh
    2    format (//,' ##### m102 exceeded maximum number of shells:',
     #           i5,//)
         call lnkerr('m102 too many shells')
      end if
c
      type=' '
      shname(nshell)=' '
      nprim(nshell)=0
      ncont(nshell)=0
      ptprim(nshell)=totpr+1
      ptcont(nshell)=totco+1
      scalef=1.0d+00
c
      pos=0
    3 continue
         token=ffnext(line,pos,start,end)
         if (token.eq.'eos') go to 10
         if (token.ne.'replacement') then
            write (iout,4) line
    4       format(//,' ##### in m102, do not understand line:',/,
     #             1x,a80,//)
            call lnkerr('m102 do not understand line')
         end if
         if (line(start:start+4).eq.'type=') then
            type=line(start+5:end)
         else if (line(start:start+4).eq.'name=') then
            shname(nshell)=line(start+5:end)
         else if (line(start:start+11).eq.'nprimitives=') then
            pos=start+11
            if (ffnext(line,pos,start,end).ne.'integer')
     $         call lnkerr(' nprimitives must be integer. i found:'//
     $                     line(start:end))
            nprim(nshell)=ctoi(line(start:end))
         else if (line(start:start+11).eq.'ncontracted=') then
            pos=start+11
            if (ffnext(line,pos,start,end).ne.'integer')
     $         call lnkerr(' ncontracted must be integer. i found:'//
     $                     line(start:end))
            ncont(nshell)=ctoi(line(start:end))
         else if (line(start:start+11).eq.'scalefactor=') then
            pos=start+11
            if (ffnext(line,pos,start,end).ne.'floating point')
     $         call lnkerr(' scalefactor must be real. i found:'//
     $                     line(start:end))
            scalef=ctofp(line(start:end))
         else if (line(start:start+5).eq.'ncore=') then
            pos=start+5
            if (ffnext(line,pos,start,end).ne.'integer')
     $         call lnkerr(' ncore must be integer. i found:'//
     $                     line(start:end))
            ncore=ctoi(line(start:end))
            atomz=float(atomno-ncore)
         else
            write (iout,5) line
    5       format (//,' ##### error in m102 on basis input:',/,1x,a80)
            line=' '
            line(start:start)='^'
            write (iout,6) line
    6       format (1x,a80,//)
         end if
      go to 3
c
c     ----- sort out names, number of contraction coefficients expected
c
   10 continue
      if (shname(nshell).eq.' ') shname(nshell)=type
      do 11 numtyp=1,ntypes
         if (type.eq.typnam(numtyp)) go to 13
   11 continue
      write (iout,12) type
   12 format (//' ##### error with shell type in m102:',a17,//)
c
   13 continue
      typeno(nshell)=numtyp
      if (ncont(nshell).eq.0) then
         if (nxtlin(line,namfil,lineno).eq.'eof')
     $      call lnkerr(' unexpected eof while reading data'
     $                 //' in m102')
         pos=0
   14    continue
         token=ffnext(line,pos,start,end)
         if (token.eq.'eos') go to 16
         if (token.ne.'floating point')
     $         call lnkerr(' exponents and contraction coefficients '//
     $                     'must be real. i found:'//line(start:end))
         ncont(nshell)=ncont(nshell)+1
         go to 14
   16    continue
         ncont(nshell)=ncont(nshell)-1
      else
         if (nxtlin(line,namfil,lineno).eq.'eof')
     $      call lnkerr(' unexpected eof while reading data'
     $                 //' in m102')
      end if
c
      if (mod(ncont(nshell),nctype(numtyp)).ne.0) then
         write (iout,17) ncont(nshell),type
   17    format (//,' ##### in m102, number of contraction ',
     #           'coefficients ',i3,' not compatible with ',
     #           'shell types ',a8,//)
         call lnkerr('m102 bad num. ccoefs for shell type')
      end if
      ncont(nshell)=ncont(nshell)/nctype(numtyp)
c
c     ----- read in the exponents and contraction coefficients -----
c
      n=0
   18 continue
            pos=0
            token=ffnext(line,pos,start,end)
            if (token.ne.'floating point') then
               if (nprim(nshell).gt.0.and.n.ne.nprim(nshell)) then
                  write (iout,29) n,nprim(nshell),type
   29             format (//,' ##### error m102 with number of ',
     #                    'primitives:',2i3,a10,//)
                  call lnkerr('m102 not as many prims. as expected')
               end if
               nprim(nshell)=n
               totpr=totpr+n
               if (totpr.gt.mxpr) then
                  write (iout,53) mxpr
   53             format (//,' ##### m102 getbas exceeded maximum ',
     #                    'number of primitives:',i6,//)
                  call lnkerr('m102 getbas: too many prims')
               end if
c
               go to 100
            end if
            ex(ptprim(nshell)+n)=ctofp(line(start:end))*scalef*scalef
            n=n+1
            nc=0
   19    continue
            token=ffnext(line,pos,start,end)
            if (token.ne.'floating point') then
               write (iout,20) line
   20          format (//,' ##### m102 error with input line:',/,1x,a80,
     #                //)
               call lnkerr('m102 error with input line')
            end if
            nc=nc+1
            cf(totco+nc)=ctofp(line(start:end))
         if (nc.lt.ncont(nshell)*nctype(numtyp)) go to 19
         totco=totco+nc
         if (totco.gt.mxcf) then
            write (iout,54) mxcf
   54       format (//' ##### m102.getbas exceeded the maximum ',
     #              'number of contraction coefficients:',i6,//)
            call lnkerr('m102 too many contraction coefficients')
         end if
         if (nxtlin(line,namfil,lineno).eq.'eof')
     $      call lnkerr(' unexpected eof while reading data'
     $                 //' in m102')
         nprim(nshell)=n
         go to 18
c
c     ----- finally ended reading input because hit $ card -----
c
   30 continue
c
c
      return
      end
