*deck @(#)getspe.f	5.1  11/6/94
      subroutine getspe(spc,noinp)
c
c***********************************************************************
c     read in a code telling what type of special orbitals to handle,  *
c     or read in the matrix itself under the 'mat' option.             *
c***********************************************************************
c
      implicit integer (a-z)
      integer numint
      character*1 multrf,valenc
      character*3 codes,words*18,gvbw*18,ossw*18
      character*8 spcode
      logical noinp
c
      common /tapes/  out,errout,input,drttap
      common /dimens/ nbf,nsym,norbs,nrowsp,nrows4p,nrows,nrows4
     #,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     #,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc,nvref
     #,               nijvir
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     #,               levocc,spec,sspesh,val
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     #,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     #,               valvir,opensh,multi,speshl
      common /drtchr/ codes(9),words(9),multrf,valenc
c
      dimension spc(nspc),gvb(16),oss(16)
c
      data gvb/ 0,0,0,0, 0,1,1,1, 0,1,1,1, 0,1,1,2 /
      data oss/ 0,0,0,1, 0,1,0,1, 0,0,1,1, 1,1,1,2 /
      data gvbw/'gvb orbital       '/
      data ossw/'open shell singlet'/
      save gvb,oss,gvbw,ossw
c
      if (noinp) then
         spcode='gvb'
      else
         read (input,1,end=99) spcode
      end if
    1 format (a3)
      if (spcode.ne.'gvb') go to 6
      if (spec.eq.2) go to 3
      write (errout,2) spec
    2 format (//,' calling for gvb reference with',i3,' special '
     #,          'orbitals: impossible',//)
      call lnkerr('incorrect number of special orbitals for gvb')
    3 continue
      do 4 i=1,16
         spc(i)=gvb(i)
    4 continue
      words(9)=gvbw
      nespec=2
      na=na+1
      write (out,16)
   16 format (//,t30,'** gvb two-reference ci **')
      return
c
    6 continue
      if (spcode.ne.'oss') go to 11
      if (spec.eq.2) go to 8
      write (errout,7) spec
    7 format (//,' calling for an open-shell singlet with',i3,' special'
     #,          ' orbitals: impossible',//)
      call lnkerr('incorrect number of special orbitals for open-shell '
     #            //'singlet')
    8 continue
      do 9 i=1,16
         spc(i)=oss(i)
    9 continue
      words(9)=ossw
      nespec=2
      na=na+1
      ns=xor(ns,sspesh)
      write (out,17)
   17 format (//,t30,'** open shell singlet reference **')
      return
c
   11 continue
      if (spcode.ne.'mat') go to 90
      read (input,12,end=99) nespec,b,s
   12 format (3i5)
      read (input,13,end=99) (spc(i),i=1,nspc)
   13 format (4(4i1,1x))
      write (out,14) nespec,b,s,(spc(i),i=1,nspc)
   14 format (//,t25,'** using input matrix for special orbitals **'
     #,       //,t30,'number of electrons:',i2
     #,        /,t30,'     spin times two:',i2
     #,        /,t30,'           symmetry:',i2
     #,       //,(t30,4(4i1,1x)))
      na=na+(nespec-b)/2
      nb=nb+b
      ns=xor(ns,s-1)
      return
c
   90 continue
      write (errout,91) spcode
   91 format (//,' no code available yet for special code of ',a3,//)
      call lnkerr('no code available for this special code')
   99 continue
      write (errout,100)
  100 format (//,' could not read in all information for',i3,' special '
     #,          'orbitals',/,' please check your input deck.',//)
      call lnkerr('input deck for special orbitals is incomplete')
c
      end
