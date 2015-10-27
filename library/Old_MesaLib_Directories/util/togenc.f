*deck @(#)togenc.f	5.1  11/6/94
      subroutine togenc(numex,numcf,ntypes,natoms,ptprim,ptcont,
     $                  nprim,ncont,ex,cf,mxpr,mxcf,nshell,typeno,
     $                  anprim,ancont,aptpr,aptco,aex,acf,
     $                  atom,iout,nctype)
c***begin prologue     togenc
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           basis set
c***author             saxe, paul (lanl)
c***source
c***purpose            rearranges basis set bookkeeping arrays into a form
c                      appropriate for general contraction schemes.
c***description
c                      call togenc(numex,numcf,ntypes,natoms,ptprim,ptcont,
c                                  nprim,ncont,ex,cf,mxpr,mxcf,nshell,typeno,
c                                  anprim,ancont,aptpr,aptco,aex,acf,atom,
c                                  iout,nctype)
c
c
c***references
c***routines called    lnkerr(mdutil),getcf(util)
c***end prologue       togenc
      implicit integer (a-z)
c
      real*8 ex(mxpr),cf(mxcf),aex(*),acf(*),small
      integer ptprim(natoms,ntypes),ptcont(natoms,ntypes)
      integer nprim(natoms,ntypes),ncont(natoms,ntypes)
      integer typeno(nshell),anprim(nshell),ancont(nshell)
      integer aptpr(nshell),aptco(nshell),nctype(ntypes)
c
      data small /1.0d-07/
      save small
c
c     ----- start timing -----
c
c
c     ----- find all shells of one type -----
c
      aptex=0
      aptcf=0
      do 100 type=1,ntypes
         ptprim(atom,type)=numex+1
         ptcont(atom,type)=numcf+1
         ncont(atom,type)=0
         nprim(atom,type)=0
         do 10 shell=1,nshell
            if (typeno(shell).eq.type) then
               do 5 apr=aptpr(shell),aptpr(shell)+anprim(shell)-1
                  do 1 pr=ptprim(atom,type),ptprim(atom,type)+
     #                                      nprim(atom,type)-1
                     if (abs(ex(pr)-aex(apr)).lt.small) go to 5
    1             continue
                  nprim(atom,type)=nprim(atom,type)+1
                  numex=numex+1
                  if (numex.gt.mxpr) then
                     write (iout,2) mxpr
    2                format (//,' ##### l031 genbas: too many ',
     #                       'primitives total!!! ',i6,//)
                     call lnkerr('l301 genbas: too many primitives')
                  end if
                  ex(numex)=aex(apr)
    5          continue
               ncont(atom,type)=ncont(atom,type)+ancont(shell)
            end if
   10    continue
         numcf=numcf+ncont(atom,type)*nprim(atom,type)*nctype(type)
c
c        ----- and form the contraction coefficient matrix -----
c
         if (nprim(atom,type).le.0) go to 100
c
         call getcf(ex(ptprim(atom,type)),cf(ptcont(atom,type)),
     #              nprim(atom,type),ncont(atom,type),aex,acf,
     #              aptpr,aptco,anprim,ancont,typeno,ntypes,type,
     #              nshell,small,iout,nctype(type))
  100 continue
c
c     ----- stop timing -----
c
c
c
      return
      end
