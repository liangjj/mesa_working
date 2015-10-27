*deck @(#)getcf.f	5.1  11/6/94
      subroutine getcf(ex,cf,nprim,ncont,aex,acf,aptpr,aptco,anprim,
     $                 ancont,typeno,ntypes,type,nshell,small,iout,
     $                 nctype)
c***begin prologue     getcf
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)  .
c***keywords           one-electron, integrals
c***author             saxe, paul (lanl).
c***source
c***purpose            associates proper contraction coefficients with
c                      exponents.
c***description
c                      call getcf(ex,cf,nprim,ncont,aex,acf,aptpr,aptco,
c                                 anprim,ancont,typeno,ntypes,type,nshell,
c                                 small,iout,nctype)
c
c
c***references
c***routines called    rzero(math),abs
c***end prologue       getcf
      implicit integer (a-z)
c
      real*8 ex(nprim),cf(nprim,ncont,nctype),aex(*),acf(*),small
      integer aptpr(nshell),aptco(nshell),anprim(nshell)
      integer ancont(nshell),typeno(nshell)
c
c     ----- start timing -----
c
c
c     ----- loop through shells of appropriate type, putting
c           coefficients with the appropriate exponents -----
c
      call rzero(cf,nprim*ncont*nctype)
c
      nco=0
      do 20 shell=1,nshell
         if (typeno(shell).eq.type) then
            pt=aptco(shell)
            do 15 pr=1,nprim
               do 11 apr=aptpr(shell),aptpr(shell)+anprim(shell)-1
                  if (abs(ex(pr)-aex(apr)).le.small) go to 12
   11          continue
               go to 15
c
   12          continue
               do 14 c=1,ancont(shell)
                  do 13 nc=1,nctype
                     cf(pr,nco+c,nc)=acf(pt)
                     pt=pt+1
   13             continue
   14          continue
   15       continue
            nco=nco+ancont(shell)
         end if
   20 continue
c
c
c     ----- stop timing -----
c
c
c
      return
      end
