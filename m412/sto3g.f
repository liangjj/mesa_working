*deck @(#)sto3g.f	5.1  11/6/94
      subroutine sto3g(atomno,atomz,typeno,nprim,ncont,
     $                 ptprim,ptcont,ex,cf,typnam,ncart,
     $                 nctype,nshell,maxsh,maxpr,maxcf,ntypes)
c***begin prologue     sto3g
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           minimum basis, sto-3g, exponents,
c                      contraction coefficients
c***author             martin, richard (lanl)
c***source             @(#)sto3g.f	5.1   11/6/94
c***purpose            retrieves sto-3g exponents, contraction coefficients,
c                      and scale factors.
c***description
c     call sto3g(atomno,atomz,typeno,nprim,ncont,
c                ptprim,ptcont,ex,cf,typnam,ncart,
c                nctype,nshell,maxsh,maxpr,maxcf,ntypes)
c***references         (none)
c***routines called    coefs(m412), scalef(m412)
c***end prologue       sto3g
      implicit integer(a-z)
      real*8 ex(maxpr),cf(maxcf),atomz,pt1
      character*(*) typnam(ntypes)
      dimension typeno(maxsh),nprim(maxsh),ncont(maxsh),ptprim(maxsh)
      dimension ptcont(maxsh),nctype(ntypes),ncart(ntypes)
      dimension minan(19),minfrz(19)
      parameter (lensto=3,pt1=0.1d00)
      data minan/0,2,2,10,10,18,18,20,36,36,38,54,54,56,56,86,86,88,88/
      data minfrz/0,2,2,10,10,28,28,18,46,46,36,78,78,68,54,110,110,
     $            100,86/
      save minan,minfrz
c
c
      totco=0
      totpr=0
      nshell=0
      nfroz=int(float(atomno)-atomz+pt1)
      do 100 j=1,19
         if(atomno.gt.minan(j).and.nfroz.le.minfrz(j)) then
            nshell=nshell+1
            ptprim(nshell)=totpr+1
            ptcont(nshell)=totco+1
            nprim(nshell)=lensto
            ncont(nshell)=1
            call coefs(ex(ptprim(nshell)),cf(ptcont(nshell)),
     $                 j,typeno(nshell),typnam,ntypes)
            call scalef(ex(ptprim(nshell)),lensto,atomno,j)
            totpr=totpr+lensto
            totco=totco+ncont(nshell)*lensto
         endif
  100 continue
c
c
c
      return
      end
