*deck @(#)reordc.f	1.1  11/30/90
      subroutine reordc(c,tc,nco,nao,nbf)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)reordc.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8 (a-h,o-z)
      dimension c(*),tc(*)
c
      common /io/ inp,iout
c
c
c  reorder the orbitals so the active orbitals
c  preceed the inactive orbitals
c
c-3      write(iout,*)' entry reordc '
      if(nco.eq.0.or.nao.eq.0) return
c
      n1=nao*nbf
      n2=nco*nbf
      n3=n1+n2
c
c
      call scopy(n1,c(n2+1),1,tc,1)
c
      call scopy(n2,c,1,tc(n1+1),1)
c
      call scopy(n3,tc,1,c,1)
c
      return
      end
