*deck derivs.f
c***begin prologue     derivs
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           functions, inverse
c***author             schneider, barry(nsf)
c***source             
c***purpose            calculate lattice representation of derivative 
c***                   operator.
c***
c***
c***references
c
c***routines called
c***end prologue
      subroutine derivs(dmat,ddmat,dbfn,ddbfn,a,work,ipvt,npt,n,
     1                  prnti,prntd)
      implicit integer (a-z)
      dimension dmat(npt,npt), ddmat(npt,npt), dbfn(npt,n), ddbfn(npt,n)
      dimension a(npt,npt), work(npt), ipvt(n), det(2)
      real*8 dmat, ddmat, dbfn, ddbfn, a, work, det
      character*80 title
      logical prnti, prntd
      common /io/ inp, iout
      call sgefa(a,npt,npt,ipvt,info)
      if (info.ne.0) then
          call lnkerr('error in inverse collocation matrix')
      else          
          call sgedi(a,npt,npt,ipvt,det,work,1)
      endif
      if (prnti) then
          title='inverse of collocation matrix'
          call prntrm(title,a,npt,npt,npt,npt,iout)          
      endif
      call ebc(dmat,dbfn,a,npt,npt,npt)
      call ebc(ddmat,ddbfn,a,npt,npt,npt)
      if (prntd) then
          title='first derivative matrix'
          call prntrm(title,dmat,npt,npt,npt,npt,iout)
          title='second derivative matrix'
          call prntrm(title,ddmat,npt,npt,npt,npt,iout)                    
      endif
      return
      end
