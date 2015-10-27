*deck solvep.f
c***begin prologue     solvep
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           lattice, grid, matrix
c***author             schneider, barry(nsf)
c***source             
c***purpose            formal solution of p-space equations
c***references
c
c***routines called
c***end prologue      solvep
      subroutine solvep(hampp,hampq,rhsp,ipvt,np,nq,n,prnt)
      implicit integer (a-z)
      dimension hampp(n,np), hampq(n,nq), rhsp(*), ipvt(*)
      real*8 hampp, hampq, rhsp
      logical prnt
      character*80 title
      common /io/ inp, iout
      call sgefa(hampp,n,np,ipvt,info)
      if (info.ne.0) then
          call lnkerr('error in p-space solution')
      endif          
      call sgesl(hampp,n,np,ipvt,rhsp,0)
      call sgesl(hampp,n,np,ipvt,hampq,0)
      if (prnt) then
          title='solution for p-space right hand side'
          call prntrm(title,rhsp,np,1,np,1,iout)
          title='solution for matrix coupling p to q-space'
          call prntrm(title,hampq,np,nq,n,nq,iout)
      endif
      return
      end
