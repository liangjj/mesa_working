*deck matslv.f
c***begin prologue     matslv
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           lattice, grid, matrix
c***author             schneider, barry(nsf)
c***source             
c***purpose            enforce boundary conditions and solve mixed collocation
c***                   numerov representation of 1-d schroedinger equation
c***                               y''  + f(x) y = g.
c***
c***
c***references
c
c***routines called
c***end prologue      matslv
      subroutine matslv(bandy,rhs,work,energy,nq,order,bw,prnth,prnts)
      implicit integer (a-z)
      dimension bandy(nq,-bw:bw), rhs(nq), work(nq)
      real*8 bandy, work, rhs, energy, amplt, pshft
      character*80 title
      logical prnth, prnts
      common /io/ inp, iout
      if (prnth) then
          title='matrix'
          call prntrm(title,bandy,nq,order,nq,order,iout)
          title='rhs'
          call prntrm(title,rhs,nq,1,nq,1,iout)
      endif                          
          call sgtsl(nq-1,bandy(1,-1),bandy(1,0),bandy(1,1),rhs(1),info)
          if (info.ne.0) then
              call lnkerr('error in solution')
          endif
          rhs(nq) = ( work(4) - work(1)*rhs(nq-2) -
     1                         work(2)*rhs(nq-1) )/work(3)
          if (prnts) then
              title='solution'
              call prntrm(title,rhs,nq,1,nq,1,iout)
          endif                          
          amplt=rhs(nq)/work(5)
          write(iout,*) '     the ampltude of the irregular solution = '
     1                    ,amplt
          pshft=atan(amplt)
          write(iout,*) '     the phase shift = ',-pshft
      return
      end

