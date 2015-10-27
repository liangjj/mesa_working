*deck numerv.f
c***begin prologue     numerv
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           lattice, grid, matrix
c***author             schneider, barry(nsf)
c***source             
c***purpose            calculate numerov representation of 
c***                          y''  + f(x) y = g.
c***
c***
c***references
c
c***routines called
c***end prologue      numerv
      subroutine numerv(f,g,rhs,bandy,bandg,work,energy,r,rlast,
     1                 ipvt,l,n,order,bw,fixed,prnth,prnts)
      implicit integer (a-z)
      dimension f(n), g(n), rhs(n), ipvt(n), r(n)
      dimension  work(n), bandy(n,-bw:bw), bandg(n,-bw:bw)
      real*8 f, g, rhs, j, y, rlast, bandy, bandg, work, stp
      real*8 jp, jpp, yp, ypp, energy, rr, k, r, logder
      real*8 amplt, pshft
      character*80 title
      logical prnth, prnts, fixed
      common /io/ inp, iout
      if (fixed) then
          stp=r(2)-r(1)
          write(iout,1) stp
      else
          stp=r(n)-r(n-1)
          write(iout,2) stp
      endif
      if (order.eq.3) then
          call num3pt(bandy,bandg,r,f,n,fixed)
          do 10 i=2,n-1
             rhs(i) = ( bandg(i,-1)*g(i-1) + bandg(i,0)*g(i) + 
     1                                bandg(i,1)*g(i+1) )
 10       continue
c         assume g vanishes outside boundary.
          rhs(n) = ( bandg(n,-1)*g(n-1) + bandg(n,0)*g(n) )
c     apply the boundary conditions by expressing the value of the
c     wavefunction at the last point in terms of known quantities
c     and the wavefunction at previous points.
          k=sqrt(energy)
          rr=rlast*k
          call rbes('ricatti-bessel',l,rr,j,jp,jpp,y,yp,ypp)
          logder=k*(yp/y)
          call bndryn(bandy,rhs,f,g,logder,stp,work,order,bw,n,n)
          if (prnth) then
              title='matrix'
              call prntrm(title,bandy,n,order,n,order,iout)
              title='rhs'
              call prntrm(title,rhs,n,1,n,1,iout)
          endif                          
          call sgtsl(n-2,bandy(2,-1),bandy(2,0),bandy(2,1),rhs(2),info)
          if (info.ne.0) then
              call lnkerr('error in solution')
          endif
          rhs(n) = ( work(4) - work(1)*rhs(n-2) -
     1                         work(2)*rhs(n-1) )/work(3)
          if (prnth) then
              title='solution'
              call prntrm(title,rhs,n,1,n,1,iout)
          endif                          
          amplt=rhs(n)/y
          write(iout,*) '     the ampltude of the irregular solution = '
     1                    ,amplt
          pshft=atan(amplt)
          write(iout,*) '     the phase shift = ',-pshft
      elseif(order.eq.5) then
          call lnkerr('order five numerov not yet implimented')
      else
          call lnkerr('order numerov incorrect')          
      endif
      return
 1    format(/,1x,'this is a fixed step size calculation.',/,1x,
     1            'the step size = ',e15.8)
 2    format(/,1x,'this is not a fixed step size calculation.',/,1x,
     1            'the last step size = 'e15.8)
      end

