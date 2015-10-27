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
      subroutine numerv(f,g,rhs,band,work,stp,energy,r,rlast,
     1                 ipvt,l,n,order,bw,prnth,prnts)
      implicit integer (a-z)
      dimension f(n), g(n), rhs(n), ipvt(n), r(n)
      dimension  work(n), band(n,-bw:bw)
      real*8 f, g, rhs, j, y, rlast, band, work, stp
      real*8 jp, jpp, yp, ypp, energy, rr, k, r, logder
      real*8 amplt, pshft
      character*80 title
      logical prnth, prnts
      common /io/ inp, iout
      if (order.eq.3) then
          call num3pt(band,f,stp,n)
          do 10 i=2,n-1
             rhs(i) = stp*stp*( g(i-1) + 10.d0*g(i) +g(i+1) )
 10       continue
c         assume g vanishes outside boundary.
          rhs(n) = stp*stp*( g(n-1) + 10.d0*g(n) )
c     apply the boundary conditions by expressing the value of the
c     wavefunction at the last point in terms of known quantities
c     and the wavefunction at previous points.
          k=sqrt(energy)
          rr=rlast*k
          call rbes('ricatti-bessel',l,rr,j,jp,jpp,y,yp,ypp)
          logder=k*(yp/y)
          call bndryn(band,rhs,f,g,logder,stp,work,order,bw,n)
          if (prnth) then
              title='matrix'
              call prntrm(title,band,n,order,n,order,iout)
              title='rhs'
              call prntrm(title,rhs,n,1,n,1,iout)
          endif                          
          call sgtsl(n-2,band(2,-1),band(2,0),band(2,1),rhs(2),info)
          if (info.ne.0) then
              call lnkerr('error in solution')
          endif
          rhs(n) = ( work(4) - work(1)*rhs(n-2) -
     1                         work(2)*rhs(n-1) )/work(3)
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
      end
