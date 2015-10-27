*deck mkgrd.f
c***begin prologue     mkgrd
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           collocation, scattering
c***author             schneider, barry (nsf)
c***source             mesa
c***purpose            grid generation 
c***description        generate grid using quadrature rules
c***references         
c
c***routines called 
c***                   
      subroutine mkgrd(type,r,ru,wt,rbeg,rend,nreg,npts,ptcnt,
     1                 nrmax,ngmax,prnt)
      implicit integer (a-z)
      real*8  r, ru, wt, rbeg, rend, wtsum
      character*(*) type
      logical prnt
      dimension r(*), ru(*), wt(*), rbeg(nrmax), rend(nrmax)
      dimension npts(nrmax)
      common /io/ inp, iout
      if (type.eq.'newton-cotes'.or.type.eq.'default') then
          ptcnt=0
          wtcnt=0
          do 10 i=1,nreg
             call necote(rbeg(i),rend(i),r(ptcnt+1),wt(wtcnt+1),
     1                   npts(i),.false.)
             if (prnt) then
                 call nwtprn(r(ptcnt+1),wt(wtcnt+1),npts(i))
             endif
             ptcnt=ptcnt+npts(i)
             wtcnt=wtcnt+npts(i)*(npts(i)-1)
 10       continue
          cnt=0
          ptcnt=0
          do 20 i=1,nreg
             call copy(r(cnt+1),ru(ptcnt+1),npts(i))
             cnt=cnt+npts(i)
             ptcnt=ptcnt+npts(i)-1
 20       continue   
          ptcnt=ptcnt+1
      else
          if (prnt) then
              write(iout,1)
          endif
          ptcnt=0
          wtsum=0.d0
          do 30 i=1,nreg-1
             if (npts(i).ge.4) then
                 do 40 pt=1,npts(i)
                    ptcnt=ptcnt+1
                    call lgndrx (npts(i),pt,wt(ptcnt),r(ptcnt))
                    wt(ptcnt)=wt(ptcnt)*(rend(i)-rbeg(i))
                    wtsum=wtsum+wt(ptcnt)
                    r(ptcnt)=r(ptcnt)*(rend(i)-rbeg(i))+rbeg(i)
                    if (prnt) then
                        write(iout,2) ptcnt, r(ptcnt), wt(ptcnt)
                    endif                
 40              continue
             else
                 call lnkerr('bad call to lgndrx')
             endif         
 30       continue
          ptcnt=ptcnt+1
          wt(ptcnt)=0.d0
          r(ptcnt)=rend(nreg)
          if (prnt) then
              write(iout,2) ptcnt, r(ptcnt), wt(ptcnt)
          endif                
          write(iout,3) wtsum
          call copy(r,ru,ptcnt)
      endif
      write(iout,4) ptcnt
      return
    1 format(/,5x,'   point  ',6x,'    value    ',10x,'    weight    ')
    2 format(8x,i3,8x,e15.8,8x,e15.8)
    3 format(/,'     sum of integration weights = ',e15.8)
    4 format(/,'the unique number of grid points = ',i5)    
      end

