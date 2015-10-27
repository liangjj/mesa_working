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
      subroutine mkgrd(r,ru,wt,wtsum,rbeg,rend,stp,nreg,npts,nrmax,prnt)
      implicit integer (a-z)
      real*8  r, ru, wt, wtsum, rbeg, rend, delreg, del, stp, sum
      logical prnt, fixed
      dimension r(*), ru(*), wt(*), rbeg(nrmax), rend(nrmax)
      dimension npts(nrmax), wtsum(*)
      common /io/ inp, iout
      delreg=rend(1)-rbeg(1)
      fixed=.true.
      ptcnt=0
      wtcnt=0
      do 10 i=1,nreg
         del=rend(i)-rbeg(i)
         if(del.ne.delreg.and.npts(i).ne.npts(1)) then
            fixed=.false.
         endif
         call necote(rbeg(i),rend(i),r(ptcnt+1),wt(wtcnt+1),
     1               npts(i),.false.)
         call sumncw(wt(wtcnt+1),wtsum(ptcnt+1),npts(i))
         if (prnt) then
             call nwtprn(r(ptcnt+1),wt(wtcnt+1),npts(i))
         endif
         ptcnt=ptcnt+npts(i)
         wtcnt=wtcnt+npts(i)*(npts(i)-1)
 10   continue
      cnt=0
      ptcnt=0
      do 20 i=1,nreg
         call copy(r(cnt+1),ru(ptcnt+1),npts(i))
         cnt=cnt+npts(i)
         ptcnt=ptcnt+npts(i)-1
 20   continue   
      ptcnt=ptcnt+1
      cnt=0
      last=0
      call copy(wtsum,wt,npts(1))
      cnt=cnt+npts(1)
      last=last+npts(1)
      do 30 i=2,nreg
         wt(last)=wt(last)+wtsum(cnt+1)
         call copy(wtsum(cnt+2),wt(last+1),npts(i)-1)
         last=last+npts(i)-1
         cnt=cnt+npts(i)
 30   continue   
      sum=0.d0
      do 40 i=1,ptcnt
         sum=sum+wt(i)
 40   continue   
      write(iout,1) ptcnt
      if(fixed) then
         write(iout,2)
         stp=ru(2)-ru(1)
      else
         write(iout,3)
      endif
c      write(iout,4) (wt(i),i=1,ptcnt)
      write(iout,5) sum
      return
    1 format(/,'the unique number of grid points = ',i5)    
 2    format(/,'this is a fixed step size calculation')
 3    format(/,'this is not a fixed step size calculation')
 4    format(/,'the total integration weights',/,(5e15.8))
 5    format(/,'the sum of the integration weights = ',e15.8)
      end

