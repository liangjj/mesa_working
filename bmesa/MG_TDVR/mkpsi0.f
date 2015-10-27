*deck mkpsi0.f
c***begin prologue     mkpsi0
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            calculate the spatial part of the initial wavepacket.
c***references         
c
c***routines called    
c***end prologue       mkpsi0
      subroutine mkpsi0(u1,u2,u3,eig1,eig2,eig3,p1,p2,p3,q1,q2,q3,
     1                  wt1,wt2,wt3,psi0,energy,t0,n,nd,
     2                  state,dim,coord,tim,i0stat,prnt)
      implicit integer (a-z)
      real*8 u1, u2, u3, eig1, eig2, eig3, p1, p2, p3
      real*8 q1, q2, q3, wt1, wt2, wt3
      real*8 psi0, energy, t0 
      character*(*) coord, i0stat
      logical prnt
      character*80 title
      dimension nd(3)
      dimension u1(nd(1),nd(1)), u2(nd(2),nd(2)), u3(nd(3),nd(3))
      dimension eig1(nd(1)), eig2(nd(2)), eig3(nd(3))
      dimension p1(nd(1),nd(1)), p2(nd(2),nd(2)), p3(nd(3),nd(3))
      dimension q1(nd(1)), q2(nd(2)), q3(nd(3))
      dimension wt1(nd(1)), wt2(nd(2)), wt3(nd(3))
      dimension psi0(n,2)
      common/io/inp, iout
      pointer(p,ind)
      if(tim.eq.1) then
         if(i0stat.eq.'one') then
            do 10 i=1,n
c               psi0(i,1)=1.d0
c               psi0(i,2)=0.d0
               psi0(i,1) = cos(.5d0*t0*t0)
               psi0(i,2) = - sin(.5d0*t0*t0)
 10         continue   
         elseif(i0stat.eq.'state-vector') then
            need=wptoin(dim*n)
            if(dim.gt.1) then
               call memory(need,p,ngot,'ind',0)
            endif
            call vect0(u1,u2,u3,eig1,eig2,eig3,psi0,energy,
     1                 ind,nd,dim,n,state)
            do 20 i=1,n
               psi0(i,2) = - psi0(i,1)*sin(energy*t0)
               psi0(i,1) = psi0(i,1)*cos(energy*t0)
 20         continue   
            if(dim.gt.1) then
               call memory(-ngot,p,idum,'ind',idum)
            endif
         elseif(i0stat.eq.'gaussian-pulse') then
            call gpaket(u1,u2,u3,p1,p2,p3,q1,q2,q3,wt1,wt2,wt3,
     1                  psi0,n,nd,dim,coord)
         else
           call lnkerr('error in initial state')
         endif 
      else
         call iosys ('read real solution from bec',n*2,psi0,0,' ')
      endif
      if(prnt) then
         title='initial state'
         call prntrm(title,psi0,n,2,n,2,iout)
      endif
      return
      end       
