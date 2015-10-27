*deck newqq.f
c***begin prologue     newqq
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           lattice, grid, matrix
c***author             schneider, barry(nsf)
c***source             
c***purpose            effective matrices in q-space 
c***
c***
c***references
c
c***routines called
c***end prologue      newqq
      subroutine newqq(hamqq,hampq,rhsp,hamqp,rhsq,bandy,np,nq,ntot,
     1                 bw,prnt)
      implicit integer (a-z)
      dimension hamqq(ntot,nq), hampq(ntot,nq), rhsp(np), hamqp(ntot,np)
      dimension rhsq(nq), bandy(nq,-bw:bw)
      real*8 hamqq, hampq, rhsp, hamqp, rhsq, bandy
      logical prnt
      character*80 title
      common /io/ inp, iout
      rhsq(1)=rhsq(1)-hamqp(1,np)*rhsp(np)
      hamqq(1,1)=hamqq(1,1)-hamqp(1,np)*hampq(np,1)
      if(prnt) then
         title='q-space effective rhs'
         call prntrm(title,rhsq,nq,1,nq,1,iout)
         title='q-space effective hamiltonian'
         call prntrm(title,hamqq,nq,nq,ntot,nq,iout)
      endif         
      do 10 i=1,nq
         bandy(i,0)=hamqq(i,i)
   10 continue
      do 20 i=1,nq-1
         bandy(i,1)=hamqq(i,i+1)
   20 continue
      do 30 i=2,nq
         bandy(i,-1)=hamqq(i,i-1)
   30 continue                      
      return
      end
