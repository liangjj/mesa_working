*deck pltwfn.f
c***begin prologue     pltwfn
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           plot, time-dependent, wavefunction
c***author             schneider, barry (nsf)
c***source             
c***purpose            plot or massage the complex time-dependent 
c***                   for analysis wavefunction.
c***                   
c***references         
c
c***routines called    
c***end prologue       pltwfn
      subroutine pltwfn(psi,psi0,q1,q2,q3,t,p1,p2,p3,pt,u1,u2,u3,
     1                  eig1,eig2,eig3,t0,energy,dim,n,m,nc,
     2                  nd,nt,type,proj,pnch)
      implicit integer (a-z)
      real*8 q1, q2, q3, t, p1, p2, p3, pt
      real*8 u1, u2, u3, eig1, eig2, eig3
      real*8 psi0, t0, psi, energy
      character*80 title
      character*(*) type
      logical proj, pnch
      dimension nd(dim)
      dimension q1(nd(1)), q2(nd(2)), q3(nd(3)), t(nt)
      dimension p1(nd(1),nd(1)), p2(nd(2),nd(2)), p3(nd(3),nd(3))
      dimension pt(nt,nt)
      dimension u1(nd(1),nd(1)), u2(nd(2),nd(2)), u3(nd(3),nd(3))
      dimension eig1(nd(1)), eig2(nd(2)), eig3(nd(3))
      dimension psi0(m,nc,2), psi(n,nc,2)
      common/io/inp, iout
      common/punch/pun(4)
c
c     if nc.ne.1 this routine is totally irrelevant
c
      if(nc.ne.1) then
         return
      endif
      if(type.eq.'one') then
c
c        special case of no spatial part
c
         call nospac(psi,psi0,pt,t,t0,n,m,nt,nc,pnch)
      elseif(type.eq.'state-vector') then
c
c        special case of eigenfunction driver
c
         call statvc(psi,psi0,q1,q2,q3,t,t0,pt,energy,dim,
     1               nd,n,m,nt,nc,pnch)
      elseif(type.eq.'general') then
c
c        analyze the wavefunction by projection onto eigenstates
c        as well as grid values.
c
         call anlyze(psi,psi0,q1,q2,q3,t,p1,p2,p3,pt,u1,u2,u3,
     1               t0,dim,nd,n,m,nt,nc,proj,pnch)
      else
         call lnkerr('error in type specification')
      endif
      return
      end       





