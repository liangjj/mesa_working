*deck ehamxt.f
c***begin prologue     ehamxt
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             timeprp
c***purpose            construct hamiltonian in dvr representation explicitly.
c***                   
c***                                                 2     2 
c***description          d                        [ d     r  ]
c                      i -  psi(r,t) + hbar*omega [ - 2 - -  ] psi(r,t)
c                        dt                       [ dr    4  ]
c                                    - v(r,t) psi(r,t) = v(r,t) psi0(r,t)
c***references         
c
c***routines called    
c***end prologue       ehamxt
      subroutine ehamxt(hr,ht,eigr,eigt,psi0,ham,rhs,hbar,omegat,omega,
     1                  ind,nply,dim,prnh)
      implicit integer (a-z)
      real*8 hr, eigr, eigt, psi0, hbar, omegat, omega
      complex*16 ht, ham, rhs, tmp
      character*80 title
      logical prnh
      dimension nply(2)
      dimension eigt(nply(1)), eigr(nply(2)), ham(dim,dim) 
      dimension ht(nply(1),nply(1)), hr(nply(2),nply(2)), psi0(nply(2))
      dimension rhs(dim), ind(dim,3)
      common/io/inp, iout 
      write(iout,1)
c     zero the hamiltonian matrix
      call czero(ham,dim*dim)
      call vscale(hr,hr,hbar*omegat,nply(2)*nply(2))
      do 10 i=1,dim
         i1=ind(i,1)
         j1=ind(i,2)
         indi = ind(i,3)
         do 20 j=1,dim
            i2=ind(j,1)
            j2=ind(j,2)
            indj = ind(j,3)
            if(j1.eq.j2) then
               ham(indi,indj) = ham(indi,indj) + ht(i1,i2)
               if(i1.eq.i2) then
                  ham(indi,indj) = ham(indi,indj) - hr(j1,j1)
               endif
            else
               if(i1.eq.i2) then
                  ham(indi,indj) = ham(indi,indj) - hr(j1,j2)
               endif
            endif
 20      continue
 10   continue
c     now add in the potential energy contribution and calculate the rhs.
      do 30 i=1,dim
         i1=ind(i,1)
         j1=ind(i,2)
         indi = ind(i,3)
         tmp=eigr(j1)*cos(omega*eigt(i1))   
         ham(indi,indi)=ham(indi,indi) - tmp
         rhs(indi) = tmp*psi0(j1)       
 30   continue
      if(prnh) then
         title='space-time hamiltonian:dvr representation'   
         call prntcm(title,ham,dim,dim,dim,dim,iout)
         title='driving term:dvr representation'
         call prntcm(title,rhs,dim,1,dim,1,iout)
      endif
 1    format(/,1x,'construction of complex time-dependent hamiltonian')
      return
      end       

