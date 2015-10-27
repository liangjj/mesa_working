*deck mkpsi0.f
c***begin prologue     mkpsi0
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            calculate the spatial part of the initial wavepacket.
c***references         
c                      n = space dimension  
c***routines called    
c***end prologue       mkpsi0
      subroutine mkpsi0(u1,u2,u3,eig1,eig2,eig3,q1,q2,q3,
     1                  wt1,wt2,wt3,psi0,energy,t0,n,nd,nc,
     2                  dim,coord,tim,i0stat,card,prn)
      implicit integer (a-z)
      real*8 u1, u2, u3, eig1, eig2, eig3
      real*8 q1, q2, q3, wt1, wt2, wt3
      real*8 psi0, energy, t0, f1, f2, fpkey, sigma, beta, x_0 
      real*8 norm, expnt
      character*(*) coord, card, i0stat
      character*2 itoc
      logical prn, dollar, prnton
      character*80 title, chrkey, typpak
      dimension nd(3)
      dimension u1(nd(1),nd(1)), u2(nd(2),nd(2)), u3(nd(3),nd(3))
      dimension eig1(nd(1)), eig2(nd(2)), eig3(nd(3))
      dimension q1(nd(1)), q2(nd(2)), q3(nd(3))
      dimension psi0(n,nc,2), sigma(3), beta(3), x_0(3)
      common/io/inp, iout
      save norm
      pointer(p,ind)
      len=n*nc
      if(tim.eq.1) then
         if( dollar('$initial-state',card,title,inp) ) then  
             i0stat=chrkey(card,'driver','state-vector',' ')
             prnton=logkey(card,'print=on',.false.,' ')
             write(iout,1) t0, i0stat
             if(i0stat.eq.'gaussian-pulse') then

             endif  
         endif  
         if(i0stat.eq.'one') then
            f1 = cos(.5d0*t0*t0)
            f2 = - sin(.5d0*t0*t0)
            call vfill(psi0(1,1,1),f1,len)
            call vfill(psi0(1,1,2),f2,len)
         elseif(i0stat.eq.'state-vector') then
            state=intkey(card,'initial-state',0,' ')
            need=dim*n
            if(dim.gt.1) then
               call getmem(need,p,ngot,'ind',0)
            endif
            call vect0(u1,u2,u3,eig1,eig2,eig3,
     1                 psi0,energy,ind,nd,dim,n,nc,state+1)
            f1 = cos(energy*t0)
            f2 = - sin(energy*t0)
            call vscale(psi0(1,1,1),psi0(1,1,1),f1,n*nc)
            call vscale(psi0(1,1,2),psi0(1,1,2),f2,n*nc)
            if(dim.gt.1) then
               call getmem(-ngot,p,idum,'ind',idum)
            endif
         elseif(i0stat.eq.'gaussian-pulse') then
            call fparr(card,'sigma',sigma,3,' ')
            call fparr(card,'x_0',x_0,3,' ')
            call fparr(card,'beta',beta,3,' ')
            call nrpaket(norm,q1,q2,q3,wt1,wt2,wt3,sigma,x_0,beta,
     1                   nd,dim,coord)     
            call gpaket(q1,q2,q3,norm,psi0,sigma,x_0,beta,n,nd,
     1                  dim,typke,coord,prnton)     
         elseif(i0stat.eq.'radial-packet') then
            expnt=fpkey(card,'exponential',1.d0,' ')
            expow=intkey(card,'power',1,' ')
            typpak=chrkey(card,'type-radial-packet','exponential',' ')
            call radpakt(psi0,q1,expnt,exppow,typpak,nd(1))
         elseif(i0stat.eq.'superpose') then
            call sppose(u1,u2,u3,psi0,n,nc,nd,dim,coord)
         else
           call lnkerr('error in initial state')
         endif 
      else
         write(iout,2) t0
         call iosys ('read real solution from bec',n*nc*2,psi0,0,' ')
      endif
      if(prn) then
         title='initial state'
         call prntrm(title,psi0,n*nc,2,n*nc,2,iout)
      endif
      return
 1    format(/,1x,'initial state construction at first time = ',e15.8,
     1       /,1x,'driver                                   = ',a24)
 2    format(/,1x,'initial state from input file at time = ',e15.8)
      end       

