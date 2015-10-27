c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Cpsi0}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck cpsi0.f
***begin prologue     cpsi0
***date written       960723   (yymmdd)
***revision date      yymmdd   (yymmdd)
***keywords           time development
***author             schneider, barry (nsf)
***source             
***purpose            calculate the spatial part of the initial wavepacket.
***references         
***                   n = space dimension  
***routines called    
***end prologue       cpsi0
      subroutine cpsi0(u1,u2,u3,eig1,eig2,eig3,q1,q2,q3,
     1                 wt1,wt2,wt3,pq1,pq2,pq3,psi0,energy,
     2                 t0,n,nd,dim,coord,tim,i0stat,
     3                 typke,imtime,card,prn)
      implicit integer (a-z)
      real*8 u1, u2, u3, eig1, eig2, eig3
      real*8 q1, q2, q3, wt1, wt2, wt3, pq1, pq2, pq3
      real*8 energy, t0, f1, f2, fpkey, sigma, beta, x_0 
      real*8 norm, st, nfac
      complex*16 psi0, eye, f, cdotc
      character*(*) coord, typke, card, i0stat
      character*2 itoc
      logical prn, dollar, prnton, imtime
      character*80 title, chrkey
      dimension nd(3)
      dimension u1(nd(1),nd(1)), u2(nd(2),nd(2)), u3(nd(3),nd(3))
      dimension eig1(nd(1)), eig2(nd(2)), eig3(nd(3))
      dimension q1(nd(1)), q2(nd(2)), q3(nd(3))
      dimension wt1(nd(1)), wt2(nd(2)), wt3(nd(3))
      dimension pq1(nd(1),nd(1)), pq2(nd(2),nd(2)), pq3(nd(3),nd(3))
      dimension psi0(n), sigma(3), beta(3), x_0(3)
      save norm
      common/io/inp, iout
      data eye / (0.d0,1.d0) /
      data sigma / 3*1.d0 /
      data x_0 / 3*0.d0 /
      data beta / 3*1.d0 /
      pointer(p,st(1)), (p,ist(1))
      if(tim.eq.1) then
         if( dollar('$initial-state',card,title,inp) ) then  
             i0stat=chrkey(card,'driver','state-vector',' ')
             state=intkey(card,'initial-state',0,' ')
             call fparr(card,'sigma',sigma,3,' ')  
             call fparr(card,'x_0',x_0,3,' ')
             call fparr(card,'beta',beta,3,' ')  
             prnton=logkey(card,'print=on',.false.,' ')
             write(iout,1) t0, i0stat
             if(i0stat.eq.'gaussian-pulse') then 
                call nrpaket(norm,q1,q2,q3,wt1,wt2,wt3,sigma,x_0,beta,
     1                       nd,dim,coord)
             endif
         endif  

c        The allowable right hand sides are; one, a state vector of the
c        unperturbed Hamiltonian, a specified superposition of states
c        of the unperturbed Hamiltonian or a Gaussian pulse with a velocity
c        component.  The routine could easily be modified to do more.

         if(i0stat.eq.'one') then
            f = exp (-eye*t0*t0*.5d0)
            call cvfill(psi0,f,n)
         elseif(i0stat.eq.'state-vector') then

c        If the right hand side is a state vector, we need some scratch
c        memory to sort the eigenstates by energy.  Go get it.

            ind=1
            eig=iadtwp(ind+dim*n)
            need=wpadti(eig+n)
            if(dim.gt.1) then
               call getmem(need,p,ngot,'state',0)
            endif
            call cvect0(u1,u2,u3,eig1,eig2,eig3,
     1                  psi0,st(eig),energy,ist(ind),nd,
     2                  dim,n,state+1)
            f = exp(-eye*energy*t0)
            call cvscal(psi0,psi0,f,n)
            if(dim.gt.1) then
               call getmem(-ngot,p,idum,'state',idum)
            endif
         elseif(i0stat.eq.'gaussian-pulse') then
            call gpaket(q1,q2,q3,wt1,wt2,wt3,pq1,pq2,pq3,norm,psi0,
     1                  sigma,x_0,beta,n,nd,dim,typke,coord,prnton)
         elseif(i0stat.eq.'superpose') then
            call sppose(u1,u2,u3,psi0,n,nd,dim,coord)
         else
           call lnkerr('error in initial state')
         endif 
      else
         write(iout,2) t0
         call iosys ('read real solution from bec',n*2,psi0,0,' ')
      endif
      if(prn) then
         title='initial state'
         call prntcm(title,psi0,n,1,n,1,iout)
      endif
      if(tim.ne.1) then
         return
      else
         if(imtime) then
            nfac=1.d0/ sqrt ( cdotc(n,psi0,1,psi0,1) )
            call cvscal(psi0,psi0,nfac,n)
         endif
      endif 
      return
 1    format(/,1x,'initial state construction at first time = ',e15.8,
     1       /,1x,'driver                                   = ',a24)
 2    format(/,1x,'initial state from input file at time = ',e15.8)
      end       
