*deck approx.f
c***begin prologue     approx
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             timeprp
c***purpose            starting approximation to gp equation
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       approx
      subroutine approx(pndvr,obj,eigr,ham0,ham,eig,tmat,psi,work,
     2                  norm0,gamma,n,nobj,type,drctv,prnt)
      implicit integer (a-z)
      real*8 pndvr, obj, eigr, ham0, ham, eig, psi
      real*8 gamma, norm0
      real*8 work, dum, tmat
      character*(*) type
      character*80 title
      logical drctv, prnt
      dimension pndvr(n,0:n-1), obj(nobj,*), eigr(n), ham0(n,n)
      dimension ham(n,n), eig(n), psi(n), work(*), tmat(n,n)
      common/io/inp, iout
      write(iout,1)
      iter=1
      if(type.eq.'fock matrix') then
         call smul(work,work,norm0,n)
         call htree(ham0,ham,work,pndvr,psi,eigr,norm0,gamma,n)
         call mkobj(obj(1,1),ham,dum,dum,iter,n,type,.false.)
         if(prnt) then
            title='initial fock matrix'
            call prntrm(title,obj(1,1),nobj,1,nobj,1,iout)
         endif
         call tred2(n,n,ham,eig,psi,tmat)
         call tql2(n,n,eig,psi,tmat,ierr)
         call copy(tmat(1,1),work,n)
      elseif(type.eq.'wavefunction') then
         call smul(work,work,norm0,n)
         if(prnt) then
            title='initial wavefunction'
            call prntrm(title,work,nobj,1,nobj,1,iout)
         endif
         call mkobj(obj(1,1),dum,work,work,iter,n,type,.false.)
         call htree(ham0,ham,work,pndvr,psi,eigr,norm0,gamma,n)
         call tred2(n,n,ham,eig,psi,tmat)
         call tql2(n,n,eig,psi,tmat,ierr)
         call copy(tmat(1,1),work,n)   
      elseif(type.eq.'hartree') then
         call smul(work,work,norm0,n)
         call htree(ham0,ham,work,pndvr,psi,eigr,norm0,gamma,n)
         if(prnt) then
            title='initial non-linear term'
            call prntrm(title,psi,nobj,1,nobj,1,iout)
         endif
         call mkobj(obj(1,1),dum,psi,psi,iter,n,type,.false.)
         call tred2(n,n,ham,eig,psi,tmat)
         call tql2(n,n,eig,psi,tmat,ierr)
         call copy(tmat(1,1),work,n)   
      else
         call lnkerr('quit error in object type')
      endif                                                    
      return 
   1  format(/,5x,'initialization')        
      end       
