*deck set3x.f
c***begin prologue     set3x
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form three-dimensional hamiltonian in dvr representation
c***                   explicitly.  
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       set3x
      subroutine set3x(eig1,eig2,eig3,ham,eig,vec,t1,t2,t3,work,ind,
     1                 ipvt,nply,npts,dim,nroots,pottyp,prnh,prnv)
      implicit integer (a-z)
      real*8 eig1, eig2, eig3
      real*8 ham, t1, t2, t3, work, eig, vec, tmp
      character*80 title
      character*(*) pottyp
      logical prnh, prnv
      dimension eig1(nply), eig2(nply), eig3(nply), work(dim,*)
      dimension ham(dim,dim), eig(dim), t1(nply,nply), t2(nply,nply)
      dimension t3(nply,nply), ind(dim,4), ipvt(dim), vec(dim,nroots)
      common/io/inp, iout 
c     zero the hamiltonian matrix
      call rzero(ham,dim*dim)
      do 10 i=1,dim
         i1=ind(i,1)
         j1=ind(i,2)
         k1=ind(i,3)
         indi=ind(i,4)
         do 20 j=1,i
            i2=ind(j,1)
            j2=ind(j,2)
            k2=ind(j,3)
            indj=ind(j,4)          
            if(k1.eq.k2) then
               if(j1.eq.j2) then
                  ham(indi,indj) = ham(indi,indj) + t1(i1,i2)
                  if(i1.eq.i2) then
                     ham(indi,indj) = ham(indi,indj) + t2(j1,j1) 
     1                                               + t3(k1,k1)
                  endif
               else
                  if(i1.eq.i2) then
                     ham(indi,indj) = ham(indi,indj) + t2(j1,j2)
                  endif
               endif
            else
               if(i1.eq.i2.and.j1.eq.j2) then
                  ham(indi,indj) = ham(indi,indj) + t3(k1,k2)
               endif
            endif
            ham(indj,indi) = ham(indi,indj)
 20      continue   
 10   continue
c     now add in the potential energy contribution.
      if(pottyp.eq.'exponential') then
         do 30 i=1,dim
            i1=ind(i,1)
            j1=ind(i,2)
            k1=ind(i,3)
            indi=ind(i,4)
            tmp=sqrt ( eig1(i1)*eig1(i1) + eig2(j1)*eig2(j1) 
     1                                   + eig3(k1)*eig3(k1) )
            ham(indi,indi)=ham(indi,indi) + exp(-tmp)
 30      continue
      elseif(pottyp.eq.'one') then   
         do 40 i=1,dim
            i1=ind(i,1)
            j1=ind(i,2)
            k1=ind(i,3)
            indi=ind(i,4)
            ham(indi,indi)=ham(indi,indi) - 1.d0
 40      continue
      elseif(pottyp.eq.'coulomb') then
         do 50 i=1,dim
            i1=ind(i,1)
            j1=ind(i,2)
            k1=ind(i,3)
            indi=ind(i,4)
            tmp=eig1(i1)*eig1(i1) + eig2(j1)*eig2(j1)+ eig3(j1)*eig3(j1)
            ham(indi,indi) = ham(indi,indi) - 1.d0/sqrt(tmp)
 50      continue
      endif
      if(prnh) then
         title='3d-hamiltonian:dvr representation'   
         call prntrm(title,ham,dim,dim,dim,dim,iout)
      endif
      nn=nroots
      if(nroots.gt.dim) then
         nn=dim
      endif
      if(nn.eq.dim) then
         call tred2(dim,dim,ham,eig,work,ham)
         call tql2(dim,dim,eig,work,ham,ierr)
      else
          call tred1(dim,dim,ham,work(1,1),work(1,2),work(1,3))
          call imtqlv(dim,work(1,1),work(1,2),work(1,3),eig,
     1                ipvt,ierr,work(1,4))
          if(ierr.ne.0) then
             call lnkerr('error in diagonalization routine')
          endif
          call tinvit(dim,dim,work(1,1),work(1,2),work(1,3),nn,
     1                eig,ipvt,vec,ierr,work(1,4),work(1,5),
     2                work(1,6),work(1,7),work(1,8))             
          if(ierr.ne.0) then
             call lnkerr('error in diagonalization routine')
          endif
          call trbak1(dim,dim,ham,work(1,2),nn,vec)
          call copy(vec,ham,dim*nn)
      endif                   
      title='3d-eigenvalues:dvr representation'
      call prntrm(title,eig,nn,1,nn,1,iout)
      if(prnv) then
         title='3d-eigenvectors'   
         call prntrm(title,ham,dim,nn,dim,nn,iout)
      endif
      return
      end       

