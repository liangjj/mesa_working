*deck set2x.f
c***begin prologue     set2x
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form two-dimensional hamiltonian in dvr representation
c***                   explicitly. the dynamical coordinates are (q1,q2) and
c***                   the third coordinate, q3, is constant.  
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       set2x
      subroutine set2x(eig1,eig2,q3,ham,eig,vec,t1,t2,work,
     1                 ind,ipvt,nply,npts,dim,nroots,pottyp,prnh,
     2                 prnv,title)
      implicit integer (a-z)
      real*8 eig1, eig2, q3
      real*8 ham, t1, t2, work, eig, vec, tmp1, tmp2
      character*(*) title
      character*(*) pottyp
      logical prnh, prnv
      dimension eig1(nply), eig2(nply), ham(dim,dim), eig(dim) 
      dimension t1(nply,nply), t2(nply,nply), work(dim,*)
      dimension vec(dim,nroots), ind(dim,3), ipvt(dim)
      common/io/inp, iout 
      write(iout,1) title
c     zero the hamiltonian matrix
      call rzero(ham,dim*dim)
      do 10 i=1,dim
         i1=ind(i,1)
         j1=ind(i,2)
         indi = ind(i,3)
         do 20 j=1,i
            i2=ind(j,1)
            j2=ind(j,2)
            indj = ind(j,3)
            if(j1.eq.j2) then
               ham(indi,indj) = ham(indi,indj) + t1(i1,i2)
               if(i1.eq.i2) then
                  ham(indi,indj) = ham(indi,indj) + t2(j1,j1)
               endif
            else
               if(i1.eq.i2) then
                  ham(indi,indj) = ham(indi,indj) + t2(j1,j2)
               endif
            endif
            ham(indj,indi)=ham(indi,indj)
 20      continue
 10   continue
c     now add in the potential energy contribution.
c     the third coordinate is held constant at q3 on this particular face.
      if(pottyp.eq.'exponential') then
         tmp1=q3*q3
         do 30 i=1,dim
            i1=ind(i,1)
            j1=ind(i,2)
            indi = ind(i,3)
            tmp2=sqrt( eig1(i1)*eig1(i1) + eig2(j1)*eig2(j1) +tmp1 )
            ham(indi,indi)=ham(indi,indi) + exp(-tmp2)          
 30      continue
      elseif(pottyp.eq.'one') then   
         do 40 i=1,dim
            i1=ind(i,1)
            j1=ind(i,2)
            indi = ind(i,3)
            ham(indi,indi) = ham(indi,indi) - 1.d0
 40      continue
      elseif(pottyp.eq.'coulomb') then
         tmp1=q3*q3
         do 50 i=1,dim
            i1=ind(i,1)
            j1=ind(i,2)
            indi = ind(i,3)
            tmp2=eig1(i1)*eig1(i1) + eig2(j1)*eig2(j1) + tmp1
            ham(indi,indi) = ham(indi,indi) - 1.d0/sqrt(tmp2)
 50      continue
      endif
      if(prnh) then
         title='2d-hamiltonian:dvr representation'   
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
      title='2d-eigenvalues:dvr representation'
      call prntrm(title,eig,nn,1,nn,1,iout)
      if (prnv) then
          title='2d-eigenvectors:dvr representation'
          call prntrm(title,ham,dim,nn,dim,nn,iout)
      endif
 1    format(a80)
      return
      end       

