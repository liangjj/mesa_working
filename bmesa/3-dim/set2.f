*deck set2.f
c***begin prologue     set2
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form two-dimensional hamiltonian in fbr representation
c***                   explicitly. the dynamical coordinates are (q1,q2) and
c***                   the third coordinate, q3, is constant.  
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       set2
      subroutine set2(p1,p2,q1,wts1,q2,wts2,q3,rho1,rho2,v,vint,ham,
     1                eig,vec,t1,t2,work,ipvt,nply,npts,tri,dim,
     2                nroots,pottyp,prnt,title)
      implicit integer (a-z)
      real*8 p1, p2, q1, wts1, q2, wts2, q3, rho1, rho2
      real*8 v, ham, vec, t1, t2, vint, work, eig, tmp1, tmp2, tmp3
      real*8 matel
      character*(*) title
      character*(*) pottyp
      logical prnt
      dimension p1(npts,0:nply-1), p2(npts,0:nply-1)
      dimension q1(npts), wts1(npts), q2(npts), wts2(npts)
      dimension rho1(npts,tri), rho2(npts,tri)
      dimension v(npts,npts), ham(dim,dim), t1(nply,nply), t2(nply,nply)
      dimension vint(npts,tri), eig(dim), work(dim,*)
      dimension vec(dim,nroots), ipvt(dim)
      common/io/inp, iout
      write(iout,1) title
c     calculate the two density matrices.  they are not identical because
c     the basis sets are not the same in each dimension.
      count=0
      do 10 i=1,nply
         do 20 j=1,i
            count=count + 1
            do 30 k=1,npts
               rho1(k,count) = p1(k,i-1)*p1(k,j-1)*wts1(k)
 30         continue
 20      continue
 10   continue
      count=0           
      do 40 i=1,nply
         do 50 j=1,i
            count=count + 1
            do 60 k=1,npts
               rho2(k,count) = p2(k,i-1)*p2(k,j-1)*wts2(k)
 60         continue
 50      continue
 40   continue
c     calculate the value of the potential on the grid.
c     the third coordinate is held constant at q3 on this particular face. 
      call rzero(v,npts*npts)
      if(pottyp.eq.'exponential') then
         tmp1=q3*q3
         do 100 k1=1,npts
            tmp2 = tmp1 + q1(k1)*q1(k1)
            do 110 k2=1,npts
               tmp3 = sqrt( tmp2 +q2(k2)*q2(k2) )
               v(k1,k2) = exp(-tmp3)
 110        continue
 100     continue   
      elseif(pottyp.eq.'one') then   
         do 200 k1=1,npts
            do 210 k2=1,npts
               v(k1,k2) = -1.d0
 210        continue
 200     continue
      elseif(pottyp.eq.'coulomb') then
         tmp1=q3*q3
         do 300 k1=1,npts
            tmp2=q1(k1)*q1(k1) + tmp1
            do 310 k2=1,npts
               tmp3 = sqrt( tmp2 + q2(k2)*q2(k2) )
               v(k1,k2) = -1.d0/tmp3
 310        continue
 300     continue
      endif
      call ebc(vint,v,rho2,npts,npts,tri)
      write(iout,*) vint
      call rzero(ham,dim*dim)
      cntik=0
      do 500 i=1,nply
         do 510 k=1,i
            cntik=cntik+1
            cntjl=0
            do 600 j=1,nply
               do 610 l=1,j
                  cntjl=cntjl+1
                  matel=0.d0
                  do 620 k1=1,npts
                     matel=matel+rho1(k2,cntik)*vint(k2,cntjl)
 620              continue
                  bigi=nply*(i-1) +j
                  bigj=nply*(k-1) +l
                  ham(bigi,bigj)=matel+ham(bigi,bigj)
                  bigi=nply*(k-1)+j
                  bigj=nply*(i-1)+l
                  ham(bigi,bigj)=matel+ham(bigi,bigj)
                  bigi=nply*(i-1) +l
                  bigj=nply*(k-1) +j
                  ham(bigi,bigj)=matel+ham(bigi,bigj)
                  bigi=nply*(k-1)+l
                  bigj=nply*(i-1)+j
                  ham(bigi,bigj)=matel+ham(bigi,bigj)
 610           continue
 600        continue
 510     continue
 500  continue
c     put in the q1-coordinate part of the kinetic energy.
c     its diagonal in (q2)
      ii1=0
      do 700 i1=1,nply
         do 710 i2=1,nply
            ii2=0
            do 720 j1=1,nply
               bigi=ii1 + j1
               bigj=ii2 + j1
               ham(bigi,bigj) = ham(bigi,bigj) + t1(i1,i2)
 720        continue
            ii2=ii2+nply
 710     continue
         ii1=ii1+nply
 700  continue
c     put in the q2-coordinate part of the kinetic energy.
c     its diagonal in (q1)
      ii1=0
      do 800 i1=1,nply
         do 810 j1=1,nply
            do 820 j2=1,nply
               bigi=ii1 + j1
               bigj=ii1 + j2
               ham(bigi,bigj)=ham(bigi,bigj) + t2(j1,j2)
 820        continue
 810     continue
         ii1=ii1+nply
 800  continue
      if(prnt) then
         title='2d-hamiltonian:finite basis representation'   
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
      endif                              
      title='2d-eigenvalues:finite basis representation '
      call prntrm(title,eig,dim,1,dim,1,iout)
      return
 1    format(a80)
      end       

