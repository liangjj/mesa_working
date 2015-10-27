*deck hamtot.f
c***begin prologue     hamtot
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            time-dependent hamiltonian not including
c***                   non-linear term
c***                   
c***description        total hamiltonian in space and time dimension
c***                   contructed explicitly.  
c***references         
c
c***routines called    
c***end prologue       hamtot
      subroutine hamtot(ham,h1,h2,h3,ht,v,ind,nd,nt,nc,n,m,dim,prn)
      implicit integer (a-z)
      real*8 h1, h2, h3, v, ham
      real*8 ht, fspace, ftime
      logical prn
      character*80 title
      dimension nd(dim)
      dimension ham(m,m), ht(nt,nt), v(n,nc,nc), ind(*)
      dimension h1(nd(1),nd(1)), h2(nd(2),nd(2)), h3(nd(3),nd(3))
      common/io/inp, iout      
      call rzero(ham,m*m)
      lcol=n*nc
      if(dim.eq.1) then
         call setd2(ind,nd(1),nt,nc)
         call hamfl1(ham(1,1),ham(1,lcol+1),ham(lcol+1,1),
     1               ham(lcol+1,lcol+1),h1,ht,v,ind,nd(1),nt,
     2               n,nc,m,prn)
      elseif(dim.eq.2) then
         call setd3(ind,nd(2),nd(1),nt,nc)
         call hamfl2(ham(1,1),ham(1,lcol+1),ham(lcol+1,1),
     1               ham(lcol+1,lcol+1),h1,h2,ht,v,ind,nd,nt,
     2               n,nc,m,prn)
      elseif(dim.eq.3) then
         call setd4(ind,nd(3),nd(2),nd(1),nt,nc)
         call hamfl3(ham(1,1),ham(1,lcol+1),ham(lcol+1,1),
     1               ham(lcol+1,lcol+1),h1,h2,h3,ht,v,ind,nd,nt,
     2               n,nc,m,prn)
      else
         call lnkerr('error in dimension')
      endif  
      if (prn) then
          title='full space and time-dependent hamiltonian'
          call prntrm(title,ham,m,m,m,m,iout)
      endif     
      return
      end       
