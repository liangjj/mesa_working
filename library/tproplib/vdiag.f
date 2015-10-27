*deck vdiag.f
c***begin prologue     vdiag
c***date written       970797   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           3-d schroedinger equation
c***author             schneider, barry (nsf)
c***source             trap3d
c***purpose            diagonal matrix element of v in h0 representation
c***                   
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vdiag
      subroutine vdiag(v,vt,u1,u2,u3,t1,t2,n,n1,n2,n3,numdim,mattyp)
      implicit integer (a-z)
      real*8 v, vt, u1, u2, u3, t1, t2, dum
      character*(*) mattyp
      character*80 title
      dimension v(n), vt(n), u1(n1,n1), u2(n2,n2), u3(n3,n3)
      dimension t1(*), t2(*)
      common/io/ inp, iout
      if(mattyp.eq.'complex'.or.
     1              mattyp.eq.'real-unsymmetric') then
         call lnkerr('not yet implimented')
      elseif(mattyp.eq.'real-symmetric') then
         if(numdim.eq.1) then
            title='v in DVR'
            call prntrm(title,v,n,1,n,1,iout)
            title='square of transformation matrix'
            call prntrm(title,u1,n1,n1,n1,n1,iout)
            call ab3d(u1,dum,dum,v,vt,dum,dum,n1,1,1,numdim)
         elseif(numdim.eq.2) then
            call ab3d(u2,u1,dum,v,vt,t1,dum,n2,n1,1,numdim)
         elseif(numdim.eq.3) then
            call ab3d(u3,u2,u1,v,vt,t2,t1,n3,n2,n1,numdim)                 
         else
            call lnkerr('dimension incorrect')
         endif
      else
         call lnkerr('improper matrix description')
      endif                  
      title='diagonal matrix elements of potential'
      call prntrm(title,vt,n,1,n,1,iout)
      return
      end       
