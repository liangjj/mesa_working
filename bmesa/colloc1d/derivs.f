*deck derivs.f
c***begin prologue     derivs
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           functions, inverse
c***author             schneider, barry(nsf)
c***source             
c***purpose            calculate lattice representation of derivative 
c***                   operator.
c***
c***
c***references
c
c***routines called
c***end prologue
      subroutine derivs(dmat,ddmat,dbfn,ddbfn,a,t,s,e,u,v,work,npt,n,
     1                  prnti,prntd,dir)
      implicit integer (a-z)
      dimension dmat(npt,npt), ddmat(npt,npt), dbfn(npt,n), ddbfn(npt,n)
      dimension a(npt,n), s(n+1), e(n), u(npt,n), vv(n,n), work(npt)
      dimension t(n,npt)
      real*8 dmat, ddmat, dbfn, ddbfn, a, t, s, e, u, vv, work, tol
      parameter ( tol=1.d-08 )
      character*80 title
      character*(*) dir
      logical prnti, prntd
      common /io/ inp, iout
      write(iout,*) 'calculation type =',dir
      if(dir.eq.'inverse') then
         call sgefa(a,n,n,t,info)
         call sgedi(a,n,n,t,e,work,1)
         call copy(a,t,n*n)
      else
         write(iout,1) tol
         call copy(a,t,npt*n)
         call ssvdc(a,npt,npt,n,s,e,u,npt,vv,n,work,21,info)
         call copy(t,a,npt*n)
         if (info.eq.0) then
             write(iout,2)
             write(iout,3) (s(i),i=1,n)
             do 10 i=1,n
                if( abs(s(i)).le.tol ) then
                    write(iout,4) i, s(i)
                    s(i)=0.d0
                else
                    s(i)=1.d0/s(i)
                endif
   10        continue                                   
         else       
             call lnkerr('error in singular value decomposition')
         endif
         do 20 i=1,n
            do 30 j=1,n
               vv(i,j)=vv(i,j)*s(j)
   30       continue
   20    continue
         call ebct(t,vv,u,n,n,npt)
      endif
      if (prnti) then
          title='pseudoinverse of collocation matrix'
          call prntrm(title,t,n,npt,n,npt,iout)          
      endif
      call ebc(dmat,dbfn,t,npt,n,npt)
      call ebc(ddmat,ddbfn,t,npt,n,npt)
      if (prntd) then
          title='first derivative matrix'
          call prntrm(title,dmat,npt,npt,npt,npt,iout)
          title='second derivative matrix'
          call prntrm(title,ddmat,npt,npt,npt,npt,iout)                    
      endif
      return
    1 format(/,1x,'eigenvalue tolerance = ',e15.8)   
    2 format(/,1x,'all singular values computed correctly. they are:')
    3 format( /,5e15.8 )      
    4 format(1x,'eigenvalue = ',i3,1x,'value = ',e15.8,1x,
     1          'is below tol')   
      end


