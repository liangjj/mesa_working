*deck trials.f
c***begin prologue     trials
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           iterative eigenvalue
c***author             schneider, barry (nsf)
c***source             
c***purpose            iterative linear system solver specialized 
c***                   to time-DVR using the Davidson algorithm.
c***                   
c***description        calculate an orthonormal set of trial vectors
c***
c***                   m = nz*ny*nx*nt*2
c***                   nxyzt = nz*ny*nx*nt
c***                   nxyz = nz*ny*nx
c***references         
c
c***routines called    
c***end prologue       trials
      subroutine trials(trial,rhs,thresh,nx,ny,nz,nt,ntrial,
     1                  m,n,dim,prn)
      implicit integer (a-z)
      real*8 trial, rhs
      real*8 thresh
      character*800 card
      character*3 ans
      character*80 cpass, chrkey, typetr
      logical dollar, prn
      dimension trial(m,n), rhs(m,*)
      common /io/ inp, iout
      if ( dollar('$trials',card,cpass,inp) ) then
           typetr=chrkey(card,'type-of-trial-vectors','unit',' ')
           ntrial=intkey(card,'number-of-trial-vectors',m,' ')
           ntrial=min(ntrial,n)
      endif 
      if(typetr.eq.'unit') then
         call rzero(trial,m*ntrial)
         add=max(m/2,1)
         do 10 i=1,ntrial
            trial(i,i)=1.d0
            trial(i+add,i)=1.d0
 10      continue
      elseif(typetr.eq.'guess-vectors') then
         call iosys('does "guess vectors" exist on ham',0,0,0,ans)
         if(ans.eq.'yes') then
            call iosys('read real "guess vectors" from ham',m*ntrial,
     1                  trial,0,' ')
         else
            call lnkerr('trial vectors do not exist')
         endif
      else
         call lnkerr('no option for other type trial vectors yet')
      endif    
      call gschmt(trial,thresh,m,1,ntrial,nout,.true.,.false.)
      ntrial=nout
      if(prn) then
         cpass='trial vectors for '//typetr//' trials'  
         call prntrm(cpass,trial,m,ntrial,m,ntrial,iout)
      endif	 
      write(iout,1) ntrial, typetr
 1    format(/,10x,'number of trial vectors = ',i4,/,10x,
     1             'type of trial vectors   = ',a24)
      return
      end       

