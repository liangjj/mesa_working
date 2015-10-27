c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Fdiff}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck fdiff.f 
***begin prologue     fdiff
***date written       960718   (yymmdd)
***revision date               (yymmdd)
***keywords           finite difference
***author             schneider, b. i.(nsf)
***source             shared
***purpose            form finite difference approximations to hamiltonian
***description        the pointer phamil is returned which points to the 
***                   starting location of the array containing these 
***                   quantities.  
***                   for the spatial
***                   variables(size) the order is;
***                   hamiltonian(nphy*nphy)
***                   potential(nphy)
***                   eigenvectors(nphy*nphy)
***                   eigenvalues(nphy*nphy) 
***                   coordinates(n) 
***references       
***routines called     
***end prologue       fdiff
      subroutine fdiff(pham,coord,nphy,ndiff)
c
      implicit integer (a-z)
      parameter ( maxreg=500 )
      character*(*) coord
      character*2 itoc, ichr
      character*24 key, typpot
      character*80 chrkey, cpass, prnkey
      character*800 card
      logical dollar, logkey
      logical prn, reuse, addv
      real*8 edge, fpkey, d, del, ham, dscale 
#ifdef DECPOINTER
      integer*8 pham 
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 pham 
#endif SGIPOINTER
      dimension edge(2), d(4)
      dimension prn(7), prnkey(7)
      common/io/inp, iout      
      data prnkey / 'global-points','global-weights','potential',
     1              'hamiltonian','eigenvalues','eigenvectors','all' /
      pointer (pham,ham(1))
      len=length(coord)
      if ( dollar('$h0('//coord(1:len)//')',card,cpass,inp) ) then
           do 10 i=1,6
              prnkey(i)='prnt=fd='//prnkey(i)
 10        continue
           prnkey(7)=chrkey(card,'print',prnkey(7),' ')   
           if(prnkey(7).eq.'all') then
              call setprn(prn,6)
           else
              call setlog(prn,prnkey,card,6)
           endif
           call fparr(card,'region-boundary',edge,2,' ')
           ndiff=intkey(card,'order-of-finite-difference-formula',3,' ')
           addv=logkey(card,'add-potential',.false.,' ')

c          Read in number of points

           nstep=intkey(card,'number-of-steps',3,' ')
           del=(edge(2)-edge(1))/nstep

c          The actual number of unspecified points is nphy.  The value of
c          the solution at the first and last point are fixed to be zero.

           nphy=nstep-1
           n=nstep+1
           write(iout,1) ndiff, n, del
           write(iout,2) (edge(i),i=1,2)
           key='$v0('//coord(1:len)//')'
      endif
      do 20 i=1,4
         d(i)=0.d0
 20   continue   
      row=ndiff/2+1
      h=1
      v=h+row*nphy
      q=v+nphy
      wt=q+n
      need=wpadti(wt+n)
      call getmem(need,pham,ngot,'ham',0)
      call points(ham(q),ham(wt),edge(1),del,n,ndiff,prn(1)) 
      call vone(key,ham(q+1),ham(v),dscale,typpot,nphy,
     1          .true.,prn(3))
      if(ndiff.eq.3) then
         call band3(ham(h),ham(v),d,del,addv,nphy,prn(4))
      elseif(ndiff.eq.5) then
         call band5(ham(h),ham(v),d,del,addv,row,nphy,prn(4))
      elseif(ndiff.eq.7) then
         call band7(ham(h),ham(v),d,del,addv,row,nphy,prn(4))
      else
         call lnkerr('error in fd type')
      endif 
      return
 1    format(/,5x,'order of finite difference formula     = ',i2,
     1       /,5x,'number of points(including end points) = ',i5,
     2       /,5x,'step size                              = ',e15.8)
 2    format(/,15x,'region edges = ',/,5(e15.5,1x))
      end

