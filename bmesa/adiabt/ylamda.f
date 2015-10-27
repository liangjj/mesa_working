*deck ylamda.f
c***begin prologue     ylamda
c***date written       001230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            compute adiabatic hyperspherical functions.
c***                   
c***references         
c
c***routines called    
c***end prologue       ylamda
      subroutine ylamda(pang,pham,n,start,prn)
      implicit integer (a-z)
      integer*8 pang, pham, pad
      real*8 angle, ham, ad, rad, fpkey
      real*8 alpha, beta, endpts, scr
      real*8 a, b, d, f
      logical dollar, nfix
      character*80 cpass, chrkey, tquad
      character*2 mstr, itoc
      character*24 str
      character*16 fptoc, key
      character*32 type 
      character*1600 card
      dimension a(36), b(36)
      dimension ngot(3), prn(*), endpts(2), nfix(2)
      common/io/inp, iout
      data nfix / .true.,.true. /
      pointer (pang,angle(1))     
      pointer (pham,ham(1))
      pointer (pad,ad(1))     
      pointer (prad,rad(1))
      pointer (pscr,scr(1))
      open(unit=90,file='curves',err=99,form='formatted',
     1             status='unknown')
c
      key='$v(hyperangle)'
      if ( dollar('$rho',card,cpass,inp) ) then
           nrad=intkey(card,'number-of-rho-points',1,' ')
           tquad=chrkey(card,'type-grid','legendre',' ')
           kpts=intkey(card,'number-fixed-endpoints',2,' ')
           endpts(1)=0.d0
           endpts(2)=0.d0
           if(kpts.ne.0) then
              call fparr(card,'region-boundaries',endpts,2,' ')
           endif
           alpha=0.d0
           beta=0.d0
           if(tquad.eq.'jacobi'.or.tquad.eq.'laguerre') then
              alpha=fpkey(card,'alpha',alpha,' ')
              beta=fpkey(card,'beta',beta,' ')
           endif
           need=wpadti(nrad+1)
           call memory(need,prad,ngot(1),'rad',0)
           need=wpadti(1+2*nrad)
           sc1=1
           sc2=sc1+nrad
           call memory(need,pscr,junk,'scr',0)
           if(nrad.gt.2) then 
              call getqpt(rad,scr(sc1),endpts(1),endpts(2),tquad,
     1                    'before',scr(sc2),nfix,nrad,nrad,1,.false.)
           else
              rad(1)=endpts(1)
              rad(2)=endpts(2)
           endif
           call memory(-junk,pscr,idum,'scr',idum)
      endif
      if(prn(1)) then
         cpass='radial grid'
         call prntrm(cpass,rad,nrad,1,nrad,1,iout)
      endif
      call vparm(type,a,b,card,cpass,key,2)
      had=1
      eigad=had+n*n 
      vad=eigad+n
      store=vad+n
      need=wpadti(store+(n+1)*nrad)
      call memory(need,pad,ngot(2),'adiabatic',0)      
      call copy(rad,ad(store),nrad)
      x=1
      y=x+n
      space=y+n 
      need=wpadti(space+5*n)
      call memory(need,pscr,junk,'scr',0)
      do 10 i=1,nrad
         call rzero(ad(vad),n)
         call vadiab(rad(i),angle(start),scr(x),scr(y),ad(vad),type,
     1               a,b,n,2)
c
c        set up the hamiltonian at this rho.
c      
         call setham(ad(had),ham,ad(vad),rad(i),n)
c
c                 diagonalize
c
         call dsyev('v','l',n,ad(had),n,ad(eigad),scr(space),
     1                              5*n,info)
         call vscale(ad(eigad),ad(eigad),1.d0/(rad(i)*rad(i)),n)
         if(prn(2)) then
            cpass='adiabatic eigenvalues at rho = '//fptoc(rad(i))
            call prntrm(cpass,ad(eigad),n,1,n,1,iout)
         endif
         if(prn(3)) then
            cpass='adiabatic eigenfunctions at rho = '//fptoc(rad(i))
            call prntrm(cpass,ad(had),n,n,n,n,iout)
         endif
         call fileig(ad(eigad),ad(store+nrad),i,n,nrad)
 10   continue   
      call wrtfil(ad(store),n,nrad)
      call memory(-ngot(1),prad,idum,'rad',idum)
      call memory(-ngot(2),pad,idum,'adiabatic',idum)      
      call memory(-junk,pscr,idum,'scr',idum)
      close(90)
      return
 99   call lnkerr('quit')
      return
 1    format(/,5x,'interaction potential type = ',a32)
      end       
