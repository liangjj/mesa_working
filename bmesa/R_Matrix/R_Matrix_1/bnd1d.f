*deck bnd1d.f
c***begin prologue     bnd1d
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            search for bound state poles using r-matrix matching.
c***                   
c***references         
c
c***routines called    
c***end prologue       scat1
      subroutine bnd1d(pham,e0,e1,rbox,cnverg,n,nroots,step,itmax,prn)
      implicit integer (a-z)
      integer*8 pham, prmat 
      real*8 rmat, ham, e0, e1, cnverg, rbox, dele, kappa
      real*8 dfob, drmat, rtest, big
      logical prn
      character*80 title
      common/io/inp, iout
      data big / 1.d+06 /
      pointer (pham,ham(1))
      pointer (prmat,rmat(1))
c
c
c     set the pointers to the r-matrix amplitudes
c
      nc=1
      no=nc
      dele=(e1-e0)/step
      ene=1
      rmt=ene+step
      zval=rmt+nc*nc*step
      root=zval+nroots+1
      fob=root+nroots
      need=wpadti(fob+nc*nc*step)
      call getmem(need,prmat,ngot,'rmat',0)         
c
c     set the other pointers
c
      hmat=1
      v=hmat+n*n
      eigvec=v+n
      gamma=eigvec+n*n
      gamma=gamma+n
      eigval=gamma+n
c
      write(iout,1)
      write(iout,2) nroots, e0, e1, dele   
      call mke(rmat(ene),e0,dele,step)
      pr=rmt
      pene=ene
      pf=fob
c
c     the objective is to find a zero of the function
c
c             F(E) = 1. + kappa * Rmat
c                    kappa = sqrt(-2.d0*E)
c
      count=0
      do 10 ns=1,step
	 if(rmat(pene).le.0) then
	    call rzero(rmat(pr),nc*nc)
            call conrmat(rmat(pr),rmat(pr),ham(eigval),ham(gamma),
     1                   ham(gamma),rmat(pene),nc,nc,nc,n,
     2                   .false.,.false.)
            count=count+1
            kappa=sqrt(-2.d0*rmat(pene))
            call fobj(rmat(pf),dfob,rmat(pr),drmat,kappa,nc,.false.)
            pr=pr+nc*nc
            pf=pf+nc*nc
            pene=pene+1
         endif
 10   continue
      call intrv(rmat(fob),rmat(zval),e0,dele,count,nroots,
     1           nrtfnd,.true.)
      if(prn) then
         write(iout,3)
      endif
      lroot=root
      locz=zval
      do 20 nr=1,nrtfnd
         call nwtrap(rmat(locz),rmat(locz+1),rmat(lroot),
     1               rtest,cnverg,itmax,iter,ham(eigval),ham(gamma),n)
c
c        check if this is a real zero
c
         if(abs(rtest).ge.big) then
	    write(iout,4) iter, rmat(lroot), rtest
	 else   
            write (iout,5) iter, rmat(lroot), rtest
         endif
         lroot=lroot+1
         locz=locz+1
   20    continue
      call getmem(-ngot,prmat,idum,'rmat',idum)
      return      
    1 format(/,5x,'searching on r-matrix')
    2 format(/,5x,'search for ',i3,' roots in range [',e15.8,',',
     1            e15.8' ]',/,5x,'function evaluation every ',
     2            e15.8,' bohr')
    3 format(/,1x,'no. newton-raphson iterations',
     1         4x,'value of root',4x,'object function')
    4 format (11x,i3,18x,e15.8,3x,e15.8,'  not a root')
    5 format (11x,i3,18x,e15.8,3x,e15.8)
      end       


















