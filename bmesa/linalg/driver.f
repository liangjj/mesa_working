c $Header$
*deck driver.f
c***begin prologue     driver
c***date written       930119   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           driver, link 6201
c***author             schneider, barry (lanl)
c***source             m6203
c***purpose            driver for solution of inhomogeneous wave equation
c***
c***description        the inhomogeneous wave equation is solved by
c***                   expansion of the solution in a set of analytic
c***                   functions. the set consists of cartesian gaussians
c***                   and a set of outgoing waves which have been 
c***                   regularized at the origin using the rescigno/orel
c***                   inhomogeneous equation technique. the objective in
c***                   this routine is to generate an iterate for the
c***                   variation-iteration method. the approach is to
c***                   use the co-ordinate space representation of the
c***                   previous vector to perform the product of vloc
c***                   on that vector efficiently. the sum over channels
c***                   is then carried out, which reduces the equation
c***                   to the form desired.
c***                   
c***
c***                   vloc(pts,ntri) = local potential matrix in 
c***                                    co-ordinate space. it is stored
c***                                    region by region. in each 
c***                                    region the points are nr*nth*nph.
c***                                    the triangular label is a
c***                                    channel label (c,c') c> c'
c***                                                          =   
c***                                    only a single region need be in
c***                                    storage at once.
c***
c***                   psi(pts,chan)  = a scattering wavefunction.
c***                                    from a previous iteration.
c***                                    also stored region by region.
c***                                    the full vector is required.
c***
c***                   rhs(pts,chan)   = sum over channels of vloc*psi
c***                                     only a single region need be in
c***                                     storage at once.
c***
c***                   phifn(nph,m)    = sin(m*phi)/cos(m*phi) angular
c***                                     functions. storage is again region
c***                                     by region. both functions appear
c***                                     for each m from m = 0 to mmax    
c***
c***                   plm(nth,m:lmax) = legendre functions p(l,m).
c***                                     storage is again compact and 
c***                                     region by region
c***
c***                   tphi(nph)       = temporary storage for the phifn
c***                                     needed for a given channel and region
c***                                     to vectorize the legendre
c***                                     decomposition
c***
c***                   tplm(nth,*)     = temporary storage for the p(l,m)
c***                                     for a given m to vectorize the
c***                                     legendre decomposition.
c***
c***                   scr             = scratch array nr*nth for the
c***                                     region under consideration to 
c***                                     hold a temporary matrix of the
c***                                     decomposition.
c***
c***                   flm             = radial expansion coefficients
c***                                     for rhs for all of the channels
c***                                     of a region. storage is packed
c***                                     region by region.   
c***references         
c***routines called
c***end prologue       driver
      subroutine driver (vloc,psi,rhs,phifn,plm,flm,greg,gireg,wphi,
     1                   wthet,tphi,tplm,scr,lval,mval,nl,nm,mpnt,lpnt,
     2                   ns,nr,nth,nph,nchan,mmax,maxl,maxm,cordir,
     3                   nrhs,prnt)
      implicit integer(a-z)
      real*8 psi, vloc, rhs, phifn, plm, flm, tphi, tplm, scr
      real*8 greg, gireg, wphi, wthet, wr
      logical prnt
      character*5 fcall
      character*(*) cordir
      dimension psi(*), vloc(*), rhs(*), nr(ns), nth(ns), nph(ns)
      dimension phifn(*), plm(*), tphi(*), tplm(*), mpnt(ns)
      dimension scr(*), flm(*), energy(nchan), wph(*), wthet(*)
      dimension lpnt(0:mmax,ns), nl(maxm,nchan), nm(nchan)
      dimension lval(maxl,maxm,nchan), mval(maxm,nchan)
      dimension greg(*), gireg(*), cordir(4)
      save fcall
      data fcall / 'vlen' /
      if (fcall.eq.'vlen') then
          len=0
          do 10 is=1,ns
             len=max(len,nr(is)*nth(is)*nph(is)*nchan)
   10     continue
          fcall='no'  
      endif
      ntri=nchan*(nchan+1)/2
      locv=1
      locpsi=1
      loclm=1
      locphi=1
      locthe=1
      if (cordir(1).ne.'in-core') then
          call iosys ('rewind "local potential matrix" '//
     1                'on lamdat read-and-write',0,0,0,' ')
      endif
      do 20 is=1,ns
c     begin calculation for points in region is 
         if (prnt) then
             write(iout,*) '     processing shell = ',is
         endif
         npts=nr(is)*nth(is)*nph(is)
c        read in potential matrix for this region.
         if (cordir(1).ne.'in-core') then
             call iosys ('read real "local potential matrix" from '//
     1                   'lamdat without rewinding',npts*ntri,
     2                    vloc(locv),0,' ')
         endif 
c
c        form the sum over channels d of vloc(c,d)*psi(d) for this region.
c        psi is needed for all physical channels to do this and we
c        produce a product called rhs for all channels. we assume
c        that this can be done region by region. i will worry about non-local
c        potentials later. if they are separable there will be no problem.  
c
         call rzero(rhs,npts*nchan)  
         call vmake(vloc(locv),psi(locpsi),rhs,npts,nchan,ntri)
c        we now have psi for all channels for this region. up the counter
         locpsi=locpsi+npts*nchan
c
         if (cordir(1).eq.'in-core') then
             locv=locv+npts*ntri
         endif
         locm=mpnt(is)
c
c        the matrix elements are accumulated region by region as we
c        go so we have to store them for all the channels.
c
c        the basis functions consist of the radial free waves multiplied
c        by spherical harmonics and the bound orbitals. the assumption is that
c        we can hold in main storage all the spherical harmonics needed for a
c        given region. for the bound 
c
         locrhs=1
         do 30 chan=1,nchan
            do 40 mu=1,nm(chan)
               m=abs(mval(mu,chan))
c           get the phi functions for this m and store in tphi
c           they are scaled by the integration weights.
               call gtherm(phifn(locm),tphi,wphi(locphi),nph(is),
     1                     mval(mu,chan))
c           get the l values for this channel and store in tplm            
c           they are scaled by the integration weights.
               call gtherl(plm(lpnt(m,is)),tplm,wthet(locthe),nth(is),
     1                     lval(1,mu,chan),nl(mu,chan),m,lmax)
c           we may now calculate the projection using matrix algebra
               if (prnt) then
                   write(iout,1000) chan, mval(mu,chan),
     1                                 (lval(k,mu,chan),k=1,nl(mu,chan))
               endif
               call matel(
               loclm=loclm+nl(mu,chan)*nr(is)
   40       continue
c           up index for next channel
            locrhs=locrhs+npts
   30    continue
         locphi=locphi+nph(is)
         locthe=locthe+nth(is)
   20 continue
      return
 1000 format (/,5x,'channel = ',i3,1x,'m value = ',i3,
     1        /,5x,'l values:',/,(10x,10(i3,1x)))
      end
