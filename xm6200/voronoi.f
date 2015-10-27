*deck %W%  %G%
c***begin prologue     voronoi
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           voronoi, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            calculate voronoi weights
c***references         becke papers in jcp on lda and poisson equation.
c
c***routines called
c***end prologue       voronoi
      subroutine voronoi (eta,center,znuc,grid,wt,yuk,rpt,rij,d,p,alpha,
     1                    cutoff,ncen,ncplus,ngrid,nwts,niter,nr,nrmax,
     1                    numshl,noscat,nonsep,prnt)
      implicit integer (a-z)
      character*3 itoc, chra
      character*30 rtyp, str
      logical prnt, noscat, nonsep
      real*8 eta, center, znuc, grid, wt, rij, d, p, norm, fac, diff
      real*8 s, pi, yuk, yukint, aprvol, yukext, alpha, rad, cutoff
      real*8 fac1, aij, rpt
      dimension center(3,ncplus), eta(ncplus), znuc(ncplus), grid(3,*)
      dimension rij(ncplus,ncplus), d(ncplus), p(*), ngrid(ncplus)
      dimension wt(*), nwts(ncplus), nr(numshl), rpt(nrmax,numshl)
      dimension yuk(*) 
      common/io/inp, iout
      data pi / 3.141592653589793d0  /
      if (ncplus.gt.ncen) then
          write(iout,1) alpha, niter
      endif    
c                  
c            calculate internuclear distances      
      do 10 atomi=1,ncplus
         do 20 atomj=1,atomi
            rij(atomi,atomj)=sqrt((center(1,atomi)-center(1,atomj))*
     1                            (center(1,atomi)-center(1,atomj)) +
     2                            (center(2,atomi)-center(2,atomj))*
     3                            (center(2,atomi)-center(2,atomj)) +
     4                            (center(3,atomi)-center(3,atomj))*
     5                            (center(3,atomi)-center(3,atomj)))
            rij(atomj,atomi)=rij(atomi,atomj)
   20    continue
   10 continue
c
c             calculate the weights for all the centers
c                  including the scattering center
c
      locgr=0
      locwt=0
      yukint=0.d0
      aprvol=0.d0   
      do 30 atomi=1,ncplus
         if (noscat) then
             chra=itoc(atomi)
             len=cskipb(chra,' ')
             str='atom-'//chra(1:len)
         else
             if (atomi.ne.ncplus) then
                 chra=itoc(atomi)
                 len=cskipb(chra,' ')
                 str='atom-'//chra(1:len)
             else
                 str='scattering center'
             endif
         endif                                          
         call iosys ('read integer "number of shells '//
     1                str//'" from lamdat',1,nshell,0,' ')
         call iosys ('read integer "number radial points per '//
     1               'shell '//str//'" from lamdat',nshell,nr,0,' ')
         call iosys('read character "radial quadrature type '//
     1               str//'" from lamdat',-1,0,0,rtyp)
         if (.not.nonsep) then   
             call iosys('read integer "theta quadrature order '//
     1                   str//'" from lamdat',1,nthet,0,' ')  
             call iosys('read integer "phi quadrature order '//
     1                   str//'" from lamdat',1,nphi,0,' ')
             nang=nthet*nphi
         else
             call iosys ('read integer "number lebedev points '//
     1                    str//'" from lamdat',1,nang,0,' ') 
             nphi=npts
             nthet=npts
         endif              
         call rdpt(rpt,nshell,nr,numshl,nrmax,str,atomi)           
         keepg=locgr
         do 40 pt=1,ngrid(atomi)
            locgr=locgr+1
            do 50 atomj=1,ncplus
               d(atomj)=sqrt((grid(1,locgr)-center(1,atomj))*
     1                       (grid(1,locgr)-center(1,atomj)) +
     2                       (grid(2,locgr)-center(2,atomj))*
     3                       (grid(2,locgr)-center(2,atomj)) +
     4                       (grid(3,locgr)-center(3,atomj))*
     5                       (grid(3,locgr)-center(3,atomj))) 
   50       continue
c   
c           initialize the sum of the weights as that from the
c                         scattering center
c 
            if (ncen.eq.ncplus) then
                norm=0.d0
            else       
                rad=sqrt( grid(1,locgr)*grid(1,locgr) +
     1                    grid(2,locgr)*grid(2,locgr) +
     2                    grid(3,locgr)*grid(3,locgr) )
                norm=0.d0
                if (rad.ge.cutoff) then
                    norm=1.d0
                endif
            endif    
            p(pt)=norm
c         
c           note that these inner loops only go over the true atomic
c           centers so that atomi can never be equal to atomj
c            
            do 60 atomj=1,ncen
               fac=1.d0
               do 70 atomk=1,atomj-1
                  diff=(d(atomj)-d(atomk))/rij(atomj,atomk)
                  diff=diff+aij(znuc(atomj),znuc(atomk))*
     1                         (1.d0-diff*diff) 
                  fac1=diff
                  do 80 itr=1,niter
                     s=(1.5d0-.5d0*fac1*fac1)*fac1
                     fac1=s
   80             continue
                  fac=fac*.5d0*(1.d0-s)
   70          continue                
               do 90 atomk=atomj+1,ncen
                  diff=(d(atomj)-d(atomk))/rij(atomj,atomk)
                  diff=diff+aij(znuc(atomj),znuc(atomk))*
     1                         (1.d0-diff*diff)  
                  fac1=diff
                  do 100 itr=1,niter
                     s=(1.5d0-.5d0*fac1*fac1)*fac1
                     fac1=s
  100             continue
                  fac=fac*.5d0*(1.d0-s)
   90          continue
               norm=norm+fac
               if (atomi.eq.atomj) then
                   p(pt)=fac
               endif
   60       continue
            p(pt)=p(pt)/norm
   40    continue
         if (prnt) then
             if (atomi.le.ncen) then
                 write(iout,2) atomi
             else
                 write(iout,4)
             endif        
             write(iout,3) (p(ii),ii=1,ngrid(atomi))
         endif
         locgr=keepg
         pcnt=0
         keepwt=locwt
         call rzero(yuk,ngrid(atomi))
         locyuk=1
         do 200 i=1,nshell
            call yukawa(yuk(locyuk),grid(1,locgr+1),eta,center,
     1                  nr(i),nthet,nphi,ncen,nang,nonsep)
            call mkvwt(p(pcnt+1),wt(locwt+1),nr(i),nthet,
     1                 nphi,pout,wtout,rtyp,nang,nonsep)
            call mkyunt(rpt(1,i),wt(locwt+1),yuk(locyuk),aprvol,
     1                   yukint,nr(i),nthet,nphi,rtyp,nang,nonsep)
            locgr=locgr+nr(i)*nang
            locyuk=locyuk+nr(i)*nang
            pcnt=pcnt+pout
            locwt=locwt+wtout
  200    continue
         call iosys ('write real "atomic weights '//str//
     1               '" to lamdat',nwts(atomi),wt(keepwt+1),0,' ')
   30 continue
      yukext=0.d0
      do 400 atomi=1,ncen
         yukext=yukext+4.d+00*pi/eta(atomi)**2
  400 continue
      write(iout,5) yukext
      write(iout,6) yukint
      return
    1 format(/,5x,'scattering grid cutoff parameter ',e15.8,
     1       /,5x,'number of voronoi iterations ',i2)
    2 format(/,5x,'voronoi weights for atom ',i3)
    3 format( (/,5x,4(e15.8,1x) ) )
    4 format(/,5x,'voronoi weights for scattering center')
    5 format(/,5x,'exact     value of yukawa integral = ',e15.8)
    6 format(/,5x,'numerical value of yukawa integral = ',e15.8)
      end

