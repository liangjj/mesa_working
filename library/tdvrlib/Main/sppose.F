*deck sppose.f
c***begin prologue     sppose
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            calculate zero time gaussian wavepacket 
c***                   
c***references         
c
c***routines called    
c***end prologue       sppose
      subroutine sppose(u01,u02,u03,psi0,n,nc,nd,dim,coord)
      implicit integer (a-z)
      real*8 u01, u02, u03
      real*8 psi0
      real*8 x, y, z
      character*(*) coord
      character*80 title
      character*320 card
      logical dollar
      dimension nd(3)
      dimension u01(nd(1),*), u02(nd(2),*), u03(nd(3),*)
      dimension psi0(n,nc,2)
      common/io/inp, iout
      pointer (px,x(1)), (px,ix(1))
      pointer (py,y(1)), (py,iy(1))
      pointer (pz,z(1)), (pz,iz(1))
      call rzero(psi0,n*nc*2)
      if ( dollar('$states',card,title,inp) ) then
         if(coord.eq.'cartesian') then
            ntot=1
            if(dim.ge.1) then
               nxs=intkey(card,'number-of-x-superposed-states',1,' ')
               ntot=ntot*nxs
               cx=1
               phix=cx+nxs
               lx=wpadti(phix+nd(1))
               nwx=lx+nxs
               call getmem(nwx,px,ngx,'x',0)
               call intarr(card,'x-state-list',ix(lx),nxs,' ')
               call fparr(card,'x-state-coefficients',x(cx),nxs,' ')
               call mkphi(x(phix),u01,x(cx),ix(lx),nxs,nd(1),'x')
            endif
            if(dim.ge.2) then
               nys=intkey(card,'number-of-y-superposed-states',1,' ')
               ntot=ntot*nys
               cy=1
               phiy=cy+nys
               ly=wpadti(phiy+nd(2))
               nwy=ly+nys
               call getmem(nwy,py,ngy,'y',0)
               call intarr(card,'y-state-list',iy(ly),nys,' ')
               call fparr(card,'y-state-coefficients',y(cy),nys,' ')
               call mkphi(y(phiy),u02,y(cy),iy(ly),nys,nd(2),'y')
            endif
            if(dim.eq.3) then
               nzs=intkey(card,'number-of-z-superposed-states',1,' ')
               ntot=ntot*nzs
               cz=1
               phiz=cz+nzs
               lz=wpadti(phiz+nd(3))
               nwz=lz+nzs
               call getmem(nwz,pz,ngz,'z',0)
               call intarr(card,'z-state-list',iz(lz),nzs,' ')
               call fparr(card,'z-state-coefficients',z(cz),nzs,' ')
               call mkphi(z(phiz),u03,z(cz),iz(lz),nzs,nd(3),'z')
            endif
         else
            nrs=intkey(card,'number-of-r-superposed-states',1,' ')
            ntot=nrs
            cr=1
            phir=cr+nrs
            phix=phir
            lr=wpadti(phir+nd(1))
            nwr=lr+nrs            
            nwx=nwr 
            call getmem(nwx,px,ngx,'r',0)
            call intarr(card,'r-state-list',ix(lr),nrs,' ')
            call fparr(card,'r-state-coefficients',x(nrs),nrs,' ')
            call mkphi(x(phir),u01,x(cr),ix(lr),nrs,nd(1),'r')
         endif
      endif
      if(dim.eq.1) then
         call fil1(psi0,x(phix),nd(1),nc)
         call getmem(-nwx,px,idum,'x',idum)   
      endif
      if(dim.eq.2) then
         call fil2(psi0,x(phix),y(phiy),nd(1),nd(2),nc)
         call getmem(-nwx,px,idum,'x',idum)   
         call getmem(-nwy,py,idum,'y',idum)   
      endif
      if(dim.eq.3) then
         call fil3(psi0,x(phix),y(phiy),z(phiz),nd(1),nd(2),nd(3),nc)
         call getmem(-nwx,px,idum,'x',idum)   
         call getmem(-nwy,py,idum,'y',idum)   
         call getmem(-nwz,pz,idum,'z',idum)   
      endif
      return
      end       
