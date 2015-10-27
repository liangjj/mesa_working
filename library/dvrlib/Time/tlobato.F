*deck tlobato.f 
c***begin prologue     tlobato
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           lobatto functions
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       tlobato
      subroutine tlobato(pham,edge,qtyp,n,m,prn)
c
      implicit integer (a-z)
      character*(*) qtyp
      real*8 ham, scr
      real*8 edge 
#ifdef DECPOINTER
      integer*8 pham
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 pham
#endif SGIPOINTER
      logical prn
      dimension edge(2), prn(*)
      common/io/inp, iout      
      pointer (pscr,scr(1))
      pointer (pham,ham(1))
      m=n-1
      hin=1
      tin=hin+n*n
      twtin=tin+n
      ptin=twtin+n
      dptin=ptin+n*n
      ddptin=dptin+n*n
      need=wpadti(ddptin+n*n)
      call getmem(need,pscr,njunk,'scr',0)
      hout=1
      t0out=hout+m*m
      tout=t0out+1
      twtout=tout+m
      ptout=twtout+m
      dptout=ptout+m*m
      ddptout=dptout+m*m
      need=wpadti(ddptout+m*m)
      call getmem(need,pham,nword,'ham',0)
c
c     calculate the functions and their derivatives.
c
      call twpt(scr(tin),scr(twtin),scr(ptin),scr(dptin),scr(ddptin),
     1          qtyp,edge,n,n,prn(1))
      call nrmlze(scr(ptin),scr(dptin),scr(twtin),n)
      call tplmat(scr(ptin),scr(dptin),scr(twtin),scr(hin),n)
      call tphy(scr(tin),scr(twtin),scr(ptin),scr(dptin),scr(ddptin),
     1          scr(hin),ham(t0out),ham(tout),ham(twtout),
     2          ham(ptout),ham(dptout),ham(ddptout),ham(hout),
     3          n,m,prn(2))
      return
      end
















