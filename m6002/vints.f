*deck @(#)vints.f	1.1 9/7/91
c***begin prologue     vints
c***date written       xxxxxx   (yymmdd)
c***revision date      890409   (yymmdd)
c***keywords           m6002, link 6003, v integral
c***authors            rescigno, tom (llnl), mccurdy, bill (ohio state)
c***                   schneider, barry (lanl)
c***source             m6002
c***purpose            numerical direct potential integrals
c***description        calculates static potential on a numerical grid
c***                   using technology of nuclear attraction integrals
c***                   thus we calculate the integral of 1/abs(r-r') 
c***                   times two primitive basis functions with r 
c***                   treated as a nuclear center.
c***references       
c
c***routines called    
c***end prologue       vints
      subroutine vints(valint,rhocn,roexst,grid,g,pc,fvec,facr,arg,
     1                 pointr,reschg,diag,sumpot,ncen,npnts,nsts,
     2                 nstri,ncntri,nptmx)
      implicit integer (a-z)
      parameter (dimpr=300 , dimcen=10)
      real *8 valint, rhocn, grid, g, pc, fvec, facr, arg, alf, cont
      real *8 charge, rloc, pi, piterm, pitern, acrcy, scale, pa
      real *8 t, t1, p1, p2, p3, exfac, pre, pcsq, nuccat, dist, sumpot
      real *8 reschg
      character *(*) diag, roexst
      dimension valint(nptmx,nstri), g(nptmx,7,3), pa(2,3)
      dimension grid(4,nptmx), pc(nptmx,3), fvec(nptmx,7)
      dimension rhocn(ncntri,nstri), facr(nstri), arg(nptmx)
      dimension diag(nstri), sumpot(nstri), roexst(nstri)
      dimension pointr(nptmx,2)
      common /aosi/ npr, ncon, nxyzc(dimpr,4), nprc(dimpr), 
     1              sncon, smllst(dimpr), junk(dimpr)
      common /aosr/ alf(dimpr), cont(dimpr)
      common /rloc/ charge(dimcen), rloc(3,dimcen)
      common /nmbrs /  pi, piterm, pitern, acrcy, scale, icanon
      common /io/ inp, iout
c----------------------------------------------------------------------c
c               loop over contracted functions                         c
c----------------------------------------------------------------------c
      frsto=0
      lsto=0
      do 10 cono=1,ncon
         frsto=frsto+1
         lsto=lsto+nprc(cono)
         if (smllst(cono).ne.0) then
             frsti=0
             lsti=0
             do 20 coni=1,cono
                cntcon=cono*(cono-1)/2+coni
                frsti=frsti+1
                lsti=lsti+nprc(coni)
                if (smllst(coni).ne.0) then
                    do 30 prmo=frsto,lsto
                       ceno=nxyzc(prmo,4)
                        do 40 prmi=frsti,lsti
                          ceni=nxyzc(prmi,4)
                          mx=nxyzc(prmo,1)+nxyzc(prmi,1)+1
                          my=nxyzc(prmo,2)+nxyzc(prmi,2)+1
                          mz=nxyzc(prmo,3)+nxyzc(prmi,3)+1
                          t1=alf(prmo)+alf(prmi)
                          t=1.d0/t1
                          p1=(alf(prmo)*rloc(1,ceno)
     1                       +alf(prmi)*rloc(1,ceni))*t
                          p2=(alf(prmo)*rloc(2,ceno)
     1                       +alf(prmi)*rloc(2,ceni))*t
                          p3=(alf(prmo)*rloc(3,ceno)
     1                       +alf(prmi)*rloc(3,ceni))*t
                          exfac=(rloc(1,ceno)-rloc(1,ceni))*
     1                          (rloc(1,ceno)-rloc(1,ceni))+
     2                          (rloc(2,ceno)-rloc(2,ceni))*
     3                          (rloc(2,ceno)-rloc(2,ceni))+
     4                          (rloc(3,ceno)-rloc(3,ceni))*
     5                          (rloc(3,ceno)-rloc(3,ceni))
                          exfac=exfac*t*alf(prmo)*alf(prmi)
                          exfac=exp(-exfac)
                          pre=2.d+00*pi*t*cont(prmo)*cont(prmi)*exfac
                          do 50 ist=1,nstri
                             facr(ist)=rhocn(cntcon,ist)*pre
   50                     continue
c----------------------------------------------------------------------c
c                 now loop over grid points                            c
c----------------------------------------------------------------------c 
                          do 60 grpt=1,npnts
                             pc(grpt,1)=p1-grid(1,grpt)
                             pc(grpt,2)=p2-grid(2,grpt)
                             pc(grpt,3)=p3-grid(3,grpt)
                             pcsq=pc(grpt,1)*pc(grpt,1)
     1                            +pc(grpt,2)*pc(grpt,2)
     2                            +pc(grpt,3)*pc(grpt,3)
                             arg(grpt)=t1*pcsq
   60                     continue
                          call izero(pointr,2*npnts)
                          nlt0=0
                          nge0=0
                          do 500 grpt=1,npnts
                             switch=arg(grpt)-34.9d+00
                             if (switch.lt.0.d+00) then
                                 nlt0=nlt0+1
                                 pointr(nlt0,1)=grpt
                             else
                                 nge0=nge0+1
                                 pointr(nge0,2)=grpt
                             endif
  500                     continue
c----------------------------------------------------------------------c
c                compute the f functions                               c
c----------------------------------------------------------------------c
                          maxtyp=mx+my+mz-2
                          if(maxtyp.eq.1) then
                             call stuff0(npnts,nlt0,nge0,arg,
     1                                   pointr,fvec)
                          elseif (maxtyp.eq.2) then
                             call stuff1(npnts,nlt0,nge0,arg,
     1                                   pointr,fvec)
                          elseif (maxtyp.eq.3) then
                             call stuff2(npnts,nlt0,nge0,arg,
     1                                   pointr,fvec)
                          elseif (maxtyp.eq.4) then
                             call stuff3(npnts,nlt0,nge0,arg,
     1                                   pointr,fvec)
                          elseif (maxtyp.eq.5) then
                             call stuff4(npnts,nlt0,nge0,arg,
     1                                   pointr,fvec)
                          elseif (maxtyp.eq.6) then
                             call stuff5(npnts,nlt0,nge0,arg,
     1                                   pointr,fvec)
                          elseif (maxtyp.eq.7) then
                             call stuff6(npnts,nlt0,nge0,arg,
     1                                   pointr,fvec)
                          else
                             call lnkerr('error in call to stuff'//
     1                                   'routine')
                          endif
                          pa(1,1)=p1-rloc(1,ceno)
                          pa(2,1)=p1-rloc(1,ceni)
                          pa(1,2)=p2-rloc(2,ceno)
                          pa(2,2)=p2-rloc(2,ceni)
                          pa(1,3)=p3-rloc(3,ceno)
                          pa(2,3)=p3-rloc(3,ceni)
                          call  gfunct(nxyzc(prmo,1),nxyzc(prmi,1),
     1                                 pa(1,1),pa(2,1),pc(1,1),t,g,1,
     2                                 npnts)
                          call  gfunct(nxyzc(prmo,2),nxyzc(prmi,2),
     1                                 pa(1,2),pa(2,2),pc(1,2),t,g,2,
     2                                 npnts)
                          call  gfunct(nxyzc(prmo,3),nxyzc(prmi,3),
     1                                 pa(1,3),pa(2,3),pc(1,3),t,g,3,
     2                                 npnts)
                          do 70 ix=1,mx
                             do 80 jy=1,my
                                do 90 kz=1,mz
                                   mxyz=ix+jy+kz-2
                                   do 100 ist=1,nstri
                                      if (roexst(ist).eq.'yes') then
                                          do 200 grpt=1,npnts
                                             valint(grpt,ist)=
     1                                              valint(grpt,ist)+
     2                                                   g(grpt,ix,1)*
     3                                                   g(grpt,jy,2)*
     4                                                   g(grpt,kz,3)*
     5                                                fvec(grpt,mxyz)*
     6                                                facr(ist)
  200                                     continue
                                      endif  
  100                              continue
   90                           continue
   80                        continue
   70                     continue
   40                  continue
   30               continue
                endif
                frsti=lsti 
   20        continue
         endif
         frsto=lsto
   10 continue
      ist=0
      do 300 is=1,nsts
         do 305 js=1,is
            ist=ist+1
            if (roexst(ist).eq.'yes') then
                if (is.eq.js) then
c----------------------------------------------------------------------c
c                   add in nuclear contribution                        c
c----------------------------------------------------------------------c
                    do 310 icen=1,ncen
                       do 320 grpt=1,npnts
                          dist=(rloc(1,icen)-grid(1,grpt))**2 +
     1                         (rloc(2,icen)-grid(2,grpt))**2 +
     2                         (rloc(3,icen)-grid(3,grpt))**2
                          nuccat= charge(icen)/sqrt(dist)
                          valint(grpt,ist)=valint(grpt,ist) - nuccat
  320                  continue
  310               continue
                    if (reschg.ne.0.d+00) then
                        do 330 grpt=1,npnts      
                           dist=grid(1,grpt)**2 +
     1                               grid(2,grpt)**2 +
     2                                    grid(3,grpt)**2
                           nuccat= reschg/sqrt(dist)
                           valint(grpt,ist)=valint(grpt,ist) + nuccat
  330                   continue
                    endif
                endif
            endif
c----------------------------------------------------------------------c
c                calculate test integral                               c
c----------------------------------------------------------------------c
            do 350 grpt=1,npnts
               sumpot(ist)=sumpot(ist)+valint(grpt,ist)*grid(4,grpt)
  350       continue
  305    continue
  300 continue
      return
      end 
