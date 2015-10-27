*deck basis.f 
c***begin prologue     basis
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           approximate inverse
c***author             schneider, b. i.(nsf)
c***source             basis
c***purpose            lobatto shape functions.
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       basis
      program basis
c
      implicit integer (a-z)
      parameter ( maxgrd=20, maxreg=50 )
      character*4096 ops
      character*1 qtyp
      character*2 itoc, idim, ireg, igr, jgr
      character*3 notim
      character*8 type
      character*56 phrse
      character*80 cpass, title, chrkey, prtit
      character*24 typfn
      character*800 card
      character*128 filbec
      logical dollar, logkey, prnply, prnwpt, prnexp
      logical check
      real*8 pi, fpkey, xl, xr, ptmp
      integer*8 px, pxply, py, pyply
      dimension ngot(2),  ngotreg(2,maxreg), npt(maxgrd)
      dimension xloc(maxgrd), pxloc(maxgrd,maxgrd)
      dimension xl(maxreg), xr(maxreg), g2use(maxreg)
      dimension py(maxreg), pyply(maxreg)
      dimension yloc(maxgrd,maxreg), pyloc(maxgrd,maxgrd,maxreg)
      dimension nptreg(maxreg)
      dimension phrse(4)
      common/io/inp, iout      
      data pi/3.1415926535897932384d0/
c
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
c
c         set dimensionality of problem and program options 
c

      nodim=intkey(ops,'number-of-space-dimensions',1,' ')
      notim=chrkey(ops,'m6290=time-dimension','no',' ')
      if(notim.ne.'no') then
         nodim=nodim+1
      endif
      type=chrkey(ops,'open-bec','unknown',' ')
      call iosys ('read character "bec filename" from rwf',
     1             -1,0,0,filbec)
      call iosys ('open bec as '//type,0,0,0,filbec)
      call iosys('read integer "no. grids" from bec',1,ngrid,0,' ')
      call iosys('read integer "no. pts per grid" from bec',
     1            ngrid,npt,0,' ')
c
c     read in the standard polynomials
c
      call rdpoly(px,pxply,xloc,pxloc,ngrid,npt,ngot)
      call iosys('write integer "no. space dimensions" to bec',1,
     1            nosdim,0,' ')
      call iosys('write character "time dimension" to bec',
     1            0,0,0,notim)
      check=logkey(ops,'m6290=check-orthogonality',.false.,' ')
      prall=logkey(ops,'print=m6290=all',.false.,' ')
      prnwpt=logkey(ops,'print=m6290=points',.false.,' ')
      prnply=logkey(ops,'print=m6290=polynomials',.false.,' ')
      prnexp=logkey(ops,'print=m6290=expansion-information',.false.,' ')
      if(prall) then
         prnwpt=.true.
         prnply=.true.
      endif
      write(iout,1)
      do 10 dim=1,nodim
         idim=itoc(dim)
         qtyp=chrkey(ops,'coordinate-'//idim,'x',' ')
         call iosys('write character "coord. type dim. =  '//ic//
     1              '" to bec',0,0,0,qtyp)
         if ( dollar('$grids-coordinate-'//qtyp,card,cpass,inp) ) then
              write(iout,2) qtyp
              fl=intkey(card,'left-boundary-condition',0,' ')
              fr=intkey(card,'right-boundary-condition',0,' ')
              call iosys('write integer "coord. type = '//qtyp//
     1                   ' left bc." to bec',1,fl,0,' ')
              call iosys('write integer "coord. type = '//qtyp//
     1                   ' right bc." to bec',1,fr,0,' ')
              nreg=intkey(card,'number-of-regions',1,' ') 
              call fparr(card,'left-endpoints',xl,nreg,' ')
              call fparr(card,'right-endpoints',xr,nreg,' ')
              call iosys('write integer "coord. type = '//qtyp
     1                   //' no. regions" to bec',1,nreg,0,' ')
              call iosys('write real "coord. type = '//qtyp
     1                   //' left endpoints" to bec',nreg,xl,0,' ')
              call iosys('write real "coord. type = '//qtyp
     1                   //' right endpoints" to bec',nreg,xr,0,' ')
              do 20 reg=1,nreg
                 ir=itoc(reg)
                 ngr=intkey(card,'number-of-grids-region-//ir,1,' ')
                 call intarr(card,'grid-list-region-//ir,g2use,ngr,' ')
c
c                compute the points, weights and polynomials
c                for this region.
c
c                go from the standard polynomials
c                to the ones required on this region.
c
                 call region(px,pxply,py(reg),pyply(reg),
     1                       xloc,pxloc,yloc(1,reg),pyloc(1,1,reg),
     2                       xl(reg),xr(reg),g2use,ngr,ngrid,maxgrd,
     3                       npt,nptreg,ngotreg(1,reg))
 20           continue   
              do 30 reg=1,nreg
c
c                drop any points and polynomials along with their
c                components which are not needed because
c                of boundary conditions and combine the polynomials at
c                the edge points to cement regions i and i+1
c
c
c                drop any points and functions required by 
c                boundary conditions if this is region 1 or region nreg
c
                 if(reg.eq.1) then
                    call drppt(pptwt(reg),fl,ngrid(reg),
     1                         npt(1,reg),nfun)
                    call drpply(ppoly(reg),pptmp,fl,ngrid(reg),
     1                          npt(1,reg),nfun,'first')
c
c                   the array npt is changed in drpply if any functions
c                   and components are dropped.
c
                 endif
                 if(reg.eq.nreg) then
                    call drpply(ppoly(reg),pptmp,fr,ngrid(reg),
     1                          npt(1,reg),nfun,'last')
                 endif
                 phrse(1)='"coord. type = '//qtyp//' region = '//ireg
                 call pt2dsk(pptwt(reg),ngrid(reg),npt(1,nreg),
     1                       phrse(1),prnwpt)
                 call ply2dsk(ppoly(reg),ngrid(reg),npt(1,nreg),
     1                        phrse(1),prnply)
 30           continue   
              call memory(-ptwds,pptmp,idum,'temp',idum)
                 call prepfn(ppoly(l,k),ppoint(k),npt(k),npt(l))
                     call iosys('write real '//title//' to bec '//
     1                          'without rewinding',npt(k)*npt(j),
     2                           z(p(k,j)),0,' ')
                     call iosys('write real '//title//' to bec '//
     1                          'without rewinding',npt(k)*npt(j),
     2                           z(dp(k,j)),0,' ')
                     call iosys('write real '//title//' to bec '//
     1                          'without rewinding',npt(k)*npt(j),
     2                           z(ddp(k,j)),0,' ')
                     call iosys('write real '//title//' to bec '//
     1                          'without rewinding',npt(k)*npt(j),
     2                           z(pn(k,j)),0,' ')
                     call iosys('write real '//title//' to bec '//
     1                          'without rewinding',npt(k)*npt(j),
     2                           z(dpn(k,j)),0,' ')
                     call iosys('write real '//title//' to bec '//
     1                          'without rewinding',npt(k)*npt(j),
     2                           z(ddpn(k,j)),0,' ')
                     if(npt(k).eq.npt(j)) then
                        if(check) then
                           call chk(z(p(k,j)),z(wt(k)),scr(sc1),
     1                              scr(sc2),npt(j),npt(k))
                           prtit='overlap matrix dvr polynomials'
                           call prntrm(prtit,scr(sc2),npt(j),npt(j),
     1                                 npt(j),npt(j),iout)
                           call chk(z(pn(k,j)),z(wtc(k)),scr(sc1),
     1                              scr(sc2),nmax(j),nmax(k))
                           prtit='modified overlap matrix dvr '//
     1                           'polynomials'
                           call prntrm(prtit,scr(sc2),nmax(j),nmax(j),
     1                                 nmax(j),nmax(j),iout)
                        endif                 
                     endif
 40               continue
 30            continue   
               call iosys('rewind '//title//' on bec read-and-write',
     1                     0,0,0,' ')
               call iosys('endfile '//title//' on bec',0,0,0,' ')
               call iosys('write integer "mod. no. points '//
     1                    'for '//phrse//'" to bec',nsubg,nmax,
     2                     0,' ')
               title='"'//qtyp(i)//' dimension pointer file"'
               titphr(2)=qtyp(i)//' pointer'
               call iosys('write character "'//titphr(2)//'" to bec',
     1                     0,0,0,title)
               call iosys('create integer '//title//' on bec',-1,
     1                     0,0,' ')
               call iosys('write integer '//title//' to bec without '//
     1                    'rewinding',nsubg,q,0,' ')
               call iosys('write integer '//title//' to bec without '//
     1                    'rewinding',nsubg,wt,0,' ')
               call iosys('write integer '//title//' to bec without '//
     1                    'rewinding',nsubg,eigc,0,' ')
               call iosys('write integer '//title//' to bec without '//
     1                    'rewinding',nsubg,wtc,0,' ')
               call iosys('write integer '//title//' to bec without '//
     1                    'rewinding',maxgrd*maxgrd,p,0,' ')
               call iosys('write integer '//title//' to bec without '//
     1                    'rewinding',maxgrd*maxgrd,dp,0,' ')
               call iosys('write integer '//title//' to bec without '//
     1                    'rewinding',maxgrd*maxgrd,ddp,0,' ')
               call iosys('write integer '//title//' to bec without '//
     1                    'rewinding',maxgrd*maxgrd,pn,0,' ')
               call iosys('write integer '//title//' to bec without '//
     1                    'rewinding',maxgrd*maxgrd,dpn,0,' ')
               call iosys('write integer '//title//' to bec without '//
     1                    'rewinding',maxgrd*maxgrd,ddpn,0,' ')
               call iosys('rewind '//title//' on bec read-and-write',
     1                     0,0,0,' ')
               call iosys('endfile '//title//' on bec',0,0,0,' ')
               if( dollar('$space-function',card,cpass,inp) ) then
                   typfn=chrkey(card,'function-type','exponential',' ')
                   call expand(z,p,pn,q,wt,eigc,wtc,phrse,titphr,
     1                         npt,nmax,typfn,maxgrd,maxd,
     2                         modfed,prnexp)
               endif 
               call memory(-ngot1,p1,idum,'space',idum)
               call memory(-ngot3,p3,idum,'scratch',idum)
               if(chebyd) then
                  call memory(-ngot4,p4,idum,'trig',idum)
               endif
            endif
 10      continue   
      endif
c
      if(notim.eq.'yes') then
         nfix(1)=.true.
         nfix(2)=.true.
         do 100 i=1,ntreg
            ic=itoc(i)
            phrse='t-'//ic
            if ( dollar('$grids-time-region-'//ic,card,cpass,inp) ) then
                 timwd=1
                 maxd=0
                 call tdat(card,qdtyp,tl,tr,ntsubg,tnpt,i,maxd)
                 call gwadd(t,twt,teigc,twtc,pt,dpt,ddpt,
     1                      ptn,dptn,ddptn,tnpt,ntsubg,maxgrd,timwd)
                 words=wpadti(timwd)
                 call memory(words,p2,ngot2,'time',0)
                 sc1=1
                 sc2=sc1+maxd*maxd
                 scrwds=wpadti(sc2+maxd*maxd)
                 call memory(scrwds,p3,ngot3,'scratch',0)
                 if(chebyd) then
                    si=1
                    ci=si+maxd*maxd
                    sj=ci+maxd*maxd
                    cj=sj+maxd*maxd
                    words = words + 4*maxd*maxd
                    need=wpadti(words)
                    call memory(need,p4,ngot4,'trig',0)
                 endif
                 call iosys('write character "q-typ for '//
     1                      phrse//'" to bec',0,0,0,qdtyp)
                 call iosys('write integer "no. subgrids for '//
     1                      phrse//'" to bec',1,ntsubg,0,' ')
                 call iosys('write real "left grid value for '//
     1                      phrse//'" to bec',1,tl,0,' ')
                 call iosys('write real "right grid value for '//
     1                      phrse//'" to bec',1,tr,0,' ')
                 call iosys('write integer "no. points for '//
     1                      phrse//'" to bec',
     2                      ntsubg,tnpt,0,' ')
                 title='"t dimension grid file-'//ic//'"'
                 titphr(1)='t title region '//ic
                 call iosys('write character "'//titphr(1)//'" to bec',
     1                       0,0,0,title)
                 call iosys('create real '//title//' on bec',-1,0,0,' ')
                 chebon=0
                 if(qdtyp.eq.'chebyshev-1'.or.
     1              qdtyp.eq.'chebyshev-2') then
                    chebon=1
                 endif
                 do 110 j=1,ntsubg         
                    jc=itoc(j)
                    if(chebon.ne.0) then
                       call chebpt(y(t(j)),y(twt(j)),tl,tr,
     1                             tnpt(j),.true.)
                    else
                       call getqpt(y(t(j)),y(twt(j)),tl,tr,
     1                             qdtyp,'before',scr(sc1),nfix,
     2                             tnpt(j),tnpt(j),1,.false.)
                    endif
                    call iosys('write real '//title//' to bec '//
     1                         'without rewinding',tnpt(j),y(t(j)),
     2                          0,' ')
                    call iosys('write real '//title//' to bec '//
     1                         'without rewinding',tnpt(j),y(twt(j)),
     2                          0,' ')
                    call droppt(y(t(j)),y(twt(j)),y(teigc(j)),
     1                          y(twtc(j)),nfix,0,1,tnpt(j),tnmax(j))
                    call iosys('write real '//title//' to bec '//
     1                         'without rewinding',tnpt(j),
     2                          y(teigc(j)),0,' ')
                    call iosys('write real '//title//' to bec '//
     1                         'without rewinding',tnpt(j),
     2                          y(twtc(j)),0,' ')
                    if(prnwpt) then
                       prtit='points for time = '//ic//' subgrid = '//jc
                       call prntrm(prtit,y(t(j)),tnpt(j),1,tnpt(j),
     1                             1,iout)
                       prtit='weights for time = '//ic//
     1                       ' subgrid = '//jc
                       call prntrm(prtit,y(twt(j)),tnpt(j),1,tnpt(j),
     1                             1,iout)
                       prtit='modified points for time = '//ic//
     1                       'subgrid = '//jc
                       call prntrm(prtit,y(teigc(j)),tnmax(j),
     1                             1,tnmax(j),1,iout)
                       prtit='modified weights for time = '//ic//
     1                       'subgrid = '//jc
                       call prntrm(prtit,y(twtc(j)),tnmax(j),1,
     1                             tnmax(j),1,iout)
                    endif 
 110             continue   
                 do 120 j=1,ntsubg
                    if(chebon.ne.0) then
                       call boxfn(trig(si),trig(ci),y(t(j)),tl,tr,
     1                            tnpt(j))
                    endif
                    do 130 k=1,ntsubg
                       if(chebon.ne.0) then            
                          call cheby(y(pt(k,j)),y(dpt(k,j)),
     1                               y(ddpt(k,j)),y(t(j)),y(twt(j)),
     2                               y(t(k)),tl,tr,trig(si),
     2                               trig(ci),trig(sj),
     3                               trig(cj),tnpt(j),tnpt(k),
     4                               tnpt(j),.true.)
                       else
                          call lgngr(y(pt(k,j)),y(dpt(k,j)),
     1                               y(ddpt(k,j)),y(t(j)),y(t(k)),
     2                               tnpt(j),tnpt(k),.false.)
                          call prepfn(y(pt(k,j)),y(dpt(k,j)),
     1                                y(ddpt(k,j)),y(twt(j)),
     2                                tnpt(j),tnpt(k))
                       endif
                       call dropfn(y(pt(k,j)),y(dpt(k,j)),y(ddpt(k,j)),
     1                             y(ptn(k,j)),y(dptn(k,j)),
     2                             y(ddptn(k,j)),nfix,0,1,
     3                             tnpt(k),tnpt(j),tnmax(k),tnmax(j))
                       if(prnply) then
                          prtit='lagrange polynomials'
                          call prntrm(prtit,y(pt(k,j)),tnpt(k),tnpt(j),
     1                                tnpt(k),tnpt(j),iout)
                          prtit='first derivative of lagrange '//
     1                          'polynomials'
                          call prntrm(prtit,y(dpt(k,j)),tnpt(k),tnpt(j),
     1                              tnpt(k),tnpt(j),iout)
                          prtit='second derivative of lagrange '//
     1                          'polynomials'
                          call prntrm(prtit,y(ddpt(k,j)),
     1                                tnpt(k),tnpt(j),
     2                                tnpt(k),tnpt(j),iout)
                          prtit='contracted lagrange polynomials'
                          call prntrm(prtit,y(ptn(k,j)),
     1                                tnmax(k),tnmax(j),
     1                                tnmax(k),tnmax(j),iout)
                          prtit='first derivative of contracted '//
     1                          'lagrange polynomials'
                          call prntrm(prtit,y(dptn(k,j)),
     1                                tnmax(k),tnmax(j),
     1                                tnmax(k),tnmax(j),iout)
                          prtit='second derivative of contracted '//
     1                          'lagrange polynomials'
                          call prntrm(prtit,y(ddptn(k,j)),
     1                                tnmax(k),tnmax(j),
     1                                tnmax(k),tnmax(j),iout)
                       endif
                       call iosys('write real '//title//' to bec '//
     1                            'without rewinding',tnpt(k)*tnpt(j),
     2                             y(pt(k,j)),0,' ')
                       call iosys('write real '//title//' to bec '//
     1                            'without rewinding',tnpt(k)*tnpt(j),
     2                             y(dpt(k,j)),0,' ')
                       call iosys('write real '//title//' to bec '//
     1                            'without rewinding',tnpt(k)*tnpt(j),
     2                             y(ddp(k,j)),0,' ')
                       call iosys('write real '//title//' to bec '//
     1                            'without rewinding',tnpt(k)*tnpt(j),
     2                             y(ptn(k,j)),0,' ')
                       call iosys('write real '//title//' to bec '//
     1                            'without rewinding',tnpt(k)*tnpt(j),
     2                             y(dptn(k,j)),0,' ')
                       call iosys('write real '//title//' to bec '//
     1                            'without rewinding',tnpt(k)*tnpt(j),
     2                             y(ddptn(k,j)),0,' ')
 130                continue
 120             continue   
                 call iosys('rewind '//title//' on bec read-and-write',
     1                       0,0,0,' ')
                 call iosys('endfile '//title//' on bec',0,0,0,' ')
                 call iosys('write integer "mod. no. points for '
     1                      //phrse//'" to bec',ntsubg,tnmax,0,' ')
                 title='"t dimension pointer file-'//ic//'"'
                 titphr(2)='t pointer region '//ic
                 call iosys('write character "'//titphr(2)//'" to bec',
     1                       0,0,0,title)
                 call iosys('create integer '//title//' on bec',-1,
     1                       0,0,' ')
                 call iosys('write integer '//title//' to bec '//
     1                      'without rewinding',ntsubg,t,0,' ')
                 call iosys('write integer '//title//' to bec '//
     1                      'without rewinding',ntsubg,twt,0,' ')
                 call iosys('write integer '//title//' to bec '//
     1                      'without rewinding',ntsubg,teigc,0,' ')
                 call iosys('write integer '//title//' to bec '//
     1                      'without rewinding',ntsubg,twtc,0,' ')
                 call iosys('write integer '//title//' to bec '//
     1                      'without rewinding',maxgrd*maxgrd,pt,0,' ')
                 call iosys('write integer '//title//' to bec '//
     1                      'without rewinding',maxgrd*maxgrd,dpt,0,' ')
                 call iosys('write integer '//title//' to bec '//
     1                      'without rewinding',maxgrd*maxgrd,ddpt,
     2                       0,' ')
                 call iosys('write integer '//title//' to bec '//
     1                      'without rewinding',maxgrd*maxgrd,ptn,
     2                       0,' ')
                 call iosys('write integer '//title//' to bec '//
     1                      'without rewinding',maxgrd*maxgrd,dptn,
     2                       0,' ')
                 call iosys('write integer '//title//' to bec '//
     1                      'without rewinding',maxgrd*maxgrd,ddptn,
     2                       0,' ')
                 call iosys('rewind '//title//' on bec read-and-write',
     1                       0,0,0,' ')
                 call iosys('endfile '//title//' on bec',0,0,0,' ')
                 if( dollar('$time-function',card,cpass,inp) ) then
                     typfn=chrkey(card,'function-type',
     1                            'exponential',' ')
                     call expand(y,pt,ptn,t,twt,teigc,twtc,
     1                           phrse,titphr,tnpt,tnmax,
     2                           typfn,maxgrd,maxd,modfed,prnexp)
                 endif 
                 call memory(-ngot2,p2,idum,'time',idum)
                 call memory(-ngot3,p3,idum,'scratch',idum)
                 if(chebyd) then
                    call memory(-ngot4,p4,idum,'trig',idum)
                 endif
            endif
 100     continue   
      endif
      call chainx(0)               
      stop
 1    format(/,20x,'grid and basis generation code')      
 2    format(/,5x,'data coordinate type = ',a24)
      end
