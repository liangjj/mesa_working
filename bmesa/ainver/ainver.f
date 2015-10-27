*deck ainver.f 
c***begin prologue     ainver
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           approximate inverse
c***author             schneider, b. i.(nsf)
c***source             ainver
c***purpose            approximate inverse of hamiltonian
c***                   using dvr.
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       ainver
      program ainver
c
      implicit integer (a-z)
      parameter ( dimmax=20) 
      character*4096 ops
      character*1 itoc, ic
      character*8 type
      character*80 cpass, title, chrkey
      character*24 coord, comp, qtyp, typfn
      character*32 qdtyp, drctv
      character*800 card
      character*128 filbec
      logical dollar, logkey, prnply, prnwpt
      logical check, nfix
      logical toau, test, useau, dospl
      real*8 z, y, scr, trig
      real*8 fpkey, xl, xr  
      real*8 pi
      real*8 dleft, dright
      dimension nmax(dimmax), npt(dimmax)
      dimension nfix(2,dimmax), xl(dimmax), xr(dimmax)
      dimension fl(dimmax), fr(dimmax), dfl(dimmax), dfr(dimmax)
      dimension qtyp(dimmax), qdtyp(dimmax)
      dimension q(dimmax), wt(dimmax)
      dimension eigc(dimmax), wtc(dimmax), iloc(dimmax,dimmax)
      common/io/inp, iout      
      pointer (p1,z(1)), (p1,iz(1))
      pointer (p2,scr(1))
      pointer (p3,y(1)), (p3,iy(1))
      pointer (p4,trig(1))
      data pi/3.1415926535897932384d0/
c
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
c
c         set spatial dimensionality of problem and program options 
c
      ngrid=intkey(ops,'number-of-grids',1,' ')
      coord=chrkey(ops,'coordinate-system','cartesian',' ')
      type=chrkey(ops,'open-bec','unknown',' ')
      check=logkey(ops,'m6290=check-orthogonality',.false.,' ')
      prall=logkey(ops,'print=m6290=all',.false.,' ')
      prnwpt=logkey(ops,'print=m6290=points',.false.,' ')
      prnply=logkey(ops,'print=m6290=polynomials',.false.,' ')
      if(prall) then
         prnwpt=.true.
         prnply=.true.
      endif
      write(iout,1) ngrid
c
      call iosys ('read character "bec filename" from rwf',
     1             -1,0,0,filbec)
      call iosys ('open bec as '//type,0,0,0,filbec)
c
c               generate needed grids
c
      ioff=1
      maxd=0
      do 1000 i=1,ngrid
         ic=itoc(i)
         if ( dollar('$grid-'//ic,card,cpass,inp) ) then
              call gdat(card,qtyp(i),qdtyp(i),xl(i),xr(i),
     1                  fl(i),fr(i),dfl(i),dfr(i),npt(i),
     2                  nmax(i),nfix(1,i),i)
              call gwadd(q(i),wt(i),eigc(i),wtc(i),npt(i),ioff)
         endif
         maxd=max(maxd,npt(i))
 1000 continue
      need=wpadti(ioff)
      call memory(need,p1,ngot1,'grid',0)      
      sc1=1
      sc2=sc1+maxd*maxd
      scrwds=wpadti(sc2+maxd*maxd)
      call memory(scrwds,p2,ngot2,'scratch',0)
      chebon=0
      do 2000 i=1,ngrid         
         if(qdtyp(i).eq.'chebyshev-1'.or.qdtyp(i).eq.'chebyshev-2') then
            call chebpt(z(q(i)),z(wt(i)),xl(i),xr(i),npt(i),.true.)
            chebon=1
         else
            call getqpt(z(q(i)),z(wt(i)),xl(i),xr(i),
     1                  qdtyp(i),'before',scr(sc1),nfix(1,i),
     2                  npt(i),npt(i),1,.false.)
         endif
         if(prnwpt) then
            title='points for grid = '//itoc(i)
            call prntrm(title,z(q(i)),npt(i),1,npt(i),1,iout)
            title='weights for grid = '//itoc(i)
            call prntrm(title,z(wt(i)),npt(i),1,npt(i),1,iout)
         endif
         call droppt(z(q(i)),z(wt(i)),z(eigc(i)),z(wtc(i)),
     1               nfix(1,i),fl(i),fr(i),npt(i),nmax(i),
     2               nmax(i),iout)
         if(prnwpt) then         
            title='modified points for grid = '//itoc(i)
            call prntrm(title,z(eigc(i)),nmax(i),1,nmax(i),1,iout)
            title='modified weights for grid = '//itoc(i)
            call prntrm(title,z(wtc(i)),nmax(i),1,nmax(i),1,iout)
         endif 
 2000 continue
c
c     calculate DVR functions on all required grids
c
      cntmem=1
      do 3000 i=1,ngrid
         do 4000 j=1,ngrid 
            cntmem=cntmem+6*npt(j)*npt(j)
 4000    continue   
 3000 continue   
      words=cntmem 
      need=wpadti(words)
      call memory(need,p3,ngot3,'dvr-fun',0)
      if(chebon.ne.0) then
         si=1
         ci=si+maxd*maxd
         sj=ci+maxd*maxd
         cj=sj+maxd*maxd
         words = words + 4*maxd*maxd
         need=wpadti(words)
         call memory(need,p4,ngot4,'trig',0)
      endif
      begin=1
      do 5000 i=1,ngrid
         if(qdtyp(i).eq.'chebyshev-1'.or.qdtyp(i).eq.'chebyshev-2') then
            call boxfn(trig(si),trig(ci),z(q(i)),xl(i),xr(i),npt(i))
         endif       
         do 6000 j=1,ngrid
            write(iout,2) i, j
            iloc(j,i)=begin
            pn=iloc(j,i)
            dpn=pn+nmax(i)*nmax(j)
            ddpn=dpn+nmax(i)*nmax(j)
            p=ddpn+nmax(i)*nmax(j)
            dp=p+npt(i)*npt(j)
            ddp=dp+npt(i)*npt(j)
            if(qdtyp(j).eq.'chebyshev-1'.or.
     1         qdtyp(j).eq.'chebyshev-2') then            
               call cheby(y(p),y(dp),y(ddp),z(q(i)),z(wt(i)),z(q(j)),
     1                    xl(j),xr(j),trig(si),trig(ci),trig(sj),
     2                    trig(cj),npt(i),npt(j),npt(i),.true.)
            else
               call lgngr(y(p),y(dp),y(ddp),z(q(i)),z(q(j)),
     1                    npt(i),npt(j),.false.)
               call prepfn(y(p),y(dp),y(ddp),z(wt(i)),npt(i),npt(j))
            endif
            if(prnply) then
               title='lagrange polynomials'
               call prntrm(title,y(p),npt(j),npt(i),
     1                     npt(j),npt(i),iout)
               title='first derivative of lagrange '//
     1               'polynomials'
               call prntrm(title,y(dp),npt(j),npt(i),
     1                     npt(j),npt(i),iout)
               title='second derivative of lagrange '//
     1               'polynomials'
               call prntrm(title,y(ddp),npt(j),npt(i),
     1                     npt(j),npt(i),iout)
            endif
            if(npt(j).ge.npt(i)) then
               if(check) then
                  call chk(y(p),z(wt(j)),scr(sc1),scr(sc2),
     1                     npt(i),npt(j))
                  title='overlap matrix dvr polynomials'
                  call prntrm(title,scr(sc2),npt(i),npt(i),
     1                        npt(i),npt(i),iout)
               endif                 
            endif
            call dropfn(y(p),y(dp),y(ddp),y(pn),y(dpn),y(ddpn),
     1                  nfix(1,i),fl(i),fr(i),
     2                  npt(j),npt(i),nmax(j),nmax(i))
            if(prnply) then
               title='contracted lagrange polynomials'
               call prntrm(title,y(pn),nmax(j),nmax(i),
     1                     nmax(j),nmax(i),iout)
               title='first derivative of contracted lagrange '//
     1               'polynomials'
               call prntrm(title,y(dpn),nmax(j),nmax(i),
     1                     nmax(j),nmax(i),iout)
               title='second derivative of contracted lagrange '//
     1               'polynomials'
               call prntrm(title,y(ddpn),nmax(j),nmax(i),
     1                     nmax(j),nmax(i),iout)
            endif
            if(nmax(j).ge.nmax(i)) then      
               if(check) then
                  call chk(y(pn),z(wtc(j)),scr(sc1),scr(sc2),
     1                     nmax(i),nmax(j))
                  title='modified overlap matrix dvr polynomials'
                  call prntrm(title,scr(sc2),nmax(i),nmax(i),
     1                        nmax(i),nmax(i),iout)
               endif
            endif
            begin=ddp+npt(i)*npt(j)
 6000    continue
 5000 continue   
      if( dollar('$function',card,cpass,inp) ) then
          typfn=chrkey(card,'function-type','exponential',' ')
      endif 
      call expand(y(1),z(1),iloc,typfn,ngrid,eigc,wtc,nmax,
     1            maxd,dimmax)
      call memory(-ngot1,p1,idum,'grid',idum)
      call memory(-ngot2,p2,idum,'scratch',idum)
      call memory(-ngot3,p3,idum,'dvr-fun',idum)
      if(chebon.ne.0) then
         call memory(-ngot4,p4,idum,'trig',idum)
      endif
      call chainx(0)               
      stop
 1    format(/,20x,'grid and basis generation code',//,20x,
     1             'number of grids = ',i2)      
 2    format(/,5x,'computing basis set = ',i2,' on grid = ',i2)
 3    format(/,15x,'code options for spatial dimension = ',i1,/,5x,
     1             'print polynomials                  = ',l1,/,5x,
     2             'print points/weights               = ',l1,/,5x,
     3             'print dvr information              = ',l1,/,5x, 
     4             'print hamiltonian information      = ',l1,/,5x, 
     5             'check orthonormality               = ',l1)
      end
