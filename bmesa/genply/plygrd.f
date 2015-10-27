*deck plygrd.f 
c***begin prologue     plygrd
c***date written       970810   (yymmdd)
c***revision date               (yymmdd)
c***keywords           plygrd
c***                   
c***author             schneider, b. i.(nsf)
c***source             plygrd
c***purpose            polynomial generation on multi-grids
c
c***description         
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       plygrd
      program plygrd
c
      implicit integer (a-z)
      parameter ( ngmx=20 )
      character*4096 ops
      character*8 prtflg
      character*2 itoc
      character*80 cpass, title, chrkey
      character*1600 card
      character*8 qtyp
      character*32 qdtyp
      logical dollar, logkey, prnply, prnwpt 
      logical prall, check, fix, nfix
      logical toau, useau
      real*8 z, y, fpkey, xl, xr, sfac
      real*8 pi
      real*8 massau, lenau, timau
      dimension npt(ngmx), nmax(ngmx), dim(ngmx)
      dimension q(ngmx), wt(ngmx), eigc(ngmx), wtc(ngmx)
      dimension fix(ngmx), nfix(2,ngmx)
      dimension qtyp(ngmx), qdtyp(ngmx), sfac(ngmx)
      dimension xl(ngmx), xr(ngmx), fl(ngmx), fr(ngmx)
      dimension dfl(ngmx), dfr(ngmx)
      common/io/inp, iout      
      pointer (p1,z),(p1,iz)
      pointer (p2,y)
c
c          mass in Kilograms length in meters
      data massau, lenau, timau / 9.109558d-31, 5.291771d-11, 
     1                            2.418884d-17 /
      call drum
      write(iout,*)
      write(iout,*) '    generate polynomials       '
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
c
c                       set print options
c
      prall=logkey(ops,'print=m7030=all',.false.,' ')
      if(prall) then
         prnply=.true.
         prnwpt=.true. 
      else         
         prnply=logkey(ops,'print=m7030=polynomials',.false.,' ')
         prnwpt=logkey(ops,'print=m7030=points/weights',
     1                    .false.,' ')
      endif
      check=logkey(ops,'m7030=check-orthogonality',.false.,' ')
      toau=logkey(ops,'to-atomic-units',.false.,' ')
      useau=logkey(ops,'use-atomic-units',.false.,' ')
      write(iout,1) prnply, prnwpt, check
c
c              set various program options
c
c               spatial basis set information
c       
      nrun=intkey(ops,'number-of-grids',1,' ')
      call iosys('write integer "number of grids" to rwf',1,nrun,0,' ')
      bigdim=0
      ioff=1
      do 1000 i=1,nrun
         if(dollar('$basis-'//itoc(i),card,cpass,inp)) then
            call bdat(.false.,card,qtyp(i),qdtyp(i),xl(i),xr(i),
     1                fl(i),fr(i),dfl(i),dfr(i),npt(i),
     2                nmax(i),nfix(1,i),1)
            title=chrkey(card,'run-title','run-number-'//itoc(i),' ')
            call iosys ('write character "coordinate type grid = '
     1                   //itoc(i)//'" to rwf',0,0,0,qtyp(i))
            call iosys('write real "left end point grid = '
     1                 //itoc(i)//'" to rwf',1,xl(i),0,' ')
            call iosys('write real "right end point grid = '
     1                 //itoc(i)//'" to rwf',1,xr(i),0,' ')
            call iosys('write integer "left derivative grid = '//
     1                  itoc(i)//'" to rwf',1,dfl(i),0,' ')
            call iosys('write integer "right derivative grid = '//
     1                  itoc(i)//'" to rwf',1,dfr(i),0,' ')
            call iosys('write integer "left boundary condition '//
     1                 'grid = '//itoc(i)//'" to rwf',1,fl(i),0,' ')
            call iosys('write integer "right boundary condition '//
     1                 'grid = '//itoc(i)//'" to rwf',1,fr(i),0,' ')
            call iosys('write integer "number of points '//
     1                 'grid = '//itoc(i)//'" to rwf',1,npt(i),
     2                  0,' ')
            call iosys('write integer "number of functions '//
     1                 'grid = '//itoc(i)//'" to rwf',1,nmax(i),
     2                  0,' ')
         endif            
         if(toau) then
c
c           converts everything to atomic units
c
            write(iout,*) 'converting to atomic units'
            xl(i)=xl(i)/lenau
            xr(i)=xr(i)/lenau
         endif
         dim(i)=max(nmax(i),npt(i))
         bigdim=max(bigdim,dim(i))
         q(i)=ioff
         wt(i)=q(i)+dim(i)
         eigc(i)=wt(i)+dim(i)
         wtc(i)=eigc(i)+dim(i)
         words=wtc(i)+dim(i)
         ioff=words
 1000 continue    
      need=wpadti(words)     
      call memory(need,p1,ngot1,'grid',0)
      sc=1
      need=wpadti(sc+bigdim*bigdim)
      call memory(need,p2,ngot2,'scratch',0)
      do 2000 i=1,nrun              
         write(iout,*) q(i), wt(i), xl(i), xr(i), qdtyp(i), sc, 
     1                 nfix(1,1), nfix(2,1), npt(i), bigdim
         call getqpt(z(q(i)),z(wt(i)),xl(i),xr(i),qdtyp(i),'before',
     1               y(sc),nfix(1,i),npt(i),npt(i),1,.true.)
         title='initial points for grid = '//itoc(i)
         call prntrm(title,z(q(i)),npt(i),1,npt(i),1,iout)
         title='initial weights for grid = '//itoc(i)
         call prntrm(title,z(wt(i)),npt(i),1,npt(i),1,iout)
         call droppt(z(q(i)),z(wt(i)),z(eigc(i)),z(wtc(i)),
     1               npt(i),nmax(i))
         if(prnwpt) then
            title='final points for grid = '//itoc(i)
            call prntrm(title,z(eigc(i)),nmax(i),1,nmax(i),
     1                  1,iout)
            title='final weights for grid = '//itoc(i)
            call prntrm(title,z(wtc(i)),nmax(i),1,nmax(i),
     1                  1,iout)
         endif
         call iosys('write real "points grid = '//itoc(i)
     1              //'" to rwf',nmax(i),z(eigc(i)),0,' ')
         call iosys('write real "weights grid = '//itoc(i)
     1              //'" to rwf',nmax(i),z(wtc(i)),0,' ')     
 2000 continue   
      call memory(-ngot2,p2,idum,'scratch',idum)
      call lnkerr('quit')
c
c           now generate the functions and derivatives on all the grids.
c
      p=1
      dp=p+bigdim*bigdim
      ddp=dp+bigdim*bigdim
      pn=ddp+bigdim*bigdim
      dpn=pn+bigdim*bigdim
      ddpn=dpn+bigdim*bigdim
      sc1=ddpn+bigdim*bigdim
      sc2=work1+bigdim*bigdim
      need=wpadti(sc2+bigdim*bigdim)
      call memory(need,p2,ngot2,'grid',0)
      do 3000 i=1,nrun
         write(iout,4) i
         do 4000 j=1,nrun
            write(iout,5) j
c
c           calculate the set of functions labeled by i which are
c           the pivot points at the grid labeled by j
c 
c            call lgngr(y(p),y(dp),y(ddp),z(q(i)),z(q(j)),
c     1                 npt(i),npt(j),.false.)
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
c
c           drop required functions and points
c
            call dropfn(y(p),y(dp),y(ddp),y(pn),y(dpn),y(ddpn),
     1                  nfix(1,i),fl(i),fr(i),npt(j),npt(i),
     2                  nmax(j),nmax(i))
            if(prnply) then
               title='contracted lagrange polynomials'
               call prntrm(title,y(p),nmax(j),nmax(i),
     1                     nmax(j),nmax(i),iout)
               title='first derivative of contracted lagrange '//
     1               'polynomials'
               call prntrm(title,y(dp),nmax(j),nmax(i),
     1                     nmax(j),nmax(i),iout)
               title='second derivative of contracted lagrange '//
     1               'polynomials'
               call prntrm(title,y(ddp),nmax(j),nmax(i),
     1                     nmax(j),nmax(i),iout)
            endif     
c            if(i.eq.j) then
               if(check) then
                  call chk(y(p),z(wtc(j)),y(work1),y(work2),
     1                  nmax(i),nmax(j))
                  title='overlap matrix dvr polynomials'
                  call prntrm(title,y(work2),nmax(i),nmax(i),
     1                        nmax(i),nmax(i),iout)
c               endif
            endif
            call iosys('write real "lgrbf = '//itoc(i)//
     1                 ' grid = '//itoc(j)//'" to rwf',
     2                  nmax(i)*nmax(j),y(p),0,' ')
            call iosys('write real "dlgrbf = '//itoc(i)//
     1                 ' grid = '//itoc(j)//'" to rwf',
     2                  nmax(i)*nmax(j),y(dp),0,' ')
            call iosys('write real "ddlgrbf = '//itoc(i)//
     1                 ' grid = '//itoc(j)//'" to rwf',
     2                  nmax(i)*nmax(j),y(ddp),0,' ')
 4000    continue
 3000 continue   
      call memory(-ngot2,p2,idum,'grid',idum)
      call chainx(0)               
      stop
 1    format(/,15x,'code options',/,
     1       /,5x,'print q1 polynomials                   = ',l1,
     2       /,5x,'print q1 points/weights                = ',l1,
     3       /,5x,'check q1 orthonormality                = ',l1)
 2    format(/,15x,a60,/,/,5x,     
     1             'order of polynomials           = ',i3,/,5x,
     2             'number of integration points   = ',i3,/,5x,
     3             'left boundary                  = ',e15.8,/,5x,
     4             'right boundary                 = ',e15.8,/,5x,
     5             'derivative at left boundary    = ',e15.8,/,5x,
     6             'derivative at right boundary   = ',e15.8,/,5x,
     7             'left-boundary-condition        = ',i2,/,5x,
     8             'right-boundary-condition       = ',i2)
 3    format(/,5x,'fixed point   = ',e15.8)
 4    format(/,10x,'generating data for polynomial set = ',i4)
 5    format(/,20x,'grid = ',i4)
 6    format(/,20x,'using tridiagonal form')
      end






