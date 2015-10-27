*deck rmtxp.f 
c***begin prologue     rmtxp
c***date written       941019   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6236, link 6236, propagation
c***author             schneider, b. i.(nsf)
c***source             m6236
c***purpose            driver for r-matrix propagation:inelastic scattering
c***description        propagation of a set of coupled second order
c***                   differential equations is accomplished using 
c***                   the r-matrix propagation approach of light-walker.
c***references         papers by light-walker j. chem. phys.
c
c***routines called    iosys, util and mdutil
c***end prologue       m6236
c                           vmat and tsave can reside in the same core
c                           tp needs to have its own space
c                           r4 needs to have its own space
      subroutine rmtxp (vmat,vcij,range,tev,ec,enrg,ipvt,tp,tsave,pdiag,
     1                  r4,space,elxxx,n,itri,nen,prntr,rstart,
     2                  rfinal,drstrt,stpmin,stpmax,beta,
     3                  prntv,npv,prnt4,npr,debug,type)
      implicit real*8 (a-h,o-z)
      common /io/ inp, iout
c
      dimension vmat(n,n), space(n,n,3), enrg(n), ec(n)
      dimension ipvt(n), tp(n,n), tsave(n,n)
      dimension pdiag(n,4), r4(n,n), rmtx1q(5), nmtx1q(2)
      dimension tev(n,n), vcij(itri)
      logical last, okr, okw, prntr, debug
      logical prntv, prnt4, okpv, okpr
      character*80 title
      character*16 fptoc
      character*5 itoc
      character*(*) type
      equivalence(rmtx1q(5),nmtx1q(1))
c
c     begin execution
c
      ienc=nen
c
c        correct the r-matrix - place in r4 form
c
      call vscale(r4,r4,-1.d0,n*n)  
      title='r-matrix at r = '//fptoc(rstart)
      call prntrm (title,r4,n,n,n,n,iout)
      call potntl (rstart,vmat,ec,vcij,range,type,n,itri)
      title='coupling matrix at r = '//fptoc(rstart)
      call prntrm (title,vmat,n,n,n,n,iout)
      call vscale(vmat,vmat,2.d0,n*n)
c     diagonalize potential and transform r-matrix to adiabatic basis
c                           at first point.
      call tred2(n,n,vmat,enrg,space(1,1,2),vmat)
      call tql2(n,n,enrg,space(1,1,2),vmat,ierr)
      call ebc(tsave,r4,vmat,n,n,n)
      call ebtc(r4,vmat,tsave,n,n,n)
c     copy eigenvectors to tev       
      call copy(vmat,tev,n*n)
      ien=ienc
      etot=2.d0*elxxx
      okw=(ien.eq.1)
      okr=ien.ne.1
      if (okr.or.okw) then
          call iosys ('rewind all on scratch read-and-write',0,0,0,' ')
      endif          
      if (okr) then
          okpr=.false.
      endif
c     this section of code decides where to diagonalize the coupling matrix
c     and then does the diagonalization and stores the information on disk.
      if (okr) go to 80
c
c     prepare for looping through the potential
c
           rl=rstart
           dr=drstrt
           nstep=1
           r=0.d0
           or=r
           couple=0.d0
           e1=0.d0
           step=drstrt
c
c     begin main loop
c
   25 continue
           rr=min(rfinal,rl+dr)
           dr=rr-rl
           r=.5d0*(rr+rl)
           nstep=nstep+1
           last=rr.eq.rfinal
           okpv=(mod(nstep-1,npv).eq.0.and.prntv).or.last
           okpr=(mod(nstep-1,npr).eq.0.and.prnt4).and.npr.ne.999
           okpr=okpr.or.last
c
c     evaluate coupling potential, diagonalize
c
           call potntl (r,vmat,ec,vcij,range,type,n,itri)
c
           if (okpv) then
               title='coupling matrix at step  = '//itoc(nstep)//
     1               ' r = '//fptoc(r)
               call prntrm (title,vmat,n,n,n,n,iout)
           endif
           call vscale(vmat,vmat,2.d0,n*n)           
           call tred2(n,n,vmat,enrg,space(1,1,2),vmat)
           call tql2(n,n,enrg,space(1,1,2),vmat,ierr)
           call copy(vmat,tev,n*n)
           if (debug) then
               title='eigenvectors at step  = '//itoc(nstep)//
     1               ' r = '//fptoc(r)               
               call prntrm (title,vmat,n,n,n,n,iout)
               title='tsave at step  = '//itoc(nstep)//
     1               ' r = '//fptoc(r)    
               call prntrm (title,tsave,n,n,n,n,iout)
           endif
           oe1=e1
           e1=0d0
           do 40 i=1,n
              e1=e1+abs(enrg(i))
   40      continue
           e1=e1/n
c
c     construct adiabatic overlap coupling matrix
c
           if (nstep.gt.1) then
                call ebtc(tsave,vmat,vmat,n,n,n) 
           endif
c           call copy(vmat,tp,n*n)    
           if (okw) then
               nmtx1q(1)=nstep
               nmtx1q(2)=last
               rmtx1q(1)=rl
               rmtx1q(2)=r
               rmtx1q(3)=rr
               rmtx1q(4)=dr
               call iosys ('write real vdiag to scratch without '//
     1                     'rewinding',5,rmtx1q,0,' ')
               call iosys ('write real vdiag to scratch without '//
     1                      'rewinding',n,enrg,0,' ')
               call iosys ('write real vdiag to scratch without '//
     1                     'rewinding',n*n,tsave,0,' ')
           endif
           couple=0d0
           if (nstep.ne.1) then
               do 50 i=1,n
                  couple=couple+abs(tsave(i,i))
   50          continue
               couple=1.d0-couple/n
           endif               
c
c     print potential and coupling information if desired
c
           nprnt=min(n,7)
           if (debug.and.ienc.eq.1) then
               write (iout,1) nstep,rl,r,rr,dr,couple,
     1                        (enrg(i),i=1,nprnt)
           endif     
           if (nstep.eq.1) go to 80
c
c     next step size code
c
               vprime=((oe1-e1)/(or-r))**2
               step=stpmax
               if (vprime.ne.0.d+00) then
                   step=(beta/vprime)**.166666666667d0
                   step=min(step,stpmax)
                   step=max(step,stpmin)
               endif
   80 continue
c
c     propagation code
c
      if (okr) then
          call iosys ('read real vdiag from scratch without rewinding',
     1                 5,rmtx1q,0,' ')
          call iosys ('read real vdiag from scratch without rewinding',
     1                 n,enrg,0,' ')
          call iosys ('read real vdiag from scratch without rewinding',
     1                 n*n,tsave,0,' ')
      endif
      call wdcalc (enrg,-dr,pdiag,etot,n)
      call rprop1 (r4,pdiag,tsave,ipvt,nstep,n,space)
      if (okpr) then
          title='r-matrix at step = '//itoc(nstep)//
     1          ' r = '//fptoc(r)      
          call prntrm (title,r4,n,n,n,n,iout)
      endif
      if (last) go to 90
          or=r
          odr=dr
          dr=step
          rl=rr
          if (okr) go to 80
              go to 25
   90 continue
      if (prnt4.and.npr.eq.999) then
          title='r-matrix at step = '//itoc(nstep)//
     1          ' r = '//fptoc(r)
          call prntrm (title,r4,n,n,n,n,iout)
      endif
c
c     end of r-matrix propagation -- apply boundary onditions
c
      write (iout,2) nstep, r
c
      if (ienc.eq.1) then
           call iosys ('write real tev to scratch',n*n,tev,0,' ')
      else
           call iosys ('read real tev from scratch',n*n,tev,0,' ')
      endif
      call ebc(tsave,r4,vmat,n,n,n)
      call ebtc(r4,vmat,tsave,n,n,n)
      call vscale(r4,r4,-1.d0,n*n)  
      title='r-matrix at last step'
      call prntrm (title,r4,n,n,n,n,iout)
      return
    1 format(/1x,'step = ',i5,1x,'r left = ',e9.3,1x,'r center = ',e9.3,
     1        1x,'r right = ',e9.3,/,1x,'step size = ',e9.3,
     2        1x,'step coupling =',e9.3,/,15x,'lowest eigenvalues',
     3                                        (/,5e9.3))         
    2 format (/,1x,'final number of steps = ',i5,
     1          1x,'final step = ',e15.6)      
      end
