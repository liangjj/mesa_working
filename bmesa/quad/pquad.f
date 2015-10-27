      subroutine pquad(z,a)
c***begin prologue     pquad
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c
c
c***keywords
c***author             schneider, barry  (nsf)
c***source             %W% %G% 
c***purpose
c                      test integration schemes
c***description
c***references
c
c***routines called
c
c***end prologue       pquad
c
      implicit integer (a-z)
c
      parameter ( mxquad=9 )
      integer a(1)
      real*8 z(1), rdel, rshel, rleft, int0f, int0b, pi
      dimension rshel(100), nrad(100), rdel(100)
      character*8 cpass
      character*20 type, chrkey, typef
      character*400 card
      character*80 title
      logical posinp, defint
      common /io/ inp, iout
      data pi / 3.141592653589793d0 /
      pi=pi-1.d-12
      if ( posinp('$quad',cpass) ) then
           call cardin(card)      
      endif
      type=chrkey(card,'quadrature-type','newton-cotes',' ')
      typef=chrkey(card,'function-type','exponential',' ')
      defint=logkey(card,'only-definite-integral',.false.,' ')
      mmax=intkey(card,'maximum-m-value',4,' ')
      nshell=intkey(card,'number-of-shells',1,' ')
      if (typef.ne.'trigonometric') then
          call fparr(card,'shell-radii',rshel,nshell,' ')
          call intarr(card,'number-of-points-per-shell',nrad,nshell,' ')
      else
          nrad(1)=intkey(card,'number-of-points-in-shell',3,' ')
          rshel(1)=2.d0*pi/nshell
          do 55 i=2,nshell
             rshel(i)=rshel(i-1)+rshel(1)
             nrad(i)=nrad(1)
   55     continue
      endif                 
      if (type.eq.'newton-cotes') then
          npts=0
          nints=0
          do 10 i=1,nshell
             npts=npts+nrad(i)
             nints=nints+nrad(i)-1
   10     continue
          pt=1
          fun=pt+npts
          if(typef.eq.'exponential') then
             wt=fun+npts
          elseif (typef.eq.'trigonometric') then
             wt=fun+npts*(mmax+1)*2
          else
             call lnkerr('error in function type')
          endif          
          sumf=wt+mxquad*mxquad*nshell
          sumb=sumf+nints
          intrig=sumb+nints
          words=intrig+(mmax+1)*2
          locpt=pt
          locfun=fun
          locwt=wt
          rleft=0.d0
          write (iout,1) nshell
          write (iout,2) (rshel(i),i=1,nshell)
          write (iout,3) (nrad(i),i=1,nshell)
          do 20 i=1,nshell
c         call mkgrd(z(locpt),rleft,rshel(i),rdel(i),nrad(i))
             call necote(rleft,rshel(i),z(locpt),z(locwt),nrad(i),
     1                   .false.)
             call mkfun(z(locpt),z(locfun),nrad(i),mmax,typef)
             rleft=rshel(i)
             locpt=locpt+nrad(i)
             if(typef.eq.'exponential') then
                locfun=locfun+nrad(i) 
             else
                locfun=locfun+nrad(i)*2*(mmax+1)
             endif
             locwt=locwt+nrad(i)*(nrad(i)-1)
   20     continue
          if (.not.defint) then
              call rzero(z(sumf),nints)
              call rzero(z(sumb),nints)   
              locfun=fun
              locwt=wt
              locsmf=sumf
              int0f=0.d0
              int0b=0.d0
              do 30 i=1,nshell
                 call forwrd(z(locfun),z(locwt),int0f,z(locsmf),nrad(i))
                 locfun=locfun+nrad(i)
                 locwt=locwt+nrad(i)*(nrad(i)-1)
                 locsmf=locsmf+nrad(i)-1
   30         continue
              write( iout,4) (z(i),i=sumf,sumf+nints-1)      
c      
              locsmb=sumb+nints
              last=locsmb-1
              do 40 i=nshell,1,-1
                 locfun=locfun-nrad(i)
                 locwt=locwt-nrad(i)*(nrad(i)-1)
                 locsmb=locsmb-nrad(i)+1
                 call bkwrd(z(locfun),z(locwt),int0b,z(locsmb),nrad(i))
   40         continue
              write( iout,5) (z(i),i=locsmb,last)
          else
              call rzero(z(intrig),(mmax+1)*2)
c              call testnt(z(intrig),mmax)
c              call lnkerr('quit')
              locfun=fun
              locwt=wt
              do 50 i=1,nshell
                 call fulint(z(locfun),z(locwt),z(intrig),nrad(i),
     1                       mmax,typef)
                 locwt=locwt+nrad(i)*(nrad(i)-1)
                 if (typef.eq.'exponential') then
                     locfun=locfun+nrad(i)
                 else    
                     locfun=locfun+nrad(i)*2*(mmax+1)
                 endif    
   50         continue
              if (typef.eq.'exponential') then
                  write(iout,6) z(intrig)
              else
                  title='trigonometric sine and cosine integrals'
                  call prntfm(title,z(intrig),mmax+1,2,mmax+1,2,iout)
              endif        
          endif
      else
          call cheby(nrad(1),mmax)
      endif
    1 format (/,5x,'number of radial shells:',1x,i3)
    2 format (/,5x,'shell radii:',5(/,20x,4(e15.8,1x)))
    3 format (/,5x,'shell quadrature orders:',10(/,20x,10(i4,1x)))
    4 format (/,5x,'values of forward integration:',5(/,10x,
     1          4(e15.8,1x)))
    5 format (/,5x,'values of backward integration:',5(/,10x,
     1          4(e15.8,1x)))
    6 format(/,1x,'exponential integral',1x,e15.8)             
      call chainx(0)
c
c
      stop
      end


