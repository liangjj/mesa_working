*deck modbes
c***begin prologue     modbes
c***date written       890801   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6004, link 6004, modified bessel, spline
c***author             schneider, barry (lanl)
c***source             m6004
c***purpose            modified spherical bessel funcions
c***description        the analog of ricatti-bessel functions, the
c***                   modified spherical bessel functions multiplied
c***                   by x are computed by backward recursion for the
c***                   the bessel like function and forward recursion for
c***                   neumann like function. this routine is for the 
c***                   exponentially decaying or growing free wave functions.
c***references       
c
c***routines called
c***end prologue       modbes
      subroutine modbes (x,mj,dmj,my,dmy,wron,scr,np,lmax,ltop,dir,prnt)
      implicit integer (a-z)
      parameter (acc=30 )
      real*8 x, mj, dmj, my, dmy, scr, tmp, wron
      character*(*) dir
      logical prnt
      dimension x(np), scr(np,3), mj(np,0:ltop), dmj(np,0:ltop)
      dimension my(np,0:top), dmy(np,0:ltop)
      common /io/ inp, iout
c----------------------------------------------------------------------c
c            estimate starting l                                       c
c----------------------------------------------------------------------c
      lfinal=min(ltop,lmax+1)     
      tmp=lmax*acc
      tmp=sqrt(tmp) 
      atmp=tmp
      strtl=lmax+atmp
      strtl=max(strtl,lfinal)
      if (strtl.gt.ltop) then
          call lnkerr('starting l bigger than ltop')
      endif
c----------------------------------------------------------------------c
c       calculate first two members to start upward recursion          c
c----------------------------------------------------------------------c
      call vmexp(scr(1,2),x,np)
      call copy(scr(1,2),my(1,0),np)
      do 10 i=1,np
         mj(i,0)=.5d0*(1.d0/scr(i,2) - scr(i,2) )
   10 continue
      call copy(mj(1,0),scr(1,1),np)
      if (dir.eq.'derivatives') then
          do 15 i=1,np
             dmj(i,0)=.5d0*(1.d0/scr(i,2) + scr(i,2) )
   15     continue
          call vneg(dmy(1,0),my(1,0),np)
      endif
      if (lmax.ge.1) then
          call vinv(scr(1,3),x,np)
          do 20 i=1,np
             mj(i,1)=-mj(i,0)*scr(i,3)+.5d0*(1.d0/scr(i,2)+scr(i,2))
             my(i,1)=(1.d+00+scr(i,3))*scr(i,2) 
   20     continue
          if (dir.eq.'derivatives') then
              do 25 i=1,np
                 dmj(i,1)=mj(i,0)-mj(i,1)*scr(i,3)
                 dmy(i,1)=-my(i,0)-my(i,1)*scr(i,3)
   25         continue
          endif 
c----------------------------------------------------------------------c
c              calculate my by upward recursion                         c
c----------------------------------------------------------------------c
          do 30 l=1,lfinal-1
             ll1=l+l+1
             lm1=l-1
             lp1=l+1
             do 40 i=1,np
                my(i,lp1)=my(i,lm1)+ll1*scr(i,3)*my(i,l)
   40        continue
   30     continue
c----------------------------------------------------------------------c
c             calculate mj by downward recursion                       c
c----------------------------------------------------------------------c
          call rzero(mj(1,strtl),np)
          onelss=strtl-1
          call vfill(mj(1,onelss),1.d0,np) 
          do 60 l=onelss-1,0,-1
             ll3=l+l+3
             lp1=l+1
             lp2=l+2
             do 70 i=1,np
                mj(i,l)=ll3*mj(i,lp1)*scr(i,3)+mj(i,lp2)
   70        continue
   60     continue
c----------------------------------------------------------------------c
c                normalize the mj                                       c
c----------------------------------------------------------------------c
          call vdiv(scr(1,1),scr(1,1),mj(1,0),np)
          do 90 l=0,lfinal
             call vmul (mj(1,l),mj(1,l),scr(1,1),np)
   90     continue
c---------------------------------------------------------------------c
c             finish calculation by getting derivatives               c
c---------------------------------------------------------------------c
          if (dir.eq.'derivatives') then
              do 200 l=0,lmax
                 lp1=l+1
                 do 300 i=1,np
                    dmj(i,l)=lp1*mj(i,l)*scr(i,3)+mj(i,lp1)            
                    dmy(i,l)=lp1*my(i,l)*scr(i,3)-my(i,lp1) 
  300           continue
  200        continue
          endif
      endif
      if (prnt) then
          do 5000 l=0,lmax
             write (iout,900) l
             write (iout,1000) (mj(i,l),i=1,np)
             write (iout,1010) l
             write (iout,1000) (dmj(i,l),i=1,np)
             write (iout,1500) l
             write (iout,1000) (my(i,l),i=1,np)
             write (iout,1600) l
             write (iout,1000) (dmy(i,l),i=1,np)
 5000     continue
      endif
  900 format(/,5x,'regular modified function l = ',1x,i3)
 1000 format( (/,5x,5(e15.8,1x) ) )
 1010 format(/,5x,'derivative regular modified function l = ',1x,i3)
 1500 format(/,5x,'irregular modified function l = ',1x,i3)
 1600 format(/,5x,'derivative irregular modified function l = ',1x,i3)
      return
      end

