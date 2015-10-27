         program dum
         implicit integer(a-z)
         parameter (imax=2,jmax=1,mhi=16)
         integer done(0:imax,0:jmax,0:mhi)
         integer needed(0:imax,0:jmax,0:mhi)
         do 60 m=0,mhi
            done(0,0,m)=0
            needed(0,0,m)=0
   60    continue
         do 50 m=0,imax+jmax
            done(0,0,m)=1
   50    continue
         do 100 i=1,imax
         do 100 j=1,jmax
         do 100 m=0,mhi
            done(i,j,m)=0
            needed(i,j,m)=0
  100    continue
         coord=1
         iatom=1
         jatom=1
         half=0.5
c
c
         lj=0
         mmax=imax+jmax
         do 180 li=0,imax-1
         do 170 m=0,mmax-li-lj-1
            done(li+1,lj,m)=1
c           call vwxy(xyz(1,li+1,lj,coord,m),xyz(1,li+1,lj,coord,m),
c    $                expon,xyz(1,li,lj,coord,m+1),+1,n)
            write(6,*) li+1,lj,m, 'needs:',li,lj,m+1
            needed(li,lj,m+1)=1
            if((needed(li,lj,m+1).eq.1).and.(done(li,lj,m+1).eq.0)) then
               write(6,*) 'conflict',li,lj,m
            endif
            if(li.gt.0) then
c              scalri=float(i)*half
c              call smul(t,ainv(1,1),scalri,n)
c              call vwxy(xyz(1,li+1,lj,coord,m),xyz(1,li+1,lj,coord,m),
c    $                   t,xyz(1,li-1,lj,coord,m),+1,n)
c              call vmul(t,t,ainv(nv,4),n)
c              call vwxy(xyz(1,li+1,lj,coord,m),xyz(1,li+1,lj,coord,m),
c    $                   t,xyz(1,li-1,lj,coord,m+1),-1,n)
            write(6,*) li+1,lj,m,'needs:',li-1,lj,m,' and',li-1,lj,m+1
            needed (li-1,lj,m)=1
            needed (li-1,lj,m+1)=1
            if((needed(li-1,lj,m).eq.1).and.(done(li-1,lj,m).eq.0)) then
               write(6,*) 'conflict',li,lj,m
            endif
        if((needed(li-1,lj,m+1).eq.1).and.(done(li-1,lj,m+1).eq.0)) then
               write(6,*) 'conflict',li,lj,m
            endif
            endif
  170   continue
  180  continue
      write(6,*) 'finished first pass'
         li=0
         mmax=imax+jmax
         do 280 lj=0,jmax-1
         do 270 m=0,0
c        do 270 m=0,mmax-li-lj-1
            done(li,lj+1,m)=1
c           call vwxy(xyz(1,li,lj+1,coord,m),xyz(1,li,lj+1,coord,m),
c    $                expon,xyz(1,li,lj,coord,m+1),+1,n)
            write(6,*) li,lj+1,m, 'needs:',li,lj,m+1
            needed(li,lj,m+1)=1
            if((needed(li,lj,m+1).eq.1).and.(done(li,lj,m+1).eq.0)) then
               write(6,*) 'conflict',li,lj,m
            endif
            if(lj.gt.0) then
c              scalrj=float(j)*half
c              call smul(t,ainv(1,1),scalrj,n)
c              call vwxy(xyz(1,li,lj+1,coord,m),xyz(1,li,lj+1,coord,m),
c    $                   t,xyz(1,li,lj-1,coord,m),+1,n)
c              call vmul(t,t,ainv(nv,4),n)
c              call vwxy(xyz(1,li,lj+1,coord,m),xyz(1,li,lj+1,coord,m),
c    $                   t,xyz(1,li,lj-1,coord,m+1),-1,n)
            write(6,*) li,lj+1,m,'needs:',li,lj-1,m,' and',li,lj-1,m+1
            needed (li,lj-1,m)=1
            needed (li,lj-1,m+1)=1
            if((needed(li,lj-1,m).eq.1).and.(done(li,lj-1,m).eq.0)) then
               write(6,*) 'conflict',li,lj,m
            endif
        if((needed(li,lj-1,m+1).eq.1).and.(done(li,lj-1,m+1).eq.0)) then
               write(6,*) 'conflict',li,lj,m
            endif
            endif
  270   continue
  280  continue
      write(6,*) 'finished second pass'
c
c        ----- bump the angular momentum on center b -----
c              this generates all (a,b) integrals
c        call ssub(expon,xyza(1,coord),c(coord,jatom),n)
         do 190 lj=0,jmax-1
            do 175 li=1,imax
c           do 174 m=0,mmax-li-lj-1
            do 174 m=0,mmax-li-lj-1 
               done(li,lj+1,m)=1
c    $                   expon,xyz(1,li,lj,coord,m+1),+1,n)
            write(6,*) li,lj+1,m,'needs:',li,lj,m+1
            needed(li,lj,m+1)=1
            if((needed(li,lj,m+1).eq.1).and.(done(li,lj,m+1).eq.0)) then
               write(6,*) 'conflict',li,lj,m
            endif
               if(lj.gt.0) then
c                 scalrj=float(j)*half
c                 call smul(t,ainv(1,2),scalrj,n)
c                 call vwxy(xyz(1,li,lj+1,coord,m),
c    $                      xyz(1,li,lj+1,coord,m),
c    $                      t,xyz(1,li,lj-1,coord,m),+1,n)
c                 call vmul(t,t,ainv(nv,5),n)
c                 call vwxy(xyz(1,li,lj+1,coord,m),
c    $                      xyz(1,li,lj+1,coord,m),
c    $                      t,xyz(1,li,lj-1,coord,m+1),-1,n)
         write(6,*) li,lj+1,m,'needs:',li,lj-1,m,'and',li,lj-1,m+1
         needed(li,lj-1,m)=1
         needed(li,lj-1,m+1)=1
            if((needed(li,lj-1,m).eq.1).and.(done(li,lj-1,m).eq.0)) then
               write(6,*) 'conflict',li,lj,m
            endif
        if((needed(li,lj-1,m+1).eq.1).and.(done(li,lj-1,m+1).eq.0)) then
               write(6,*) 'conflict',li,lj,m
            endif
               endif
               if(li.gt.0) then
c                 scalri=float(i)*half
c                 call smul(t,ainv(1,3),scalri,n)
c                 call vwxy(xyz(1,li,lj+1,coord,m),
c    $                      xyz(1,li,lj+1,coord,m),
c    $                      t,xyz(1,li-1,lj,coord,m+1),+1,n)
          write(6,*) li,lj+1,m, 'needs:',li-1,lj,m+1
         needed(li-1,lj,m+1)=1
       if((needed(li-1,lj,m+1).eq.1).and.(done(li-1,lj,m+1).eq.0)) then
               write(6,*) 'conflict',li,lj,m
            endif
               endif
  174       continue
  175    continue
  190    continue
  200 continue
         do 300 i=0,imax
         do 300 j=0,jmax
         do 300 m=0,mhi
            write(6,*) i,j,m,'need done',needed(i,j,m),done(i,j,m)
c           if (done(i,j,m).and.not.needed(i,j,m)) then
c              write(6,*) 'overkll',i,j,m
c           endif
  300    continue
      stop
      end
