         subroutine loops(imax,jmax)
         implicit integer(a-z)
         parameter (imx=16,jmx=16,mhi=16)
         integer done(0:imx,0:jmx,0:mhi)
         integer needed(0:imx,0:jmx,0:mhi)
c
         common/io/inp,iout
c
c        this routine runs through the loops which the two-electron
c        recursion routine will take, and checks to see that everything
c        which must be calculated is done before it is needed.
c        a value of 0 implies it is not done, while a value of 1 implies
c        it has been done.
         do 60 m=0,mhi
            done(0,0,m)=0
            needed(0,0,m)=0
   60    continue
         do 50 m=0,imax+jmax
            done(0,0,m)=1
   50    continue
         do 100 i=1,imx
         do 100 j=1,jmx
         do 100 m=0,mhi
            done(i,j,m)=0
            needed(i,j,m)=0
  100    continue
c
c
         lj=0
         mmax=imax+jmax
         do 180 li=0,imax-1
         do 170 m=0,mmax-li-lj-1
            done(li+1,lj,m)=1
c           write(iout,*) li+1,lj,m, 'needs:',li,lj,m+1
            needed(li,lj,m+1)=1
            if((needed(li,lj,m+1).eq.1).and.(done(li,lj,m+1).eq.0)) then
               write(iout,*) 'conflict',li,lj,m
               call lnkerr('conflict in loops')
            endif
            if(li.gt.0) then
c              write(iout,*) li+1,lj,m,'needs:',li-1,lj,m,
c    $                              ' and',li-1,lj,m+1
               needed (li-1,lj,m)=1
               needed (li-1,lj,m+1)=1
               if((needed(li-1,lj,m).eq.1)
     $             .and.(done(li-1,lj,m).eq.0)) then
                  write(iout,*) 'conflict',li,lj,m
                  call lnkerr('conflict in loops')
               endif
               if((needed(li-1,lj,m+1).eq.1)
     $             .and.(done(li-1,lj,m+1).eq.0)) then
                  write(iout,*) 'conflict',li,lj,m
                  call lnkerr('conflict in loops')
               endif
            endif
  170   continue
  180  continue
c
c        ----- bump the angular momentum on center b -----
c              this generates all (a,b) integrals
         do 190 lj=0,jmax-1
            do 175 li=0,imax
               do 174 m=0,mmax-li-lj-1
                  done(li,lj+1,m)=1
c                 write(iout,*) li,lj+1,m,'needs:',li,lj,m+1
                  needed(li,lj,m+1)=1
                  if((needed(li,lj,m+1).eq.1)
     $               .and.(done(li,lj,m+1).eq.0)) then
                     write(iout,*) 'conflict',li,lj,m
                     call lnkerr('conflict in loops')
                  endif
                  if(lj.gt.0) then
c                    write(iout,*) li,lj+1,m,'needs:',li,lj-1,m,
c    $                                    'and',li,lj-1,m+1
                     needed(li,lj-1,m)=1
                     needed(li,lj-1,m+1)=1
                     if((needed(li,lj-1,m).eq.1)
     $                  .and.(done(li,lj-1,m).eq.0)) then
                        write(iout,*) 'conflict',li,lj,m
                        call lnkerr('conflict in loops')
                     endif
                     if((needed(li,lj-1,m+1).eq.1)
     $                  .and.(done(li,lj-1,m+1).eq.0)) then
                        write(iout,*) 'conflict',li,lj,m
                        call lnkerr('conflict in loops')
                     endif
                  endif
                  if(li.gt.0) then
c                    write(iout,*) li,lj+1,m, 'needs:',li-1,lj,m+1
                     needed(li-1,lj,m+1)=1
                     if((needed(li-1,lj,m+1).eq.1)
     $                  .and.(done(li-1,lj,m+1).eq.0)) then
                        write(iout,*) 'conflict',li,lj,m
                        call lnkerr('conflict in loops')
                  endif
               endif
  174       continue
  175    continue
  190    continue
  200 continue
c        do 300 i=0,imax
c        do 300 j=0,jmax
c        do 300 m=0,mhi
c           write(iout,*) i,j,m,'need done',needed(i,j,m),done(i,j,m)
c           if (done(i,j,m).and.not.needed(i,j,m)) then
c              write(iout,*) 'overkll',i,j,m
c           endif
c 300    continue
c
c
      write(iout,*) 'loops looks ok'
      return
      end
