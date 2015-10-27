*deck linslv
      subroutine linslv (oper,status,a1,a2,a3,a4,v1,v2,v3,wt,ia,n,m,nmax
     1 ,nvec,iwrit,nowgt)
c***begin prologue     linslv
c***date written       861117   (yymmdd)
c***revision date      881223   (yymmdd)
c***                   program cleaned up and documented. now well tested
c***                   for many problems.
c***keywords           solve linear equations , linear equations
c***author             schneider, barry (lanl)
c***source             mylib
c***purpose            solve one or more sets of linear equations
c***description        linslv uses an iteration-variation method to
c***                   solve a set of linear algebraic equations of
c***                   large dimension by building an increasing vector
c***                   space based on the iterates of the born sequence
c***                   the set of equations is assumed to be of the form
c***                   (i + a) * x  =  y
c***                   the matrix a is used to form the iterates by
c***                   matrix multiplication on a previously computed
c***                   vector. the number of right hand sides in y may
c***                   be larger than one. the vector space is built by
c***                   using y as the initial set of vectors.
c
c***references         schneider, b. i. and collins, l. a. phys rev
c
c***routines called    sgefa(clams), sgeslv(clams), sgemm(clams), snorm
c***                   saxpy(clams)
      implicit integer(a-z)
      real *8 v1(m,*), v2(m,*), v3(n,*)
      real *8 a1(n,*), a2(n,*), a3(m,*), a4(m,*)
      real *8 wt(n)
      real *8 thresh, consln, delmax
      real *8 anorm, snorm, diff, sum, zero, one
      dimension ia(*)
      character *(*) oper, status
      logical nowgt
      save thresh, consln, mxiter, iout, nrhs, nold, ntrial
      save cycle
      parameter (zero=0.d0,one=1.d0)
c
      if (oper.eq.'prepare for solution') then
c
c----------------------------------------------------------------------c
c         prepare for solution and save local copies of important      c
c                          variables                                   c
c----------------------------------------------------------------------c
          thresh=a3(1,1)
          consln=a4(1,1)
          mxiter=ia(1)
          iout=m
          nrhs=nmax
          ntrial=nrhs
          cycle=0
          write (iout,290) mxiter,thresh,consln,iout,nrhs
c
c----------------------------------------------------------------------c
c              copy into initial set of vectors                        c
c----------------------------------------------------------------------c
c
          do 10 i=1,n
             do 11 j=1,ntrial
                a1(i,j)=a2(i,j)
   11        continue
   10     continue
c
          nvec=0
          if (iwrit.gt.0) then
              write (iout,300) (v1(i,1),i=1,n)
              do 20 i=1,ntrial
                 write (iout,310) i,(a1(j,i),j=1,n)
   20         continue
          endif
c
      elseif (oper.eq.'generate vectors') then
c
c----------------------------------------------------------------------c
c              schmidt orthonormalization of vectors                   c
c----------------------------------------------------------------------c
              nold=nvec
              nstart=nvec+1
              nfin=nvec+ntrial
              ntrial=0
              do 60 i=nstart,nfin
                 anorm=snorm(a1(1,i),a1(1,i),wt,n,nowgt)
                 if (anorm.ne.0.d+00) then
                     anorm=1.d+00/sqrt(anorm)
                 endif
                 do 30 j=1,n
                    a1(j,i)=anorm*a1(j,i)
   30            continue
                 do 40 j=1,i
                    if (i.ne.j) then
                        sum=-snorm(a1(1,j),a1(1,i),wt,n,nowgt)
                        call saxpy (n,sum,a1(1,j),1,a1(1,i),1)
                    endif
   40            continue
                 anorm=snorm(a1(1,i),a1(1,i),wt,n,nowgt)
                 if (anorm.gt.thresh) then
                     nvec=nvec+1
                     ntrial=ntrial+1
                     anorm=1.d+00/sqrt(anorm)
                     do 50 j=1,n
                        a1(j,nvec)=anorm*a1(j,i)
   50                continue
                 endif
   60         continue
   70         status='ok'
              if (ntrial.eq.0) status='none left'
              ia(1)=ntrial
              if (iwrit.gt.0) then
                  do 80 i=nstart,nvec
                     do 85 j=1,i
                        v1(i,j)=snorm(a1(1,i),a1(1,j),wt,n,nowgt)
   85                continue
   80             continue
                  do 90 i=nstart,nvec
                     write (iout,280) i,(v1(i,j),j=1,i)
   90             continue
                  do 100 i=1,nvec
                     write (iout,320) i,(a1(j,i),j=1,n)
  100             continue
              endif
c----------------------------------------------------------------------c
c          here is where user should generate iterates by              c
c                    calling appropriate routine                       c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c             now solve small set of linear equations                  c
c----------------------------------------------------------------------c
      elseif (oper.eq.'solve equations') then
c----------------------------------------------------------------------c
c             make copies of old matrix and rhs                        c
c----------------------------------------------------------------------c
              nstart=nold+1
              if (nold.eq.0) go to 140
              do 110 i=1,nold
                 do 111 j=1,nold
                    a3(i,j)=a4(i,j)
  111            continue
  110         continue
              do 120 i=1,nold
                 do 121 j=1,nrhs
                    v1(i,j)=v2(i,j)
  121            continue
  120         continue
c----------------------------------------------------------------------c
c          calculate new matrix elements of small matrix               c
c----------------------------------------------------------------------c
              do 130 i=1,nold
                 do 131 j=nstart,nvec
                    a3(i,j)=snorm(a1(1,i),a2(1,j),wt,n,nowgt)
                    a3(j,i)=snorm(a1(1,j),a2(1,i),wt,n,nowgt)
                    a4(i,j)=a3(i,j)
                    a4(j,i)=a3(j,i)
  131            continue
  130         continue
  140         do 150 i=nstart,nvec
                 do 151 j=1,nrhs
                    v1(i,j)=snorm(a1(1,i),v3(1,j),wt,n,nowgt)
                    v2(i,j)=v1(i,j)
  151            continue
  150         continue
              do 160 i=nstart,nvec
                 do 161 j=nstart,nvec
                    a3(i,j)=snorm(a1(1,i),a2(1,j),wt,n,nowgt)
  161            continue
  160         continue
              do 170 i=nstart,nvec
                 a3(i,i)=1.d+00+a3(i,i)
  170         continue
              do 180 i=nstart,nvec
                 do 181 j=nstart,nvec
                    a4(i,j)=a3(i,j)
  181            continue
  180         continue
              if (iwrit.gt.0) then
                  do 190 i=1,nvec
                     write (iout,330) i,(a3(j,i),j=1,nvec)
  190             continue
                  do 200 i=1,nrhs
                     write (iout,340) i,(v1(j,i),j=1,nvec)
  200             continue
              endif
c----------------------------------------------------------------------c
c                solve small set of equations                          c
c----------------------------------------------------------------------c
             call sgefa (a3,m,nvec,ia,info)
             do 210 i=1,nrhs
                call sgesl (a3,m,nvec,ia,v1(1,i),0)
  210        continue
c----------------------------------------------------------------------c
c                   test for convergence                               c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c                   calculate residuals                                c
c----------------------------------------------------------------------c
             cycle=cycle+1
             if (iwrit.lt.0) write (iout,360) cycle,nvec
             delmax=0.d+00
             do 240 irhs=1,nrhs
                diff=0.d+00
                do 230 i=1,n
                   sum=0.d+00
                   do 220 j=1,nvec
                      sum=sum+v1(j,irhs)*(a1(i,j)+a2(i,j))
  220              continue
                   sum=sum-v3(i,irhs)
                   diff=diff+sum*sum
  230           continue
                diff=sqrt(diff)
                if (iwrit.lt.0) write (iout,380) irhs,diff
                delmax=max(diff,delmax)
  240        continue
             if (delmax.gt.consln) then
                 mm=nvec
                 do 260 i=nstart,nvec
                    mm=mm+1
                    do 261 j=1,n
                       a1(j,mm)=a2(j,i)
  261               continue
  260            continue
c----------------------------------------------------------------------c
c                continue iterations or quit if converged              c
c----------------------------------------------------------------------c
                 status='continue iterations'
             else
                 status='converged'
                 call sgemm('n','n',n,nrhs,nvec,one,a1,n,v1,m,zero,v3,n)
c                call sgmm (n,nvec,nrhs,a1,n,v1,m,v3,n,0,1)
                 if (iwrit.gt.0) then
                     do 270 i=1,nrhs
                        write (iout,350) i,(v3(j,i),j=1,n)
  270                continue
                 endif
             endif
      elseif (oper.eq.'final best solution') then
              call sgemm('n','n',n,nrhs,nvec,one,a1,n,v1,m,zero,v3,n)
c             call sgmm (n,nvec,nrhs,a1,n,v1,m,v3,n,0,1)
              status='hit iteration limit'
      elseif (oper.eq.'cleanup') then
              thresh=0.0d+00
              consln=0.0d+00
              mxiter=0
              iout=0
              nvec=0
              nrhs=0
              nold=0
      else
      write (iout,370)
      endif
c
      return
c
  280 format (/,1x,'ovlp',1x,i3,(2x,5d15.8))
  290 format (//,5x,'initialize linslv:',/,30x,'max. iterations       =
     1 ',2x,i4,/,30x,'overlap thresh        = ',2x,d15.8,/,30x,'converge
     2nce criterion = ',2x,d15.8,/,30x,'write file            =  ',2x,i2
     3 ,/,30x,'no. right hand sides  = ',2x,i4)
  300 format ((/,1x,'diagonal',1x,5d15.8))
  310 format (/,1x,'initial vector',1x,i3,(1x,5d15.8))
  320 format (/,1x,'vec',1x,i3,(2x,5d15.8))
  330 format (/,1x,'small eqns',1x,i3,(2x,5d15.8))
  340 format (/,1x,'small rhs',1x,i3,(2x,5d15.8))
  350 format (/,1x,'converged soln',1x,i3,(2x,5d15.8))
  360 format (/,1x,25('-'),'cycle ',i3,' using ',i4,' vectors ',10('-'),
     1 /,t3,'solution',t30,'convergence')
  370 format (//,15x,'*****  error in call to linslv  *****')
  380 format (t4,i4,t26,d15.8)
      end
