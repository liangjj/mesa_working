*deck hamslv
c***begin prologue     hamslv
c***date written       861117   (yymmdd)
c***revision date      890806   (yymmdd)
c***                   separated matrix formation and solution in two
c***                   calls.   
c***keywords           solve linear equations , linear equations
c***author             schneider, barry (lanl)
c***source             mylib
c***purpose            solve one or more sets of linear equations
c***description        hamslv uses an iteration-variation method to
c***                   solve a set of linear algebraic equations of
c***                   large dimension by building an increasing vector
c***                   space based on an initial set of input guesses.
c***                   primarily for use in forming optical potentials.
c***                   the coefficient matrix is used to form the
c***                   iterates by matrix multiplication on previous
c***                   vectors. the number of right hand sides may
c***                   be larger than one.
c
c***references         schneider, b. i.
c
c***routines called    sgefa(clams), sgeslv(clams), sgemm(clams), saxpy(clams)
      subroutine hamslv (oper,status,a1,a2,a3,a4,v1,v2,v3,ia,n,m,nmax
     1 ,nvec,iwrit)
      implicit integer(a-z)
      real *8 v1(m,*), v2(m,*), v3(n,*)
      real *8 a1(n,*), a2(n,*), a3(m,*), a4(m,*)
      real *8 thresh, converg, delmax
      real *8 anorm, sdot, diff, sum, zero, one
      dimension ia(*)
      character *(*) oper, status
      save thresh, converg, mxiter, iout, nrhs, nold, ntrial
      save cycle
      parameter (zero=0.d0,one=1.d0)
c
      if (oper.eq.'prepare for solution') then
c
c      ------ prepare matrix for iteration-variation sequence ------
c
c     ------ save local copies of important variables ------
c
         thresh=a3(1,1)
         converg=a4(1,1)
         mxiter=ia(1)
         iout=m
         nrhs=nmax
         ntrial=nvec
         cycle=0
         write (iout,290) mxiter,thresh,converg,iout,nrhs,ntrial
c
c
c     ------   copy into first ntrial vectors  -----
c
         do 20 i=1,n
            do 21 j=1,ntrial
               a1(i,j)=a2(i,j)
   21       continue
   20    continue
         if (iwrit.gt.0) then
             do 30 i=1,ntrial
                write (iout,310) i,(a1(j,i),j=1,n)
   30        continue
          endif
      elseif (oper.eq.'prepare for new energy') then
         thresh=a3(1,1)
         converg=a4(1,1)
         mxiter=ia(1)
         iout=m
         nrhs=nmax
         ntrial=nvec
         cycle=0
         write (iout,290) mxiter,thresh,converg,iout,nrhs,ntrial
c
      elseif (oper.eq.'generate vectors') then
c
c     ----- schmidt orthogonalize to old vectors and then -----
c     -----                 normalize        -----
         nold=nvec
         nstart=nvec+1
         nfin=nvec+ntrial
         ntrial=0
         icnt=nvec
         do 70 i=nstart,nfin
            anorm=sdot(n,a1(1,i),1,a1(1,i),1)
            if (anorm.gt.thresh) then
                anorm=1.d+00/sqrt(anorm)
            endif
            do 40 j=1,n
               a1(j,i)=anorm*a1(j,i)
   40       continue
            do 50 j=1,i
               if (i.ne.j) then
                   sum=-sdot(n,a1(1,j),1,a1(1,i),1)
                   call saxpy (n,sum,a1(1,j),1,a1(1,i),1)
               endif
   50       continue
            anorm=sdot(n,a1(1,i),1,a1(1,i),1)
            if (anorm.gt.thresh) then
                icnt=icnt+1
                ntrial=ntrial+1
                anorm=1.d+00/sqrt(anorm)
                do 60 j=1,n
                   a1(j,icnt)=anorm*a1(j,i)
   60           continue
            endif
   70    continue
   80    status='ok'
         if (ntrial.eq.0) status='none left'
         ia(1)=ntrial
         if (iwrit.gt.0) then
             do 90 i=nstart,nvec
                do 91 j=1,i
                   v1(i,j)=sdot(n,a1(1,i),1,a1(1,j),1)
   91           continue
   90        continue
             do 100 i=nstart,nvec
                write (iout,280) i,(v1(i,j),j=1,i)
  100        continue
             do 110 i=1,nvec
                write (iout,320) i,(a1(j,i),j=1,n)
  110        continue
         endif
c
c     ----------  in the main code the user should  ----------
c     ------ generate iterates corresponding to new vectors ------
c     ------ by operating with matrix with zero diagonals ------
c     ------ on the previous iterates and then divide by energy ------
c     ------ shifted inverse of the diagonal ------
c
      elseif (oper.eq.'form small matrix') then
c
c     ------ copy old small matrix and right hand side into ------
c     ------ new matrix and right hand side. ------
         nstart=nold+1
         if (nold.eq.0) go to 150
c
c     ------   copy old matrix and rhs  -----
c
         do 120 i=1,nold
            do 121 j=1,nold
               a3(i,j)=a4(i,j)
  121       continue
  120    continue
         do 130 i=1,nold
            do 131 j=1,nrhs
               v1(i,j)=v2(i,j)
  131       continue
  130    continue
c
c     ------ calculate interaction with old vectors -----
c
         do 140 i=1,nold
            do 141 j=nstart,nvec
               a3(i,j)=-sdot(n,a1(1,i),1,a2(1,j),1)
               a3(j,i)=-sdot(n,a1(1,j),1,a2(1,i),1)
               a4(i,j)=a3(i,j)
               a4(j,i)=a3(j,i)
  141       continue
  140    continue
c
c     ------ calculate interactions in new block -----
c
  150    do 160 i=nstart,nvec
            do 161 j=1,nrhs
               v1(i,j)=sdot(n,a1(1,i),1,v3(1,j),1)
               v2(i,j)=v1(i,j)
  161       continue
  160    continue
         do 170 i=nstart,nvec
            do 171 j=nstart,nvec
               a3(i,j)=-sdot(n,a1(1,i),1,a2(1,j),1)
  171       continue
  170    continue
         do 180 i=nstart,nvec
            a3(i,i)=1.d+00+a3(i,i)
  180    continue
         do 190 i=nstart,nvec
            do 191 j=nstart,nvec
               a4(i,j)=a3(i,j)
  191       continue
  190    continue
         if (iwrit.gt.0) then
             do 200 i=1,nvec
                write (iout,330) i,(a3(j,i),j=1,nvec)
  200        continue
             do 210 i=1,nrhs
                write (iout,340) i,(v1(j,i),j=1,nvec)
  210        continue
         endif
      elseif (oper.eq.'solve small set of equations') then
         call sgefav (a3,m,nvec,ia,info)
         do 220 i=1,nrhs
            call sgeslv (a3,m,nvec,ia,v1(1,i),0)
  220    continue
c
c     ------          test for convergence      ------
c
c     -----    calculate residuals   -----
c
         cycle=cycle+1
         if (iwrit.lt.0) write (iout,360) cycle, nvec
         delmax=0.d+00
         do 250 irhs=1,nrhs
            diff=0.d+00
            do 240 i=1,n
               sum=0.d+00
               do 230 j=1,nvec
                  sum=sum+v1(j,irhs)*(a1(i,j)-a2(i,j))
  230          continue
               sum=sum-v3(i,irhs)
               diff=diff+sum*sum
  240       continue
            diff=sqrt(diff)
            if (iwrit.lt.0) write (iout,380) irhs,diff
            delmax=max(diff,delmax)
  250    continue
         if (delmax.gt.converg) then
            mm=nvec
            do 260 i=nstart,nvec
               mm=mm+1
               do 261 j=1,n
                  a1(j,mm)=a2(j,i)
  261          continue
  260       continue
            status='continue iterations'
         else
            status='converged'
            call sgemm('n','n',n,nrhs,nvec,one,a1,n,v1,m,zero,v3,n)
c            call sgmm (n,nvec,nrhs,a1,n,v1,m,v3,n,0,1)
            if (iwrit.gt.0) then
                do 270 i=1,nrhs
                   write (iout,350) i,(v3(j,i),j=1,n)
  270           continue
            endif
         endif
      elseif (oper.eq.'replace new vectors') then
          call sgemm('n','n',n,nrhs,nvec,one,a1,n,v1,m,zero,v3,n)
c         call sgmm (n,nvec,nrhs,a1,n,v1,m,v3,n,0,1)
         mm=nvec
            do 400 i=1,nrhs
               mm=mm+1
               do 410 j=1,n
                  a1(j,mm)=v3(j,i)
  410          continue
  400       continue
      elseif (oper.eq.'cleanup') then
         thresh=0.0d+00
         converg=0.0d+00
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
  290 format (//,5x,'initialize hamslv:',/,30x,'max. iterations       =
     1 ',2x,i4,/,30x,'overlap thresh        = ',2x,d15.8,/,30x,'converge
     2nce criterion = ',2x,d15.8,/,30x,'write file            =  ',2x,i2
     3 ,/,30x,'no. right hand sides  = ',2x,i4,/,30x,'no. trial vectors
     4    = ',2x,i4)
  300 format ((/,1x,'diagonal',1x,5d15.8))
  310 format (/,1x,'initial vector',1x,i3,(1x,5d15.8))
  320 format (/,1x,'vec',1x,i3,(2x,5d15.8))
  330 format (/,1x,'small eqns',1x,i3,(2x,5d15.8))
  340 format (/,1x,'small rhs',1x,i3,(2x,5d15.8))
  350 format (/,1x,'converged soln',1x,i3,(2x,5d15.8))
  360 format (/,1x,25('-'),'cycle ',i3,' using ',i4,' vectors ',10('-'),
     1 /,t3,'solution',t30,'convergence')
  370 format (//,15x,'*****  error in call to hamslv  *****')
  380 format (t4,i4,t26,d15.8)
      end
