c $Header: linslv.f,v 1.2 92/12/12 09:34:18 bis Exp $
*deck linslv
      subroutine linslv (oper,status,a1,a2,a3,a4,v1,v2,v3,wt,ia,n,m,nmax
     1 ,nvec,iwrit,nowgt,dim)
c***begin prologue     linslv
c***date written       861117   (yymmdd)
c***revision date      921206   (yymmdd)
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
c***                                   a * x  =  y
c***                                       or
c***                             ( 1 + a ) * x = y
c***                             (   -   )
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
      real *8 v1(m,*), v2(m,*), v3(dim,*)
      real *8 a1(dim,*), a2(dim,*), a3(m,*), a4(m,*)
      real *8 wt(n)
      real *8 thresh, consln, delmax, sdot
      real *8 anorm, snorm, diff, sum, zero, one, mone
      dimension ia(*)
      character *(*) oper, status
      character*30 slndir 
      character*80 title
      logical nowgt
      save thresh, consln, mxiter, iout, nrhs, nold, ntrial
      save cycle, slndir

      parameter (zero=0.d0,one=1.d0,mone=-1.d0)
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
          slndir=status
          write (iout,290) mxiter,thresh,consln,iout,nrhs,slndir
c
c----------------------------------------------------------------------c
c              copy into initial set of vectors                        c
c----------------------------------------------------------------------c
c
          call copy(a2,a1,dim*ntrial)
c
          nvec=0
          if (iwrit.gt.0) then
              title='initial vectors'
              call prntfm(title,a1,dim,ntrial,dim,ntrial,iout)
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
                 call sscal(dim,anorm,a1(1,i),1)
                 do 40 j=1,i
                    if (i.ne.j) then
                        sum=-snorm(a1(1,j),a1(1,i),wt,n,nowgt)
                        call saxpy (dim,sum,a1(1,j),1,a1(1,i),1)
                    endif
   40            continue
                 anorm=snorm(a1(1,i),a1(1,i),wt,n,nowgt)
                 if (anorm.gt.thresh) then
                     nvec=nvec+1
                     ntrial=ntrial+1
                     anorm=1.d+00/sqrt(anorm)
                     do 50 j=1,dim
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
                        v1(j,i)=v1(i,j)
   85                continue
   80             continue
                  nel=nvec-nstart+1
                  title='overlap matrix'
                  call prntfm(title,v1(nstart,1),nel,nvec,m,nvec,iout)
                  title='vectors'
                  call prntfm(title,a1,dim,nvec,dim,nvec,iout)
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
              nel=nvec-nstart+1
              if (nold.eq.0) go to 140
              do 110 i=1,nold
                 call copy(a4(1,i),a3(1,i),nold)
  110         continue
              do 120 i=1,nrhs
                 call copy(v2(1,i),v1(1,i),nold)
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
                    a4(i,j)=a3(i,j)
  161            continue
  160         continue
              if (slndir.eq.'one-plus-matrix') then
                  call rzero(a3,m*nvec)
                  do 180 i=1,nvec
                     a3(i,i)=1.d+00
  180             continue
                  call madd(one,a3,m,one,a4,m,a3,m,nvec,nvec)
              endif
              if (slndir.eq.'one-minus-matrix') then
                  call rzero(a3,m*nvec)
                  do 175 i=1,nvec
                     a3(i,i)=1.d+00
  175             continue
                  call madd(one,a3,m,mone,a4,m,a3,m,nvec,nvec)
              endif
              if (iwrit.gt.0) then
                  title='small equations'
                  call prntfm(title,a3,nvec,nvec,m,nvec,iout)
                  title='small rhs'
                  call prntfm(title,v1,nvec,nrhs,m,nvec,iout)
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
             if (slndir.eq.'matrix') then
                 do 240 irhs=1,nrhs
                    diff=0.d+00
                    do 230 i=1,n
                       sum=sdot(nvec,a2(i,1),dim,v1(1,irhs),1)
     1                          - v3(i,irhs)
                       diff=diff+sum*sum
  230               continue
                    diff=sqrt(diff)
                    if (iwrit.lt.0) write (iout,380) irhs,diff
                    delmax=max(diff,delmax)
  240            continue
             elseif (slndir.eq.'one-plus-matrix') then
                 do 245 irhs=1,nrhs
                    diff=0.d+00
                    do 235 i=1,n
                       sum=sdot(nvec,a2(i,1),dim,v1(1,irhs),1) +
     1                     sdot(nvec,a1(i,1),dim,v1(1,irhs),1) - 
     2                     v3(i,irhs)
                       diff=diff+sum*sum
  235               continue
                    diff=sqrt(diff)
                    if (iwrit.lt.0) write (iout,380) irhs,diff
                    delmax=max(diff,delmax)
  245            continue
             elseif (slndir.eq.'one-minus-matrix') then
                 do 250 irhs=1,nrhs
                    diff=0.d+00
                    do 255 i=1,n
                       sum=sdot(nvec,a1(i,1),dim,v1(1,irhs),1) -
     1                     sdot(nvec,a2(i,1),dim,v1(1,irhs),1) - 
     2                     v3(i,irhs)
                       diff=diff+sum*sum
  255               continue
                    diff=sqrt(diff)
                    if (iwrit.lt.0) write (iout,380) irhs,diff
                    delmax=max(diff,delmax)
  250            continue
             else
                 call lnkerr('error in linslv for solution call')
             endif
             if (delmax.gt.consln) then
                 mm=nvec
                 do 260 i=nstart,nvec
                    mm=mm+1
                    call copy(a2(1,i),a1(1,mm),dim)
  260            continue
c----------------------------------------------------------------------c
c                continue iterations or quit if converged              c
c----------------------------------------------------------------------c
                 status='continue iterations'
             else
                 status='converged'
                 call sgemm('n','n',dim,nrhs,nvec,one,a1,dim,v1,m,
     1                       zero,v3,dim)
c                call sgmm (n,nvec,nrhs,a1,n,v1,m,v3,n,0,1)
                 if (iwrit.gt.0) then
                     title='converged solution'
                     call prntfm(title,v3,dim,nrhs,dim,nrhs,iout)  
                     iwrit=0 
                 endif
             endif
      elseif (oper.eq.'final best solution') then
              call sgemm('n','n',dim,nrhs,nvec,one,a1,dim,v1,m,zero,
     1                                             v3,dim)
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
              slndir=' '
      else
      write (iout,370)
      endif
c
      return
c
  290 format (//,6x,'initialize linslv:',/,30x,'max. iterations       =
     1 ',2x,i4,/,30x,'overlap thresh        = ',2x,d15.8,/,30x,'converge
     2nce criterion = ',2x,d15.8,/,30x,'write file            =  ',2x,i2
     3 ,/,30x,'no. right hand sides  = ',2x,i4,/,30x,'form of matrix   
     4    = ',a30)
  360 format (/,1x,25('-'),'cycle ',i3,' using ',i4,' vectors ',10('-'),
     1 /,t3,'solution',t30,'convergence')
  370 format (//,15x,'*****  error in call to linslv  *****')
  380 format (t4,i4,t26,d15.8)
      end
