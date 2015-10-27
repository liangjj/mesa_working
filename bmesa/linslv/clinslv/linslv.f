*deck linslv
      subroutine linslv (oper,status,a1,a2,a3,a4,v1,v2,v3,avdum,ia,n,
     1                   m,nmax,nvec,prdir,schmdt)
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
c***                              a * x = b
c***                                 or 
c***                             (i + a) * x  =  y
c***                             (  -  )
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
      complex*16 v1(m,*), v2(m,*), v3(n,*)
      complex*16 a1(n,*), a2(n,*), a3(m,*), a4(m,*)
      real*8  ovtol, cnverg, delmax, diff, avdum(7)
      complex*16 vv, cdotu, cdotc, sum
      complex*16 zeroc, onec, monec
      dimension ia(*), prdir(*), prnt(8)
      character *(*) oper, status
      character*32 slndir
      character*80 title
      logical schmdt, prdir, prnt
      save ovtol, cnverg, mxiter, iout, nrhs, nold, ntrial
      save cycle, slndir, prnt
      parameter (zeroc=(0.d0,0.d0), onec=(1.d0,0.d0), 
     1           monec=(-1.d0,0.d0))
c
      if (oper.eq.'initialize') then
c
c----------------------------------------------------------------------c
c         prepare for solution and save local copies of important      c
c                          variables                                   c
c----------------------------------------------------------------------c
c
c         threshold for convergence of residual.
c
          cnverg=avdum(1)
c
c         threshold for acceptance of a new vector as being linearly
c         independent of the other vectors.
c
          ovtol=avdum(1)
c
c         maximum number of iterations or vectors allowed.
c
          mxiter=ia(1)
c
c         the output unit.
c
          iout=m
c
c         the number of right hand sides.
c
          nrhs=nmax
c
c         the number of initial trial vectors.  this is taken as the number
c         right hand sides.
c
          ntrial=min(nrhs,n)
c
c         set the cycle equal to zero.
c
          cycle=0
c
c         initialize the form of the equations
c
          slndir=status
c            
c
c         print directives
c
          prnt(1)=prdir(1)          
          prnt(2)=prdir(2)          
          prnt(3)=prdir(3)          
          prnt(4)=prdir(4)          
          prnt(5)=prdir(5)          
          prnt(6)=prdir(6)          
          prnt(7)=prdir(7)          
          prnt(8)=prdir(8)          
c
          write (iout,1) mxiter, cnverg, ovtol, iout, nrhs, slndir
c
c----------------------------------------------------------------------c
c              form initial set of vectors from zeroth order solution  c
c----------------------------------------------------------------------c
c
          call cc2opy(a2,a1,n*ntrial)
c
c         initialize the number of vectors as zero
c
          nvec=0
c
c
      elseif (oper.eq.'generate basis functions') then
c
c----------------------------------------------------------------------c
c              schmidt orthonormalization of vectors                   c
c----------------------------------------------------------------------c
c           orthonormalize the input vectors, twice, if desired to 
c           ensure linear independence.
c
              nold=nvec
              nstart=nold+1
              nfin=nold+ntrial
              if (prnt(1).or.prnt(2)) then
                  title='initial vectors'
                  call prntcm(title,a1(1,nstart),n,ntrial,n,
     1                        ntrial,iout)
              endif
              call cschmt(a1,ovtol,n,nstart,nfin,nout,schmdt)
              if(nfin.eq.0) then
                 title='final vectors'
                 call prntcm(title,a1(1,nstart),n,ntrial,n,
     1                       ntrial,iout)
                 call lnkerr('number vectors is zero')
              endif
              nvec=nfin
              ntrial=nout
              if (prnt(1).or.prnt(2)) then
                  title='final vectors'
                  call prntcm(title,a1(1,nstart),n,ntrial,n,
     1                        ntrial,iout)
              endif
c
c             test if there are any vectors left for further calculation.
c
              status='ok'
              if (ntrial.eq.0) then
                  status='none left'
              endif    
              ia(1)=ntrial
c
c             print information if requested.
c
              if (prnt(1).or.prnt(3)) then
                  title='partial overlap matrix'
                  write(iout,2) title
                  do 50 i=nstart,nvec
                     do 60 j=1,i
                        vv=cdotc(n,a1(1,j),1,a1(1,i),1)
                        write(iout,3) i, j, vv
   60                continue
   50             continue
              endif
c----------------------------------------------------------------------c
c          here is where user should generate iterates by              c
c                    calling appropriate routine                       c
c          the actual operation is done in the calling routine. 
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c             now solve small set of linear equations                  c
c----------------------------------------------------------------------c
      elseif (oper.eq.'solve equations') then
c----------------------------------------------------------------------c
c             make copies of old matrix and rhs                        c
c----------------------------------------------------------------------c
              nstart=nold+1
              if (nold.ne.0) then
                  do 70 i=1,nold
                     do 80 j=1,nold
                        a3(i,j)=a4(i,j)
   80                continue
   70             continue
                  do 90 i=1,nold
                     do 100 j=1,nrhs
                        v1(i,j)=v2(i,j)
  100                continue
   90             continue
c----------------------------------------------------------------------c
c          calculate new matrix elements of small matrix               c
c----------------------------------------------------------------------c
                  do 110 i=1,nold
                     do 120 j=nstart,nvec
                        a3(i,j)=cdotc(n,a1(1,i),1,a2(1,j),1)
                        a3(j,i)=cdotc(n,a1(1,j),1,a2(1,i),1)
                        a4(i,j)=a3(i,j)
                        a4(j,i)=a3(j,i)
  120                continue
  110             continue
              endif
              do 130 i=nstart,nvec
                 do 140 j=1,nrhs
                    v1(i,j)=cdotc(n,a1(1,i),1,v3(1,j),1)
                    v2(i,j)=v1(i,j)
  140            continue
  130         continue
              do 150 i=nstart,nvec
                 do 160 j=nstart,nvec
                    a3(i,j)=cdotc(n,a1(1,i),1,a2(1,j),1)
                    a4(i,j)=a3(i,j)
  160            continue
  150         continue
              call czero(a3,m*nvec)
              do 170 i=1,nvec
                 a3(i,i)=onec
  170         continue
              if(slndir.eq.'one-plus-matrix') then
                 call cmadd(onec,a3,m,onec,a4,m,a3,m,nvec,nvec)
              endif
              if(slndir.eq.'one-minus-matrix') then
                 call cmadd(onec,a3,m,monec,a4,m,a3,m,nvec,nvec)
              endif                 
c
c             print the desired information, if requested.
c
              if (prnt(1).or.prnt(4)) then
                  title='small matrix'
                  call prntcm(title,a3,nvec,nvec,m,nvec,iout) 
                  title='small rhs'
                  call prntcm(title,v1,nvec,nrhs,m,nrhs,iout)
              endif
c----------------------------------------------------------------------c
c                solve small set of equations                          c
c----------------------------------------------------------------------c
             call cgefa (a3,m,nvec,ia,info)
             do 180 i=1,nrhs
                call cgesl (a3,m,nvec,ia,v1(1,i),0)
  180        continue
c----------------------------------------------------------------------c
c                   test for convergence                               c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c                   calculate residuals                                c
c----------------------------------------------------------------------c
             cycle=cycle+1
             if (prnt(1).or.prnt(5)) then
                 write (iout,4) cycle,nvec
             endif   
             delmax=0.d+00
             if(slndir.eq.'matrix') then
                do 190 irhs=1,nrhs
                   diff=0.d+00
                   do 200 i=1,n
                      sum=cdotu(nvec,a2(i,1),n,v1(1,irhs),1)
     1                               -v3(i,irhs)
                      diff=diff+sum*conjg(sum)
  200              continue
                   diff=sqrt(diff)
                   if (prnt(1).or.prnt(5)) then
                       write (iout,5) irhs,diff
                    endif
                    delmax=max(diff,delmax)
  190           continue     
             elseif (slndir.eq.'one-plus-matrix') then
                 do 210 irhs=1,nrhs
                    diff=0.d+00
                    do 220 i=1,n
                       sum=cdotu(nvec,a1(i,1),n,v1(1,irhs),1) +
     1                     cdotu(nvec,a2(i,1),n,v1(1,irhs),1) - 
     2                                  v3(i,irhs)
                       diff=diff+sum*conjg(sum)
  220               continue
                    diff=sqrt(diff)
                    if (prnt(1).or.prnt(5)) then
                        write (iout,5) irhs,diff
                    endif
                    delmax=max(diff,delmax)
  210            continue
             elseif (slndir.eq.'one-minus-matrix') then
                 do 230 irhs=1,nrhs
                    diff=0.d+00
                    do 240 i=1,n
                       sum=cdotu(nvec,a1(i,1),n,v1(1,irhs),1) -
     1                     cdotu(nvec,a2(i,1),n,v1(1,irhs),1) - 
     2                                  v3(i,irhs)
                       diff=diff+sum*conjg(sum)
  240               continue
                    diff=sqrt(diff)
                    if (prnt(1).or.prnt(5)) then
                        write (iout,5) irhs,diff
                    endif
                    delmax=max(diff,delmax)
  230            continue
             else
                 write(iout,*) slndir
                 call lnkerr('error in linslv for solution call')
             endif    
             if (delmax.gt.cnverg) then
c----------------------------------------------------------------------c
c                continue iterations or quit if converged              c
c----------------------------------------------------------------------c
                 ntrial=nvec-nstart+1
                 ntest=nvec+ntrial
                 if(ntest.gt.mxiter) then
                    ntrial=mxiter-nvec
                 endif
                 nmax=ntrial
                 status='continue iterations'
             else
                 status='converged'
                 call cgemm('n','n',n,nrhs,nvec,onec,a1,n,v1,m,zeroc,
     1                                                    v3,n)
                 if (prnt(1).or.prnt(7)) then
                     title='converged solutions'
                     call prntcm(title,v3,n,nrhs,n,nrhs,iout)  
                 endif
             endif
      elseif(oper.eq.'new trials') then
             call cc2opy(a2,a1,n*ntrial)
             if(prnt(1).or.prnt(6)) then
                title='new input vectors'
                call prntcm(title,a1,n,ntrial,n,ntrial,iout)
             endif   
      elseif (oper.eq.'final best solution') then
              call cgemm('n','n',n,nrhs,nvec,onec,a1,n,v1,m,zeroc,
     1                    v3,n)
              if (prnt(1).or.prnt(7)) then
                  title='best solutions'
                  call prntcm(title,v3,n,nrhs,n,nrhs,iout)  
              endif
              if(status.eq.'converged') then
                 return
              else
                 status='iteration limit exceeded'
              endif
      elseif (oper.eq.'cleanup') then
              cnverg=0.0d+00
              ovtol=0.0d+00
              mxiter=0
              iout=0
              nvec=0
              nrhs=0
              nold=0
      else
              write (iout,6)
      endif
c
      return
c
    1 format(/,5x,'initialize linslv:',/,5x,
     1            '      maximum number of iterations = ',i4,/,5x,
     2            '      convergence criterion        = ',e15.8,/,5x,
     3            '      overlap threshold            = ',e15.8,/,5x,
     4            '      write file                   = ',i2,/,5x,
     5            '      number of right hand sides   = ',i4,/,5x,
     6            '      form of equations            = ',a32)
    2 format(/,5x,a32)   
    3 format(5x,'i = ',i3,2x,'j = ',i3,2x,
     1         'overlap = ', e15.8,1x,e15.8)   
    4 format (/,1x,25('-'),'cycle ',i3,' using ',i4,' vectors ',10('-'),
     1 /,t3,'solution',t30,'convergence')
    5 format (t4,i4,t26,e15.8)                      
    6 format (//,15x,'*****  error in call to linslv  *****')
      end
