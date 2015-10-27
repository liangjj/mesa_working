*deck  @(#)prwalk.f	1.3 7/30/91
      subroutine prwalk(ops,coef,walk,brkdn,arc,weight,nimprt,nrows,
     #                  orbtbf,norbs,ref,alpha,beta,bfsym,symnum,
     #                  nsym,nbf,symlab,iout,srow)
c
c***begin prologue prwalk
c***date written   851003   (yymmdd)
c***revision date  yymmdd   (yymmdd)
c***keywords  configuration printing, print, breakdown
c
c***author  saxe, paul,    (lanl)
c***purpose  to print the most important configurations in a ci vector.
c
c***description   prwalk will print the orbital occupancies of the most
c        important configuration, and then the excitations between
c        other important configurations and the most important
c        configuration. to do this, prwalk needs the configuration
c        numbers and some drt information.
c
c        on input:
c
c           coef     real (nimprt)
c                    the coefficients in the ci vector of the
c                    configurations to print.
c
c           walk     integer (nimprt)
c                    the configuration numbers of the configurations to
c                    print. the first is assumed to be the reference.
c
c           arc      integer (4,nrows)
c                    the arc array for the drt.
c
c           weight   integer (4,nrows)
c                    the weights of arcs in the drt.
c
c           nimprt   integer
c                    the number of configurations to print.
c
c           nrows    integer
c                    the number of rows in the drt.
c
c           orbtbf   integer (norbs)
c                    array giving number of scf orbital corresponding
c                    to each ci orbital.
c
c           norbs    integer
c                    number of orbitals in the ci.
c
c           bfsym    integer (nbf)
c                    symmetries of the scf orbitals, counting from 0.
c
c           nsym     integer
c                    number of irreducible representations of the
c                    point group actually used in the ci.
c
c           nbf      integer
c                    the number of scf orbitals.
c
c           symlab   character*3 (nsym)
c                    the labels for the irreducible representations.
c
c           iout     integer
c                    the output unit number or name.
c
c        scratch:
c
c           brkdn    integer (norbs)
c           ref      integer (norbs)
c           alpha    integer (norbs)
c           beta     integer (norbs)
c           symnum   integer (nbf)
c
c***references
c
c***routines called  rmvnb (char), lnkerr (mdutil)
c***end prologue  prwalk
c
      implicit integer (a-z)
c
      character*(*) ops
      real*8 coef(nimprt)
      integer walk(nimprt),brkdn(norbs),arc(4,nrows),weight(4,nrows)
      integer orbtbf(norbs),ref(norbs),alpha(norbs),beta(norbs)
      integer bfsym(nbf),symnum(nbf)
      integer alp(4,4),bet(4,4)
      character*3 lets(4),symlab(0:*)*3
      real*8 sum
      logical full
      logical logkey
c
      data alp / 0, 1, 0, 1,
     #          -1, 0,-1, 0,
     #           0, 1, 0, 1,
     #          -1, 0,-1, 0/
c
      data bet / 0, 0, 1, 1,
     #           0, 0, 1, 1,
     #          -1,-1, 0, 0,
     #          -1,-1, 0, 0/
c
      data lets /'  ','/ ','\\ ','x '/
c
      full=logkey(ops,'print=walks=full',.false.,' ')
c
c     ----- set up symnum: the number of each bf in a symmetry block -
c
      do 200 i=1,nsym
         pt=0
         do 199 j=1,nbf
            if (bfsym(j)+1.eq.i) then
               pt=pt+1
               symnum(j)=pt
            end if
  199    continue
  200 continue
c
      sum=0.0d+00
      do 100 imprt=1,nimprt
         left=walk(imprt)-1
         row=srow
         do 3 orb=norbs,1,-1
            do 1 case=4,1,-1
               if (weight(case,row).le.left
     #             .and.arc(case,row).gt.0) go to 2
    1       continue
            call lnkerr('19')
    2       continue
            left=left-weight(case,row)
            brkdn(orb)=case
            row=arc(case,row)
    3    continue
c
c        ----- excitations from the reference -----
c
         sum=sum+coef(imprt)**2
         if (imprt.eq.1.or.full) then
            write (iout,4) imprt,coef(imprt),sqrt(sum),walk(imprt)
            pos=32
            i=norbs
  110       continue
               start=i
               min=symnum(orbtbf(i))
               sym=bfsym(orbtbf(i))
               occ=brkdn(i)
  120          continue
                  i=i-1
                  if (i.lt.1) go to 125
               if (occ.eq.brkdn(i).and.sym.eq.bfsym(orbtbf(i)))
     #                go to 120
c
  125          continue
               if (i.eq.start-1) then
                  write (iout,130) min,symlab(sym),lets(occ)
  130             format (2x,i4,a3,a3)
               else
                  write (iout,140) min,symnum(orbtbf(i+1)),
     #                               symlab(sym),lets(occ)
  140             format (2x,i4,'-',i4,a3,a3)
               end if
            if (i.ge.1) go to 110
c
c
            do 10 i=1,norbs
               ref(i)=brkdn(i)
   10       continue
         else
            write (iout,4) imprt,coef(imprt),sqrt(sum),walk(imprt)
    4       format (1x,i3,2f9.4,i7)
            pos=32
            do 20 i=1,norbs
               alpha(i)=alp(ref(i),brkdn(i))
               beta(i)=bet(ref(i),brkdn(i))
   20       continue
c
c           ----- check for ab-->ab excitations -----
c
            do 30 i=norbs,1,-1
               if (alpha(i).ne.1.or.beta(i).ne.1) go to 30
               do 21 j=norbs,1,-1
                  if (alpha(j).eq.-1.and.beta(j).eq.-1) go to 22
   21          continue
               go to 30
   22          continue
               write (iout,23) symnum(orbtbf(i)),
     #                        symlab(bfsym(orbtbf(i))),
     #                      symnum(orbtbf(j)),symlab(bfsym(orbtbf(j)))
   23          format (5x,i3,a3,'x-->',i3,a3,'x')
               alpha(i)=0
               alpha(j)=0
               beta(i)=0
               beta(j)=0
   30       continue
c
c           ----- ab --> a, b excitations -----
c
            do 40 i=norbs,1,-1
               if (alpha(i).ne.1.or.beta(i).ne.1) go to 40
               do 31 j=norbs,1,-1
                  if (alpha(j).eq.-1) go to 32
   31          continue
               call lnkerr('31')
   32          continue
               do 33 k=norbs,1,-1
                  if (beta(k).eq.-1) go to 34
   33          continue
               call lnkerr('33')
   34          continue
               write (iout,35) symnum(orbtbf(i)),
     #                      symlab(bfsym(orbtbf(i))),
     #                      symnum(orbtbf(j)),symlab(bfsym(orbtbf(j))),
     #                      symnum(orbtbf(k)),symlab(bfsym(orbtbf(k)))
   35          format (5x,i3,a3,'x-->',i3,a3,'/  ',i3,a3,'\ ')
               alpha(i)=0
               alpha(j)=0
               beta(i)=0
               beta(k)=0
   40       continue
c
c           ----- the rest -----
c
            do 50 i=norbs,1,-1
               if (alpha(i).ne.1) go to 50
               do 41 j=norbs,1,-1
                  if (alpha(j).eq.-1) go to 42
   41          continue
               call lnkerr('41')
   42          continue
               write (iout,43) symnum(orbtbf(i)),
     #                         symlab(bfsym(orbtbf(i))),
     #                      symnum(orbtbf(j)),symlab(bfsym(orbtbf(j)))
   43          format(i3,a3,'/-->',i3,a3,'/')
               alpha(i)=0
               alpha(j)=0
   50       continue
            do 60 i=norbs,1,-1
               if (beta(i).ne.1) go to 60
               do 51 j=norbs,1,-1
                  if (beta(j).eq.-1) go to 52
   51          continue
               call lnkerr('51')
   52          continue
               write (iout,53) symnum(orbtbf(i)),
     #                         symlab(bfsym(orbtbf(i))),
     #                      symnum(orbtbf(j)),symlab(bfsym(orbtbf(j)))
   53          format(5x,i3,a3,'\-->',i3,a3,'\ ')
               beta(i)=0
               beta(j)=0
   60       continue
         end if
c
c
         if (left.ne.0) call lnkerr('20')
  100 continue
c
c
      return
      end
