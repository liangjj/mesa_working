*deck @(#)genrtr.f	5.1  11/6/94
      subroutine genrtr(gengrp,naxis,gen,nirrep,ngen,gamma,nop,lambda,
     #                  genpt,lengen,char,t,lirrep)
c
      implicit integer (a-z)
c     this routine creates the generators of the various point groups,
c     and from that the irreducible representation matrices.
c     (in the cases of higher symmetry, the representation matrices
c     are generated explicitly without the use of generators)
c     the character table is also returned. note that the generators
c     or rotations correspond to a clockwise rotation about the z-axis.
c
      character*(*) gengrp
      real*8 gen(lengen,ngen),gamma(lengen,nop),char(nirrep,nop)
      real*8 t(3,3,nop)
      character*1 lab1,lab2,lab3
      character*8 label,lirrep(*)
      real*8 theta,pi,trace,m
      integer naxis,lambda(nirrep),genpt(nirrep)
      logical nodd
      common/io/inp,iout
c
c     ----- fill in pointers to gamma matrices -----
c
      n=1
      do 301 irrep=1,nirrep
         genpt(irrep)=n
         n=n+lambda(irrep)**2
 301  continue
c
c
      if (naxis.gt.0) then
         pi=4.0d+00*atan(1.0)
         nodd=mod(naxis,2).ne.0
         theta=2*pi/naxis
      end if
c
c     ----- the simple generators are tabulated -----
c
      if (gengrp.eq.'c1') then
         gen(1,1)=1.0d+00
      else if (gengrp.eq.'cs'.or.gengrp.eq.'ci') then
         gen(1,1)= 1.0d+00
         gen(2,1)=-1.0d+00
      end if
c
c     the cn groups.
      if (gengrp.eq.'cn') then
c        the sole generator is cn(z).      
         do 1 k=1,nirrep
            n=genpt(k)
            label=lirrep(k)
            lab1=label(1:1)
            lab2=label(2:2)
            lab3=label(3:3)
            if(lab1.eq.'a') then
               gen(n,1)= 1.0d+00
            else if(lab1.eq.'b') then
               gen(n,1)=-1.0d+00
            else if (lab1.eq.'e') then
               if(lab2.eq.' '.or.lab2.eq.'1') then
                  m=1.0d+00
               else
                  m=float(ctoi(lab2))
               end if
               gen(n  ,1)= cos(m*theta)
               gen(n+1,1)= sin(m*theta)
               gen(n+2,1)=-sin(m*theta)
               gen(n+3,1)= cos(m*theta)
            end if
    1    continue
      end if
c
c     the dn groups.
      if (gengrp.eq.'dn') then
c        the generators are cn(z), c2'(y).      
         do 11 k=1,nirrep
            n=genpt(k)
            label=lirrep(k)
            lab1=label(1:1)
            lab2=label(2:2)
            lab3=label(3:3)
            if(lab1.eq.'a') then
               gen(n,1)= 1.0d+00
               if(lab2.eq.' '.or.lab2.eq.'1') then
                  gen(n,2)= 1.0d+00
               else if(lab2.eq.'2') then
                  gen(n,2)=-1.0d+00
               end if
            else if(lab1.eq.'b') then
c              handle the special case of d2.
               if (naxis.eq.2) then
                  if(lab2.eq.'1') then
                     gen(n,1)= 1.0d+00
                     gen(n,2)=-1.0d+00
                  else if(lab2.eq.'2') then
                     gen(n,1)=-1.0d+00
                     gen(n,2)= 1.0d+00
                  else if(lab2.eq.'3') then
                     gen(n,1)=-1.0d+00
                     gen(n,2)=-1.0d+00 
                  end if
               else
c                 the general case.
                  gen(n,1)=-1.0d+00
                  if(lab2.eq.'1') then
                     gen(n,2)= 1.0d+00
                  else if(lab2.eq.'2') then 
                     gen(n,2)=-1.0d+00
                  end if         
               end if
            else if (lab1.eq.'e') then
               if(lab2.eq.' '.or.lab2.eq.'1') then
                  m=1.0d+00
               else
                  m=float(ctoi(lab2))
               end if
               gen(n  ,1)= cos(m*theta)
               gen(n+1,1)= sin(m*theta)
               gen(n+2,1)=-sin(m*theta)
               gen(n+3,1)= cos(m*theta)
c
               gen(n  ,2)=-1.0d+00
               gen(n+1,2)= 0.0d+00
               gen(n+2,2)= 0.0d+00
               gen(n+3,2)= 1.0d+00
            end if
   11    continue
      end if
c
c     the cnv groups.
      if (gengrp.eq.'cnv') then
c        the generators are cn(z), sigma(xz).
c        the group is isomorphic with dn.
         do 21 k=1,nirrep
            n=genpt(k)
            label=lirrep(k)
            lab1=label(1:1)
            lab2=label(2:2)
            lab3=label(3:3)
            if(lab1.eq.'a') then
               gen(n,1)= 1.0d+00
               if(lab2.eq.' '.or.lab2.eq.'1') then
                  gen(n,2)= 1.0d+00
               else if(lab2.eq.'2') then
                  gen(n,2)=-1.0d+00
               end if
            else if(lab1.eq.'b') then
               gen(n,1)=-1.0d+00
               if(lab2.eq.'1') then
                  gen(n,2)= 1.0d+00
               else if(lab2.eq.'2') then
                  gen(n,2)=-1.0d+00
               end if         
            else if (lab1.eq.'e') then
               if(lab2.eq.' '.or.lab2.eq.'1') then
                  m=1.0d+00
               else
                  m=float(ctoi(lab2))
               end if
               gen(n  ,1)= cos(m*theta)
               gen(n+1,1)= sin(m*theta)
               gen(n+2,1)=-sin(m*theta)
               gen(n+3,1)= cos(m*theta)
c
               gen(n  ,2)= 1.0d+00
               gen(n+1,2)= 0.0d+00
               gen(n+2,2)= 0.0d+00
               gen(n+3,2)=-1.0d+00
            end if
   21    continue
      end if
c
c     the cnh groups.
      if (gengrp.eq.'cnh') then
c        the generators are cn(z), sigmah(xy) if n is odd,      
c                                  inversion  if n is even.
         do 31 k=1,nirrep
            n=genpt(k)
            label=lirrep(k)
            lab1=label(1:1)
            lab2=label(2:2)
            lab3=label(3:3)
            if(lab1.eq.'a') then
               gen(n,1)= 1.0d+00
               if(lab3.eq.''''.or.lab3.eq.'g') then
                  gen(n,2)= 1.0d+00
               else if(lab3.eq.'"'.or.lab3.eq.'u') then
                  gen(n,2)=-1.0d+00
               end if
            else if(lab1.eq.'b') then
               gen(n,1)=-1.0d+00
               if(lab3.eq.''''.or.lab3.eq.'g') then
                  gen(n,2)= 1.0d+00
               else if(lab3.eq.'"'.or.lab3.eq.'u') then 
                  gen(n,2)=-1.0d+00
               end if         
            else if (lab1.eq.'e') then
               if(lab2.eq.' '.or.lab2.eq.'1') then
                  m=1.0d+00
               else
                  m=float(ctoi(lab2))
               end if
               gen(n  ,1)= cos(m*theta)
               gen(n+1,1)= sin(m*theta)
               gen(n+2,1)=-sin(m*theta)
               gen(n+3,1)= cos(m*theta)
c
               if(lab3.eq.''''.or.lab3.eq.'g') then
                  gen(n  ,2)= 1.0d+00
                  gen(n+1,2)= 0.0d+00
                  gen(n+2,2)= 0.0d+00
                  gen(n+3,2)= 1.0d+00
c
               else if(lab3.eq.'"'.or.lab3.eq.'u') then
                  gen(n  ,2)=-1.0d+00
                  gen(n+1,2)= 0.0d+00
                  gen(n+2,2)= 0.0d+00
                  gen(n+3,2)=-1.0d+00
               end if
            end if
   31    continue
      end if
c
c     the dnh groups.
      if (gengrp.eq.'dnh') then
c        the generators are cn(z), c2'(y);
c                                  sigmah(xy) if n is odd,      
c                                  inversion  if n is even.
         do 41 k=1,nirrep
            n=genpt(k)
            label=lirrep(k)
            lab1=label(1:1)
            lab2=label(2:2)
            lab3=label(3:3)
            if(lab1.eq.'a') then
               gen(n,1)= 1.0d+00
               if(lab2.eq.' '.or.lab2.eq.'1') then
                  gen(n,2)= 1.0d+00
               else if(lab2.eq.'2') then
                  gen(n,2)=-1.0d+00
               endif
               if(lab3.eq.''''.or.lab3.eq.'g') then
                  gen(n,3)= 1.0d+00
               else if(lab3.eq.'"'.or.lab3.eq.'u') then
                  gen(n,3)=-1.0d+00
               end if
            else if(lab1.eq.'b') then
c              handle the special case of d2h.
               if(naxis.eq.2) then
                  if(lab2.eq.'1') then
                     gen(n,1)= 1.0d+00
                     gen(n,2)=-1.0d+00
                  else if(lab2.eq.'2') then
                     gen(n,1)=-1.0d+00
                     gen(n,2)= 1.0d+00
                  else if(lab2.eq.'3') then
                     gen(n,1)=-1.0d+00
                     gen(n,2)=-1.0d+00
                  end if
               else
c                 general case.
                  gen(n,1)=-1.0d+00
                  if(lab2.eq.'1') then
                     gen(n,2)= 1.0d+00
                  else if(lab2.eq.'2') then 
                     gen(n,2)=-1.0d+00
                  end if
               end if
               if(lab3.eq.''''.or.lab3.eq.'g') then
                  gen(n,3)= 1.0d+00
               else if(lab3.eq.'"'.or.lab3.eq.'u') then 
                  gen(n,3)=-1.0d+00
               end if         
            else if (lab1.eq.'e') then
               if(lab2.eq.' '.or.lab2.eq.'1') then
                  m=1.0d+00
               else
                  m=float(ctoi(lab2))
               end if
               gen(n  ,1)= cos(m*theta)
               gen(n+1,1)= sin(m*theta)
               gen(n+2,1)=-sin(m*theta)
               gen(n+3,1)= cos(m*theta)
c
               gen(n,  2)=-1.0d+00
               gen(n+1,2)= 0.0d+00
               gen(n+2,2)= 0.0d+00
               gen(n+3,2)= 1.0d+00
c
               if(lab3.eq.''''.or.lab3.eq.'g') then
                  gen(n  ,3)= 1.0d+00
                  gen(n+1,3)= 0.0d+00
                  gen(n+2,3)= 0.0d+00
                  gen(n+3,3)= 1.0d+00
c
               else if(lab3.eq.'"'.or.lab3.eq.'u') then
                  gen(n  ,3)=-1.0d+00
                  gen(n+1,3)= 0.0d+00
                  gen(n+2,3)= 0.0d+00
                  gen(n+3,3)=-1.0d+00
               end if
            end if
   41    continue
      end if
c
c     the dnd groups.
      if (gengrp.eq.'dnd') then
c        the generators are s2n(z), c2'(y)            if n is even,      
c                            cn(z), c2'(y),inversion  if n is odd.
         do 51 k=1,nirrep
            n=genpt(k)
            label=lirrep(k)
            lab1=label(1:1)
            lab2=label(2:2)
            lab3=label(3:3)
            if(lab1.eq.'a') then
               gen(n,1)= 1.0d+00
               if(nodd) then
                  if(lab3.eq.'g') then
                     gen(n,2)= 1.0d+00
                  else if(lab3.eq.'u') then
                     gen(n,2)=-1.0d+00
                  end if
               else
                  if(lab2.eq.'1') then
                     gen(n,2)= 1.0d+00
                  else if(lab2.eq.'2') then
                     gen(n,2)=-1.0d+00
                  end if
               end if
            else if(lab1.eq.'b') then
               gen(n,1)=-1.0d+00
               if(lab2.eq.'1') then
                  gen(n,2)= 1.0d+00
               else if(lab2.eq.'2') then 
                  gen(n,2)=-1.0d+00
               end if         
            else if (lab1.eq.'e') then
               if(lab2.eq.' '.or.lab2.eq.'1') then
                  m=1.0d+00
               else
                  m=float(ctoi(lab2))
               end if
               if(nodd) then
c                 operator is cn(z)
                  gen(n  ,1)= cos(m*theta)
                  gen(n+1,1)= sin(m*theta)
                  gen(n+2,1)=-sin(m*theta)
                  gen(n+3,1)= cos(m*theta)
c
c                 operator 2 is c2'(y)
                  gen(n  ,2)=-1.0d+00
                  gen(n+1,2)= 0.0d+00
                  gen(n+2,2)= 0.0d+00
                  gen(n+3,2)= 1.0d+00
c
c                 operator is inversion.
                  if(lab3.eq.'g') then
                     gen(n  ,3)= 1.0d+00
                     gen(n+1,3)= 0.0d+00
                     gen(n+2,3)= 0.0d+00
                     gen(n+3,3)= 1.0d+00
c
                  else if(lab3.eq.'u') then
                     gen(n  ,3)=-1.0d+00
                     gen(n+1,3)= 0.0d+00
                     gen(n+2,3)= 0.0d+00
                     gen(n+3,3)=-1.0d+00
                  end if
               else
c                 operator is s2n; rotation angle is (m/2)*theta
                  m=0.5d+00*m
                  gen(n  ,1)=-cos(m*theta)
                  gen(n+1,1)=-sin(m*theta)
                  gen(n+2,1)=+sin(m*theta)
                  gen(n+3,1)=-cos(m*theta)
c
c                 apply c2'(y).
                  gen(n  ,2)=-1.0d+00
                  gen(n+1,2)= 0.0d+00
                  gen(n+2,2)= 0.0d+00
                  gen(n+3,2)= 1.0d+00
               endif
c
            end if
   51    continue
      end if
c
c     the sn groups.
      if (gengrp.eq.'sn') then
c        the sole generator is sn(z).      
c        proper rotation by theta, followed by reflection in 
c        perpendicular plane.
         if (nodd) call lnkerr('only s2, s4, and s6 groups allowed')
         do 61 k=1,nirrep
            n=genpt(k)
            label=lirrep(k)
            lab1=label(1:1)
            lab2=label(2:2)
            lab3=label(3:3)
            if(lab1.eq.'a') then
               gen(n,1)= 1.0d+00
            else if(lab1.eq.'b') then
               gen(n,1)=-1.0d+00
            else if (lab1.eq.'e') then
               if(lab2.eq.' ') then
                  m=1.0d+00
               else
                  m=float(ctoi(lab2))
               end if
               gen(n  ,1)=-cos(m*theta)
               gen(n+1,1)=-sin(m*theta)
               gen(n+2,1)=+sin(m*theta)
               gen(n+3,1)=-cos(m*theta)
            end if
   61    continue
      end if
c
c     ----- work out the representation matrices for the simple and axial
c     groups
c
      if (gengrp.eq.'c1'.or.gengrp.eq.'cs'.or.gengrp.eq.'ci'.or.
     #     gengrp.eq.'cn'.or.gengrp.eq.'cnv'.or.gengrp.eq.'cnh'.or.
     #     gengrp.eq.'dn'.or.gengrp.eq.'dnd'.or.gengrp.eq.'dnh'.or.
     #     gengrp.eq.'sn') then
c
         if (gengrp.eq.'c1') then
            nmax=1
         else if (gengrp.eq.'cs'.or.gengrp.eq.'ci') then
            nmax=2
         else
            nmax=naxis
         end if
c
         do 100 irrep=1,nirrep
            n=genpt(irrep)
            l=lambda(irrep)
            call runit(gamma(n,1),l)
            do 20 op=2,nmax
               call ebc(gamma(n,op),gamma(n,op-1),gen(n,1),l,l,l)
 20         continue
            if (ngen.ge.2) then
               do 30 op=nmax+1,2*nmax
                  call ebc(gamma(n,op),gamma(n,op-nmax),
     #                 gen(n,2),l,l,l)
 30            continue
            end if
            if (ngen.eq.3) then
               do 40 op=2*nmax+1,4*nmax
                  call ebc(gamma(n,op),gamma(n,op-2*nmax),
     #                 gen(n,3),l,l,l)
 40            continue
            end if
 100     continue
c
      end if
c
c     --- if this is a higher symmetry group, work out the gamma
c         matrices a little differently
      if (gengrp.eq.'t'.or.gengrp.eq.'th'.or.gengrp.eq.'td'
     $    .or.gengrp.eq.'o'.or.gengrp.eq.'oh'
     $    .or.gengrp.eq.'ih') then
         if(gengrp.eq.'t'.or.gengrp.eq.'th'.or.gengrp.eq.'td') then
            call lnkerr('tetrahedral groups not yet implemented')
         else if(gengrp.eq.'o'.or.gengrp.eq.'oh') then
c           --- work out representation matrices for o ---
c               see the routine tmat for the order in which the operators
c               are generated. oh is generated from o by applying the 
c               inversion.
            do 200 irrep=1,nirrep
               n=genpt(irrep)
               l=lambda(irrep)
               label=lirrep(irrep)
               lab1=label(1:1)
               lab2=label(2:2)
               lab3=label(3:3)
c              set up the generator for inversion if oh
               if(lab3.eq.'g') then
                  call runit(gen(n,1),l)
               else if(lab3.eq.'u') then
                  call runit(gen(n,1),l)
                  call vneg(gen(n,1),gen(n,1),l*l)
               endif
c              the identity
               call runit(gamma(n,1),l)
               if(lab1.eq.'a') then
                  if(lab2.eq.'1') then
                     do 120 op=2,24
                        gamma(n,op)=1.0d+00
  120                continue
                  else if(lab2.eq.'2') then
c                    the eight c3's
                     do 130 op=2,9
                        gamma(n,op)=1.0d+00
  130                continue
c                    the 9 c4's
                     do 140 i=0,2
                        gamma(n,10+3*i)=-1.0d+00
                        gamma(n,11+3*i)= 1.0d+00
                        gamma(n,12+3*i)=-1.0d+00
  140                continue
c                    the 6 c2's
                     do 150 op=19,24
                        gamma(n,op)=-1.0d+00
  150                continue
                  endif
               else if(lab1.eq.'e') then
c
c as with t2, e is generated from the (x,y,z) transformation matrices
c
                  do 155 op=2,24
                     call oe(gamma(n,op),t(1,1,op))
 155              continue 
               else if(lab1.eq.'t') then
                  if(lab2.eq.'1') then
c                    t1 is represented by the coordinate transformations
c                    which we've already done.
                     do 170 op=2,24
                        call vmove(gamma(n,op),t(1,1,op),9)
  170                continue
                  else if(lab2.eq.'2') then
c
c t2 representations are generated from the t1 guys by a nonlinear change
c     of basis, similar to that done in trmat
c
                     do 175 op=2,24
                        call ot2(gamma(n,op),t(1,1,op))
 175                 continue 
                  endif
               endif
c              generate the last 24 via the inversion
               if(gengrp.eq.'oh') then
                  do 180 op=1,24
                     call ebc(gamma(n,op+24),gamma(n,op),gen(n,1),l,l,l)
  180             continue
               endif
  200       continue
         endif
      endif
                     
c
c     ------ sum the representations down to character tables ------
c
      do 250 irrep=1,nirrep
         n=genpt(irrep)
         l=lambda(irrep)
c        remove the blank entries that characterize irrep labels, i.e. "a g"
         call rmvnb(lirrep(irrep),lirrep(irrep))
         do 240 op=1,nop
            char(irrep,op)=trace(gamma(n,op),l)
 240     continue
 250  continue
c
c
      return
      end
