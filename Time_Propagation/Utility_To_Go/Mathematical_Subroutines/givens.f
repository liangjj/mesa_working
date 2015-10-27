*deck @(#)givens.f	5.1  11/6/94
      subroutine givens(nx,nrootx,njx,a,b,root,vect)
c***begin prologue     givens
c***date written       670901  (yymmdd)
c***revision date      920427  (yymmdd)  
c   april 27, 1992     rlm at lanl
c                      replacing constants with parameter statements
c   august 1, 1972     mel levy(tulane university).
c***keywords           matrix diagonalization, givens
c***author             prosser, franklin (indiana university).
c***source             @(#)givens.f	5.1   11/6/94
c***purpose            eigenvalues and eigenvectors of a real symmetric
c                      matrix by the givens method.
c***description
c                      call givens(nx,nrootx,njx,a,b,root,vect)
c                        nx       order of the matrix.
c                        nrootx   the nrootx most negative eigenvalues are
c                                 returned. if nrootx is negative, then no
c                                 vectors are returned.
c                        njx      row dimension of vect array (njx.ge.nx).
c                        a        input matrix of order nx packed upper
c                                 triangular.
c                        b        scratch array dimensioned at least (5*nx).
c                        root     output eigenvalues at least nrootx cells
c                                 long.  the nrootx most negative roots are
c                                 ordered largest first in this array.
c                        vect     output eigenvectors dimensioned at least
c                                 (njx,nrootx).  each column holds a
c                                 normalized eigenvector.  if nrootx is
c                                 negative, this is just a dummy argument.
c
c                      the arrays a and b are destroyed by the computation.
c                      the sequence of events is:
c                      1) reduce the matrix to tridiagonal form by the
c                         householder technique,
c                      2) locate the roots by the sturm sequence method,
c                      3) evaluate the vectors of the tridiagonal form,
c                      4) rotate the vectors back to those of the original
c                         array.
c
c***references
c                      j.m.ortega,"mathematics for digital computers", vol. 2,
c                      ralston and wilf,eds. wiley(1967), p.94.
c***routines called    abs,max,sqrt,sign,mod
c***end prologue       givens
c
c
cvax  extended dummy a,b,root,vect
c
c
c     real version by mel levy 8/72
c 62.3  givens  -eigenvalues and eigenvectors by the givens method.
c     by franklin prosser, indiana university.
c     september, 1967
c     calculates eigenvalues and eigenvectors of real symmetric matrix
c     stored in packed upper triangular form.
c
c     thanks are due to f. e. harris (stanford university) and h. h.
c     michels (united aircraft research laboratories) for excellent
c     work on numerical difficulties with earlier versions of this
c     program.
c
c     the parameters for the routine are...
c         nx     order of matrix
c         nrootx number of roots wanted.  the nrootx smallest (most
c                 negative) roots will be calculated.  if no vectors
c                 are wanted, make this number negative.
c         njx    row dimension of vect array.  see -vect- below.
c                 njx must be not less than nx.
c         a      matrix stored by columns in packed upper triangular
c                form, i.e. occupying nx*(nx+1)/2 consecutive
c                locations.
c         b      scratch array used by givens.  must be at least
c                 nx*5 cells.
c         root   array to hold the eigenvalues.  must be at least
c                nrootx cells long.  the nrootx smallest roots are
c                 ordered largest first in this array.
c         vect   eigenvector array.  each column will hold an
c                 eigenvector for the corresponding root.  must be
c                 dimensioned with -njx- rows and at least -nrootx-
c                 columns, unless no vectors
c                 are requested (negative nrootx).  in this latter
c                 case, the argument vect is just a dummy, and the
c                 storage is not used.
c                 the eigenvectors are normalized to unity.
c
c     the arrays a and b are destroyed by the computation.  the results
c     appear in root and vect.
c     for proper functioning of this routine, the result of a floating
c     point underflow should be a zero.
c     to convert this routine to real (e.g. on ibm 360
c     machines), be sure that all real variables and function
c     references are properly made real.
c     the value of -eta- (see below) should also be changed, to reflect
c     the increased precision.
c
c     the original reference to the givens technique is in oak ridge
c     report number ornl 1574 (physics), by wallace givens.
c     the method as presented in this program consists of four steps,
c     all modifications of the original method...
c     first, the input matrix is reduced to tridiagonal form by the
c     householder technique (j. h. wilkinson, comp. j. 3, 23 (1960)).
c     the roots are then located by the sturm sequence method (j. m.
c     ortega (see reference below).  the vectors of the tridiagonal
c     form are then evaluated (j. h. wilkinson, comp. j. 1, 90 (1958)),
c     and last the tridiagonal vectors are rotated to vectors of the
c     original array (first reference).
c     vectors for degenerate (or near-degenerate) roots are forced
c     to be orthogonal, using a method suggested by b. garbow, argonne
c     national labs (private communication, 1964).  the gram-schmidt
c     process is used for the orthogonalization.
c
c     an excellent presentation of the givens technique is found in
c     j. m. ortega-s article in -mathematics for digital computers,-
c     volume 2, ed. by ralston and wilf, wiley (1967), page 94.
c
      real*8 b(nx,5),a(1),root(nrootx),vect(njx,nrootx)
      real*8 eta,theta,del1,delta,small,delbig,theta1,toler
      real*8 rpower,rpow1,rand1,factor,anorm,alimit,sum,temp
      real*8 ak,rootl,rootx,trial,f0,aroot
      real*8 elim1,elim2
      real*8 zero,one,two,three,half,four99
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00,three=3.0d+00)
      parameter (half=0.5d+00,four99=4099.0d+00)
c
c ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c **  users please note...
c **  the following two parameters, eta and theta, should be adjusted
c **  by the user for his particular machine.
c **  eta is an indication of the precision of the floating point
c **  representation on the computer being used (roughly 10**(-m),
c **  where m is the number of decimals of precision ).
c **  theta is an indication of the range of numbers that can be
c **  expressed in the floating point representation (roughly the
c **  largest number).
c **  some recommended values follow.
c **  for ibm 7094, univac 1108, etc. (27-bit binary fraction, 8-bit
c **  binary exponent), eta=1.d-8, theta=1.d37.
c **  for control data 3600 (36-bit binary fraction, 11-bit binary
c **  exponent), eta=1.d-11, theta=1.d307.
c **  for control data 6600 (48-bit binary fraction, 11-bit binary
c **  exponent), eta=1.d-14, theta=1.d307.
c **  for ibm 360/50 and 360/65 real (56-bit hexadecimal
c **  fraction, 7-bit hexadecimal exponent), eta=1.d-16, theta=1.d75.
c **
      theta = 1.d38
      eta = 1.d-16
c ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     print 1,loc(a),loc(b),loc(root),loc(vect)
c   1 format (' givenx',4(1x,z8))
c
      del1=eta/100.0d+00
      delta=eta**2*100.0d+00
      small=eta**2/100.0d+00
      delbig=theta*delta/1000.0d+00
      theta1=1000.0d+00/theta
c     toler  is a factor used to determine if two roots are close
c     enough to be considered degenerate for purposes of orthogonali-
c     zing their vectors.  for the matrix normed to unity, if the
c     difference between two roots is less than toler, then
c     orthogonalization will occur.
      toler = eta*100.
c
c     initial value for pseudorandom number generator... (2**23)-3
      rpower=8388608.0d+00
      rpow1 = rpower/two
      rand1=rpower-three
c
      n = nx
      nroot = iabs(nrootx)
      if (nroot.eq.0) go to 1001
      if (n-1) 1001,1003,105
1003  root(1) = a(1)
      if(nrootx .gt. 0)vect(1,1)=one
      go to 1001
105   continue
c     nsize    number of elements in the packed array
      nsize = (n*(n+1))/2
      nm1 = n-1
      nm2 = n-2
c
c     scale matrix to euclidean norm of 1.  scale factor is anorm.
      factor=zero
      do 70 i=1,nsize
 70   factor=max(factor,abs(a(i)))
      if(factor .ne. zero)go to 72
c     null matrix.  fix up roots and vectors, then exit.
      do 78 i=1,nroot
           if (nrootx.lt.0) go to 78
           do 77 j=1,n
 77        vect(j,i)=zero
           vect(i,i)=one
 78   root(i)=zero
      go to 1001
c
 72   anorm=zero
      j = 1
      k = 1
      do 80 i=1,nsize
           if (i.ne.j) go to 81
           anorm=anorm+(a(i)/factor)**2/two
           k = k+1
           j = j+k
           go to 80
81         anorm = anorm + (a(i)/factor)**2
80    continue
      anorm=sqrt(anorm*two)*factor
      do 91 i=1,nsize
91    a(i) = a(i)/anorm
      alimit=one
c
c     tridia section.
c     tridiagonalization of symmetric matrix
      id = 0
      ia = 1
      if (nm2.eq.0) go to 201
      do 200  j=1,nm2
c     j       counts row  of a-matrix to be diagonalized
c     ia      start of non-codiagonal elements in the row
c     id      index of codiagonal element on row being codiagonalized.
           ia = ia+j+2
           id = id + j + 1
           jp2 = j+2
c     sum squares of non-codiagonal elements in row j
           ii = ia
           sum=zero
           do 100 i=jp2,n
                sum = sum + a(ii)**2
100        ii = ii + i
           temp = a(id)
           if (sum.gt.small) go to 110
c     no transformation necessary if all the non-codiagonal
c     elements are tiny.
           b(j,1) = temp
           a(id)=zero
           go to 200
c     now complete the sum of off-diagonal squares
 110       sum=sqrt(sum+temp**2)
c     new codiagonal element
           b(j,1)=-sign(sum,temp)
c     first non-zero element of this w-vector
           b(j+1,2)=sqrt((one+abs(temp)/sum)/two)
c     form rest of the w-vector elements
           temp=sign(half/(b(j+1,2)*sum),temp)
           ii = ia
           do 130 i=jp2,n
                b(i,2) = a(ii)*temp
130        ii = ii + i
c     form p-vector and scalar.  p-vector = a-matrix*w-vector.
c     scalar = w-vector*p-vector.
           ak=zero
c     ic      location of next diagonal element
           ic = id + 1
           j1 = j + 1
           do 190  i=j1,n
                jj = ic
                temp=zero
                do 180  ii=j1,n
c     i       runs over the non-zero p-elements
c     ii      runs over elements of w-vector
                     temp = temp + b(ii,2)*a(jj)
c     change incrementing mode at the diagonal elements.
                     if (ii.lt.i) go to 210
                     jj = jj + ii
                     go to 180
210                  jj = jj + 1
180             continue
c     build up the k-scalar (ak)
                ak = ak + temp*b(i,2)
                b(i,1) = temp
c     move ic to top of next a-matrix -row-
190        ic = ic + i
c     form the q-vector
           do 150  i=j1,n
150        b(i,1) = b(i,1) - ak*b(i,2)
c     transform the rest of the a-matrix
c     jj      start-1 of the rest of the a-matrix
           jj = id
c     move w-vector into the old a-matrix locations to save space
c     i       runs over the significant elements of the w-vector
           do 160  i=j1,n
                a(jj) = b(i,2)
                do 170  ii=j1,i
                     jj = jj + 1
 170            a(jj)=a(jj)-two*(b(i,1)*b(ii,2)+b(i,2)*b(ii,1))
160        jj = jj + j
200   continue
c     move last codiagonal element out into its proper place
201   continue
      b(nm1,1) = a(nsize-1)
      a(nsize-1)=zero
c
c     sturm section.
c     sturm sequence iteration to obtain roots of tridiagonal form.
c     move diagonal elements into second n elements of b-vector.
c     this is a more convenient indexing position.
c     also, put square of codiagonal elements in third n elements.
      jump=1
      do 320 j=1,n
           b(j,2)=a(jump)
           b(j,3) = b(j,1)**2
320   jump = jump+j+1
      do 310 i=1,nroot
310   root(i) = +alimit
      rootl = -alimit
c     isolate the roots.  the nroot lowest roots are found, lowest first
      do 330 i=1,nroot
c     find current -best- upper bound
           rootx = +alimit
           do 335 j=i,nroot
 335       rootx=min(rootx,root(j))
           root(i) = rootx
c     get improved trial root
 500       trial=(rootl+root(i))*0.5d+00
           if (trial.eq.rootl.or.trial.eq.root(i)) go to 330
c     form sturm sequence ratios, using ortega-s algorithm (modified).
c     nomtch is the number of roots less than the trial value.
           nomtch=n
           j=1
360        f0 = b(j,2) - trial
370        continue
           if(abs(f0) .lt. theta1)go to 380
           if(f0 .ge. zero)nomtch=nomtch-1
           j = j + 1
           if (j.gt.n) go to 390
c     since matrix is normed to unity, magnitude of b(j,3) is less than
c     one, and since f0 is greater than theta1, overflow is not possible
c     at the division step.
           f0 = b(j,2) - trial - b(j-1,3)/f0
           go to 370
380        j = j + 2
           nomtch = nomtch - 1
           if (j.le.n) go to 360
390        continue
c     fix new bounds on roots
           if (nomtch.ge.i) go to 540
           rootl = trial
           go to 500
540        root(i) = trial
           nom = min(nroot,nomtch)
           root(nom) = trial
           go to 500
330   continue
c     reverse the order of the eigenvalues, since custom dictates
c     -largest first-.  this section may be removed if desired without
c     affecting the remainder of the routine.
c     nrt = nroot/2
c     do 10 i=1,nrt
c     save = root(i)
c     nmip1 = nroot - i + 1
cc    root(i) = root(nmip1)
c10   root(nmip1) = save
c
c     trivec section.
c     eigenvectors of codiagonal form
c807  continue
c     quit now if no vectors were requested.
      if (nrootx.lt.0) go to 1002
c     initialize vector array.
      do 15 i=1,n
           do 15 j=1,nroot
 15   vect(i,j)=one
      do 700 i=1,nroot
           aroot = root(i)
c     orthogonalize if roots are close.
           if (i.eq.1) go to 710
c     the absolute value in the next test is to assure that the trivec
c     section is independent of the order of the eigenvalues.
           if(abs(root(i-1)-aroot) .lt. toler)go to 720
710        ia = -1
720        ia = ia + 1
           elim1 = a(1) - aroot
           elim2 = b(1,1)
           jump = 1
           do 750  j=1,nm1
                jump = jump+j+1
c     get the correct pivot equation for this step.
                if(abs(elim1) .le. abs(b(j,1)))go to 760
c     first (elim1) equation is the pivot this time.  case 1.
                b(j,2) = elim1
                b(j,3) = elim2
                b(j,4)=zero
                temp = b(j,1)/elim1
                elim1 = a(jump) - aroot - temp*elim2
                elim2 = b(j+1,1)
                go to 755
c     second equation is the pivot this time.  case 2.
760             b(j,2) = b(j,1)
                b(j,3) = a(jump) - aroot
                b(j,4) = b(j+1,1)
                temp=one
                if(abs(b(j,1)) .gt. theta1)temp=elim1/b(j,1)
                elim1 = elim2 - temp*b(j,3)
                elim2 = -temp*b(j+1,1)
c     save factor for the second iteration.
755             b(j,5) = temp
750        continue
           b(n,2) = elim1
           b(n,3)=zero
           b(n,4)=zero
           b(nm1,4)=zero
           iter = 1
           if (ia.ne.0) go to 801
c     back substitute to get this vector.
790        l = n + 1
           do 780 j=1,n
                l = l - 1
786             continue
                elim1=vect(l,i)-vect(l+1,i)*b(l,3)-vect(l+2,i)*b(l,4)
c     if overflow is conceivable, scale the vector down.
c     this approach is used to avoid machine-dependent and system-
c     dependent calls to overflow routines.
                if(abs(elim1) .gt. delbig)go to 782
                temp = b(l,2)
                if(abs(b(l,2)) .lt. delta)temp=delta
                vect(l,i) = elim1/temp
                go to 780
c     vector is too big.  scale it down.
782             do 784 k=1,n
784             vect(k,i) = vect(k,i)/delbig
                go to 786
780        continue
           go to (820,800), iter
c     second iteration.  (both iterations for repeated-root vectors).
820        iter = iter + 1
890        elim1 = vect(1,i)
           do 830 j=1,nm1
                if (b(j,2).eq.b(j,1)) go to 840
c     case one.
                vect(j,i) = elim1
                elim1 = vect(j+1,i) - elim1*b(j,5)
                go to 830
c     case two.
840             vect(j,i) = vect(j+1,i)
                elim1 = elim1 - vect(j+1,i)*temp
830        continue
           vect(n,i) = elim1
           go to 790
c     produce a random vector
801        continue
           do 802 j=1,n
c     generate pseudorandom numbers with uniform distribution in (-1,1).
c     this random number scheme is of the form...
c     rand1 = amod((2**12+3)*rand1,2**23)
c     it has a period of 2**21 numbers.
                rand1=mod(four99*rand1,rpower)
 802       vect(j,i)=rand1/rpow1-one
           go to 790
c
c     orthogonalize this repeated-root vector to others with this root.
800        if (ia.eq.0) go to 885
           do 860 j1=1,ia
                k = i - j1
                temp=zero
                do 870 j=1,n
870             temp = temp + vect(j,i)*vect(j,k)
                do 880 j=1,n
880             vect(j,i) = vect(j,i) - temp*vect(j,k)
860        continue
885        go to (890,900), iter
c     normalize the vector
 900       elim1=zero
           do 904 j=1,n
 904       elim1=max(abs(vect(j,i)),elim1)
           temp=zero
           do 910 j=1,n
                elim2=vect(j,i)/elim1
910        temp = temp + elim2**2
           temp=one/(sqrt(temp)*elim1)
           do 920 j=1,n
                vect(j,i) = vect(j,i)*temp
                if(abs(vect(j,i)) .lt. del1)vect(j,i)=zero
920        continue
700   continue
c
c     simvec section.
c     rotate codiagonal vectors into vectors of original array
c     loop over all the transformation vectors
      if (nm2.eq.0) go to 1002
      jump = nsize - (n+1)
      im = nm1
      do 950  i=1,nm2
           j1 = jump
c     move a transformation vector out into better indexing position.
           do 955  j=im,n
                b(j,2) = a(j1)
955        j1 = j1 + j
c     modify all requested vectors.
           do 960  k=1,nroot
                temp=zero
c     form scalar product of transformation vector with eigenvector
                do 970  j=im,n
970             temp = temp + b(j,2)*vect(j,k)
                temp = temp + temp
                do 980  j=im,n
980             vect(j,k) = vect(j,k) - temp*b(j,2)
960        continue
           jump = jump - im
950   im = im - 1
1002  continue
c     restore roots to their proper size.
      do 95 i=1,nroot
95    root(i) = root(i)*anorm
1001  return
      end
