*deck lzit.f
      subroutine lzit(nd,nn,a,b,x,iter,wantx,eiga,eigb)                  183.
c                                                                        184.
c this subroutine solves the generalized eigenvalue problem              185.
c              a x  = lambda b x                                         186.
c where a is a complex upper hessenberg matrix of order nn and b is      187.
c a complex upper triangular matrix of order nn                          188.
c                                                                        189.
c                                                                        190.
c input parameters                                                       191.
c                                                                        192.
c                                                                        193.
c nd     row dimension of the matrices a,b,x                             194.
c                                                                        195.
c nn     order of the problem                                            196.
c                                                                        197.
c a      an nn x nn upper hessenberg complex matrix                      198.
c                                                                        199.
c b      an nn x nn upper triangular complex matrix                      200.
c                                                                        201.
c                                                                        202.
c x      contains transformations to obtain eigenvectors of              203.
c        original system                                                 204.
c                                                                        205.
c wantx  logical variable which should be set to .true. if eigenvectors  206.
c        are wanted. otherwise it should be set to false                 207.
c                                                                        208.
c                                                                        209.
c     continue                                                           210.
c                                                                        211.
c output parameters                                                      212.
c                                                                        213.
c                                                                        214.
c x      the ith column contains the ith eigenvector if eigenvectors are 215.
c        requested                                                       216.
c                                                                        217.
c iter   an integer array of length nn whose ith entry contains the      218.
c        number of iterations needed to find the ith eigenvalue          219.
c        for any i if iter(i) =-1 then after 30 iterations there         220.
c        has not been a sufficient decrease in the last                  221.
c        subdiagonal elementof a to continue iterating.                  222.
c                                                                        223.
c eiga   an nn array containing the diagonal of a                        224.
c                                                                        225.
c eigb   an nn array containing the diagonal of b                        226.
c                                                                        227.
c the ith eigenvalue can be found by dividing eiga(i) by eigb(i)         228.
c watch out for eigb(i) being zero                                       229.
c                                                                        230.
c                                                                        231.
c problems with this subroutine should be directed to                    232.
c         linda kaufman                                                  233.
c         serra house, computer science department                       234.
c         stanford university                                            235.
c                                                                        236.
c                                                                        237.
      complex*16 a(nd,nd),b(nd,nd),eiga(nd),eigb(nd),x(nd,nd)            238.
      complex*16 s,w,y,z                                                 239.
      integer  iter(nd)                                                  240.
      complex*16 annm1,alfm,betm,d,sl,den,num,anm1m1                     241.
      real*8 epsa,epsb,ss,r,anorm,bnorm,ani,bni,c                          242.
      real*8 d0,d1,d2,e0,e1,rabs                                         244.
      logical wantx                                                      245.
      n=nn                                                               246.
c                                                                        247.
c   compute the machine precision times the norm of a and b              248.
c                                                                        249.
      anorm = 0.d0                                                         250.
      bnorm = 0.d0                                                         251.
      do   5 i=1,n                                                       252.
         ani = 0.d0                                                        253.
         if (i.ne.1) ani = rabs(a(i,i-1))                                254.
         bni = 0.d0                                                        255.
         do 3   j=i,n                                                    256.
            ani = ani + rabs(a(i,j))                                     257.
            bni = bni + rabs(b(i,j))                                     258.
    3    continue                                                        259.
         if (ani.gt.anorm) anorm = ani                                   260.
         if (bni.gt.bnorm) bnorm = bni                                   261.
    5 continue                                                           262.
      if (anorm.eq.0.d0) anorm = 1.d0                                      263.
      if (bnorm.eq.0.d0) bnorm = 1.d0                                      264.
      epsb=bnorm                                                         265.
      epsa=anorm                                                         266.
    6 epsa=epsa/2.d0                                                     267.
      epsb=epsb/2.d0                                                     268.
      c=anorm+epsa                                                       269.
      if (c.gt.anorm) go to 6                                            270.
c                                                                        271.
      if (n.le.1) go to 100                                              272.
   10 its=0                                                              273.
      nm1=n-1                                                            274.
c                                                                        275.
c check for negligible subdiagonal elements                              276.
c                                                                        277.
   11 do 12 lb=2,n                                                       278.
         l=n+2-lb                                                        279.
         ss=rabs(a(l-1,l-1))+rabs(a(l,l))                                280.
         r=ss+rabs(a(l,l-1))                                             281.
         if(r.eq.ss) go to 13                                            282.
   12 continue                                                           283.
      l=1                                                                284.
   13 if(l.eq.n) go to 100                                               285.
      if (its.lt.30) go to 20                                            286.
      iter(n)=-1                                                         287.
      if (rabs(a(n,nm1)).gt.0.8d0*rabs(annm1)) return                    288.
   20 if(its.eq.10.or.its.eq.20) go to 28                                289.
c                                                                        290.
c compute shift as eigenvalue of lower 2 by 2                            291.
c                                                                        292.
      annm1=a(n,nm1)                                                     293.
      anm1m1=a(nm1,nm1)                                                  294.
      s=a(n,n)*b(nm1,nm1)-(annm1)*b(nm1,n)                               295.
      w=annm1*b(n,n)*(a(nm1,n)*b(nm1,nm1)-                               296.
     i b(nm1,n)*anm1m1)                                                  297.
      y=(anm1m1*b(n,n)-s)/2.d0                                             298.
      z=sqrt(y*y+w)                                                     299.
      if (rabs(z).eq.0.d0) go to 26                                      300.
      d0=y/z                                                             301.
      if(d0.lt.0.d0) z=-z                                                302.
   26 den=(y+z)*b(nm1,nm1)*b(n,n)                                        303.
      if (rabs(den).eq.0.d0) den=epsa                                    304.
      num=(y+z)*s-w                                                      305.
      go to 30                                                           306.
c                                                                        307.
c ad-hoc shift                                                           308.
c                                                                        309.
   28 num=dcmplx(rabs(a(n,n-1)),rabs(a(n-1,n-2)))                        310.
      den=(1.d0,0.d0)                                                    311.
c                                                                        312.
c check for 2 consecutive small subdiagonal elements                     313.
c                                                                        314.
   30 if(n.eq.l+1) go to 35                                              315.
      d2=rabs(a(n-1,n-1))                                                316.
      e1=rabs(a(n,n-1))                                                  317.
      d1=rabs(a(n,n))                                                    318.
      nl=n-(l+1)                                                         319.
      do 34 mb=1,nl                                                      320.
         m=n-mb                                                          321.
         e0=e1                                                           322.
         e1=rabs(a(m,m-1))                                               323.
         d0=d1                                                           324.
         d1=d2                                                           325.
         d2=rabs(a(m-1,m-1))                                             326.
         d0=(d0+d1+d2)*rabs(a(m,m)*den-b(m,m)*num)                       327.
         e0=e0*e1*rabs(den)+d0                                           328.
         if(e0.eq.d0) go to 36                                           329.
   34 continue                                                           330.
   35 m=l                                                                331.
   36 continue                                                           332.
      its=its+1                                                          333.
      w=a(m,m)*den-b(m,m)*num                                            334.
      z=a(m+1,m)*den                                                     335.
c                                                                        336.
c find l and m and set a=lam and b=lbm                                   337.
c                                                                        338.
       np1=n+1                                                           339.
       lor1=l                                                            340.
       nnorn=n                                                           341.
       if (.not.wantx) go to 42                                          342.
           lor1=1                                                        343.
           nnorn=nn                                                      344.
   42  do 90 i=m,nm1                                                     345.
           j=i+1                                                         346.
c                                                                        347.
c find row transformations to restore a to                               348.
c upper hessenberg form. apply transformations                           349.
c to a and b                                                             350.
c                                                                        351.
           if (i.eq.m) go to 50                                          352.
           w=a(i,i-1)                                                    353.
           z=a(j,i-1)                                                    354.
           if(rabs(z).eq.0.d0) go to 11                                  355.
   50      if (rabs(w).ge.rabs(z)) go to 60                              356.
c                                                                        357.
c must pivot                                                             358.
c                                                                        359.
              do 55 k=i,nnorn                                            360.
                 y=a(i,k)                                                361.
                 a(i,k)=a(j,k)                                           362.
                 a(j,k)=y                                                363.
                 y=b(i,k)                                                364.
                 b(i,k)=b(j,k)                                           365.
                 b(j,k)=y                                                366.
   55         continue                                                   367.
              if (i.gt.m) a(i,i-1)=a(j,i-1)                              368.
              if(rabs(w).eq.0.d0) go to 65                               369.
              y=w/z                                                      370.
              go to 62                                                   371.
   60      y=z/w                                                         372.
   62      do 64 k=i,nnorn                                               373.
              a(j,k)=a(j,k)-y*a(i,k)                                     374.
              b(j,k)=b(j,k)-y*b(i,k)                                     375.
   64      continue                                                      376.
           if (i.gt.m) a(j,i-1)=(0.d0,0.d0)                              377.
c                                                                        378.
c perform transformations from right to restore b to                     379.
c   trianglular form                                                     380.
c apply transformations to a                                             381.
c                                                                        382.
   65      if (rabs(b(j,i)).eq.0.d0) go to 11                            383.
           if (rabs(b(j,j)).ge.rabs(b(j,i))) go to 81                    384.
c                                                                        385.
c must pivot columns                                                     386.
c                                                                        387.
              do 70 k=lor1,j                                             388.
                 y=a(k,j)                                                389.
                 a(k,j)=a(k,i)                                           390.
                 a(k,i)=y                                                391.
                 y=b(k,j)                                                392.
                 b(k,j)=b(k,i)                                           393.
                 b(k,i)=y                                                394.
   70         continue                                                   395.
              if (i.eq.nm1) go to 75                                     396.
                 y=a(j+1,j)                                              397.
                 a(j+1,j)=a(j+1,i)                                       398.
                 a(j+1,i)=y                                              399.
   75         if(.not.wantx) go to 80                                    400.
              do 78 k=1,nn                                               401.
                 y=x(k,j)                                                402.
                 x(k,j)=x(k,i)                                           403.
                 x(k,i)=y                                                404.
   78         continue                                                   405.
   80         if (rabs(b(j,i)).eq.0.d0) go to 90                         406.
   81      z=b(j,i)/b(j,j)                                               407.
           do 82 k=lor1,j                                                408.
              a(k,i)=a(k,i)-z*a(k,j)                                     409.
              b(k,i)=b(k,i)-z*b(k,j)                                     410.
   82      continue                                                      411.
           b(j,i)=(0.d0,0.d0)                                            412.
           if (i.lt.nm1) a(i+2,i)=a(i+2,i)-z*a(i+2,j)                    413.
           if(.not.wantx) go to 90                                       414.
           do 85 k=1,nn                                                  415.
              x(k,i)=x(k,i)-z*x(k,j)                                     416.
   85      continue                                                      417.
   90 continue                                                           418.
      go to 11                                                           419.
c                                                                        420.
  100 continue                                                           421.
       eiga(n)=a(n,n)                                                    422.
       eigb(n)=b(n,n)                                                    423.
       if (n.eq.1) go to 110                                             424.
       iter(n)=its                                                       425.
       n=nm1                                                             426.
       if (n.gt.1) go to  10                                             427.
       iter(1)=0                                                         428.
       go to 100                                                         429.
c                                                                        430.
c      find eigenvectors using b for intermediate storage                431.
c                                                                        432.
  110   if(.not.wantx) return                                            433.
       m=nn                                                              434.
  115  continue                                                          435.
       alfm=a(m,m)                                                       436.
       betm=b(m,m)                                                       437.
       b(m,m)=(1.d0,0.d0)                                                438.
         l = m-1                                                         439.
         if (l.eq.0) go to 140                                           440.
  120    continue                                                        441.
            l1 = l+1                                                     442.
            sl = 0.d0                                                      443.
            do 130 j=l1,m                                                444.
               sl = sl + (betm*a(l,j)-alfm*b(l,j))*b(j,m)                445.
  130       continue                                                     446.
            d = betm*a(l,l)-alfm*b(l,l)                                  447.
            if (rabs(d).eq.0.d0) d = (epsa+epsb)/2.d0                      448.
            b(l,m) = -sl/d                                               449.
            l = l-1                                                      450.
  140    if (l.gt.0) go to 120                                           451.
         m=m-1                                                           452.
         if (m.gt.0) go to 115                                           453.
c                                                                        454.
c  transform to original coordinate system                               455.
c                                                                        456.
      m = nn                                                             457.
  200 continue                                                           458.
         do 220 i=1,nn                                                   459.
            s = 0.d0                                                       460.
            do 210 j=1,m                                                 461.
               s = s + x(i,j)*b(j,m)                                     462.
  210       continue                                                     463.
            x(i,m) = s                                                   464.
  220    continue                                                        465.
         m = m-1                                                         466.
      if (m.gt.0) go to 200                                              467.
c                                                                        468.
c  normalize so that largest component = 1.                              469.
c                                                                        470.
      m = nn                                                             471.
  230 continue                                                           472.
         ss = 0.d0                                                         473.
         do 235 i=1,nn                                                   474.
            r =  rabs(x(i,m))                                            475.
            if (r.lt.ss) go to 235                                       476.
            ss = r                                                       477.
            d = x(i,m)                                                   478.
  235    continue                                                        479.
         if (ss.eq.0.d0) go to 245                                       480.
         do 240 i=1,nn                                                   481.
            x(i,m) = x(i,m)/d                                            482.
  240    continue                                                        483.
  245 m = m-1                                                            484.
      if (m.gt.0) go to 230                                              485.
      return                                                             486.
      end                                                                487.
