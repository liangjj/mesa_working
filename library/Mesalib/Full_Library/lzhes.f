*deck lzhes.f
      subroutine lzhes(nd,n,a,b,x,wantx)                                  15.
c                                                                         16.
c    this subroutine reduces the complex matrix a to upper                17.
c    hessenberg form and reduces the complex matrix b to                  18.
c    triangular form                                                      19.
c                                                                         20.
c    input parameters:                                                    21.
c                                                                         22.
c    nd  the row dimension of the matrices a,b,x                          23.
c                                                                         24.
c    n   the order of the problem                                         25.
c                                                                         26.
c    a   a complex matrix                                                 27.
c                                                                         28.
c    b   a complex matrix                                                 29.
c                                                                         30.
c    wantx a logical variable which is set to .true. if                   31.
c          the eigenvectors are wanted. otherwise it should               32.
c        be set to .false.                                                33.
c                                                                         34.
c      continue                                                           35.
c    output parameters:                                                   36.
c                                                                         37.
c    a  on output a is an upper hessenberg matrix, the                    38.
c       original matrix has been destroyed                                39.
c                                                                         40.
c    b  an upper triangular matrix, the original matrix                   41.
c       has been destroyed                                                42.
c                                                                         43.
c    x  contains the transformations needed to compute                    44.
c       the eigenvectors of the original system                           45.
c                                                                         46.
c                                                                         47.
c     problems with this subroutine should be directed to:                48.
c                                                                         49.
c         linda kaufman                                                   50.
c         serra house                                                     51.
c         computer science department                                     52.
c         stanford university                                             53.
c                                                                         54.
      complex*16 y,a(nd,nd),b(nd,nd),x(nd,nd)
      real*8 c,rabs,d                                                     56.
      logical wantx                                                       57.
      nm1=n-1                                                             58.
c                                                                         59.
c     reduce b to triangular form using elementary transformations        60.
c                                                                         61.
      do 30 i=1,nm1                                                       62.
      d=0.0d0                                                             63.
           ip1=i+1                                                        64.
           do 10 k=ip1,n                                                  65.
                c=rabs(b(k,i))                                            66.
                if (c.le.d) go to 10                                      67.
                     d=c                                                  68.
                     ii=k                                                 69.
   10      continue                                                       70.
           if (d.eq.0.d0) go to 30                                        71.
           if (d.le.rabs(b(i,i))) go to 15                                72.
c                                                                         73.
c must pivot                                                              74.
c                                                                         75.
              do 11 j=1,n                                                 76.
                   y=a(i,j)                                               77.
                   a(i,j)=a(ii,j)                                         78.
                   a(ii,j)=y                                              79.
   11         continue                                                    80.
              do 12 j=i,n                                                 81.
                   y=b(i,j)                                               82.
                   b(i,j)=b(ii,j)                                         83.
                   b(ii,j)=y                                              84.
   12         continue                                                    85.
   15      do 20 j=ip1,n                                                  86.
                y=b(j,i)/b(i,i)                                           87.
                if (rabs(y).eq.0.d0) go to 20                             88.
                do 18 k=1,n                                               89.
                     a(j,k)=a(j,k)-y*a(i,k)                               90.
   18           continue                                                  91.
                do 19 k=ip1,n                                             92.
                     b(j,k)=b(j,k)-y*b(i,k)                               93.
   19           continue                                                  94.
   20      continue                                                       95.
           b(ip1,i)=(0.d0,0.d0)                                           96.
   30 continue                                                            97.
c                                                                         98.
c initialize x                                                            99.
c                                                                        100.
      if (.not.wantx) go to 40                                           101.
      do 38 i=1,n                                                        102.
           do 37 j=1,n                                                   103.
                x(i,j)=(0.d0,0.d0)                                       104.
   37      continue                                                      105.
           x(i,i)=(1.d0,0.0d0)                                           106.
   38 continue                                                           107.
c                                                                        108.
c     reduce a to upper hessenberg form                                  109.
c                                                                        110.
   40 nm2=n-2                                                            111.
      if (nm2.lt.1) go to 100                                            112.
      do 90 j=1,nm2                                                      113.
           jm2=nm1-j                                                     114.
           jp1=j+1                                                       115.
               do 80 ii=1,jm2                                            116.
                i=n+1-ii                                                 117.
                im1=i-1                                                  118.
                imj=i-j                                                  119.
                if (rabs(a(i,j)).le.rabs(a(im1,j))) go to 50             120.
c                                                                        121.
c       must pivot rows                                                  122.
c                                                                        123.
                   do 45 k=j,n                                           124.
                        y=a(i,k)                                         125.
                        a(i,k)=a(im1,k)                                  126.
                        a(im1,k)=y                                       127.
   45              continue                                              128.
                   do 46 k=im1,n                                         129.
                        y=b(i,k)                                         130.
                        b(i,k)=b(im1,k)                                  131.
                        b(im1,k)=y                                       132.
   46              continue                                              133.
   50           if (rabs(a(i,j)).eq.0.d0) go to 58                       134.
                y=a(i,j)/a(im1,j)                                        135.
                do 52 k=jp1,n                                            136.
                     a(i,k)=a(i,k)-y*a(im1,k)                            137.
   52           continue                                                 138.
                do 54 k=im1,n                                            139.
                     b(i,k)=b(i,k)-y*b(im1,k)                            140.
   54           continue                                                 141.
c                                                                        142.
c     transformation from the right                                      143.
c                                                                        144.
   58           if (rabs(b(i,im1)).le.rabs(b(i,i))) go to 70             145.
c                                                                        146.
c       must pivot columns                                               147.
c                                                                        148.
                   do 60 k=1,i                                           149.
                        y=b(k,i)                                         150.
                        b(k,i)=b(k,im1)                                  151.
                        b(k,im1)=y                                       152.
   60              continue                                              153.
                   do 64 k=1,n                                           154.
                        y=a(k,i)                                         155.
                        a(k,i)=a(k,im1)                                  156.
                        a(k,im1)=y                                       157.
   64              continue                                              158.
                   if (.not.wantx) go to 70                              159.
                   do 68 k=imj,n                                         160.
                        y=x(k,i)                                         161.
                        x(k,i)=x(k,im1)                                  162.
                        x(k,im1)=y                                       163.
   68              continue                                              164.
   70           if (rabs(b(i,im1)).eq.0.d0) go to 80                     165.
                y=b(i,im1)/b(i,i)                                        166.
                do 72 k=1,im1                                            167.
                     b(k,im1)=b(k,im1)-y*b(k,i)                          168.
   72           continue                                                 169.
                b(i,im1)=(0.d0,0.d0)                                     170.
                do 74 k=1,n                                              171.
                     a(k,im1)=a(k,im1)-y*a(k,i)                          172.
   74           continue                                                 173.
                if (.not.wantx) go to 80                                 174.
                do 76 k=imj,n                                            175.
                     x(k,im1)=x(k,im1)-y*x(k,i)                          176.
   76           continue                                                 177.
   80      continue                                                      178.
           a(jp1+1,j)=(0.d0,0.d0)                                        179.
   90 continue                                                           180.
  100 return                                                             181.
      end                                                                182.

