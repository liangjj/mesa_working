*deck @(#)lupdat.f	5.1  11/6/94
      block data lupdat
cbhl      common /tapes/itape2,itape5,itape6,itape8,itap12,itap03,itap04
cbhl     *,             itape3,itap05,itap06
      common /tables/ jsegnr(22),jsegpt(22),iarcmn(228),iarcsb(228)
     *,itrk(228),jcond(228),kcond(228),nxtseg(228),jsegpx(3)
cbhl      data itape2,itape3,itape5,itape6,itape8,itap12/2,3,5,6,8,12/
cbhl      data itap03,itap04,itap05,itap06/93,94,95,94/
      data jsegnr/16,34,52,63,75,92,102,118,128,137,148,155,162,172,
     a 179,186,193,200,207,214,221,228/
      data itrk  / 1,  3,  1,  3,  1,  2,  9,  1,  1,  7
     a ,           2,  9,  1,  7, 10,  9,  0,  4,  4,  3
     b ,           0,  2,  4,  9,  0, 10,  0, 10,  4,  0
     c ,           9,  3,  0, 10,  0,  4,  4,  9,  0, 10
     d ,           0, 10,  4,  3,  0,  2,  4,  0,  3,  9
     e ,           0, 10,  0,  0,  0,  0,  0,  0,  0,  0
     f ,           0,  0,  0,  0,  0,  0,  0,  0,  1,  0
     g ,           0,  0,  0,  0,  0,  0,  5,  5,  0,  0
     h ,           0,  3,  5,  0,  3,  0,  0,  5,  0,  0
     i ,           0,  0,  0,  0,  0,  1,  1,  0,  0,  0
     j ,           1,  0,  0,  0,  0,  1,  1,  0,  0,  0
     k ,           1,  0,  1,  0,  0,  1,  1,  0,  0,  0
     l ,           1,  0,  1,  0,  0,  0,  1,  0,  0,  0
     m ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     n ,           0,  0,  0,  0,  1,  0,  0,  0,  0,  0
     o ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     p ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     q ,           0,  0,  0,  0,  0,  0,  0,  8,  0,  0
     r ,           0,  0,  0,  0,  8,  0,  0,  0,  0,  0
     s , 0,6,6*0,6,29*0/
c    s ,           0,  6,  0,  0,  0,  0,  0,  0,  6,  0
c    t ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
c    u ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
c    v ,           0,  0,  0,  0,  0,  0,  0,  0/
      data jcond /-1,  1, -1,  1,  1,  1,  1,  1, -1,  1
     a ,           1,  1, -1,  1,  1,  1, -1,  1,  1,  1
     b ,          -1,  1,  1,  1, -1,  1, -1,  1,  1,  1
     c ,           1,  1, -1,  1, -1,  1,  1,  1, -1,  1
     d ,          -1,  1,  1,  1, -1,  1,  1,  1,  1,  1
     e ,          -1,  1,  0,  0,  0,  0,  0,  0,  0,  0
     f ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     g ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     h ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     i ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     j ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     k ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     l ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     m ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     n ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     o ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     p ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     q ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     r ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     s , 38*0/
c    s ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
c    t ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
c    u ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
c    v ,           0,  0,  0,  0,  0,  0,  0,  0/
      data kcond / 0,  1,  0,  1,  0,  1,  0,  0,  0,  1
     a ,           1,  0,  0,  1,  1,  0,  0,  0,  0,  0
     b ,           0,  1,  0,  0,  0,  1,  0,  1,  0,  1
     c ,           0,  0,  0,  1,  0,  0,  0,  0,  0,  1
     d ,           0,  1,  0,  0,  0,  1,  0,  1,  0,  0
     e ,           0,  1,  0,  1,  1,  1,  0,  1,  1,  0
     f ,           1,  1,  0,  0,  1,  1,  1,  0,  0,  1
     g ,           1,  0,  1,  1,  0,  0,  1,  1,  1,  1
     h ,           0,  0,  1,  1,  0,  1,  0,  1,  1,  1
     i ,           1,  0,  0,  1,  0,  1,  1,  0,  0,  1
     j ,           1,  0,  0,  1,  1,  1,  1,  0,  0,  1
     k ,           1,  0,  1,  0,  1,  1,  1,  0,  0,  1
     l ,           1,  0,  1,  0,  1,  0,  1,  0,  0,  1
     m ,           0,  1,  0,  1,  1,  1,  0,  0,  1,  0
     n ,           0,  1,  0,  0,  1,  1,  1,  0,  0,  0
     o ,           1,  0,  0,  1,  0,  0,  1,  0,  0,  0
     p ,           1,  0,  0,  1,  0,  0,  1,  0,  0,  1
     q ,           1,  0,  0,  0,  0,  0,  0,  0,  0,  0
     r ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     s , 38*0/
c    s ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
c    t ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
c    u ,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
c    v ,           0,  0,  0,  0,  0,  0,  0,  0/
      data nxtseg/ 3, 18,  2, 17, 10,  0,  5,  7,  2, 15
     a ,           0,  5,  3, 16,  0,  4,  2, 11, 12,  7
     b ,           2, 21, 12,  6,  3, 22,  2, 21, 11, 20
     c ,           6,  7,  2, 21,  3, 13, 11,  6,  3, 22
     d ,           2, 21, 11,  9,  3, 22, 13, 19,  9,  6
     e ,           3, 22,  4, 22, 21,  0,  4, 21,  0,  4
     f ,          22,  0,  4,  5, 22, 21,  0,  5,  7, 21
     g ,           0,  5, 22,  0,  5,  6, 22, 21, 20,  0
     h ,           6,  7, 21, 19,  9,  0,  6, 22, 19, 20
     i ,           0,  6,  7, 21,  7, 20,  0,  8,  7, 21
     j ,          20,  7,  8, 22, 21, 20,  0,  8,  7, 21
     k ,          19,  9,  0,  8, 22, 19, 20,  8,  9, 22
     l ,          19,  9,  0,  8, 22,  9, 19,  9, 10, 21
     m ,          10, 22, 10,  0, 22, 21, 10, 11, 21, 11
     n ,          12, 22, 13, 11,  0, 22, 21, 11, 12, 12
     o ,          21, 14, 12, 21, 12, 13, 22, 13, 14, 13
     p ,          22, 13, 14, 21, 14, 12, 22, 13, 14, 22
     q ,21,14,15,15, 0,16,15, 0,15,16, 0,16,15,16, 0,16,17,17, 0,18
     r ,17, 0,17,18, 0,18,17,18, 0,18,19, 0,19,20, 0,19,19,20, 0,20
     s ,19,20, 0,20,21,21, 0,22,21, 0,21,22, 0,22,21,22, 0,22/
c    q ,          21, 14, 15, 15,  0, 16, 15,  0, 15, 16
c    r ,           0, 16, 15, 16,  0, 16, 17, 17,  0, 18
c    s ,          17,  0, 17, 18,  0, 18, 17, 18,  0, 18
c    t ,          19,  0, 19, 20,  0, 19, 19, 20,  0, 20
c    u ,          19, 20,  0, 20, 21, 21,  0, 22, 21,  0
c    v ,          21, 22,  0, 22, 21, 22,  0, 22/
      data iarcmn/ 2,  2,  3,  3,  4,  2,  2,  3,  4,  4
     a ,           3,  3,  4,  4,  4,  4,  1,  2,  3,  1
     b ,           2,  2,  4,  1,  2,  2,  3,  3,  4,  1
     c ,           2,  3,  4,  4,  1,  2,  3,  1,  2,  2
     d ,           3,  3,  4,  1,  3,  3,  4,  1,  2,  3
     e ,           4,  4,  1,  2,  3,  2,  2,  4,  3,  3
     f ,           4,  4,  4,  1,  2,  3,  2,  2,  3,  4
     g ,           3,  3,  4,  4,  4,  1,  2,  3,  1,  2
     h ,           2,  3,  4,  1,  2,  3,  3,  4,  2,  3
     i ,           4,  4,  1,  2,  2,  1,  2,  2,  3,  4
     j ,           2,  4,  1,  2,  3,  1,  2,  2,  3,  4
     k ,           1,  2,  3,  3,  4,  2,  3,  4,  1,  3
     l ,           1,  2,  3,  3,  4,  3,  3,  4,  1,  1
     m ,           2,  1,  3,  1,  2,  3,  4,  1,  1,  2
     n ,           3,  1,  2,  3,  1,  2,  3,  4,  1,  2
     o ,           1,  2,  3,  2,  4,  1,  1,  2,  3,  3
     p ,           3,  4,  1,  1,  2,  3,  1,  2,  3,  2
     q ,3,4,1,2,1,2,3,2,4,1, 1,2,3,3,3,4,1,2,1,2
     r ,3,2,4,1,1,2,3,3,3,4, 1,3,2,3,4,3,4,1,2,2
     s ,2,3,4,4,1,2,1,2,3,2, 4,1,1,2,3,3,3,4/
c    q ,           3,  4,  1,  2,  1,  2,  3,  2,  4,  1
c    r ,           1,  2,  3,  3,  3,  4,  1,  2,  1,  2
c    s ,           3,  2,  4,  1,  1,  2,  3,  3,  3,  4
c    t ,           1,  3,  2,  3,  4,  3,  4,  1,  2,  2
c    u ,           2,  3,  4,  4,  1,  2,  1,  2,  3,  2
c    v ,           4,  1,  1,  2,  3,  3,  3,  4/
      data iarcsb/ 1,  1,  1,  1,  1,  2,  2,  2,  2,  2
     a ,           3,  3,  3,  3,  4,  4,  1,  1,  1,  2
     b ,           2,  2,  2,  3,  3,  3,  3,  3,  3,  4
     c ,           4,  4,  4,  4,  1,  1,  1,  2,  2,  2
     d ,           2,  2,  2,  3,  3,  3,  3,  4,  4,  4
     e ,           4,  4,  1,  1,  1,  2,  2,  2,  3,  3
     f ,           3,  4,  4,  1,  1,  1,  2,  2,  2,  2
     g ,           3,  3,  3,  4,  4,  1,  1,  1,  2,  2
     h ,           2,  2,  2,  3,  3,  3,  3,  3,  4,  4
     i ,           4,  4,  1,  1,  2,  3,  3,  3,  3,  3
     j ,           4,  4,  1,  1,  1,  2,  2,  2,  2,  2
     k ,           3,  3,  3,  3,  3,  4,  4,  4,  1,  1
     l ,           2,  2,  2,  2,  2,  3,  4,  4,  1,  2
     m ,           2,  3,  3,  4,  4,  4,  4,  1,  2,  2
     n ,           2,  3,  3,  3,  4,  4,  4,  4,  1,  2
     o ,           3,  3,  3,  4,  4,  1,  2,  2,  2,  3
     p ,           4,  4,  1,  2,  2,  2,  3,  3,  3,  4
     q ,4,4,1,2,3,3,3,4,4,1, 2,2,2,3,4,4,1,2,3,3
     r ,3,4,4,1,2,2,2,3,4,4, 1,1,2,2,2,3,4,1,1,2
     s ,3,3,3,4,1,2,3,3,3,4, 4,1,2,2,2,3,4,4/
c    q ,           4,  4,  1,  2,  3,  3,  3,  4,  4,  1
c    r ,           2,  2,  2,  3,  4,  4,  1,  2,  3,  3
c    s ,           3,  4,  4,  1,  2,  2,  2,  3,  4,  4
c    t ,           1,  1,  2,  2,  2,  3,  4,  1,  1,  2
c    u ,           3,  3,  3,  4,  1,  2,  3,  3,  3,  4
c    v ,           4,  1,  2,  2,  2,  3,  4,  4/
      data jsegpx/12,29,47/
      end
