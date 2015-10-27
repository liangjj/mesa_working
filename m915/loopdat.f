*deck @(#)loopdat.f	1.1  11/30/90
      block data loopdat
c
c      common /tables/ jsegnr(22),jsegpt(22),iarcmn(228),iarcsb(228)
c     *,itrk(228),jcond(228),kcond(228),nxtseg(228),jsegpx(3)
      common /tables/ jsegnr(22),jsegpt(22),iarcmn(150),iarcmq(78),
     #                iarcsb(150),iarcsq(78),itrk(150),itrkq(78),
     #                jcond(150),jcondq(78),kcond(150),kcondq(78),
     #                nxtseg(150),nxtseq(78),jsegpx(3)
c
c
      data jsegnr/16,34,52,63,75,92,102,118,128,137,148,155,162,172,
     a179,186,193,200,207,214,221,228/
      data itrk  / 1,  3,  1,  3,  1,  2,  9,  1,  1,  7
     a,           2,  9,  1,  7, 10,  9,  0,  4,  4,  3
     b,           0,  2,  4,  9,  0, 10,  0, 10,  4,  0
     c,           9,  3,  0, 10,  0,  4,  4,  9,  0, 10
     d,           0, 10,  4,  3,  0,  2,  4,  0,  3,  9
     e,           0, 10,  0,  0,  0,  0,  0,  0,  0,  0
     f,           0,  0,  0,  0,  0,  0,  0,  0,  1,  0
     g,           0,  0,  0,  0,  0,  0,  5,  5,  0,  0
     h,           0,  3,  5,  0,  3,  0,  0,  5,  0,  0
     i,           0,  0,  0,  0,  0,  1,  1,  0,  0,  0
     j,           1,  0,  0,  0,  0,  1,  1,  0,  0,  0
     k,           1,  0,  1,  0,  0,  1,  1,  0,  0,  0
     l,           1,  0,  1,  0,  0,  0,  1,  0,  0,  0
     m,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     n,           0,  0,  0,  0,  1,  0,  0,  0,  0,  0/
       data itrkq/0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     p,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     q,           0,  0,  0,  0,  0,  0,  0,  8,  0,  0
     r,           0,  0,  0,  0,  8,  0,  0,  0,  0,  0
     s,           0,  6,  0,  0,  0,  0,  0,  0,  6,  0
     t,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     u,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     v,           0,  0,  0,  0,  0,  0,  0,  0/
      data jcond /-1,  1, -1,  1,  1,  1,  1,  1, -1,  1
     a,           1,  1, -1,  1,  1,  1, -1,  1,  1,  1
     b,          -1,  1,  1,  1, -1,  1, -1,  1,  1,  1
     c,           1,  1, -1,  1, -1,  1,  1,  1, -1,  1
     d,          -1,  1,  1,  1, -1,  1,  1,  1,  1,  1
     e,          -1,  1,  0,  0,  0,  0,  0,  0,  0,  0
     f,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     g,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     h,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     i,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     j,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     k,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     l,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     m,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     n,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0/
      data jcondq/0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     p,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     q,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     r,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     s,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     t,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     u,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     v,           0,  0,  0,  0,  0,  0,  0,  0/
      data kcond / 0,  1,  0,  1,  0,  1,  0,  0,  0,  1
     a,           1,  0,  0,  1,  1,  0,  0,  0,  0,  0
     b,           0,  1,  0,  0,  0,  1,  0,  1,  0,  1
     c,           0,  0,  0,  1,  0,  0,  0,  0,  0,  1
     d,           0,  1,  0,  0,  0,  1,  0,  1,  0,  0
     e,           0,  1,  0,  1,  1,  1,  0,  1,  1,  0
     f,           1,  1,  0,  0,  1,  1,  1,  0,  0,  1
     g,           1,  0,  1,  1,  0,  0,  1,  1,  1,  1
     h,           0,  0,  1,  1,  0,  1,  0,  1,  1,  1
     i,           1,  0,  0,  1,  0,  1,  1,  0,  0,  1
     j,           1,  0,  0,  1,  1,  1,  1,  0,  0,  1
     k,           1,  0,  1,  0,  1,  1,  1,  0,  0,  1
     l,           1,  0,  1,  0,  1,  0,  1,  0,  0,  1
     m,           0,  1,  0,  1,  1,  1,  0,  0,  1,  0
     n,           0,  1,  0,  0,  1,  1,  1,  0,  0,  0/
      data kcondq/1,  0,  0,  1,  0,  0,  1,  0,  0,  0
     p,           1,  0,  0,  1,  0,  0,  1,  0,  0,  1
     q,           1,  0,  0,  0,  0,  0,  0,  0,  0,  0
     r,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     s,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     t,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     u,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     v,           0,  0,  0,  0,  0,  0,  0,  0/
      data nxtseg/ 3, 18,  2, 17, 10,  0,  5,  7,  2, 15
     a,           0,  5,  3, 16,  0,  4,  2, 11, 12,  7
     b,           2, 21, 12,  6,  3, 22,  2, 21, 11, 20
     c,           6,  7,  2, 21,  3, 13, 11,  6,  3, 22
     d,           2, 21, 11,  9,  3, 22, 13, 19,  9,  6
     e,           3, 22,  4, 22, 21,  0,  4, 21,  0,  4
     f,          22,  0,  4,  5, 22, 21,  0,  5,  7, 21
     g,           0,  5, 22,  0,  5,  6, 22, 21, 20,  0
     h,           6,  7, 21, 19,  9,  0,  6, 22, 19, 20
     i,           0,  6,  7, 21,  7, 20,  0,  8,  7, 21
     j,          20,  7,  8, 22, 21, 20,  0,  8,  7, 21
     k,          19,  9,  0,  8, 22, 19, 20,  8,  9, 22
     l,          19,  9,  0,  8, 22,  9, 19,  9, 10, 21
     m,          10, 22, 10,  0, 22, 21, 10, 11, 21, 11
     n,          12, 22, 13, 11,  0, 22, 21, 11, 12, 12/
      data nxtseq/21, 14, 12, 21, 12, 13, 22, 13, 14, 13
     p,          22, 13, 14, 21, 14, 12, 22, 13, 14, 22
     q,          21, 14, 15, 15,  0, 16, 15,  0, 15, 16
     r,           0, 16, 15, 16,  0, 16, 17, 17,  0, 18
     s,          17,  0, 17, 18,  0, 18, 17, 18,  0, 18
     t,          19,  0, 19, 20,  0, 19, 19, 20,  0, 20
     u,          19, 20,  0, 20, 21, 21,  0, 22, 21,  0
     v,          21, 22,  0, 22, 21, 22,  0, 22/
      data iarcmn/ 2,  2,  3,  3,  4,  2,  2,  3,  4,  4
     a,           3,  3,  4,  4,  4,  4,  1,  2,  3,  1
     b,           2,  2,  4,  1,  2,  2,  3,  3,  4,  1
     c,           2,  3,  4,  4,  1,  2,  3,  1,  2,  2
     d,           3,  3,  4,  1,  3,  3,  4,  1,  2,  3
     e,           4,  4,  1,  2,  3,  2,  2,  4,  3,  3
     f,           4,  4,  4,  1,  2,  3,  2,  2,  3,  4
     g,           3,  3,  4,  4,  4,  1,  2,  3,  1,  2
     h,           2,  3,  4,  1,  2,  3,  3,  4,  2,  3
     i,           4,  4,  1,  2,  2,  1,  2,  2,  3,  4
     j,           2,  4,  1,  2,  3,  1,  2,  2,  3,  4
     k,           1,  2,  3,  3,  4,  2,  3,  4,  1,  3
     l,           1,  2,  3,  3,  4,  3,  3,  4,  1,  1
     m,           2,  1,  3,  1,  2,  3,  4,  1,  1,  2
     n,           3,  1,  2,  3,  1,  2,  3,  4,  1,  2/
      data iarcmq/1,  2,  3,  2,  4,  1,  1,  2,  3,  3
     p,           3,  4,  1,  1,  2,  3,  1,  2,  3,  2
     q,           3,  4,  1,  2,  1,  2,  3,  2,  4,  1
     r,           1,  2,  3,  3,  3,  4,  1,  2,  1,  2
     s,           3,  2,  4,  1,  1,  2,  3,  3,  3,  4
     t,           1,  3,  2,  3,  4,  3,  4,  1,  2,  2
     u,           2,  3,  4,  4,  1,  2,  1,  2,  3,  2
     v,           4,  1,  1,  2,  3,  3,  3,  4/
      data iarcsb/ 1,  1,  1,  1,  1,  2,  2,  2,  2,  2
     a,           3,  3,  3,  3,  4,  4,  1,  1,  1,  2
     b,           2,  2,  2,  3,  3,  3,  3,  3,  3,  4
     c,           4,  4,  4,  4,  1,  1,  1,  2,  2,  2
     d,           2,  2,  2,  3,  3,  3,  3,  4,  4,  4
     e,           4,  4,  1,  1,  1,  2,  2,  2,  3,  3
     f,           3,  4,  4,  1,  1,  1,  2,  2,  2,  2
     g,           3,  3,  3,  4,  4,  1,  1,  1,  2,  2
     h,           2,  2,  2,  3,  3,  3,  3,  3,  4,  4
     i,           4,  4,  1,  1,  2,  3,  3,  3,  3,  3
     j,           4,  4,  1,  1,  1,  2,  2,  2,  2,  2
     k,           3,  3,  3,  3,  3,  4,  4,  4,  1,  1
     l,           2,  2,  2,  2,  2,  3,  4,  4,  1,  2
     m,           2,  3,  3,  4,  4,  4,  4,  1,  2,  2
     n,           2,  3,  3,  3,  4,  4,  4,  4,  1,  2/
      data iarcsq/3,  3,  3,  4,  4,  1,  2,  2,  2,  3
     p,           4,  4,  1,  2,  2,  2,  3,  3,  3,  4
     q,           4,  4,  1,  2,  3,  3,  3,  4,  4,  1
     r,           2,  2,  2,  3,  4,  4,  1,  2,  3,  3
     s,           3,  4,  4,  1,  2,  2,  2,  3,  4,  4
     t,           1,  1,  2,  2,  2,  3,  4,  1,  1,  2
     u,           3,  3,  3,  4,  1,  2,  3,  3,  3,  4
     v,           4,  1,  2,  2,  2,  3,  4,  4/
      data jsegpx/12,29,47/
      end
