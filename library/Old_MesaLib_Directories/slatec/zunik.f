*deck zunik
      subroutine zunik (zrr, zri, fnu, ikflg, ipmtr, tol, init, phir,
     +   phii, zeta1r, zeta1i, zeta2r, zeta2i, sumr, sumi, cwrkr, cwrki)
c***begin prologue  zunik
c***subsidiary
c***purpose  subsidiary to zbesi and zbesk
c***library   slatec
c***type      all (cunik-a, zunik-a)
c***author  amos, d. e., (snl)
c***description
c
c        zunik computes parameters for the uniform asymptotic
c        expansions of the i and k functions on ikflg= 1 or 2
c        respectively by
c
c        w(fnu,zr) = phi*exp(zeta)*sum
c
c        where       zeta=-zeta1 + zeta2       or
c                          zeta1 - zeta2
c
c        the first call must have init=0. subsequent calls with the
c        same zr and fnu will return the i or k function on ikflg=
c        1 or 2 with no change in init. cwrk is a complex work
c        array. ipmtr=0 computes all parameters. ipmtr=1 computes phi,
c        zeta1,zeta2.
c
c***see also  zbesi, zbesk
c***routines called  d1mach, zdiv, zlog, zsqrt
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c   930122  added external statement with zlog and zsqrt.  (rwc)
c***end prologue  zunik
c     complex cfn,con,cone,crfn,cwrk,czero,phi,s,sr,sum,t,t2,zeta1,
c    *zeta2,zn,zr
      double precision ac, c, con, conei, coner, crfni, crfnr, cwrki,
     * cwrkr, fnu, phii, phir, rfn, si, sr, sri, srr, sti, str, sumi,
     * sumr, test, ti, tol, tr, t2i, t2r, zeroi, zeror, zeta1i, zeta1r,
     * zeta2i, zeta2r, zni, znr, zri, zrr, d1mach
      integer i, idum, ikflg, init, ipmtr, j, k, l
      dimension c(120), cwrkr(16), cwrki(16), con(2)
      external zlog, zsqrt
      data zeror,zeroi,coner,conei / 0.0d0, 0.0d0, 1.0d0, 0.0d0 /
      data con(1), con(2)  /
     1 3.98942280401432678d-01,  1.25331413731550025d+00 /
      data c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9), c(10),
     1     c(11), c(12), c(13), c(14), c(15), c(16), c(17), c(18),
     2     c(19), c(20), c(21), c(22), c(23), c(24)/
     3     1.00000000000000000d+00,    -2.08333333333333333d-01,
     4     1.25000000000000000d-01,     3.34201388888888889d-01,
     5    -4.01041666666666667d-01,     7.03125000000000000d-02,
     6    -1.02581259645061728d+00,     1.84646267361111111d+00,
     7    -8.91210937500000000d-01,     7.32421875000000000d-02,
     8     4.66958442342624743d+00,    -1.12070026162229938d+01,
     9     8.78912353515625000d+00,    -2.36408691406250000d+00,
     a     1.12152099609375000d-01,    -2.82120725582002449d+01,
     b     8.46362176746007346d+01,    -9.18182415432400174d+01,
     c     4.25349987453884549d+01,    -7.36879435947963170d+00,
     d     2.27108001708984375d-01,     2.12570130039217123d+02,
     e    -7.65252468141181642d+02,     1.05999045252799988d+03/
      data c(25), c(26), c(27), c(28), c(29), c(30), c(31), c(32),
     1     c(33), c(34), c(35), c(36), c(37), c(38), c(39), c(40),
     2     c(41), c(42), c(43), c(44), c(45), c(46), c(47), c(48)/
     3    -6.99579627376132541d+02,     2.18190511744211590d+02,
     4    -2.64914304869515555d+01,     5.72501420974731445d-01,
     5    -1.91945766231840700d+03,     8.06172218173730938d+03,
     6    -1.35865500064341374d+04,     1.16553933368645332d+04,
     7    -5.30564697861340311d+03,     1.20090291321635246d+03,
     8    -1.08090919788394656d+02,     1.72772750258445740d+00,
     9     2.02042913309661486d+04,    -9.69805983886375135d+04,
     a     1.92547001232531532d+05,    -2.03400177280415534d+05,
     b     1.22200464983017460d+05,    -4.11926549688975513d+04,
     c     7.10951430248936372d+03,    -4.93915304773088012d+02,
     d     6.07404200127348304d+00,    -2.42919187900551333d+05,
     e     1.31176361466297720d+06,    -2.99801591853810675d+06/
      data c(49), c(50), c(51), c(52), c(53), c(54), c(55), c(56),
     1     c(57), c(58), c(59), c(60), c(61), c(62), c(63), c(64),
     2     c(65), c(66), c(67), c(68), c(69), c(70), c(71), c(72)/
     3     3.76327129765640400d+06,    -2.81356322658653411d+06,
     4     1.26836527332162478d+06,    -3.31645172484563578d+05,
     5     4.52187689813627263d+04,    -2.49983048181120962d+03,
     6     2.43805296995560639d+01,     3.28446985307203782d+06,
     7    -1.97068191184322269d+07,     5.09526024926646422d+07,
     8    -7.41051482115326577d+07,     6.63445122747290267d+07,
     9    -3.75671766607633513d+07,     1.32887671664218183d+07,
     a    -2.78561812808645469d+06,     3.08186404612662398d+05,
     b    -1.38860897537170405d+04,     1.10017140269246738d+02,
     c    -4.93292536645099620d+07,     3.25573074185765749d+08,
     d    -9.39462359681578403d+08,     1.55359689957058006d+09,
     e    -1.62108055210833708d+09,     1.10684281682301447d+09/
      data c(73), c(74), c(75), c(76), c(77), c(78), c(79), c(80),
     1     c(81), c(82), c(83), c(84), c(85), c(86), c(87), c(88),
     2     c(89), c(90), c(91), c(92), c(93), c(94), c(95), c(96)/
     3    -4.95889784275030309d+08,     1.42062907797533095d+08,
     4    -2.44740627257387285d+07,     2.24376817792244943d+06,
     5    -8.40054336030240853d+04,     5.51335896122020586d+02,
     6     8.14789096118312115d+08,    -5.86648149205184723d+09,
     7     1.86882075092958249d+10,    -3.46320433881587779d+10,
     8     4.12801855797539740d+10,    -3.30265997498007231d+10,
     9     1.79542137311556001d+10,    -6.56329379261928433d+09,
     a     1.55927986487925751d+09,    -2.25105661889415278d+08,
     b     1.73951075539781645d+07,    -5.49842327572288687d+05,
     c     3.03809051092238427d+03,    -1.46792612476956167d+10,
     d     1.14498237732025810d+11,    -3.99096175224466498d+11,
     e     8.19218669548577329d+11,    -1.09837515608122331d+12/
      data c(97), c(98), c(99), c(100), c(101), c(102), c(103), c(104),
     1     c(105), c(106), c(107), c(108), c(109), c(110), c(111),
     2     c(112), c(113), c(114), c(115), c(116), c(117), c(118)/
     3     1.00815810686538209d+12,    -6.45364869245376503d+11,
     4     2.87900649906150589d+11,    -8.78670721780232657d+10,
     5     1.76347306068349694d+10,    -2.16716498322379509d+09,
     6     1.43157876718888981d+08,    -3.87183344257261262d+06,
     7     1.82577554742931747d+04,     2.86464035717679043d+11,
     8    -2.40629790002850396d+12,     9.10934118523989896d+12,
     9    -2.05168994109344374d+13,     3.05651255199353206d+13,
     a    -3.16670885847851584d+13,     2.33483640445818409d+13,
     b    -1.23204913055982872d+13,     4.61272578084913197d+12,
     c    -1.19655288019618160d+12,     2.05914503232410016d+11,
     d    -2.18229277575292237d+10,     1.24700929351271032d+09/
      data c(119), c(120)/
     1    -2.91883881222208134d+07,     1.18838426256783253d+05/
c***first executable statement  zunik
      if (init.ne.0) go to 40
c-----------------------------------------------------------------------
c     initialize all variables
c-----------------------------------------------------------------------
      rfn = 1.0d0/fnu
c-----------------------------------------------------------------------
c     overflow test (zr/fnu too small)
c-----------------------------------------------------------------------
      test = d1mach(1)*1.0d+3
      ac = fnu*test
      if (abs(zrr).gt.ac .or. abs(zri).gt.ac) go to 15
      zeta1r = 2.0d0*abs(log(test))+fnu
      zeta1i = 0.0d0
      zeta2r = fnu
      zeta2i = 0.0d0
      phir = 1.0d0
      phii = 0.0d0
      return
   15 continue
      tr = zrr*rfn
      ti = zri*rfn
      sr = coner + (tr*tr-ti*ti)
      si = conei + (tr*ti+ti*tr)
      call zsqrt(sr, si, srr, sri)
      str = coner + srr
      sti = conei + sri
      call zdiv(str, sti, tr, ti, znr, zni)
      call zlog(znr, zni, str, sti, idum)
      zeta1r = fnu*str
      zeta1i = fnu*sti
      zeta2r = fnu*srr
      zeta2i = fnu*sri
      call zdiv(coner, conei, srr, sri, tr, ti)
      srr = tr*rfn
      sri = ti*rfn
      call zsqrt(srr, sri, cwrkr(16), cwrki(16))
      phir = cwrkr(16)*con(ikflg)
      phii = cwrki(16)*con(ikflg)
      if (ipmtr.ne.0) return
      call zdiv(coner, conei, sr, si, t2r, t2i)
      cwrkr(1) = coner
      cwrki(1) = conei
      crfnr = coner
      crfni = conei
      ac = 1.0d0
      l = 1
      do 20 k=2,15
        sr = zeror
        si = zeroi
        do 10 j=1,k
          l = l + 1
          str = sr*t2r - si*t2i + c(l)
          si = sr*t2i + si*t2r
          sr = str
   10   continue
        str = crfnr*srr - crfni*sri
        crfni = crfnr*sri + crfni*srr
        crfnr = str
        cwrkr(k) = crfnr*sr - crfni*si
        cwrki(k) = crfnr*si + crfni*sr
        ac = ac*rfn
        test = abs(cwrkr(k)) + abs(cwrki(k))
        if (ac.lt.tol .and. test.lt.tol) go to 30
   20 continue
      k = 15
   30 continue
      init = k
   40 continue
      if (ikflg.eq.2) go to 60
c-----------------------------------------------------------------------
c     compute sum for the i function
c-----------------------------------------------------------------------
      sr = zeror
      si = zeroi
      do 50 i=1,init
        sr = sr + cwrkr(i)
        si = si + cwrki(i)
   50 continue
      sumr = sr
      sumi = si
      phir = cwrkr(16)*con(1)
      phii = cwrki(16)*con(1)
      return
   60 continue
c-----------------------------------------------------------------------
c     compute sum for the k function
c-----------------------------------------------------------------------
      sr = zeror
      si = zeroi
      tr = coner
      do 70 i=1,init
        sr = sr + tr*cwrkr(i)
        si = si + tr*cwrki(i)
        tr = -tr
   70 continue
      sumr = sr
      sumi = si
      phir = cwrkr(16)*con(2)
      phii = cwrki(16)*con(2)
      return
      end
