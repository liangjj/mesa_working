*
*     ------------------------------------------------------------------
*               d a t a
*     ------------------------------------------------------------------
*
*       data concerning the number of configurations (ncfg), the number
*   and type of electrons in each  configuration,  as  well  as  data
*   associated with the energy expression are read and stored.
*
      SUBROUTINE data
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
*
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
      CHARACTER config*50,atom*6,term*6,ans*1,string*50
      COMMON /LABEL/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      LOGICAL fail,omit,rel,all,trace,first,strong,done,orderd
      COMMON /TEST/fail,omit,rel,all,trace
*
      INTEGER ind(nwd)
      CHARACTER*3 el,el1,el2,elcsd(18)
      CHARACTER*1 aster,w
      data aster/'*'/
*
    1 format(18(1x,a3))
    7 format(a3,f6.0,i3,i3,f3.1)
*
*  *****  read 'atom' card
*
5     write(iscw,'(a/a)') ' Enter ATOM, TERM, Z',
     :   ' Examples: O,3P,8. or Oxygen,AV,8.'
      read(5,'(a50)') string
      i = index(string,',')
      if ( i .eq. 0) then
          write(iscw,*)' ATOM, TERM, and Z must be separated by commas '
          go to 5
      end if
      atom = string(1:i-1)
      j = index(string(i+1:),',')
      if ( j .eq. 0) then
          write(iscw,*)' ATOM, TERM, and Z must be separated by commas '
          go to 5
      end if
      term = string(i+1:i+j-1)
      read(string(i+j+1:), '(f3.0)') z
*
*  *****  input COMMON closed shells
*
   2  write(iscw,*)
      write(iscw,'(a,a)')' List the closed shells in fields indicated',
     :        ' (blank line if none)'
      write(iscw,'(a)') ' ... ... ... ... ... ... ... ... etc.'
      read(5,1) (elcsd(i),i=1,18)
*
*  *****  input the configuration
*
      write(iscw,'(/a,a/a)')' Enter electrons outside closed shells ',
     : '(blank line if none)',' example: 2s(1)2p(3)'
      read(5,'(a)')  string
      call reform(string, config)
*
*      determine the number of closed shells
*
      i = 0
      ss = 0.d0
      lmax = 0
   12 if (elcsd(i+1) .ne. '   ') then
         i = i+1
         el(i) = elcsd(i)
         j = 3
         if (el(i)(1:1) .ne. ' ') j = 2
         l(i) = lval(el(i)(j:j))
         n(i) = ichar(el(i)(j-1:j-1)) - ichar('1') + 1
         ifull = 2*(2*l(i)+1)
         sum(i) = ifull
         s(i) = ss + ifull/2
         ss = ss + ifull
         meth(i) = 1
         acc(i) = 0.d0
         ind(i) = 0
         if (iuf .ne. 0)  ind(i) = -1
         if( i .lt. 18) go to 12
         stop ' too many closed shells: max = 18'
      end if
      nclosd = i
*
*  *****  determine the other electrons
*
      maxorb = nclosd
      string = config
      j = 2
      i = 0
 16   if (string(j:j+2) .ne. '   ' ) then
*
*  --------- an electron has been found; is it a new one?
*
         i = i+1
         if (i .gt. (5)) stop ' too many shells: max= (5)'
         el1 = string(j:j+2)
         k = nclosd + 1
 17      if (k .le. maxorb) then
            if ( el(k) .ne. el1 ) then
               k = k+1
               if (k .gt. (nwd)) then
                  write(iscw,'(a,i4)')' too many electrons: max =',nwd
                  go to 2
               else
                  go to 17
               end if
            end if
         else
*
*  ------------  a new electron has been found; add it to the list
*
            maxorb = k
            el(maxorb) = el1
            read(string(j+4:j+7),'(f4.0)') sum(k)
         end if
         j = j+10
         if (j .lt. (50)) go to 16
      end if
*
*  -----  the list of electrons has been determined
*
      write(iscw,19) maxorb,(el(j),j=1,maxorb)
   19 format(/' there are ',i3,' orbitals as follows:'/(1x,18(1x,a3)))
      nwf = maxorb
      if (nit .lt. 0) nit=nwf
21    write(iscw,'(/a,a)') ' orbitals to be varied: ',
     :  'all/none/=i (last i)/comma delimited list/h'
      read(5, '(a)') string
      if (string(1:1) .eq. 'h' .or. string(1:1) .eq. 'h') then
       call help(1)
       go to 21
      else if (string(1:3) .eq. 'all' .or. string(1:3) .eq. 'all') then
         nit = nwf
      else if (string(1:4).eq.'none' .or. string(1:4).eq.'none') then
         nit = 0
      else if (index(string,'=') .ne. 0) then
         j = index(string,'=')
	 jj = index(string,' ')
	 if ( jj .eq. j+2) then
	    read(string(j+1:j+1),'(i1)') nit
	 else if (jj .eq. j+3) then
	    read(string(j+1:j+2),'(i2)') nit
	 else
	    write(iscw,'(a,a)') ' nit must be specified by one or two',
     :        'digits immediately following = sign (no blanks): Renter'
	    go to 21
	 end if
      else
         nit = 0
         j = 1
22       next = index(string(j:),',')
*
*        ***  search for last electron LABEL which need not be followed
*             by a comma
*
         if (next .eq. 0 .and. string(j:j+2) .ne. '   ')
     :       next = index(string(j+1:),' ') + 1
         if (next .ge. 1) then
            if (next .eq. 4) then
               el1 = string(j:j+2)
            else if (next .eq. 3) then
               el1 = ' '//string(j:j+1)
            else
               write(iscw,*)
     :              ' electron LABELs must be separated by commas;'
               write(iscw,*)' each LABEL must contain 2 or 3 CHARACTERs'
               go to 21
            end if
            call reord(el,el1,nwf,ierr)
            if (ierr .eq. 0) then
               nit = nit + 1
               j = j + next
               if (j .lt. 72) go to 22
            else
               write(iscw,*) ' case must match as well as position of',
     :                  ' imbedded blanks'
               write(iscw,*) ' for 3rd CHARACTER of LABEL to be blank',
     :                 ' follow blank with comma'
     
               go to 21
            end if
         end if
      end if
*
      ib = nwf - nit + 1
      if (nit .ne. 0) then
23       write(iscw,'(/a)') ' default electron parameters ? (y/n/h) '
         read(5,'(a)') ans
         if ( ans .eq. 'h' .or. ans .eq. 'h' ) then
            call help(2)
            go to 23
         end if
      else
         ans = 'y'
*
      end if
      if ( ans .ne. 'y' .and. ans .ne. 'y') write(iscw,'(a,a)')
     :   ' s, ind, meth, acc for non-closed shell electrons: '
      do 20 i = nclosd+1,nwf
         if ( ans .eq. 'y' .or. ans .eq. 'y') then
            s(i) = ss + (sum(i)-1.d0)/2.d0
            ss = ss + sum(i)
            meth(i) = 1
            acc(i) = 0.d0
            ind(i) = 0
            if (iuf .ne. 0)  ind(i) = -1
         else
            write(iscw,'(a,a)') el(i),':  '
            read(5,*) s(i),ind(i),meth(i),acc(i)
         end if
         j = 2
         if (el(i)(1:1) .eq. ' ') j = 3
         l(i) = lval(el(i)(j:j))
         n(i) = ichar(el(i)(j-1:j-1)) - ichar('1') + 1
 20   continue
*
*  *****  define all orbitals in the configuration to be orthogonal
*
       do 34 i = 1,nwf
         e(i,i) = 0.d0
	 lmax = max0(lmax,l(i))
         do 33 j = 1,i-1
            e(i,j) = 0.d0
            if (l(i) .eq. l(j)) e(i,j) = 1.d-5
            e(j,i) = e(i,j)
   33       continue
   34 continue
      kmax = 2*lmax
      ib = nwf - nit + 1
      no = ns
      strong = .false.
      write (3,62) atom,term,z,(el(i),int(sum(i)),i=1,nclosd)
62    format(1h1///9x,33hhartree-fock wave functions for  ,2a6,4h z =,
     : f5.1//14x,'core =',5(1x,a3,'(',i4,')')/(20x,5(1x,a3,'(',i4,')')))
      write (3,'(5x,a15,a50)')  'configuration =',config
68    continue
      write (3,71)
71    format(//9x,10hinput data/9x,10h----- ----//13x,13hwave function,
     :   11h  procedure/17x,22hnl  sigma meth acc opt///)
      do 79 i = 1,nwf
         write (3,78) i,el(i),n(i),l(i),s(i),meth(i),acc(i),ind(i)
78       format(i8, 2x,a3,2i3,f7.1,i4,f4.1,i4)
79    continue
      omit = .not. strong
*
      call array
      call enexpr(term, done)
      if (.not. done) then
*
*  ---  case needs additional data
*
         write(iscw,85)
85       format(/' The program could not derive the energy expression'/
     :           ' Select one of the following options and enter:'/
     :           '    1  re-enter the term and configuration'/
     :           '    2  enter the deviations from eav as input'/
     :           '    3  stop'/)
         read(5,*) iselec
         go to (5,86,99) iselec
86       write(iscw,83)
83       format(/' Input data for deviations from the average energy'/
     :   ' first fk integrals, then gk integrals in indicated format'/
     :   '  cc.ccccccccccfkk(el1,el2)  - terminate each list with an *',
     :   ' in the f column')
         first = .true.
*
*  *****  read 'fk' and 'gk' cards, omitting the header if a file
*
82       read (5,84) cfg,w,kfg,el1,el2
84       format(f14.8,a1,i2,1x,a3,1x,a3)
         if ( w .ne. aster ) then
            call eptr(el,el1, ifg, *99)
            call eptr(el,el2, jfg, *99)
            call add(cfg,kfg,ifg,jfg,first)
            call add(cfg,kfg,jfg,ifg,first)
            go to 82
         else if (first) then
            first = .false.
            go to 82
         end if
      end if
*
*  *****  compute the initial array and initial radial functions
*
      call initsp
      call wavefn(ind)
*
*      ... define an order for the functions to be iterated
*
      do 90 jp = 1,nwf
         iord(jp) = jp
90    continue
91    orderd = .true.
      do 92 jp = ib,nwf-1
         n1 = n(iord(jp))
         l1 = l(iord(jp))
         n2 = n(iord(jp+1))
         l2 = l(iord(jp+1))
         if (n1.gt.n2 .or. (n1.eq.n2 .and. l1.gt.l2)) then
            itemp = iord(jp)
            iord(jp) = iord(jp+1)
            iord(jp+1) = itemp
            orderd = .false.
         end if
92    continue
      if (.not. orderd) go to 91
      return
99    stop
      end
