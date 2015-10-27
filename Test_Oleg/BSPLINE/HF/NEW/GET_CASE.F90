!=======================================================================
  SUBROUTINE get_case(z) 
!=======================================================================
!  Get information from the user about the problem to be solved
!----------------------------------------------------------------------
!
    USE atomic_state
    USE periodic_table
    IMPLICIT NONE

    CHARACTER(LEN=1) :: ans
    LOGICAL :: ok

    write(*,'(A)',ADVANCE='NO') 'Atomic Number                         : '
    read(*,*) z
    atom = atom_label(nint(z))
    
    DO
      ! obtain closed shells
      call closed_shells
  
      ! obtain configuration
      write(*,'(A)',ADVANCE='NO') 'Configuration [ Ex. 2s(1)2p(3) ]      : '
      read(*,*) config
  
      ! obtain term
      write(*,'(A)',ADVANCE='NO') 'Term [ Ex. 3P or AV ]                 : '
      read(*,*) term
  
      ! determine the list of electrons
      call find_electrons
      IF (ok) EXIT
    END DO
    
    ! determine orbitals to be varied
    DO
      write(*,'(A,I3,A)') 'There are',nwf,' orbitals as follows:'
      write(*,'(4x,18(1X,A3))') el
      write(*,'(A/A)',ADVANCE='NO') 'Indicate those to be varied', &
                                '[Ex. all/none/=i/comma delimited list/h]:'
      call orbitals_varied
      if (ok) EXIT
    END DO

    ! determine orbital parameters
    write(*,'(A)',ADVANCE='NO') 'Default orbital parameters  (y/n)    : '
    read(*,*) ans
    call orbital_parameters(nwf,el,nclosd,ans)
    
    ! get the energy expression
    call get_energy_expression
   
    CONTAINS 

    !===================================================================
      SUBROUTINE closed_shells
    !===================================================================
    !  Get information from the user about closed shells
    !------------------------------------------------------------------
    !
      IMPLICIT NONE
      INTEGER :: i,ii, core_type

      i = 0
      DO
	IF (number_of_electrons(i+1) > z) then
	  EXIT
        ELSE
          i = i + 1
        END IF
      END DO

      Write(*,*) 'Possible cores are:
      Write(*,'(3X,I3,A,2X,A)') (ii,'.', shells(ii),ii=1,i)
      Write(*,*,ADVANCE='no') 'Enter the number for your selection: '
      Read(*,*) core_type
      i = 0
      DO
	IF (shell_number(i+1) > number_or_electrons(core_type)) then
	  EXIT
        ELSE
          i = i + 1
        END IF
      END DO
      nclosed = i
     
    END SUBROUTINE closed_shells

    !===================================================================
      SUBROUTINE find_electrons
    !===================================================================
    !  Analyze the configuration and report the list of electrons
    !------------------------------------------------------------------
    !
      IMPLICIT NONE
      INTEGER :: i,is,ipl,ipr,iel
      CHARACTER(LEN=3) :: ell
      CHARACTER(LEN=80) :: sconfig

      ok = .true.
      write(*,'(A,F6.0,A,A)')'Atom         : ', z,atom
      write(*,'(A,18(1X,A2)) 'Closed Shells: ',(shell_order(1:nclosd))
      Write(*,*)             'Configuration: ', trim(config), '  ',term
      ! count the number of shells in the configuration
      i = 0
      is = 1
      DO
       ipl = index(config(is:),'(')
       IF (ipl > 0) then 
         ! there was a shell
         IF (ipl-is > 3) then
           write(*,*) 'ERROR: Electron label ', config(is:ipl-1), &
                      ' has more than 3 characters'
           ok = .false.
           RETURN
         END IF
	 ell = config(is:ipl-1)
	 ! check this is not in closed shells 
         DO iel = 1,nclosd
           IF (ell == shell_order(iel)) then
             write(*,*) 'ERROR: Electron', ell, ' appears both in ', &
                        config, 'and the closed shells'
             ok = .false.
             RETURN
           END IF
         END DO
	 ipr= index(config(ipl+1:),')')
         ! check for matching parentheses 
         IF (ipr = 0 ) then
           write(*,*) 'ERROR: Matching parentheses not found in configuration'
           ok = .false.
           RETURN
         ELSE IF (ipr-ipl .gt. 5) then
           write(*,*) 'ERROR:'                                       &
             'Occupation number should not exceed 4 characters'
           ok = .false.
           RETURN 
         END IF
	 i = i+1
	 is = ipr+1
       ELSE
	 EXIT
       END IF
      END DO
      nwf = nclosd+i
 
      ! at this point the orbital data should be correct
      call allocate_atomic_state

      ! set the list of orbitals
      el(1:nclosd) = shell_order(1:nclosd)
      el(1:nclosd) = ADJUSTR(el(1:nclosd))
      is = 1
      call reform(TRIM(ADJUSTL(config)), sconfig)
      read(sconfig,'(8(1x,a3,6x)))') el(nclosd+1:nwf)
 
    END SUBROUTINE find_electrons
      
  !===================================================================
    SUBROUTINE orbitals_varied
  !===================================================================
  !  Determine which orbitals are to be varied
  !------------------------------------------------------------------
  !
     IMPLICIT NONE
     CHARACTER(LEN=72) :: string
     
     ok = .true. 
     read(5, '(a)') string
     if (string(1:1) .eq. 'h' .or. string(1:1) .eq. 'H') then
       call help(1)
       ok = .false.
       RETURN
     else if (string(1:3) .eq. 'all' .or. string(1:3) .eq. 'ALL') then
       nit = nwf
     else if (string(1:4).eq.'none' .or. string(1:4).eq.'NONE') then
       nit = 0
     else if (index(string,'=') .ne. 0) then
       j = index(string,'=')
       jj = index(string,' ')
       if ( jj .eq. j+2) then
	 read(string(j+1:j+1),'(i1)') nit
       else if (jj .eq. j+3) then
	 read(string(j+1:j+2),'(i2)') nit
       else
	 write(iscw,'(a,a)') ' nit must be specified by one or two' &
              'digits immediately following = sign (no blanks)
	 ok = .false.
         RETURN
       end if
     else
       nit = 0
       j = 1
       DO
         next = index(string(j:),',')
         ! search for last electron LABEL which need not be followed by a comma 
         if (next .eq. 0 .and. string(j:j+2) .ne. '   ') &
             next = index(string(j+1:),' ') + 1
         if (next .ge. 1) then
            if (next .eq. 4) then
               el1 = string(j:j+2)
            else if (next .eq. 3) then
               el1 = ' '//string(j:j+1)
            else
               write(iscw,*) &
                 ' ERROR: electron LABELs must be separated by commas;'
               write(iscw,*)' each LABEL must contain 2 or 3 CHARACTERs'
               ok = .false.
               RETURN
            end if
            call reord(el,el1,nwf,ierr)
            if (ierr .eq. 0) then
               ! no error has occurred
               nit = nit + 1
               j = j + next
               if (next > 72) EXIT
            else
               write(*,*) 'ERROR: Electron ', el1, ' not found in the list'
               write(*,*) 'Case and position of imbedded blanks must match'
               ok = .false.
               RETURN
            end if
         END DO
         end if
      end if
    END SUBROUTINE orbitals_varied
      
  !===================================================================
    SUBROUTINE orbital_parameters
  !===================================================================
  !  Set orbital parameters
  !------------------------------------------------------------------
  !
     IMPLICIT NONE
     INTEGER :: i
     REAL(KIND=8) :: ss
      
     ! Set occupation number of non-core orbitals
     Read(sconfig,'(8(5x,F4.0,1X)') (sum(nclosd+1:nwf)
     ! Determine nl values and related quantities
     
     ss = 0
     lmax = 0
     DO i = 1,nwf
       l(i) = lval(el(i)(3:3)
       read(el(i)(1:2), '(I2)')  n(i)
       IF (i <= nclosd) sum(i) = 2*( 2*l(i) + 1)
       s(i) = ss + sum(i)/2
       ss = ss + sum(i)
       meth(i) = 1
       acc(i) = 0
       lmax = max(lmax,l(i))
     END DO

     ! For HF calculations we set
     strong = .false.

     kmax = 2*lmax

    END SUBROUTINE orbital_parameters
     
  !===================================================================
    SUBROUTINE get_energy_expression
  !===================================================================
  !  Obtain the data that defines the HF equations from an
  !  energy expression stored in a coef list and an ijptr array.
  !------------------------------------------------------------------
  !
     IMPLICIT NONE
     INTEGER :: i,j, ip
     REAL(KIND=8) ::C, CA, CB`

     ! Find the size of the coef array
     ip = 0
     DO i=1,nwf
       DO j= 1,nwf
	 ! direct contribution
         ip = ip +1
         if (i == j ) ip = ip + l(i)
         ! exchange contribution
         ip = ip + min(l(i),l(j))
       END DO
     END DO

     ALLOCATE ( coef(ip) )

     ip = 0
     do i = nclosd+1,nwf
       isumi = sum(i)
       dsumi = sum(i) - isumi
       do j = nclosd+1,nwf
         isumj = sum(j)
         dsumj = sum(j) - isumj
         if ( i .ne. j) then
           c = sum(j)
           if (dsumi .ne. 0.d0 .and. dsumj .ne. 0.d0) &
             c = (dsumi*(isumi+1)*isumj +a            &
                  dsumj*(isumj+1)*isumi)/sum(i)
           else
             c = sum(i) - 1.d0
             if (dsumi .ne. 0.d0)                     &
               c = (isumi*(sum(i)+dsumi-1))/sum(i)
           end if

           ijptr(i-nclosd,j-nclosd) = ip

!           ...direct contribution

           do k = 0,2*min0(l(i),l(j)),2
             ip = ip + 1
             if (k .eq. 0) then
               coef(ip) = c
             else if (i .eq. j) then
               coef(ip) = -c*ca(l(i),k)
             end if
           end dO

!           ... exchange contribution

           if (i .ne. j) then
             do k = abs(l(i)-l(j)),l(i)+l(j),2
               ip = ip + 1
               coef(ip) = -c*cb(l(i),l(j),k)
             END Do
           end if
         END Do
       END Do

       CALL enexpR

    END SUBROUTINE get_energy_expressioN
     
  !===================================================================
    SUBROUTINE enexpr
  !===================================================================
  ! Determine the deviations from the average energy for the term
  ! dependent energy expression for the following cases:
  !      i) an open p- or d-shell
  !     ii) a single electron or hole, any l
  !    iii) an s-electron and a single electron, any l
  !     iv) an s-electron and an open p- or d-shell
  !      v) an open p-shell and a single electron, any l
  !------------------------------------------------------------------
  !
      IMPLICIT NONE

      iNTEGER, DIMENSION(5) ::  sumtab(5)
      INTEGER, DIMENSION(11) :: partab,ptrtab,plval
      INTEGER, DIMENSION(54) :: ltaB
      INTEGER, DIMENSION(2)  :: nos(2)
      INTEGER :: pacval,sp,ps1,ps2
 
      ... fint, gint1, and gint2 are coefficients of polynomials
          in l, tabulated by slater,
 
 
      cHARACTER(LEN=1) :: parch(11)

      !... coefficients of f2 integrals for p(n)l(1) configurations
      INTEGER, DIMENSION(3,54) :: fint = &
        (/2,-1,0,-4,-4,3,2,5,3,2,-1,0,-4,-4,3,2,5,3,             &
         -2,1,0,4,4,-3,-2,-5,-3,-2,1,0,4,4,-3,-2,-5,-3,          &
         4,-2,0,-2,-11,6,-4,-4,15,-2,7,15,4,10,6,0,0,0,          &
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  &
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  &
         2,-1,0,-4,-4,3,2,5,3,2,-1,0,-4,-4,3,2,5,3,              &
         -4,2,0,2,11,-6,4,4,-15,2,-7,-15,-4,-10,-6,0,0,0,        &
         -2,1,0,4,4,-3,-2,-5,-3,-2,1,0,4,4,-3,-2,-5,-3/)
 
      !... coefficients of g(l-1) integrals
      INTEGER, DIMENSION(3,54) :: gint1 `= &
       (/-10,5,0,2,11,-6,2,-1,-6,14,-7,0,2,-13,6,2,-1,6,               &
         -8,4,0,-8,4,0,4,10,0,10,-5,0,10,-5,0,4,-8,0,                  &
         -8,4,0,-2,13,-6,2,11,-12,4,4,-12,4,-2,-12,0,0,0,-6,3,0,10,-5,0,&
         -6,3,0,-6,3,0,-4,8,3,0,12,3,6,9,0,-6,3,0,6,21,12,12,12,-27,    &
         12,-6,27,6,-15,-12,-6,3,0,0,6,-3,0,0,-3,6,-3,0,0,-6,3,12,6,3,  &
         -4,2,0,-4,2,0,-4,2,0,-4,2,0,14,11,-9,14,-7,-9,-4,-2,0,-4,2,0,  &
         -2,7,3,2,11,3,8,8,0,0,0,0,-2,1,0,-2,1,0,-2,1,0,-2,1,0,-2,1,0,  &
         22,13,0/)
 
      !... coefficients of g(l+1) integrals
      INTEGER, DIMENSION(3,54) :: gint2 `= &
       (/2,5,-3,2,-7,-15,-10,-25,-15,2,5,9,2,17,21,14,35,21,       &
         4,-2,-6,-8,-20,-12,-8,-20,-12,4,16,12,10,25,15,10,25,15,  &
         4,10,0,4,4,-12,2,-7,-21,-2,-17,-21,-8,-20,-12,0,0,0,      &
         -6,-15,-9,10,25,15,6,3,-3,0,-12,-9,-4,-16,-9,-6,-15,-9,   &
         -6,-15,-9,6,27,9,12,30,-9,12,12,-27,6,-9,-27,-6,-15,-9,   &
         0,0,-3,0,-6,-9,-6,-15,-9,12,18,9,0,6,9,6,15,9,            &
         -4,-10,-6,-4,-10,-6,-4,-10,-6,14,35,12,14,17,-6,-4,-10,-6,&
         8,8,0,2,-7,-6,-2,-11,-6,-4,-10,-6,-4,-10,-6,0,0,0,        &
         -2,-5,-3,-2,-5,-3,-2,-5,-3,22,31,9,-2,-5,-3,-2,-5,-3/)
 
      !.. encoded term value --   s = ltab/10
      !                       lterm = l + (ltab mod 10 - 5)
      !   example: ltab = 36 with l = 2  is 3f
 
      data ltab/36,35,34,16,15,14,46,45,44,26,25,24,27,26,25,24,
     :   23,25,55,35,37,36,35,34,33,17,16,15,14,13,36,35,34,16,15,14,
     :   46,45,44,26,25,24,27,26,25,24,23,25,36,35,34,16,15,14/
      data sumtab/1,4,7,10,11/
      data partab/2,3,1,1,4,2,2,3,1,1,2/
      data ptrtab/6,12,17,18,20,30,36,42,47,48,54/
      data plval/1,1,2,0,0,2,1,1,2,0,1/
      data parch/'p','p','d','s','s','d','p','p','d','s','p'/
 
       
      sl = term
      senor = ' '
 
      !  convert lowercase l symbol to uppercase
 
      if (sl(2:2).gt.'a' .and. sl(2:2).lt.'z')                      &
          sl(2:2) = char(ichar(sl(2:2)) + ichar('a') - ichar('a'))
 
      !  determine if fk or gk data needs to be input
 
      il = 0
      is = 0
      j = 1
      do I = nclosd+1, nwf
         if (sum(i) .ne. 4*l(i)+2 .and. sum(i) .ne. 0.d0) then
            if (j.gt.2) then
               if (sl .ne. 'av' .and. sl .ne. 'av') then
                  done=.false.
               else
                  done = .true.
               end if
               return
            endif
            nos(j) = i
            j = j + 1
            if (l(i) .eq. 0 .and. is .eq. 0) then
               is = is + 1
               iis = i
            else
               il = il + 1
               iil = i
            end if
         end if
      END Do
      if (sl .ne. 'av' .and. sl .ne. 'av' .and. is+il.ne.0) then
         done = .false.
         c = 0.d0
         if (is+il .le. 2 .and. il .le. 1) then
            if (is .eq. 0 .and. il .eq. 1) then
	      Do
               call looktm(l(iil),sl,senor,sum(iil),ip,nsl)
	       if (nsl > 1) then
		  write(iscw,*)' ambiguous term: enter seniority'
		  read (5,'(a1)') senor
	       else
		  EXIt
	       end if
	      END Do
              call dev(iil,l(iil),sum(iil),ip,done)
            else if (is .eq. 1 .and. il .eq. 1) then
               slm = sl
	       slp = sl
               slm(1:1) = char(ichar(slm(1:1)) - 1)
               slp(1:1) = char(ichar(slp(1:1)) + 1)
	       call looktm(l(iil),slm,senor,sum(iil),ipm,nslm)
	       call looktm(l(iil),slp,senor,sum(iil),ipp,nslp)
	       if (nslm+nslp .eq. 0) then
	          done = .false.
	          return
	       else if (nslm .eq. 1 .and. nslp .eq. 0) then
		  sl = slm
		  ip = ipm
	       else if (nslm .eq. 0 .and. nslp .eq. 1) then
		  sl = slp
		  ip = ipp
	       else if (nslm .eq. 1 .and. nslp .eq. 1) then
		 Do
                  write(iscw,'(a,a3,a,a3)')                          &
                       ' ambiguous l**n term: enter',slm, ' or ',slp 
	          read(5,'(a2)') sl
		  if (sl .eq. slm) then
		     ip = ipm
		     EXIt
		  else if (sl .eq. slp) then
		     ip = ipp
		     EXIt
		  else
	             write(iscw,*) ' term not allowed: re-enter'
		  end if
		 END Do
	       else 
		 Do
                  write(iscw,'(a,a)') ' ambiguous l**n parent term:', &
                         'enter term and seniority'
	          read(5,'(a2,a1)') sl, senor
	          call looktm(l(iil),sl,senor,sum(iil),ip,nsl)
	          if (nsl <> 1) then
		    write(iscw,'(a,a3,a,a3,a)') ' allowed terms are ',  &
                          slm, ' or ', slp,' plus seniority'
		  ELSE
		     EXIt
	          end if
		 END Do
	       end if
               call dev(iil,l(iil),sum(iil),ip,done)
               if (done ) then
                  sp = ichar(sl(1:1)) - ichar('0')
                  csp = (sp - 1)/2.
		  if (sl .eq. slm) then
                     c = -csp/(2*l(iil)+1)
                  else
                     c = (csp + 1)/(2*l(iil)+1)
                  end if
                  call add(c,l(iil),iis,iil,.false.)
                  call add(c,l(iil),iil,iis,.false.)
               end if
            else if (is .eq. 1 .and. il .eq. 0) then
               done = .true.
            end if
         else
            if (((l(nos(1)).eq.1).and.(sum(nos(2)).eq.1.d0)).or.    &
                ((l(nos(2)).eq.1).and.(sum(nos(1)).eq.1.d0))) then
               if (l(nos(1)).eq.1.and.sum(nos(2)).eq.1.d0) then
                  isump=sum(nos(1))
                  np = nos(1)
                  nl = nos(2)
               else
                  isump = sum(nos(2))
                  np = nos(2)
                  nl = nos(1)
               endif
               sp=ichar(sl(1:1))-ichar('0')
               lp=lval(sl(2:2))
               ps1=sp+1
               ps2=sp-1
               if (isump.eq.1) then
                  iptr1=1
               else
                  iptr1=sumtab(isump-1)+1
               end if
               iptr2=sumtab(isump)
               nomach=0
               call lookup(partab,iptr1,iptr2,ind,nomach,ps1)
               call lookup(partab,iptr1,iptr2,ind,nomach,ps2)
               psl(1:1)=char(partab(ind)+ ichar('0'))
               psl(2:2)=parch(ind)
               if (nomach.gt.1) then
                 write(iscw,*)' ambiguous parent case'
		 Do
                  write(iscw,*)' enter the sl term for p(n) subshell'
                  read (5,'(a)')psl
                  if (psl(2:2).gt.'a'.and.psl(2:2).lt.'z')           &
                         psl(2:2)=char(ichar(psl(2:2))+ichar('a')    &
                                                      -ichar('a'))  
                  ps1=ichar(psl(1:1))-ichar('0')
                  ps2=lval(psl(2:2))
                  call lookup(plval,iptr1,iptr2,ind,nomach,ps2)
                  if (.NOT.(nomach.ne.1).and.(partab(ind).ne.ps1)) EXIT
                 END Do
               end if
               if (isump.eq.1) then
                  iptr1=1
               else
                  iptr1=ptrtab(ind-1)+1
               end if
               iptr2 = ptrtab(ind)
               lv=l(nl)
               pacval=sp*10+lp-lv+5
               nomach=0
               call lookup(ltab,iptr1,iptr2,ind,nomach,pacval)
               if (nomach.ne.1) then
                  done=.false.
                  return
               endif
               val1=((fint(1,ind)*lv+fint(2,ind))*lv+fint(3,ind))    &
                        /(5.d0*(2*lv-1)*(2*lv+3))
               val2=((gint1(1,ind)*lv+gint1(2,ind))*lv+gint1(3,ind)) &
                        /(2.d0*(2*lv+1)*(2*lv-1)**2)
               val3=((gint2(1,ind)*lv+gint2(2,ind))*lv+gint2(3,ind)) &
                        /(2.d0*(2*lv+1)*(2*lv+3)**2)
 
               !  add contributions from between p-subshell and l-electron
 
               call add(val1,2,np,nl,.true.)
               call add(val1,2,nl,np,.true.)
               call add(val2,lv-1,np,nl,.false.)
               call add(val2,lv-1,nl,np,.false.)
               call add(val3,lv+1,np,nl,.false.)
               call add(val3,lv+1,nl,np,.false.)
 
               !... add deviations for p-subshell
 
               call looktm(1,psl,' ',sum(np),ip,nsl)
               call dev(np,1,sum(np),ip,done)
            else
               done = .false.
            end if
         end if
      else
         done = .true.
      end if
     END SUBROUTINE enexpr
 
  !=====================================================================
    SUBROUTINE looktm( l, sl, sen, q, ip, nsl)
  !=====================================================================
  !
  !    add the deviations to the average energy for a partially filled
  !      p- or d- shell
  ! 
  !--------------------------------------------------------------------
  !
    IMPLICIT NONE
      CHARACTER sl*2, sen*1
      INTEGER iptr(5)
      CHARACTER*3 terms(51)
      data    iptr/6,11,19,35,51/
      data        terms/'3p2','1d2','1s0','4s3','2d3','2p1',
*             .. d2 and d3 terms
     :      '3f2','3p2','1g2','1d2','1s0','4f3','4p3','2h3','2g3',
     :      '2f3','2d1','2d3','2p3',
*            ... d4 terms ...
     :       '5d4','3h4','3g4','3f2','3f4','3d4','3p2','3p4',
     :       '1i4','1g2','1g4','1f4','1d2','1d4','1s0','1s4',
*            ... d5 terms ...
     :       '6s5','4g5','4f3','4d5','4p3','2i5','2h3','2g3',
     :       '2g5','2f3','2f5','2d1','2d3','2d5','2p3','2s5'/
 
*
*  --- search for a partially unfilled p- or d-shell
*
      n = q
      if (n .gt. 2*l+1) n = 4*l+2 - n
      ip = 0
      nsl = 0
      if (n .gt. 1  .and. l .le. 2) then
         if (l .eq. 1) then
            ibegin = 1
            iend = 6
         else
            ibegin = iptr(n-1) + 1
            iend = iptr(n)
         end if
1        i = ibegin
         DO
         if (sl .eq. terms(i)(1:2)) then
            if (sen .eq. ' ' .or. sen .eq. terms(i)(3:3)) then
	       nsl = nsl + 1
	       ip = i
            end if
         end if
         i = i+1
         if (i > iend) EXIT
         END DO
      else if ( n .eq. 1 .and. sl(1:1) .eq. '2') then
	 nsl = 1
      end if
    END SUBROUTINE looktm
 
  !=======================================================================
     SUBROUTINE dev(iel, l, q, i, done) 
  !=======================================================================
  !
  !     add the deviations to the average energy for a partially filled
  !       p- or d- shell
  !----------------------------------------------------------------------
  !
      IMPLICIT NONE
      INTEGER, INTENT(IN):: iel, l, q, i
      LOGICAL, INTENT(OUT) :: done 
      INTEGER f2pp(6), f2dd(45), f4dd(45)
      data    f2pp/-3,3,12,-9,0,6/
      data    f2dd/-58,77,50,-13,140,
!             ... d3 coefficients
     :       -93,42,-12,-57,123,105,69,-12,
!             ... d4 coefficients
     :        -105,-69,-24,66,12,39,21,57,
     :        -51,30,48,84,219,111,210,138,
!             ... d5 coefficients
     :        -175,-85,23,-22,-112,-76,-58,167,
     :        23,-85,59,140,104,86,320,113/
      data    f4dd/5,-70,15,50,140,
!             ... d3 coefficients
     :        -30,-105,30,55,-45,105,-15,30,
    
!             ... d4 coefficients
     :        -105,15,-10,45,-30,-45,70,-55,
     :        75,135,20, 0,30,-15,210,-30,
!             ... d5 coefficients
     :        -175,-50,-40,-85,35,50,110,-15,
     :        -5,125,-25,140,20,-40,-100,-55/
     
      done = .true.
      n = q
      if (n .gt. 2*l+1) n = 4*l+2 - n
      if (n .gt. 1) then
         if (l .eq. 1) then
            call add(2*f2pp(i)/25.d0,2,iel,iel,.true.)
         else if (l .eq. 2) then
            i = i-6
            call add(2*f2dd(i)/441.d0,2,iel,iel,.true.)
            call add(2*f4dd(i)/441.d0,4,iel,iel,.true.)
	 else
	    done = .false.
         end if
      end if
     END SUBROUTINE dev

   END SUBROUTINE get_case
