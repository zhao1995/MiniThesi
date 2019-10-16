      subroutine arclen(ug,u,u1,u2,f,f0,xal,xau,xa,jp,id,nneq,neq,
     +           time,rlnew1)
c----------------------------------------------------------------------
c
c       Purpose: constant arc length method with original newton/raphson
c                normal plane assumption and displacement control

c       kflag   - flag for solution type
c                            mod n/r (eq.0)
c                           orig n/r (eq.2)
c                updated normal plane iteration
c                            mod n/r (eq.1)
c                           orig n/r (eq.3)
c                 displacement control
c                            mod n/r (eq.4)
c                           orig n/r (eq.5)
c                 new tangent plane
c                                    (eq.6)
c
c      Input:
c       ug(*)    - displacements of converged solution
c       u(*)     - actual displacements increments
c       u1(*)    - displacement before first iteration
c       u2(*)    - displacement for load level 1.0 (only for n/r)
c       f(*)     - load vector (only for n/r)
c       f0(*)    - load vector (only for n/r)
c       jp(*)    - Pointer array for row/columns of tangent
c                  only necessary for dasol (only for n/r)
c       rlnew    - actual load level
c       rlold    - load level of old iteration (only updated plane )
c       ds0      - constant arc length
c       alfa0    - reducing factor for n/r for prescribed displ.
c       ite      - number of iteration
c       rr       - energy at beginning of load increment
c       ndis     - equation no. for displacement control
c       u2(ndis) - displacement for single displ.control
c       c0       - current stiffness parameter value from start
c       cs1,cs2  - stiffness parameter values of previous two load steps
c       det      - current determinant of tangent stiffness
c
c      Output:
c
c----------------------------------------------------------------------
c
      USE arcl
      USE conv
      USE fdata
      USE iofile
      USE prlod
      USE stepc
      USE tdata
      implicit real*8(a-h,o-z)
c.... Declare variable types
c.... Declare arrays types
      integer jp(*),id(*)
      real*8  ug(*),u(*),u1(*),u2(*),f(*),f0(*),xal(*),xau(*),xa(*)
c.... Intrinsics
      intrinsic abs, log, max, min, sqrt
c.... External
      external ddot

      save  alfold
      data zero,one/0.0d0,1.0d0/
      ifirst=0
      if (ds0.lt.1.0d-20) ifirst=1
c
      if (time.ne.timold) then
c..... first iteration step:  G = 1.0*P
        ite    = 0
        timold = time
        rold   = rlnew1    !
        cs1o = cs1        ! save old values for restart
        cs2o = cs2        !
c.....  calculate arc length at beginning of time step
        ds = sqrt(ddot(neq,u,1,u,1))
c.......set arc length constant  (only in very first time step)
c.......sign of r positive
        if (time.le.dt .or. ifirst.eq.1) then
          r   = one
          ds0 = ds  ! initial value
          if (ifirst.eq.1) then
            ds0 = ds/time  ! to prevent large values if arcl starts later
          end if
        end if
c.......check for limit points
        call ploads(ug,u1,one,.false.,.false.,.false.) ! ww 08/04/2010 compressed!
        call pload(id,f,f0,u1,nneq,one)
        rr = ddot(neq,u,1,u1,1)
c.......safe energy value at beginning of loading
        if (time.le.dt .or. ifirst.eq.1) c0=rr
c.......current stiffness parameter of current load increment
        cs = c0/rr
        rr = sign(1.0d0,rr)
c.......check if stiffness parameter passes infinity
c.......or limit point (then change of sign!)
        if (time.gt.(dt+dt)) then
          if (cs2.gt.zero.and.cs.lt.zero) then
            dcs = cs2 - cs1
            if (dcs.lt.zero) r = rr
          else if (cs2.lt.zero.and.cs.gt.zero) then
            dcs = cs2 - cs1
            if (dcs.gt.zero) r = rr
          end if
        end if
c.......save current stiffness values
        if (time.eq.dt .or. ifirst.eq.1) cs2 = cs
        cs1 = cs2
        cs2 = cs
c.......for constant arc length method
        if (kflag.le.3.or.kflag.eq.6) then
c.......calculate reducing factor
          alfa0 = ds0/ds
c.........change direction if passing limit points
          alfa0 = alfa0*r
c.........actual load level
          rlnew1 = rlnew1 + alfa0 - 1.0d0/time
c.........save displacement vector for this load step
c.........for modified and original newton/raphson
          do i = 1,neq
            u(i)  = u(i)*alfa0
            u1(i) = u(i)
          end do
c
          write (iow,1001) time,ds,alfa0
          write (iow,1002) rlnew1*prop,r
          write (iow,1005) cs,cs1,cs2
          if(ior.lt.0.and.pfr) then
            write (*,1001) time,ds,alfa0
            write (*,1002) rlnew1*prop,r
            write (*,1005) cs,cs1,cs2
          end if
        end if
c.......displacement control (beginning of time step)
c
c.......calculate displacements for load level 1.0
        if (kflag.eq.4.or.kflag.eq.5) then
          call ploads(ug,u1,prop,.false.,.false.,.false.)
          call pload(id,f,f0,u1,nneq,prop)
          call dasol(xal,xau,xa,u1,jp,neq,energy)
c
c.......  load factor for displacement control
c
          alfa = (alfa0 - u(ndis))/u1(ndis)
c
c......   update displacements
c
          do i=1,neq
            u(i) = u(i) + alfa*u1(i)
          end do
c
c......   update load level
c
          rlnew1 = rlnew1 + alfa
          write (iow,1006) alfa*prop,rlnew1*prop
          write (iow,1005) cs,cs1,cs2
          if(ior.lt.0) then
            if(pfr) then
              write (*,1006) alfa*prop,rlnew1*prop
              write (*,1005) cs,cs1,cs2
            else
              write (*,1002) rlnew1*prop,r
            end if
          end if
        end if
c
      else
c
c       iteration starts
c
        ite = ite + 1
        if(kflag.gt.1.and.kflag.ne.4) then
c.......calculate newton displacement for load level 1.0*P
          call ploads(ug,u2,prop,.false.,.false.,.false.)
          call pload(id,f,f0,u2,nneq,prop)
          call dasol(xal,xau,xa,u2,jp,neq,energy)
        end if
c-----------------------------------------
c       Update 1. modified and original newton raphson
c-----------------------------------------
c.......calculate reducing factor
        if (kflag.eq.0) then
c         iteration on normal plane (unscaled displ. vector)
          alfa = - ddot(neq,u1,1,u,1)*alfa0/(ds0*ds0)
        else if (kflag.eq.1) then
c         iteration on updated normal plane
          if (ite.eq.1) then
            do i=1,neq
              u2(i) = u1(i)
            end do
          end if
          alfan =  ddot(neq,u,1,u2,1)
          alfad =  ddot(neq,u1,1,u2,1)
          alfa =  alfan/alfad
        else if (kflag.eq.4) then
c         single displacement control
          alfa = -u(ndis)/u1(ndis)
        else if (kflag.eq.5) then
c         single displacement control
          alfa = -u(ndis)/u2(ndis)
        else
c         arc length (unscaled displacement vector)
          alfan = ddot(neq,u1,1,u,1)
          alfad = ddot(neq,u1,1,u2,1)
          alfa = -alfan/alfad
         end if
c.......numerical damping if sign of alfa is varying
        damp = 1.d0
        if (ndamp.eq.0) then
          if (ite.eq.1) alfold=alfa
          if(alfold.eq.0) goto 100 ! added ww
          rat = alfa/alfold
          if (-1.0d0.lt.rat.and.rat.lt.0.0d0) damp = 0.5d0
100       alfa = alfa*damp
          alfold = alfa
        end if
c.......update new load level
        rlnew1 = rlnew1 + alfa
        write (iow,1003) ite,alfa*prop,rlnew1*prop
        write (iow,1004) damp
        if(ior.lt.0) then
          if(pfr) then
            write (*,1003) ite,alfa*prop,rlnew1*prop
            write (*,1004) damp
          else
            write (*,1002) rlnew1*prop,r
          end if
        end if
c.......update displacements in iteration step
        if(kflag.le.1.or.kflag.eq.4) then
          do i = 1,neq
            u(i) = u(i)*damp + u1(i)*alfa
          end do
c.......  for updated normal plane iteration
          if (kflag.eq.1) then
            do i=1,neq
              u2(i) = u1(i) + u(i)
            end do
          end if
        else
c.......  update displacements in iteration step
          do i=1,neq
            u(i) = u2(i)*alfa + u(i)*damp
          end do
c.......  only for updated normal plane iteration
          if (kflag.eq.3) then
            do i=1,neq
              u1(i) = u(i) + u1(i)
            end do
c.......  new tangent plane
          else if (kflag.eq.6) then
            do i=1,neq
              u1(i) =  u(i)
            end do
          end if
        end if
c.... calculate determinant
                     ifdeta = 0
      if(time.le.dt) ifdeta = 1
      call detkt(xa,neq,ifdeta)
                             write(iow,1008) detc,nneg
        if(ior.lt.0.and.pfr) write(*  ,1008) detc,nneg
      end if
      iti = ite
c.... formats
 1001 format(/3x,'const. arc length method - necessary values'//
     1   3x,'time = ',g7.2,' arc length = ',g12.5,' red. factor = ',
     2   g12.5)
 1002 format('   load    level = ',g12.5,' direc ',g9.2)
 1003 format(/3x,'factors during iterations'//3x,'iterat.nr = ',i5,
     1  ' reduct.factor = ',g12.5,' load level = ',g12.5)
 1004 format(3x,'damp. fact. = ',g12.5)
 1005 format(3x,'cs-param.   = ',g12.5,/,
     *       3x,'cs1         = ',g12.5,/,
     *       3x,'cs2         = ',g12.5)
 1006 format(3x,'displacement control parameters',/,
     1       3x,'reduction factor = ',g12.5,' load level = ',g12.5)
 1008 format(3x,'determinant = ',g12.5,3x,'neg. diagonals =',i3)
c
      end
c
      subroutine ploa1(time,rlnew1,prop,propq)
c----------------------------------------------------------------------
c
c      Purpose: calculate load level for constant arc length method
c               rlnew_new*p*t_new = rlnew_old*p*t_old + p*1
c
c      Input:
c      time      - actual time
c      prop      - load factor from plod
c
c      Output:
c      rlnew     - actual load multiplier
c      propq     - actual load factor
c
c----------------------------------------------------------------------
      USE arcl
      USE ext2
      implicit real*8(a-h,o-z)
c..... Declare variable types
      if(.not. arcf .and. .not. extflg) then
        rlnew1 = 1.
      end if
c.... for restart only
      if (refl) then
cww        rlnew1 =  rlold
        refl  = .false.
      end if
c.... check if new time step and set new load level
      if(time .ne. timold) then
        if (arcf) rlnew1 = rlnew1*(time-1)/time + 1.0d0/time
      end if
c.... actual load factor
      propq = prop*rlnew1
      end
c
      subroutine paddi(vk,ve,nneq,xsi,id)
c----------------------------------------------------------------------
c
c      Purpose: Add imperfection(eigenvector) to displacement vector
c               vk(i) = vk(i) + xsi * ve(j)
c
c      Input:
c       vk(nneq)   - displacement vector
c       ve(neq)    - eigenvector
c       nneq       - No of unknowns in problem
c       xsi        - factor
c       id(*)      - Equation numbers for each active dof
c
c      Output:
c       vk(nneq)   - displacement vector
c
c----------------------------------------------------------------------
      USE iofile
      implicit real*8(a-h,o-z)
      integer id(*)
      real*8  vk(*),ve(*)
      do 100 i = 1,nneq
        j = id(i)
        if (j.gt.0) vk(i) = vk(i) + xsi * ve(j)
100   continue
      end
c
      subroutine paddv(vk,ve,nneq,tau,id)
c----------------------------------------------------------------------
c
c      Purpose: Add imperfection(eigenvector) to displacement vector
c               vk(i) = vk(i) + xsi * ve(j)
c
c      Input:
c       vk(nneq)   - displacement vector
c       ve(neq)    - eigenvector
c       nneq       - No of unknowns in problem
c       tau        - factor -> xsi   = vknorm / venorm / tau
c       id(*)      - Equation numbers for each active dof
c
c      Output:
c       vk(nneq)   - displacement vector
c
c----------------------------------------------------------------------
c
      USE iofile
      implicit real*8(a-h,o-z)
c..... Declare variable types
c..... Declare array types
      integer id(*)
      real*8  vk(*),ve(*)
c..... Intrinsics
      intrinsic sqrt
c..... External
      external ddot
c
      vknorm = sqrt(ddot(nneq,vk,1,vk,1))
      venorm = 0.0d0
      do 100 i = 1,nneq
        j = id(i)
        if (j.gt.0) venorm = venorm + ve(j) * ve(j)
100   continue
      venorm = sqrt (venorm)
      xsi    = vknorm / venorm / tau
      do 110 i = 1,nneq
        j = id(i)
        if (j.gt.0) vk(i) = vk(i) + xsi * ve(j)
110   continue
      write(iow,2000) vknorm,venorm,xsi
      if(ior.lt.0) write(*,2000) vknorm,venorm,xsi
c
2000  format(/,3x,'Norm displ. vector  = ',g12.5,/,
     1         3x,'Norm eigenvector    = ',g12.5,/,
     2         3x,'Scaling factor      = ',g12.5,/)
      end
c
      subroutine dicont(id,numnp,ndf,lflag)
c----------------------------------------------------------------------
c
c      Purpose: changes arc-length for lflag .ne. 0
c               provides information for displacement control
c               reads in :
c               node, dof, nr's size of assigned displ for displ. control
c               determines:
c               eq. nr of assigned displacement
c               factors for scaled arc length control
c
c      Input:
c       id(ndf,*)    - Equation numbers for each active dof
c       numnp        - Number of nodes in mesh
c       ndf          - Number dof/node
c       lflag        -
c
c      Output:
c
c----------------------------------------------------------------------
      USE arcl
      USE iofile
      implicit real*8(a-h,o-z)
c
c.... Declare variable types
      character*1 ch
c.... Declare array types
      integer id(ndf,*)
      real*8  td(3)

      if (lflag .ne. 0) go to 100
c.... read in if numerical damping desired or not
c.kne      if(ior.lt.0) write(*,3002)
c.kne      call dinput(td,1)
c.kne      ndamp = td(1)
c.kne      write (iow,2003) ndamp
      ndamp = 0                  ! damping set to zero
c
c..... restart flag
      if (refl) go to 100
50    if (kflag.eq.4.or.kflag.eq.5) then ! disp control
        if(ior.lt.0) write(*,3001)
        call dinput(td,3)
        nodis = td(1)
        nddis = td(2)
        alfa0 = td(3)
        if(ior.lt.0) then
          if(nodis.le.0 .or. nodis.gt.numnp) go to 50
          if(nddis.le.0 .or. nddis.gt.ndf  ) go to 50
        end if
c...    calculation of eq.nr for prescribed displacement
        ndis = id(nddis,nodis)
        write (iow,3000) nodis,nddis,ndis,alfa0
        if( ndis .le. 0 ) then
          if(ior.lt.0) then
            write(*,2001)
            go to 50
          else
            write(iow,2001)
            stop 'SR DICONT'
          end if
        end if
      end if
      return
c.... for restart only
 100  continue
c.... any method (displacement control stiff.param. just for chance)
      write(iow,2005) rlnew,c0,cs1,cs2

      if (kflag.lt.4.or.kflag.eq.6) then
c....   arc length method (any)
        if(ior.lt.0) then
          assign 110 to iend
 101      write(*,2006) ds0,r
          write(*,2007)
          read (*,1000,err=101,end=900) ch
        else
          read (ior,1000,end=900) ch
        end if
 110    if(ch.eq.'n' .or. ch.eq.'N') then
          if(ior.lt.0) write(*,3003)
          call dinput(td,2)
          ds0 = td(1)
          r   = td(2) ! +-1
          if(r.ge.0.d0) r= 1.d0
          if(r.lt.0.d0) r=-1.d0
          write(iow,2006) ds0,r
        end if
      else
c....   displacement control
        if(ior.lt.0) then
          assign 120 to iend
 102      write(*,2009) nodis,nddis,alfa0
          write(*,2007)
cww          write(*,2008)
          read (*,1000,err=102) ch
        else
          read (ior,1000,end=900) ch
        end if
 120    if(ch.eq.'n' .or. ch.eq.'N') go to 50
      end if
      return
c.... eof encountered
900   call  endclr ('DICONT',ch)
      go to  iend
c
1000  format(a1)
2001  format('   Displacement control specified on restrained node')
cww2003  format('   Numerical damping = ',i3,3x,'(0 = no damping)')
2005  format('   v a l u e s  for  r e s t a r t:',/,
     * '     Current load level      = ',g12.4,/,
     * '     S t i f f n e s s  parameter values ',/,
     * '     Stiff.param first step  = ',g12.4,/,
     * '     Stiff.param 1.prev.step = ',g12.4,/,
     * '     Stiff.param 2.prev.step = ',g12.4,/)
2006  format('   Given arc length         = ',g14.7/
     *       '   Load direction           = ',f12.4)
2007  format('   Keep values (y or n): ',$)
cww2007  format('   Keep arc-length and load-direction (y or n): ',$)
c2008  format('   Keep displacement control parameters (y or n): ',$)
2009  format('   Node number      = ',i7,/,
     1     '   Ndof number      = ',i3,/,
     2     '   Prescribed disp. = ',e12.5,/)
3000  format('   s i n g l e   d i s p l a c e m e n t   control ',/,
     1      '   node nr.  ndof.nr.  equat.nr.  prescribed displ. ',/,
     2      4x,i4,7x,i2,8x,i4,6x,g12.5)
3001  format(' Input: node, dof, value-> ',$)
cww3002  format('   Input: numerical damping (0 = no damping)->',$)
3003  format(' Input: new arc-length+load direction->',$)
      end
c
      subroutine stepcntl (u,m,ct,isw,ic)
c-----------------------------------------------------------------------
c
c      Purpose: step lenght control for arcl-length method
c      (coded according to crisfield and diss. reitinger by K.Knebel jan. 96)
c
c      isw = 1 modify ds0
c      parameter:
c             itd = max. number iter. (desired) [default = 6]
c             iti = number of iter. of previous time step
c                   reduce   ds0  if  itd < iti
c                   rise     ds0  if  itd > iti
c             sp =  exponent   0.5 < sp < 1.0  [default=0.5]
c                   set fixed in ini
c
c      isw = 2 monotonie check    do restart and reduce ds0
c      parameter:
c             unorm  = norm of the iterativ displacements
c             unorm1 = unorm of previous iteration
c             cnorm  = ratio  unorm/unorm1
c             cmax   = max. cnorm (user defined) [default=1.0]
c                      0.5 < cmax < 5.0
c                      reduce ds0 if cnorm > cmax
c             alfak  = reduktion factor  (ds0 = alfak*ds0)
c             rm     = safety factor     [default=0.5]
c
c-----------------------------------------------------------------------
      USE arcl
      USE cdata
      USE endata
      USE iofile
      USE ldata
      USE stepc
      USE tdata
      implicit real*8(a-h,o-z)
      dimension u(*),ct(3,*)
      ic=0
      if(ttim.le.dt) return          ! not in first time step !!!
      go to (10,20) isw
10    continue                       ! modify arc-length
        if(ttim.eq.timold) return    ! only in predictor step
          alphas  = (1.d0*itd/(1.d0*iti))**sp
          ds0 = alphas*ds0
                     write(iow,1000) alphas,ds0
                     write(*,  1000) alphas,ds0
c        if(ior.lt.0) write(*,  1000) alphas,ds0
      return
c
20    continue          ! check for divergence
        unorm1 = unorm
        unorm  = sqrt(ddot(neq,u,1,u,1))
      if(ite.ge.2) then
        cnorm = unorm/unorm1
        if(cnorm .gt. cmax) then
          alfak = sqrt(rm*cmax/cnorm)     ! reduction factor
          ds0=alfak*ds0
c
          ct(1,lve(lv)) = ct(1,lvs(lv))     ! stop current loop
          timold = timold-dt
          rlnew  = rold - 1.d0 ! ww used?? geht das???? ist das rlnew aus /arcle/???
          ttim   = ttim - dt                   ! reset time
          aengy  = 1.d+10
          iti = itd
          ic  = 1
          cs1 = cs1o
          cs2 = cs2o
                       write(iow,1100) alfak,ds0,ttim
                       write(*,  1100) alfak,ds0,ttim
c         if(ior.lt.0) write(*,  1100) alfak,ds0,ttim
          end if
      end if
      return
1000  format( 'ds0 modified by factor =',g10.4,'new ds0=',g12.4)
1100  format( 'ds0 reset    by factor =',g10.4,'new ds0=',g12.4,/
     +        'restart time step ',f10.2)
      end
c

