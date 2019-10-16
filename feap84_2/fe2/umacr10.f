c$Id:$
      subroutine umacr10(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    27/03/2011
c       1. Correct set of is(*) for prtype = 1              15/01/2012
c       2. Change to umacr10                                25/01/2012
c       3. Increase rbuf to 24 and sbuf to 56; reset ns,    08/05/2012
c          nw & nsends; add finflg & prtype to osets10 call;
c          Correct set of temperature into sbuf; corrrect
c          set of is(i) for prtypes; Change name to 'fe2 '.
c       4. Add periodic temperature set to osetu10          09/05/2012
c          Change a_avg to v_avg; rve_rho to v_rho,
c          rve_c to v_c
c       5. Return averaged density for prtype = 2           10/05/2012
c       6. Remove 'tol' from osetdu10 routine               18/05/2012
c       7. Change TEMP1,2,3 to HILLI,G,X                    13/04/2013
c       8. Change fm(*) to gradu(*)                         18/04/2013
c          Divide shear stres by 2 for isw = 6, etc.
c       9. Export convergence flage                         06/05/2013
c      10. Move stop using -999 to preserve counters        08/05/2013
c      11. Replace 'omacr1.h' by 'elpers.h'                 21/05/2013
c      12. Use td(1) to retrieve vmxits                     15/08/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: MPI Interface for OpenMPI: Thermo-mechanical problem
c               classes for small and finite deformation problems

c      Problem Types:
c               1. Thermal
c               2. Stress: Small and finite deformation.
c               3. Thermal-mechanical (incomplete)

c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters

c      Outputs:
c         N.B.  Interprocessor communications
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'cdat2.h'
      include   'comnds.h'
      include   'counts.h'
      include   'debugs.h'
      include   'elpers.h'
      include   'endata.h'
      include   'fdata.h'
      include   'idptr.h'
      include   'iofile.h'
      include   'ldata.h'
      include   'oelmt.h'
      include   'print.h'
      include   'prlod.h'
      include   'rdata.h'
      include   'rdat1.h'
      include   'sdata.h'
      include   'setups.h'
      include   'tdata.h'
      include   'tdato.h'
      include   'umac1.h'

      include   'mpif.h'

      include   'p_int.h'
      include   'pointer.h'
      include   'comblk.h'
      include   'chdata.h'

      logical    pcomp,setval,palloc,flgu,flgh,hflg,ioflg
      logical    lflg,fflg,wflg,swapfl,startfl, exst, isopen
      character  lct*15
      real*8     ctl(3)

      character  HistRD*20, HistWR*20, SWitch*15, fext*2
      integer    i,j,k,nc,ni,nn,ierr, msg_stat(MPI_STATUS_SIZE)
      integer    nw, nwrv, nq, ns,nss, option, lenu,lenh,preu,preh
      integer    nntot,nitot,nsends, err_no
      integer    rdunit,wrunit, ostep,isw, psw
      real*4     etime, tary(2), tt
      real*8     volmr
      real*8     rbuf(24), sbuf(57), sig(6),dd(6,6), td(1)
      integer    usr_msg

      character  vvv*80
      logical    elnkfl,vmflg,vinput
      real*8     valoop,vmxits

      save

c     Set command word

      if(pcomp(uct,'ma10',4)) then       ! Usual    form
        uct   = 'fe2 '                   ! Specify 'name' as FE2
      elseif(urest.eq.1) then            ! Read  restart data

      elseif(urest.eq.2) then            ! Write restart data

      else                               ! Perform user operation

c       Start process and compute arrays for projections

        if(ntasks.le.1) then

          write(ilg,*) ' *ERROR* MPI solution for np = 2 only'
          write(iow,*) ' *ERROR* MPI solution for np = 2 only'
          call plstop()

        elseif(pcomp(lct,'star',4)) then ! Set array for stress projects

c         Set switch and maxium iterations for tangent computations

          call acheck(lzz(l),vvv,15,80,80)

          SWitch = vvv(1:15)
          vmflg  = vinput(vvv(16:30),15,td(1),1)
          vmxits = td(1)

          if(nint(vmxits).le.0) then
            vmxits = 15.0d0
          endif

          if(pfr) then
            write(iow,*) ' MPI: START - Parameter = ', SWitch(1:2)
     &                  ,' Max TANG Iters =',nint(vmxits)
          endif

c         Set number of stress/flux components

          if(debug) then
            write(  *,*) ' PRTYPE = ',prtype
            write(iow,*) ' PRTYPE = ',prtype
          endif

          if(prtype.eq.1) then         ! Thermal problem (scalar)
            ns     = ndm
            nw     = 13
            nsends = 17  ! Tangent + stress/flux + 2 + 1 (error)
          elseif(prtype.eq.2) then     ! Stress problem (small/finite)
            if(ndm.eq.1) then
              ns = 1
            elseif(ndm.eq.2) then
              ns = 4
            else
              ns = 6
            endif
            nw     = 21
            nsends = 45   ! Tangent + stress/flux + 2 + 1 (error)
          elseif(prtype.eq.3) then     ! Thermo-mechanical problem
            if(ndm.eq.1) then
              ns = 1
            elseif(ndm.eq.2) then
              ns = 4
            else
              ns = 6
            endif
            nq = ndm
            ns = ns + nq
            nw     = 24
            nsends = 44 + 12 + 1! Tangent + stress/flux + 2 + 1 (error)
          endif
          ostep = 0
          hflg  = .false.

c         Set size of tangent tensor

          nss    = ns               ! Size of G/H arrays
          nwrv   = nw

          setval = palloc(331,'HILLI',nen*ndf   , 1)    ! For ixl
          if(neq.gt.0) then
            setval = palloc(332,'HILLG',nss*neq*2 , 2)    ! g for stress
          endif
          setval = palloc(333,'HILLX',nen*ndm   , 2)    ! xs for coord
          option = nint(ctl(1))
          call psetvol(hr(np(43)),ndm,numnp) ! Volume of RVE
          volmr  = 1.d0/volm0

          HistRD = 'HistIO_1'
          HistWR = 'HistIO_2'
          fext   = '00'
          if(rank.lt.10) then
            write(fext(2:2),'(i1)') rank
          else
            write(fext(1:2),'(i2)') rank
          endif
          call addext(HistRD,fext,20,2)
          call addext(HistWR,fext,20,2)
          rdunit = 3          ! History read  unit (at start)
          wrunit = 4          ! History write unit (at start)

c         Check if file exists and destroy if true

          inquire(file = HistRD, exist = exst, opened = isopen)
          if(exst) then
            write(*,*) ' UMACR10: Delete exisiting file =',HistRD
            if(.not.isopen) then
              open (unit = rdunit, file = HistRD, status = 'old')
            endif
            close(unit = rdunit, status = 'delete')
          endif

          inquire(file = HistWR, exist = exst, opened = isopen)
          if(exst) then
            write(*,*) ' UMACR10: Delete exisiting file =',HistWR
            if(.not.isopen) then
              open (unit = wrunit, file = HistWR, status = 'old')
            endif
            close(unit = wrunit, status = 'delete')
          endif

c         Set flags for links and periodic b.c.

          elnkfl = .true.
          if(np(257).ne.0) then
            perflg = .true.
          endif

c         Open HistRD and HistWR files for history variable I/O

          nitot  = 0
          nntot  = 0
          psw    = 0
          ioflg  = .false.
          open(unit=rdunit, file=HistRD, form='unformatted',
     &         status='new')
          open(unit=wrunit, file=HistWR, form='unformatted',
     &         status='new')
          swapfl  = .false.
          startfl = .true.

c         Get length of U and H

          call pgetd('U    ',fp(1),lenu,preu,flgu)
          if(np(49).ne.0) then
            call pgetd('H    ',fp(2),lenh,preh,flgh)
          else
            lenh = 0
            flgh = .false.
          endif

c       Kill all processes

        elseif(pcomp(lct,'stop',4)) then ! Stop all processes

          sbuf(1) = -999  ! Stop indicator
          usr_msg =  12
          do i = 1,ntasks-1
            call MPI_SSend(sbuf, nw, MPI_DOUBLE_PRECISION, i, usr_msg,
     &                     MPI_COMM_WORLD, ierr)
          end do ! i
          write(*,*) ' --> MPI STOP EXECUTION'
          tt = etime(tary)
          write(iow,2000) tary
c         return
          call plstop()

c       Send Stress/Flux and Tangent Moduli to 'rank = 0' process

        elseif(pcomp(lct,'send',4)) then ! Send data to main processor

c         Output displacement and history record to history files

          if(rank.ge.1) then

            if(hflg) then
              if(startfl) then  ! Write to both files first time
                call usetio(ni,hr(np(40)),hr(np(49)),lenu,lenh,
     &                      flgh,rdunit,2)
                if(lflg) startfl = .false.
              endif
              if(wflg) then
                call usetio(ni,hr(np(40)),hr(np(49)),lenu,lenh,
     &                      flgh,wrunit,2)
                wflg = .false.
              endif
            endif

c           Compute Stress/Flux and Tangent Moduli accumulations

            if(isw.ne.12) then
              if(np(331).ne.0) then
                call osets10(mr(np(331)),mr(np(33)),mr(id31),nss, isw)
              else
                write(iow,*) ' *ERROR* Use MPI_START in slave process'
                call plstop()
              endif

c             Set send buffer

              sbuf(1)     = ni   ! Return 'ni' stress
              if(niter.gt.0. and. flncon) err_no = 1  ! No conv flag
              sbuf(nsends) = err_no

c             Store thermal flux and moduli

              if(prtype.eq.1) then

                k = 7
                do i = 1,ndm
                  sbuf(i+1) = pflux(i)*volmr
                  do j = 1,ndm
                    k = k + 1
                    sbuf(k) = pcflux(j,i)*volmr
                  end do ! j
                end do ! i

c               Return density and specific heat averages

                if(v_avg.gt.0.0d0) then
                  sbuf(5) = v_rho/v_avg
                  sbuf(6) = v_c/v_avg
                endif

c             Convert Kirchhoff stress to Cauchy stress and moduli

              elseif(prtype.eq.2 .or. prtype.eq.3) then

                k = 7
                if(finflg) then

                  call tau2sig(ptau,pctau,volmr,fdet, sig,dd, ndm,6)

c                 Place results and moduli into send buffer
c                 N.B. Stress and moduli multiplied by 'volmr'

                  do i = 1,6
                    sbuf(i+1) = sig(i)
                    do j = 1,6
                      k = k + 1
                      sbuf(k) = dd(j,i)
                    end do ! j
                  end do ! i

c               Small deformation Cauchy stress and moduli

                else

                  do i = 1,6
                    sbuf(i+1) = ptau(i)*volmr
                    do j = 1,6
                      k = k + 1
                      sbuf(k) = pctau(j,i)*volmr
                    end do ! j
                  end do ! i

                endif
                if(v_avg.gt.0.0d0) then
                  sbuf(44) = v_rho/v_avg
                else
                  sbuf(44) = 0.0d0
                endif

c               2-d case thickness stress

                if(ndm.eq.2 .and. v_avg.gt.0.0d0) then
                  sbuf(4) = sig_33/v_avg
                endif
              endif

            else
              sbuf(1) = rbuf(1) ! Return received value
            end if !isw

            usr_msg = 13
            call MPI_SSend( sbuf, nsends, MPI_DOUBLE_PRECISION, 0,
     &                      usr_msg, MPI_COMM_WORLD, ierr)
          endif

c       Get deformation gradient

        elseif(pcomp(lct,'get' ,3)) then ! Get data from Rank '0'
          nc = 0                         ! Set counter to zero
          if(rank.ge.1) then
            usr_msg = 12
            call MPI_Recv(rbuf, nw, MPI_DOUBLE_PRECISION, 0, usr_msg,
     &                    MPI_COMM_WORLD, msg_stat, ierr)

c           Set point parameter

            ni    = nint(rbuf(1))
c           rnmax = 1.0d0

c           Stop execution

            if(ni.eq.-999) then
              close(unit = rdunit,status='delete')
              close(unit = wrunit,status='delete')
              tt = etime(tary)
              write(iow,2000) tary
              call plstop()  ! Stop process
            endif

c           Set control parameters

            nstep = nint(rbuf(2))
            niter = nint(rbuf(3))
            dt    = rbuf(4)
            isw   = nint(rbuf(5))
            hflg  = rbuf(6).eq. 1.d0 ! Flag for history save: true = 1
            fflg  = rbuf(7).eq.-1.d0 ! Flag for first receive: true = 1
            lflg  = rbuf(8).eq. 1.d0 ! Flag for last receive: true = 1
            wflg  = ni.gt.0 ! .true.
            opar  = rbuf(9)          ! History variable
            call setparam(SWitch, rbuf(5), pfr)

            if(ni.eq.0) niter = 0

            if(isw.eq.3) then
              if(niter.le.0) then
                valoop = 1.d00
                j      = lv + 1    ! Loop level for tangent loop
                if(lvs(j).gt.0) then
                  ct(1,lvs(j)) = valoop
                  if(debug) then
                    write(iow,*)'SET LOOP1',j,lvs(j),ct(1,lvs(j))
                  endif
                endif
              else
                valoop = vmxits
                j      = lv + 1    ! Loop level for tangent loop
                if(lvs(j).gt.0) then
                  ct(1,lvs(j)) = valoop
                  if(debug) then
                    write(iow,*)'SET LOOP2',j,lvs(j),ct(1,lvs(j))
                  endif
                endif
              endif
            else
              valoop = 1.0d0
            endif

            if(debug) then
              write(iow,*) 'N =',ni,'TIME =',ttim,' NITER =',niter
              write(iow,*) 'ISW =',isw,' LOOP =',nint(valoop)
            endif

c           Initialize error number

            err_no = 0

c           Stress outputs (not needed)

            if(isw.eq.4) then

c           Receive macro model data

            else

c             Count number of total unit cells (depends on isw=6
c             being first call with a 'get'!)

              if(startfl) then
                nntot = nntot + 1
              else
                nitot = nntot
              endif

c             Swap and reopen files

              if(swapfl .or. ni.eq.0) then
                close(unit = rdunit, status='keep')
                close(unit = wrunit, status='keep')
                ioflg  = .true.
                open(unit=rdunit, file=HistRD, form='unformatted',
     &               status='old')
                open(unit=wrunit, file=HistWR, form='unformatted',
     &               status='old')
                swapfl  = .false.
              elseif(fflg) then  ! Start of new iteration
                rewind rdunit
                rewind wrunit
              endif

c             Time update

              if(isw.eq.12) then

c               Write solution status to log file

                if(nstep.gt.0) then
                  if(niter.gt.1) then
                    write(ilg,3001) nstep,niter,nform,ttim,dt,rel0,
     &                              rnorm,rnmax,aengy,etime(tary)
                  else
                    write(ilg,3002) nstep,niter,nform,ttim,dt,rel0,
     &                              rnorm,rnmax,etime(tary)
                  endif
                  call pflush(ilg)
                endif

                if(echo .and. ior.gt.0) then
                  write(*,3003) nstep,ttim
                endif

c               Increment time

                ttim = ttim + dt

c               Set iteration counters

                titer  = titer + niter
                niter  = 0
                taugm  = taugm + naugm
                naugm  = 0
                iaugm  = 0
                tform  = tform + nform
                nform  = 0
                iform  = 0
                intvc  = intvc + 1
                dtold  = dt

                do nn = 1,nitot

c                 Input displacement and history

                  call usetio(ni,hr(np(40)),hr(np(49)),lenu,lenh,
     &                        flgh,wrunit, 1)

c                 Zero displacement increment for time step

                  call pzero(hr(np(40)+nneq),nneq)

c                 Reset history variables and save to unit 'wrunit'

                  call reshis(mr(np(33)+nen),nen1,numel,2, 1)

                  call usetio(ni,hr(np(40)),hr(np(49)),lenu,lenh,flgh,
     &                        rdunit,2)
                end do ! nn

                if(lflg) then
                  swapfl = .true. ! Force close & reopen
                else
                  rewind rdunit
                  rewind wrunit
                endif
                wflg = .false.

                sbuf(1) =  ni

c             Solution step

              else
                if(nw.ge.nwrv) then

c                 Input displacement and history record

                  if(ioflg) then
                    call usetio(ni,hr(np(40)),hr(np(49)),lenu,lenh,
     &                          flgh,rdunit, 1)
                  else
                    call pzero(hr(np(40)),lenu)
                  endif

c                 Set gradient terms from rbuf

                  call usetrv0(rbuf(10), ns)

c                 Set edge values

                  call osetd10(mr(id31),hr(np(43)),hr(np(27)))
                  if(np(257).ne.0) then
                    if(elnkfl) then
                      elnkfl = .false.
                    endif
                    call osetu10(mr(np(257)),hr(np(43)),hr(np(40)))
                  endif

                else
                  write(iow,*) ' Insufficient data received in MPI_GET'
                  call plstop()
                endif
              endif

            endif
          endif
        else
          write(  *,*) ' *WARNING* MPI2 ',lct(1:4),' not implemented'
          write(iow,*) ' *WARNING* MPI2 ',lct(1:4),' not implemented'
        endif

      endif

c     Formats

2000  format(' *End of Solution Execution*',31x,'t=',2f9.2)

3001  format(2i6,i5,1p,1e11.3,1p,5e10.2,    0p,1f10.2)

3002  format(2i6,i5,1p,1e11.3,1p,4e10.2,10x,0p,1f10.2)

3003  format(' --> Step',i6,' Solution: Time =',1p,1e12.5)

      end

      subroutine ochkp10(ix,elnk,ndm,nel,sflg)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Check for active link in each element

c     Inputs:
c        ix(nel)     - Element nodes
c        elnk(ndm,*) - Link information at global nodes
c        ndm         - Mesh space dimension
c        nel         - Number nodes on element
c     Outputs:
c        sflg        - True if links exist
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      logical    sflg
      integer    ndm,nel, i,n
      integer    ix(nel),elnk(ndm,*)

      save

      do n = 1,nel
        if(ix(n).gt.0) then
          do i = 1,ndm
            if(elnk(i,ix(n)).gt.0) then
              sflg = .true.
              return
            endif
          end do ! i
        endif
      end do ! n

      end

      subroutine osetd10(id,x, f)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Set boundary displacement values using received values

c     Inputs:
c         id(ndf,numnp,2) - Equation/boundary codes
c         x(ndm,numnp)    - Reference coordinates

c     Outputs:
c         f(ndf,numnp,2)  - Force/Displacements at nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'                 ! numnp,numel
      include   'elpers.h'                ! xc(3)
      include   'idptr.h'                 ! id31
      include   'sdata.h'                 ! ndm,ndf

      integer    id(ndf,numnp,2)
      real*8     x(ndm,numnp), f(ndf,numnp,2)

      integer    i,j,n

      save

c     Thermal model

      if(prtype.eq.1) then
        do n = 1,numnp
          if(id(1,n,2).ne.0) then        ! Restrained dof
            f(1,n,2) = ttemp
            do i = 1,ndm
              f(1,n,2) = f(1,n,2) + gradt(i)*x(i,n)
            end do ! i
          else
            f(1,n,2) = 0.0d0
          end if
        end do ! n

c     Mechanical model

      elseif(prtype.eq.2 .or. prtype.eq.3) then
        do n = 1,numnp
          do j = 1,ndm
            f(j,n,2) = 0.0d0
            if(id(j,n,2).ne.0) then        ! Restrained dof
              do i = 1,ndm
                f(j,n,2) = f(j,n,2) + gradu(j,i)*x(i,n)
              end do ! i
            end if
          end do ! j
        end do ! n

      endif

c     Thermo-mechanical temperature set

      if(prtype.eq.3) then

c       Set temperatures

        if(id(ndm+1,n,2).ne.0) then        ! Restrained dof
          f(ndm+1,n,2) = ttemp
          do i = 1,ndm
            f(ndm+1,n,2) = f(ndm+1,n,2) + gradt(i)*(x(i,n) - xc(i))
          end do ! i
        else
          f(ndm+1,n,2) = 0.0d0
        endif

      endif

      end

      subroutine osets10(ixl,ix,id,nss,isw)
c                        331 33 31  <--- Pr no in call

c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Compute sums of boundary displacements for fine scale
c               model to obtain stress and its tangent moduli
c               (both are returned in matrix form)

c      Inputs:
c        finflg       - Finite deformation flag
c        prtype       - Problem type
c        ix(nen1,*)   - Element connection data
c        id(ndf,*,2)  - Boundary condition codes
c        nss          - Number tau stress components
c        isw          - FEAP element switch value

c      Outputs:
c        ixl(ndf,*)   - DOF indicators: -1 = no equation
c                                        0 = free dof
c                                        1 = boundary dof
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'iofile.h'
      include   'hdatam.h'
      include   'elpers.h'
      include   'sdata.h'
      include   'tdata.h'
      include   'pointer.h'
      include   'comblk.h'

      include   'counts.h'

      logical    sflg
      integer    i,j,nn, nel, nss, isw
      integer    ix(nen1,*), id(ndf,numnp,2), ixl(ndf,*)
      real*8     ht(144)

      save

c     Initialize averaged stress and tangent modulus arrays

      if(prtype.eq.1 .or. prtype.eq.3) then
        do i = 1,3
          pflux(i) = 0.0d0
        end do ! i
      endif
      if(prtype.eq.2 .or. prtype.eq.3) then
        do i = 1,6
          ptau(i) = 0.0d0
        end do ! i
      endif
      do i = 1,nss*nss
        ht(i) = 0.0d0                                ! Zero h array
      end do ! i
      if(neq.gt.0) call pzero(hr(np(332)),neq*nss*2) ! Zero g array

      do nn = 1,numel

c       Set 'ixl' array to mark dofs with fixed boundaries

        sflg = .false.
        nel  = 0
        do i = 1,nen
          if(ix(i,nn).gt.0) then
            nel    = i
            do j = 1,ndf
              if(id(j,ix(i,nn),2).ne.0) then  ! Look at boundary code
                ixl(j,i) = 1
                sflg   = .true.
              else
                ixl(j,i) = 0
              endif
            end do ! j
          else !  No node
            do j = 1,ndf
              ixl(j,i) = -1
            end do ! j
          endif
        end do ! i

c       Check for periodic case

        if(.not.sflg .and. np(257).ne.0) then
          call ochkp10(ix(1,nn),mr(np(257)),ndm,nel,sflg)
        endif

c       Get element tangent and residual: No assembly

        if(sflg) then

          hflgu  = .false.
          h3flgu = .false.

c                        U      B
          call formfe(np(40),np(26),np(26),np(26),
     &                .false.,.false.,.false.,.false.,3,nn,nn,1)

c         Project to Kirchhoff stress and tangent modulus

          if(isw.eq.3) then

c           Periodic boundary case

            if(np(257).ne.0) then

              call uprojpp10(ixl,id,ix(1,nn),mr(np(257)),
     &                       hr(np(43)),hr(np(44)),hr(np(333)),
     &                       hr(np(35)),hr(np(36)),hr(np(332)), ht,
     &                       ndm,ndf,nel,nst, neq, nss)

c           Displacement boundary case

            else

c             Arrays in pointer: XL=44, P/R=35, S=36
              call uprojpd10(ixl,id,ix(1,nn),hr(np(44)),
     &                       hr(np(333)),hr(np(35)),hr(np(36)),
     &                       hr(np(332)),ht,ndm,ndf,nel,nst, neq, nss)

            endif

c         Compute stress only

          else

            if(np(257).ne.0) then
              call uprojp10(ixl,ix(1,nn),mr(np(257)),
     &                      hr(np(43)),hr(np(44)),hr(np(333)),
     &                      hr(np(35)),ndm,ndf,nel)
            else
              call uprojd10(ixl,hr(np(44)),hr(np(333)),
     &                      hr(np(35)),ndm,ndf,nel)
            endif

          endif
        endif ! sflg
      end do ! nn

c     Use pmove to copy G(*,*,1) into G(*,*,2)

      if(isw.eq.3) then

c       Modify shear terms

        if(prtype.eq.2 .or. prtype.eq.3) then
          call utang10(hr(np(332)),ht, neq,nss)
        endif

        if(neq.gt.0) then
          call pmove(hr(np(332)),hr(np(332)+neq*nss),neq*nss)

c         Form material moduli by static condensation
c         N.B. Moduli returned in h array.

c                     h     g        g_cols
          call formhh(ht,hr(np(332)),nss,neq)

        endif

c       Store into cflux tangent

        if(prtype.eq.1 .or. prtype.eq.3) then

          nn = 0
          do i = 1,nss
            do j = 1,nss
              nn = nn + 1
              pcflux(j,i) = ht(nn)
            end do ! j
          end do ! i

        endif

c       Store into pctau tangent

        if(prtype.eq.2 .or. prtype.eq.3) then

          nn = 0
          do i = 1,nss
            do j = 1,nss
              nn = nn + 1
              pctau(j,i) = ht(nn)
            end do ! j
          end do ! i

        endif

c     Other cases

      else

c       Modify shear values

        if(prtype.eq.2 .or. prtype.eq.3) then
          do i = 4,nss
            ptau(i) = ptau(i)*0.5d0
          end do ! i
        endif

      endif

      end

      subroutine osetu10(elnk,x,u)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Set boundary displacement values from input values

c     Inputs:
c         prtype          - Problem type
c         elnk(ndm,numnp) - Linked edge indicators
c         x(ndm,numnp)    - Nodal coordinates

c     Outputs:
c         u(ndf,numnp,*)  - Displacements at nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'                 ! numnp,numel
      include   'cdat2.h'                 ! dxlnk(3)
      include   'elpers.h'                ! gradu(3,3)
      include   'sdata.h'                 ! ndm,ndf

      integer    elnk(ndm,numnp)
      real*8     x(ndm,numnp), u(ndf,numnp,*)
      real*8     ulnk(3)

      integer    i,j,n, e, it

      save

c     Compute linking incremental temperatures

      if(prtype.eq.1 .or. prtype.eq.3) then

        if(prtype.eq.1) then
          it = 1
        else
          it = ndm + 1
        endif

c       Loop over linked nodes to set temperature

        do n = 1,numnp
          ulnk(it)  = 0.0d0
          u(it,n,4) = 0.0d0
          e = elnk(it,n)
          if(e.gt.0) then        ! Node n is linked to node elnk(j,n)
            do i = 1,ndm
              ulnk(it) = ulnk(it) + gradt(i)*(x(i,n) - x(i,e))
            end do ! k
            u(it,n,4) = (u(it,e,1) - u(it,n,1)) + ulnk(it)
          endif
        end do ! n

      endif

c     Compute linking incremental displacements

      if(prtype.eq.2 .or. prtype.eq.3) then

c       Loop over linked nodes to set displacements

        do n = 1,numnp
          do j = 1,ndm
            ulnk(j)  = 0.0d0
            u(j,n,4) = 0.0d0
          end do ! j
          do j = 1,ndm
            e = elnk(j,n)
            if(e.gt.0) then        ! Node n is linked to node elnk(j,n)
              do i = 1,ndm
                ulnk(j) = ulnk(j) + gradu(j,i)*(x(i,n) - x(i,e))
              end do ! k
              u(j,n,4) = (u(j,e,1) - u(j,n,1)) + ulnk(j)
            endif
          end do ! j

        end do ! n

      endif ! prtype

      end

      subroutine uprojd10(ixl,xl,xs,p, ndm,ndf,nel)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set averaged tau stress

c      Inputs:
c        ixl(ndf,*)   - DOF indicators: -1 = no equation
c                                        0 = free dof
c                                        1 = boundary dof
c        xl(ndm,nel)  - Element nodal coordinates
c        p(ndf,nel)   - Element residual
c        ndm          - Spatial dimension of mesh
c        nel          - Number of maximum node on element

c      Working:
c        xs(ndm,nel)  - Element nodal coordinates

c      Outputs:
c        ptau(6)      - Kirchhoff stress (through common)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'elpers.h'

      integer    ndm,ndf,nel
      integer    i,ib,ir,a, it

      integer    ixl(ndf,*),isIb(3,3)
      real*8     xl(ndm,nel),xs(ndm,nel),p(ndf,nel)

      save

      data       isIb / 1, 4, 6,
     &                  4, 2, 5,
     &                  6, 5, 3/

c     Thermal problem

      if(prtype.eq.1 .or. prtype.eq.3) then
        if(prtype.eq.1) then
          it = 1
        else
          it = ndm + 1
        endif

        do ir = 1,nel
          if(ixl(it,ir).eq.1) then  ! Assemble from 'P'
            do ib = 1,ndm
              pflux(ib) = pflux(ib) - p(it,ir)*xl(ib,ir)  ! Thermal flux
            end do ! ib
          endif
        end do ! ir

      endif

c     Mechanical problem

      if(prtype.eq.2 .or. prtype.eq.3) then

c       Form current coordinates

        do ir = 1,nel
          do i = 1,ndm
            xs(i,ir) = xl(i,ir)
            if(finflg) then
              do a = 1,ndm
                xs(i,ir) = xs(i,ir) + gradu(i,a)*xl(a,ir)
              end do ! a
            endif
          end do ! i
        end do ! ir

c       Form tau

        do ir = 1,nel
          do i = 1,ndm
            if(ixl(i,ir).eq.1) then
              do ib = 1,ndm
                a      = isIb(i,ib)
                ptau(a) = ptau(a) - p(i,ir)*xs(ib,ir)  ! Stress
              end do ! ib
            endif
          end do ! i
        end do ! ir

      endif ! prtype

      end

      subroutine uprojp10(ixl,ix,elnk,x,xl,xs,p,
     &                  ndm,ndf,nel)

c-----[--.----+----.----+----.-----------------------------------------]

c      Inputs:
c        ixl(ndf,*)   - DOF indicators: -1 = no equation
c                                        0 = free dof
c                                        1 = boundary dof
c        ix(nel)      - Element connection list
c        elnk(ndm,*)  - Link indicators
c        x(ndm,*)     - Global nodal coordinate list
c        xl(ndm,nel)  - Element nodal coordinates
c        p(ndf,nel)   - Element residual
c        ndm          - Spatial dimension of mesh
c        nel          - Number of maximum node on element

c      Working:
c        xs(ndm,nel)  - Element nodal coordinates

c      Outputs:
c        ptau(6)      - Stress (through common)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'elpers.h'

      integer    ndm,ndf,nel
      integer    i,ib,ir,a, it

      integer    ixl(ndf,*),ix(nel),elnk(ndm,*), isIb(3,3), lelnk(3,8)
      real*8     x(ndm,*), xl(ndm,nel),xs(ndm,nel),p(ndf,nel)

      save

      data       isIb / 1, 4, 6,
     &                  4, 2, 5,
     &                  6, 5, 3/

c     Form linked coordinates

      do ir = 1,nel
        if(ix(ir).gt.0) then
          do i = 1,ndm
            lelnk(i,ir) = elnk(i,ix(ir))
          end do ! i
          do i = 1,ndm
            if(lelnk(i,ir).gt.0) then
              a = lelnk(i,ir)
              xl(i,ir) = xl(i,ir) - x(i,a)
            endif
          end do ! i
        else
          do i = 1,ndm
            lelnk(i,ir) = 0
          end do ! i
        endif
      end do ! ir

c     Thermal problem

      if(prtype.eq.1 .or. prtype.eq.3) then

        if(prtype.eq.1) then
          it = 1
        else
          it = ndm + 1
        endif

        do ir = 1,nel
          if(ixl(it,ir).eq.1 .or. lelnk(it,ir).gt.0) then  ! Assemble P1
            do ib = 1,ndm
              pflux(ib) = pflux(ib) - p(it,ir)*xl(ib,ir) ! flux
            end do ! ib
          endif
        end do ! ir

      endif

c     Mechanical problem

      if(prtype.eq.2 .or. prtype.eq.3) then

c       Form current coordinates

        do ir = 1,nel
          do i = 1,ndm
            xs(i,ir) = xl(i,ir)
            if(finflg) then
              do a = 1,ndm
                xs(i,ir) = xs(i,ir) + gradu(i,a)*xl(a,ir)
              end do ! a
            endif
          end do ! i
        end do ! ir

        do ir = 1,nel
          do i = 1,ndm
            if(ixl(i,ir).eq.1 .or. lelnk(i,ir).gt.0) then
              do ib = 1,ndm
                a      = isIb(i,ib)
                ptau(a) = ptau(a) - p(i,ir)*xs(ib,ir)  ! Stress
              end do ! ib
            endif
          end do ! i
        end do ! ir

      endif ! prtype

      end

      subroutine uprojpd10(ixl,id,ix,xl,xs,p,s, g, ht,
     &                     ndm,ndf,nel,nst, neq, nss)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set Coupling and diagonal arrays to compute averaged
c               TAU stress and tangent modulus arrays

c      Inputs:
c        finflg       - Finite deformation flag
c        prtype       - Problem type
c        ixl(ndf,*)   - DOF indicators: -1 = no equation
c                                        0 = free dof
c                                        1 = boundary dof
c        id(ndf,*)    - Equation numbers at nodes
c        ix(*)        - Element connection list
c        xl(ndm,nel)  - Element nodal coordinates
c        p(ndf,nel)   - Element residual
c        s(nst,nst)   - Element tangent matrix
c        ndm          - Spatial dimension of mesh
c        nel          - Number of maximum node on element
c        nst          - Dimension of tangent matrix
c        neq          - Number of active equations in mesh
c        nss          - Number of modes to project (generally = 9)

c      Working:
c        xs(ndm,nel)  - Element nodal coordinates

c      Outputs:
c        g(neq,nss,1) - Coupling matrix
c        ht(nss,nss)  - Block matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'elpers.h'

      integer    ndm,ndf,nel,nst, neq, nss
      integer    i,j,ib,jb,ir,jc,ia,a,b

      integer    ixl(ndf,*),id(ndf,*),ix(*), is(64),isIb(3,3)
      real*8     xl(ndm,nel),xs(ndm,nel),p(ndf,nel),s(nst,nst)
      real*8     g(neq,nss),ht(nss,nss)

      save

      data       isIb / 1, 4, 6,
     &                  4, 2, 5,
     &                  6, 5, 3/

c     Set assembly pointers

      if(prtype.eq.1) then
        is(1) = 1
      else
        is(1) = 0
      endif
      do i = 2,nel
        is(i) = is(i-1) + ndf
      end do ! i

c     Thermal problem

      if(prtype.eq.1) then

c       Form g_ib and h_ab arrays

        do ir = 1,nel
          if(ixl(1,ir).eq.0) then      ! Assemble G
            ia = id(1,ix(ir))
            if(ia.gt.0) then ! Equation number exists
              do jc = 1,nel
                if(ixl(1,jc).eq.1) then
                  do jb = 1,ndm
                    g(ia,jb) = g(ia,jb)
     &                       + s(is(ir),is(jc))*xl(jb,jc)
                  end do ! jb
                endif
              end do ! jc
            endif

          elseif(ixl(1,ir).eq.1) then  ! Assemble P1 and H1

            do ib = 1,ndm
              pflux(ib) = pflux(ib) - p(1,ir)*xl(ib,ir)  ! Thermal flux
              do jc = 1,nel
                if(ixl(1,jc).eq.1) then
                  do jb = 1,ndm
                    ht(ib,jb) = ht(ib,jb) + xl(ib,ir)
     &                        * s(is(ir),is(jc))*xl(jb,jc)
                  end do ! jb
                endif
              end do ! jc
            end do ! ib

          endif

        end do ! ir

c     Mechanical problem

      elseif(prtype.eq.2 .or. prtype.eq.3) then

c       Form current coordinates

        do ir = 1,nel
          do i = 1,ndm
            xs(i,ir) = xl(i,ir)
            if(finflg) then
              do j = 1,ndm
                xs(i,ir) = xs(i,ir) + gradu(i,j)*xl(j,ir)
              end do ! j
            endif
          end do ! i
        end do ! ir

c       Form g_ib and h_ab arrays

        do ir = 1,nel
          do i = 1,ndm
            if(ixl(i,ir).eq.0) then      ! Assemble G
              ia = id(i,ix(ir))
              if(ia.gt.0) then ! Equation number exists
                do jc = 1,nel
                  do j = 1,ndm
                    if(ixl(j,jc).eq.1) then
                      do jb = 1,ndm
                        b       = isIb(j,jb)
                        g(ia,b) = g(ia,b)
     &                          + s(is(ir)+i,is(jc)+j)*xs(jb,jc)
                      end do ! jb
                    endif
                  end do ! j
                end do ! jc
              endif

            elseif(ixl(i,ir).eq.1) then  ! Assemble P1 and H1

              do ib = 1,ndm
                a      = isIb(i,ib)
                ptau(a) = ptau(a) - p(i,ir)*xs(ib,ir)  ! Stress
                do jc = 1,nel
                  do j = 1,ndm
                    if(ixl(j,jc).eq.1) then
                      do jb = 1,ndm
                        b       = isIb(j,jb)
                        ht(a,b) = ht(a,b) + xs(ib,ir)
     &                          * s(is(ir)+i,is(jc)+j)*xs(jb,jc)
                      end do ! jb
                    endif
                  end do ! j
                end do ! jc
              end do ! ib

            endif
          end do ! i
        end do ! ir

      endif ! prtype

      end

      subroutine uprojpp10(ixl,id,ix,elnk,x,xl,xs,p,
     &                     s,g,ht,ndm,ndf,nel,nst,neq,nss)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set Coupling and diagonal arrays to compute averaged
c               TAU stress and tangent modulus arrays: Periodic case

c      Inputs:
c        finflg       - Finite deformaiton flag
c        prtype       - Problem type
c        ixl(ndf,*)   - DOF indicators: -1 = no equation
c                                        0 = free dof
c                                        1 = boundary dof
c        id(ndf,*)    - Equation numbers at nodes
c        ix(*)        - Element connection list
c        elnk(ndm,*)  - Linked information
c        x(ndm,*)     - Nodal coordinates
c        xl(ndm,nel)  - Element nodal coordinates
c        xs(ndm,nel)  - Element nodal coordinates
c        p(ndf,nel)   - Element residual
c        s(nst,nst)   - Element tangent matrix
c        ndm          - Spatial dimension of mesh
c        nel          - Number of maximum node on element
c        nst          - Dimension of tangent matrix
c        neq          - Number of active equations in mesh
c        nss          - Number of modes to project (generally = 9)

c      Working:
c        xs(ndm,nel)  - Element nodal coordinates

c      Outputs:
c        g(neq,nss,1) - Coupling matrix
c        h(nss,nss)   - Block matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'elpers.h'

      integer    ndm,ndf,nel,nst, neq, nss
      integer    i,j,ib,jb,ir,jc,ia,a,b

      integer    ixl(ndf,*),id(ndf,*),ix(*),elnk(ndm,*),is(64),isIb(3,3)
      real*8     x(ndm,*),xl(ndm,nel),xs(ndm,nel),p(ndf,nel),s(nst,nst)
      real*8     g(neq,nss),ht(nss,nss)

      integer     lelnk(3,8)

      save

      data       isIb / 1, 4, 6,
     &                  4, 2, 5,
     &                  6, 5, 3/

c     Set assembly pointers

      if(prtype.eq.1) then
        is(1) = 1
      else
        is(1) = 0
      endif
      do i = 2,nel
        is(i) = is(i-1) + ndf
      end do ! i

c     Modify nodal coordinates for links

      do ir = 1,nel
        if(ix(ir).gt.0) then
          do i = 1,ndm
            lelnk(i,ir) = elnk(i,ix(ir))
            if(lelnk(i,ir).gt.0) then
              xl(i,ir) = xl(i,ir) - x(i,lelnk(i,ir))
            endif
          end do ! i
        else
          do i = 1,ndm
            lelnk(i,ir) = 0
          end do ! i
        endif
      end do ! ir

c     Thermal problem

      if(prtype.eq.1) then

c       Compute thermal flux

        do ir = 1,nel
          if(ixl(1,ir).eq.1 .or. lelnk(1,ir).gt.0) then  ! Assemble P1
            do ib = 1,ndm
              pflux(ib) = pflux(ib) - p(1,ir)*xl(ib,ir)  ! Thermal flux
            end do ! ib
          endif
        end do ! ir

c       Form g_ib array

        do ir = 1,nel
          if(ixl(1,ir).eq.0) then      ! Assemble G
            ia = id(1,ix(ir))
            if(ia.gt.0) then ! Equation number exists
              do jc = 1,nel
                if(ixl(1,jc).eq.1 .or. lelnk(1,jc).gt.0) then
                  do jb = 1,ndm
                    g(ia,jb) = g(ia,jb)
     &                       + s(is(ir),is(jc))*xl(jb,jc)
                  end do ! jb
                endif
              end do ! jc
            endif
          endif
        end do ! ir

c       Form h_ab array

        do ir = 1,nel
          if(ixl(1,ir).eq.1 .or. lelnk(1,ir).gt.0) then  ! Assemble H1
            do ib = 1,ndm
              do jc = 1,nel
                if(ixl(1,jc).eq.1 .or. lelnk(1,jc).gt.0) then
                  do jb = 1,ndm
                    ht(ib,jb) = ht(ib,jb) + xl(ib,ir)
     &                        * s(is(ir),is(jc))*xl(jb,jc)
                  end do ! jb
                endif
              end do ! jc
            end do ! ib
          endif
        end do ! ir

c     Mechanical problem

      elseif(prtype.eq.2 .or. prtype.eq.3) then

c       Form current coordinates

        do ir = 1,nel
          do i = 1,ndm
            xs(i,ir) = xl(i,ir)
            if(finflg) then
              do j = 1,ndm
                xs(i,ir) = xs(i,ir) + gradu(i,j)*xl(j,ir)
              end do ! j
            endif
          end do ! i
        end do ! ir

c       Compute stress

        do ir = 1,nel
          do i = 1,ndm
            if(ixl(i,ir).eq.1 .or. lelnk(i,ir).gt.0) then
              do ib = 1,ndm
                a      = isIb(i,ib)
                ptau(a) = ptau(a) - p(i,ir)*xs(ib,ir)  ! Stress
              end do ! ib
            endif
          end do ! i
        end do ! ir

c       Form g_ib array

        do ir = 1,nel
          do i = 1,ndm
            if(ixl(i,ir).eq.0) then      ! Assemble G
              ia = id(i,ix(ir))
              if(ia.gt.0) then ! Equation number exists
                do jc = 1,nel
                  do j = 1,ndm
                    if(ixl(j,jc).eq.1 .or. lelnk(j,jc).gt.0) then
                      do jb = 1,ndm
                        b       = isIb(j,jb)
                        g(ia,b) = g(ia,b)
     &                          + s(is(ir)+i,is(jc)+j)*xs(jb,jc)
                      end do ! jb
                    endif
                  end do ! j
                end do ! jc
              endif
            endif
          end do ! i
        end do ! ir

c       Form h_ab array

        do ir = 1,nel
          do i = 1,ndm
            if(ixl(i,ir).eq.1 .or. lelnk(i,ir).gt.0) then  ! Assemble H1
              do ib = 1,ndm
                a      = isIb(i,ib)
                do jc = 1,nel
                  do j = 1,ndm
                    if(ixl(j,jc).eq.1 .or. lelnk(j,jc).gt.0) then
                      do jb = 1,ndm
                        b       = isIb(j,jb)
                        ht(a,b) = ht(a,b) + xs(ib,ir)
     &                          * s(is(ir)+i,is(jc)+j)*xs(jb,jc)
                      end do ! jb
                    endif
                  end do ! j
                end do ! jc
              end do ! ib
            endif
          end do ! i
        end do ! ir

      endif ! prtype

      end

      subroutine utang10(g,ht, neq,nss)

      implicit   none

      include   'elpers.h'

      integer    neq,nss, a,b,ia
      real*8     g(neq,nss), ht(nss,nss)

      save

c     Modify shear values

      do a = 4,nss
        ptau(a) = ptau(a)*0.5d0
        do ia = 1,neq
          g(ia,a) = g(ia,a)*0.5d0
        end do ! ia
        do b = 1,nss
          ht(a,b) = ht(a,b)*0.5d0
        end do ! b
        do b = 1,nss
          ht(b,a) = ht(b,a)*0.5d0
        end do ! b
      end do ! a

      end

      subroutine usetrv0(rbuf, ns)

c     Set temperature and/or displacement gradient.

      implicit   none

      include   'elpers.h'    ! gradu(3,3), gradt(3),

      integer    ns, i,j,k
      real*8     rbuf(*)

      save

c     Thermal problem

      if(prtype.eq.1) then

        do i = 1,3
          gradt(i) = rbuf(i)
        end do ! i
        ttemp = rbuf(4)  ! Temperature

c     Mechanical problem

      elseif(prtype.eq.2 .or. prtype.eq.3) then

        if(finflg) then
          k = 0
          do j = 1,3
            do i = 1,3
              k = k + 1
              gradu(i,j) = rbuf(k)
            end do ! i
          end do ! j
        else
          do i = 1,3
            gradu(i,i) = rbuf(i)
          end do ! i
          gradu(1,2) = 0.5d0*rbuf(4)
          gradu(2,1) = gradu(1,2)
          if(ns.gt.4) then
            gradu(2,3) = 0.5d0*rbuf(5)
            gradu(3,2) = gradu(2,3)
            gradu(3,1) = 0.5d0*rbuf(6)
            gradu(1,3) = gradu(3,1)
          endif
        endif
        fdet  = rbuf(10)  ! Det F
        ttemp = rbuf(11)  ! Temperature

      endif

c     Thermo-mechanical problem

      if(prtype.eq.3) then

        do i = 1,3
          gradt(i) = rbuf(i+11)
        end do ! i

      endif

      end
