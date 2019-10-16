c$Id:$
      subroutine phillmandel(lct,ct)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    27/03/2011
c       1. Correct set of is(*) for prtype = 1              15/01/2012
c       2. Compute fdet for finite deformation              03/10/2012
c          Set 'finflg' to true in modlfd
c       3. Remove prtype.eq.4, combine prtype 2 & 3         13/04/2013
c       4. Modify shear stress for other isw values         18/04/2013
c       5. Output sig(3) for 2-d case                       07/07/2013
c       6. 'fdet' moved to 'elpers.h'                       21/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Hill-Mandel: Thermo-mechanical problem class for small
c               and finite deformation problems

c      Inputs:
c         lct       - Command character parameters
c         ct(3)     - Command numerical parameters

c      Outputs:
c         N.B.  Interprocessor communications
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'debugs.h'
      include   'elpers.h'
      include   'idptr.h'
      include   'iofile.h'
      include   'oelmt.h'
      include   'setups.h'
      include   'print.h'
      include   'sdata.h'

      include   'pointer.h'
      include   'comblk.h'

      logical    pcomp,setval,palloc
      logical    strefl,tangfl
      character  lct*15
      real*8     ct(3)

      integer    ns,nss, option, i,j
      real*8     volmr
      real*8     sig(6),dd(6,6)

      logical    elnkfl

      save

c     Check for projection type

      if(pcomp(lct,'tang',4) .or. pcomp(lct,'    ',4)) then
        tangfl = .true.
        strefl = .true.
      elseif(pcomp(lct,'stre',4)) then
        tangfl = .false.
        strefl = .true.
      else
        tangfl = .false.
        strefl = .false.
      endif

c     Set class of problem

      if(prtype.gt.0) then

c       prtype = 1           ! Thermal
c       prtype = 2           ! Mechanical
c       prtype = 3           ! Thermo-Mechanical

c       Set number of stress/flux components

        if(prtype.eq.1) then       ! Thermal case
          ns = ndm
        elseif(prtype.eq.2) then   ! Stress on solid case
          if(ndm.eq.1) then
            ns = 1
          elseif(ndm.eq.2) then
            ns = 4
          else
            ns = 6
          endif
        elseif(prtype.eq.3) then  ! Thermo-mechanical case (not coded)
          if(ndm.eq.1) then
            ns = 1
          elseif(ndm.eq.2) then
            ns = 4
          else
            ns = 6
          endif
          ns = ns + ndm
        endif

c       Set size of tangent tensor

        nss    = ns               ! Size of G/H arrays

        setval = palloc(331,'HILLI',nen*ndf   , 1)    ! For ixl
        if(neq.gt.0) then
          setval = palloc(332,'HILLG',nss*neq*2 , 2)  ! g for stress
        endif
        setval = palloc(333,'HILLX',nen*ndm   , 2)    ! xs for coord
        option = nint(ct(1))
        call psetvol(hr(np(43)),ndm,numnp) ! Volume of RVE
        volmr  = 1.d0/volm0

c       Set flags for links and periodic b.c.

        elnkfl = .true.
        if(np(257).ne.0) then      ! Periodic case exists
          perflg = .true.
        endif

c       Compute Stress and Tangent Moduli accumulations

        if(tangfl .or. strefl) then
          call osets(mr(np(331)),mr(np(33)),mr(id31),hr(np(332)),
     &               hr(np(333)),nss,tangfl)

c         Convert Kirchhoff to Cauchy stress and moduli

          if(finflg) then
            call pdetf(gradu,fdet)
            call tau2sig(ptau,pctau,volmr,fdet, sig,dd, ndm,6)
          else
            do i = 1,nss
              sig(i) = ptau(i)*volmr
              do j = 1,nss
                dd(j,i) = pctau(j,i)*volmr
              end do ! j
            end do ! i
          endif

c         2-d thickness stress

          if(ndm.eq.2 .and. v_avg.gt.0.0d0) then
            sig(3) = sig_33/v_avg
          endif

        endif

c       Delete temp arrays

        setval = palloc(331,'HILLI',0 , 1)    ! For ixl
        if(neq.gt.0) then
          setval = palloc(332,'HILLG',0 , 2)  ! g for stress
        endif
        setval = palloc(333,'HILLX',0 , 2)    ! xs for coord

c       Output homogenized results

        if(strefl) then
          if(ior.lt.0) then
            if(prtype.eq.1) then
              write(*,2000) 'T h e r m a l    F l u x'
              write(*,2001) (sig(j),j=1,nss)
            else
              write(*,2000) 'C a u c h y    S t r e s s'
              write(*,2001) (sig(j),j=1,nss)
            endif
          endif
          if(prtype.eq.1) then
            write(iow,2000) 'T h e r m a l    F l u x'
            write(iow,2001) (sig(j),j=1,nss)
          else
            write(iow,2000) 'C a u c h y    S t r e s s'
            write(iow,2001) (sig(j),j=1,nss)
          endif
        endif

        if(tangfl) then
          if(ior.lt.0) then
            write(*,2000) 'T a n g e n t    M o d u l i'
            do i = 1,nss
              write(*,2001) (dd(i,j),j=1,nss)
            end do ! i
          endif
          write(iow,2000) 'T a n g e n t    M o d u l i'
          do i = 1,nss
            write(iow,2001) (dd(i,j),j=1,nss)
          end do ! i
        endif

      endif

c     Formats

2000  format(/5x,a)
2001  format(1p,6e12.4)

      end

      subroutine ochkp(ix,elnk,ndm,nel,sflg)

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

      subroutine osets(ixl,ix,id,g,   xs, nss,tangfl)
c                      331 33 31 332 `333 <--- Pointer no in call

c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Compute sums of boundary displacements for fine scale
c               model to obtain stress and its tangent moduli
c               (both are returned in matrix form)

c      Inputs:
c        ix(nen1,*)   - Element connection data
c        id(ndf,*,2)  - Boundary condition codes
c        nss          - Number tau stress components
c        tangfl       - Compute tangent if true

c      Outputs:
c        ixl(ndf,*)   - DOF indicators: -1 = no equation
c                                        0 = free dof
c                                        1 = boundary dof
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'elpers.h'
      include   'iofile.h'
      include   'sdata.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    sflg, tangfl
      integer    i,j,nn, nel, nss
      integer    ix(nen1,*), id(ndf,numnp,2), ixl(ndf,*)
      real*8     g(neq,nss,2),ht(144), xs(ndm,*)

      save

c     Initialize averaged stress and tangent modulus arrays

      do i = 1,6
        ptau(i) = 0.0d0
        do j = 1,6
          pctau(j,i) = 0.0d0
        end do ! j
      end do ! i
      do i = 1,nss*nss
        ht(i) = 0.0d0                      ! Zero h array
      end do ! i
      if(neq.gt.0) call pzero(g,neq*nss*2) ! Zero g array

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
          call ochkp(ix(1,nn),mr(np(257)),ndm,nel,sflg)
        endif

c       Get element tangent and residual: No assembly

        if(sflg) then
c                        U      B
          call formfe(np(40),np(26),np(26),np(26),
     &                .false.,.false.,.false.,.false.,3,nn,nn,1)

c         Project to Kirchhoff stress and tangent modulus
c         Arrays in pointer: X=43, XL=44, P/R=35, S=36

          if(tangfl) then

c           Periodic boundary case

            if(np(257).ne.0) then

              call uprojpp(ixl,id,ix(1,nn),mr(np(257)),
     &                     hr(np(43)),hr(np(44)),xs,
     &                     hr(np(35)),hr(np(36)),g, ht,
     &                     ndm,ndf,nel,nst, neq, nss)

c           Displacement boundary case

            else

              call uprojpd(ixl,id,ix(1,nn),hr(np(44)),xs,
     &                     hr(np(35)),hr(np(36)),g, ht,
     &                     ndm,ndf,nel,nst, neq, nss)

            endif

c         Compute stress only

          else

            if(np(257).ne.0) then
              call uprojp(ixl,ix(1,nn),mr(np(257)),
     &                    hr(np(43)),hr(np(44)),xs,
     &                    hr(np(35)),ndm,ndf,nel)
            else
              call uprojd(ixl,hr(np(44)),xs,
     &                    hr(np(35)),ndm,ndf,nel)
            endif

          endif
        endif ! sflg
      end do ! nn

c     Solve for tangent moduli

      if(tangfl) then

c       Modify shear terms

        if(prtype.eq.2 .or. prtype.eq.3) then

          call utangm(g,ht, neq,nss)   ! Include half factors on shears

c         Inform warning is o.k.

          if(ndm.lt.3) then
            if(ior.lt.0) then
              write(  *,2000)
            else
              write(iow,2000)
            endif
          endif
        endif

        if(neq.gt.0) then

c         Use pmove to copy G(*,*,1) into G(*,*,2)
          call pmove(g(1,1,1),g(1,1,2),neq*nss)

c         Form material moduli by static condensation
c         N.B. Moduli returned in h array.

c                     h  g g_cols
          call formhh(ht,g,nss,neq)

        endif

c       Store into tau tangent

        nn = 0
        do i = 1,nss
          do j = 1,nss
            nn = nn + 1
            pctau(j,i) = ht(nn)
          end do ! j
        end do ! i

c     Modify shear stress for other cases of isw

      else

        do i = 4,nss
          ptau(i) = ptau(i)*0.5d0
        end do ! i

      endif

c     Formats

2000  format(/'-->N.B. G-has no entries for thickness direction.',
     &        '  Results in following warning'/
     &     7x,' and no entries in third row/column of tangent and',
     &        ' stress.')
      end

      subroutine uprojd(ixl,xl,xs,p, ndm,ndf,nel)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set averaged tau stress

c      Inputs:
c        prtype       - Problem type
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
      integer    i,ib,ir,a

      integer    ixl(ndf,*),isIb(3,3)
      real*8     xl(ndm,nel),xs(ndm,nel),p(ndf,nel)

      save

      data       isIb / 1, 4, 6,
     &                  4, 2, 5,
     &                  6, 5, 3/

c     Thermal problem

      if(prtype.eq.1) then
        do ir = 1,nel
          if(ixl(1,ir).eq.1) then  ! Assemble from 'P'
            do ib = 1,ndm
              ptau(ib) = ptau(ib) - p(1,ir)*xl(ib,ir)  ! Thermal flux
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

      subroutine uprojp(ixl,ix,elnk,x,xl,xs,p,
     &                  ndm,ndf,nel)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set averaged TAU stress

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
      integer    i,ib,ir,a

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

      if(prtype.eq.1) then

        do ir = 1,nel
          if(ixl(1,ir).eq.1 .or. lelnk(1,ir).gt.0) then  ! Assemble P1
            do ib = 1,ndm
              ptau(ib) = ptau(ib) - p(1,ir)*xl(ib,ir) ! Kirchhoff stress
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

      subroutine uprojpd(ixl,id,ix,xl,xs,p,s, g, ht,
     &                    ndm,ndf,nel,nst, neq, nss)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set Coupling and diagonal arrays to compute averaged
c               TAU stress and tangent modulus arrays

c      Inputs:
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
              ptau(ib) = ptau(ib) - p(1,ir)*xl(ib,ir)  ! Thermal flux
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

            elseif(ixl(i,ir).eq.1) then  ! Assemble stress and H

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

      subroutine uprojpp(ixl,id,ix,elnk,x,xl,xs,p,
     &                   s,g,ht,ndm,ndf,nel,nst,neq,nss)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set Coupling and diagonal arrays to compute averaged
c               TAU stress and tangent modulus arrays: Periodic case

c      Inputs:
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

      integer    lelnk(3,8)

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
              ptau(ib) = ptau(ib) - p(1,ir)*xl(ib,ir)  ! Thermal flux
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

c     Mechanical small and finite deformation problem

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

      subroutine utangm(g,ht, neq,nss)

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

      subroutine pdetf(g,fdet)

c     Computation of deformation gradient determinant from displacement
c     gradient

      implicit   none

      real*8     g(3,3),fdet

      fdet = (1.d0+g(1,1))*((1.d0+g(2,2))*(1.d0+g(3,3)) - g(2,3)*g(3,2))
     &     + g(1,2)*(g(2,3)*g(3,1) - g(2,1)*(1.d0+g(3,3)))
     &     + g(1,3)*(g(2,1)*g(3,2) - (1.d0+g(2,2))*g(3,1))

      end
