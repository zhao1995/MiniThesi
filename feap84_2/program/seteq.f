c$Id:$
      subroutine seteq(id,irb,irbc,ixt,jnt,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1.  Add array 213 in call to newlnum.               02/07/2007
c       2.  Modify calls to setlagf, newlnum                13/07/2007
c       3.  Add 'point' and 'p_point.h'                     23/08/2007
c       4.  Add 'ip' (np(120)) to argument of 'pelink'      13/12/2007
c       5.  Add 'elnk' (np(257)) to argument of pelink'     25/12/2007
c       6.  Add global equation number set                  27/03/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set equation numbers for solution

c      Inputs:
c        id(ndf,*) - Boundary condition indicators
c        irbc(nrbdof,*) - Rigid body boundary restraint codes
c        ixt(*)    - Nodal rigid body number array
c        jnt(6,*)  - Joint identifier
c        ndf       - Number dof/node
c        ndm       - Spatial dimension of problem
c        numnn     - Number of nodes in mesh
c        prt       - Flag, output results if true

c      Outputs:
c        id(ndf,*) - Equation numbers for each dof.  Active dof
c                    have positive numbers for equation, fixed
c                    dof have negative numbers
c        irb(*)    - Rigid body equation numbers
c        jnt(6,*)  - Joint equations in position 6
c        nee       - Number of active equations in problem
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'idptr.h'
      include  'mxsiz.h'
      include  'part0.h'
      include  'part7.h'
      include  'pglob1.h'
      include  'pointer.h'
      include  'p_point.h'
      include  'rigid1.h'
      include  'rjoint.h'
      include  'sdata.h'
      include  'comblk.h'

      logical   prt, setlagf, setvar, palloc
      integer   n,nn,nad, i, j
      integer   id(ndf,numnp), jnt(6,numjts), ixt(numnp)
      integer   irb(nrbdof,nrbody),irbc(nrbdof,nrbody)

      save

c     If rigid body, set all rigid body nodes to fixed condition

      if(rbody .and. nrbprt.eq.npart ) then

        do n = 1,numnp

          if(ixt(n).gt.0) then

            do i = 1,ndf
              if(ndfp(i).eq.npart) then
                id(i,n) = 1
              end if
            end do ! i

          end if

        end do ! n

      end if

c     Set equation numbers for all types of points

      neq  = 0
      nad  = 0
      do n = 0,numnp-1
        nn = mr(np(89)+n)
        do i = 1,ndf
          if(ndfp(i).eq.npart) then
            j = id(i,nn)
            if(j.eq.0) then
              neq      = neq + 1
              id(i,nn) = neq
            else
              nad      = nad - 1
              id(i,nn) = nad
            endif
          endif
        end do ! i
      end do ! n

c     Adjust for master-slave nodes

      if(np(167).ne.0) then
        call setrlk(id,mr(np(167)),mr(np(100)),ndf,numnp,neq)
      endif

c     Link nodes from data

      if(lkflg .or. leflg .or. lcflg) then
        if(lkflg) then
          call plink(id,mr(np(190)),ndf,numnp,neq,prt)
        endif
        if(leflg) then
          setvar = palloc(120,'TEMP0',numnp,1)
          call pelink(id,hr(np(43)),ndm,ndf,numnp,neq,prt,mr(np(257)),
     &                mr(np(120)))
          setvar = palloc(120,'TEMP0',    0,1)
        endif
        if(lcflg) then
          call pclink(id,mr(np(79)),ndf,numnp,neq)
        endif
      endif

c     Add nodal Lagrange multiplier equations

      if(lagrfl) then
        setvar = palloc(120,'TEMP0',numnp*3,1)
        point = np(31) + ndf*numnp
        call nodlnum(mr(np(120)),mr(id31),mr( point ),mr(np( 33)),
     &               mr(np(89)),   mr(np(99)),mr(np(101)))
        setvar = palloc(120,'TEMP0',    0,1)
      endif

c     Add element Lagrange multiplier equations

      if(setlagf(mr(np(33)),mr(np(32)))) then

        if(np(211).eq.0) then
          setvar = palloc(211,'LAGRE',numel,1)
          setvar = palloc(212,'LAGRN',numnp,1)
        endif
        setvar = palloc(120,'TEMP0',numel,1)

c       Determine number of multipliers for each element

        call setlagm(mr(np(211)),mr(np(33)), mr(np(32)))

c       Adjust equation numbers

        call newlnum(mr(np(211)),mr(np(212)),mr(np(120)),mr(id31),
     &               mr(np( 33)),mr(np( 32)),mr(np(89)),   mr(np(99)),
     &               mr(np(101)))

c       Allocate storage for multipliers

        setvar = palloc(120,'TEMP0',        0,1)
        if(np(213).eq.0) then
          setvar = palloc(213,'ULAGR',ndl*numel,2)
        endif

c     Set maximum number of element multipliers to zero

      else
        ndl = 0
      endif
      ndlp(npart) = ndl

c     Check for global equations

      gneq = neq
      neq  = neq + geqnum

c     Save number of fe-equations for possible symmetry considerations

      nfeqs = neq

c     Set list of rigid body equations

      if(rbody .and. nrbprt.eq.npart ) then

c       Joint equations: Flexible/Rigid or Flexible/Flexible Coupling
c                        (Lagrange multiplier equations)

        do n = 1,numjts

          if(jnt(2,n).lt.0 .or. jnt(3,n).lt.0) then

c           Ball and Socket Joint

            if(jnt(1,n).eq.1) then
              jnt(5,n) = ndm
              jnt(6,n) = neq
              neq      = neq + ndm

c           Revolute Joint

            elseif(jnt(1,n).eq.2) then
              jnt(5,n) = ndm - 1
              jnt(6,n) = neq
              neq      = neq + ndm - 1

c           Basic Constraint Type 2

            elseif(jnt(1,n).eq.8) then
              jnt(5,n) = 1
              jnt(6,n) = neq
              neq      = neq + 1

c           Slider Joint

            elseif(jnt(1,n).eq.3) then
              jnt(5,n) = 2*(ndm - 1)
              jnt(6,n) = neq
              neq      = neq + 2*(ndm - 1)

c           Plane Joint

            elseif(jnt(1,n).eq.4) then
              jnt(5,n) = 1
              jnt(6,n) = neq
              neq      = neq + 1

c           Translation Joint

            elseif(jnt(1,n).eq.5) then
              jnt(5,n) = 2*ndm - 1
              jnt(6,n) = neq
              neq      = neq + 2*ndm - 1

c           Angle Control

            elseif(jnt(1,n).eq.6.and.hr(np(102)+9*n-1).eq.0) then
              jnt(5,n) = 1
              jnt(6,n) = neq
              neq      = neq + 1

c           Displacement Control

            elseif(jnt(1,n).eq.7.and.hr(np(102)+9*n-1).eq.0) then
              jnt(5,n) = 1
              jnt(6,n) = neq
              neq      = neq + 1
            endif

          endif ! End test for flexible/rigid

        end do ! n

        do n = 1,nrbody
          do i = 1,nrbdof
            if(irbc(i,n).eq.0) then
              neq      = neq + 1
              irb(i,n) = neq
            else
              irb(i,n) = 0
            endif
          end do ! i
        end do ! n

        neqr = neq ! Set number of equations to reduce by Gauss Elimin.

c       Joint equations: Rigid/Rigid Coupling (Lagrange multipliers)

        do n = 1,numjts

          if(jnt(2,n).ge.0 .and. jnt(3,n).ge.0) then

c           Ball and Socket Joint

            if(jnt(1,n).eq.1) then
              jnt(5,n) = ndm
              jnt(6,n) = neq
              neq      = neq + ndm

c           Revolute Joint

            elseif(jnt(1,n).eq.2) then
              jnt(5,n) = ndm - 1
              jnt(6,n) = neq
              neq      = neq + ndm - 1

c           Basic Constraint Type 2

            elseif(jnt(1,n).eq.8) then
              jnt(5,n) = 1
              jnt(6,n) = neq
              neq      = neq + 1

c           Slider Joint

            elseif(jnt(1,n).eq.3) then
              jnt(5,n) = 2*(ndm - 1)
              jnt(6,n) = neq
              neq      = neq + 2*(ndm - 1)

c           Plane Joint

            elseif(jnt(1,n).eq.4) then
              jnt(5,n) = 1
              jnt(6,n) = neq
              neq      = neq + 1

c           Translation Joint

            elseif(jnt(1,n).eq.5) then
              jnt(5,n) = 2*ndm - 1
              jnt(6,n) = neq
              neq      = neq + 2*ndm - 1

c           Angle Control

            elseif(jnt(1,n).eq.6.and.hr(np(102)+9*n-1).eq.0) then
              jnt(5,n) = 1
              jnt(6,n) = neq
              neq      = neq + 1

c           Displacement Control

            elseif(jnt(1,n).eq.7.and.hr(np(102)+9*n-1).eq.0) then
              jnt(5,n) = 1
              jnt(6,n) = neq
              neq      = neq + 1
            endif

          endif ! End test for rigid/rigid

        end do ! n

      else

c       Joint equations: Flexible/Flexible Coupling - no rigid bodies
c                        (Lagrange multiplier equations)

        do n = 1,numjts

          if(jnt(2,n).le.0 .and. jnt(3,n).le.0) then

c           Ball and Socket Joint

            if(jnt(1,n).eq.1) then
              jnt(5,n) = ndm
              jnt(6,n) = neq
              neq      = neq + ndm

c           Revolute Joint

            elseif(jnt(1,n).eq.2) then
              jnt(5,n) = ndm - 1
              jnt(6,n) = neq
              neq      = neq + ndm - 1

c           Basic Constraint Type 2

            elseif(jnt(1,n).eq.8) then
              jnt(5,n) = 1
              jnt(6,n) = neq
              neq      = neq + 1

c           Slider Joint

            elseif(jnt(1,n).eq.3) then
              jnt(5,n) = 2*(ndm - 1)
              jnt(6,n) = neq
              neq      = neq + 2*(ndm - 1)

c           Plane Joint

            elseif(jnt(1,n).eq.4) then
              jnt(5,n) = 1
              jnt(6,n) = neq
              neq      = neq + 1

c           Translation Joint

            elseif(jnt(1,n).eq.5) then
              jnt(5,n) = 2*ndm - 1
              jnt(6,n) = neq
              neq      = neq + 2*ndm - 1

c           Angle Control

            elseif(jnt(1,n).eq.6.and.hr(np(102)+9*n-1).eq.0) then
              jnt(5,n) = 1
              jnt(6,n) = neq
              neq      = neq + 1

c           Displacement Control

            elseif(jnt(1,n).eq.7.and.hr(np(102)+9*n-1).eq.0) then
              jnt(5,n) = 1
              jnt(6,n) = neq
              neq      = neq + 1
            endif

          endif ! End test for flexible/flexible

        end do ! n

        neqr = neq

      end if

      end
