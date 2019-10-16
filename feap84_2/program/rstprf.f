c$Id:$
      subroutine rstprf(jp,idl,id,ix,ie,  irb,ixt, jnt,
     &                  ndf,nen1,nen,neq,numnp,numel,nummat)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Loop over material sets to identify element      21/07/2007
c          multipliers. Add 'nummat' to argument list.
c       2. Check material number for lagrange multipliers   27/03/2008
c          Add global equation numbers
c       3. Skip counting explicit elements for profile      18/05/2012
c       4. Change ix(nen+7 to ix(nen+6: implicit/explicit   01/07/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Reset profile for elements, contact, and rigid bodies

c      Inputs:
c         id(ndf,*)     - Equation numbers for nodal dof
c         ix(nen1,*)    - Element nodal connections
c         ie(nie,*)     - Element properties
c         irb(nrbdof,*) - Rigid body equations
c         ixt(*)        - List of nodes connected to rigid bodies
c         jnt(6,*)      - Joint equations
c         ndf           - Number dof/node
c         nen1          - Dimension for ix array
c         nen           - Number nodes connected to an element
c         neq           - Number of equations active
c         numnp         - Number of nodes in mesh
c         numel         - Number of elements in mesh
c         nummat        - Number of material sets

c      Scratch:
c         idl(*)        - Store element active equations, etc.

c      Outputs:
c         jp(*)         - Row/column lengths for each equation
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdat1.h'
      include  'part0.h'
      include  'pglob1.h'
      include  'pointer.h'
      include  'rigid1.h'
      include  'rjoint.h'
      include  'comblk.h'

      logical   rbpart, rn
      integer   ndf,nen1,nen,neq,numnp,numel,nummat
      integer   i,j,ii, m,ma,ml,mm, n,nad
      integer   jp(*),idl(*),id(ndf,*),ix(nen1,*),ie(nie,*)
      integer   irb(nrbdof,nrbody), ixt(numnp), jnt(6,*)

      save

c     Zero pointer array

      do n = 1,neq
        jp(n) = 0
      end do ! n

c     Compute column heights

      do n = 1,numel

c       Test for active region for element

        rn  = (ix(nen1-1,n).ge.0) .and. (ix(nen+6,n).eq.0)
        if(rn) then

          mm  = 0
          nad = 0
          do i = 1,nen
            ii = ix(i,n)

c           Set element profile

            if(ii.gt.0) then
              do j = 1,ndf
                if( ndfp(j).eq.npart .and. id(j,ii).gt.0 ) then
                  if(mm.eq.0) mm = id(j,ii)
                  mm       = min(mm,id(j,ii))
                  nad      = nad + 1
                  idl(nad) = id(j,ii)
                end if
              end do ! j

c             Test for rigid bodies and active partition

              rbpart = rbody .and. nrbprt.eq.npart

              if( rbpart ) then

c               Add a rigid node equation to list

                if( ixt(ii).gt.0 ) then
                  do j = 1,nrbdof
                    if(irb(j,ixt(ii)) .gt. 0) then
                      nad            = nad + 1
                      idl(nad)       = irb(j,ixt(ii))
                      if(mm.eq.0) mm = idl(nad)
                      mm             = min(mm,idl(nad))
                    end if
                  end do ! j
                end if
              end if
            endif
          end do ! i

c         Add any Lagrange multiplier equations

          ma = ix(nen1,n)
          do m = 1,nummat
            if(ie(nie-2,m).eq.ma) then
              if((ie(nie-9,m).eq.0 .or. ie(nie-9,m).eq.npart) .and.
     &            ie(nie-8,m).gt.0) then
                ml = ix(nen+4,n)
                if(ml.gt.0) then
                  ml = mr(np(211)-1+ml) + ix(nen+5,n) - 1
                  do j = 1,ie(nie-8,m)
                    nad      = nad + 1
                    idl(nad) = ml + j
                    if(mm.eq.0) mm = idl(nad)
                    mm       = min(mm,idl(nad))
                  end do ! j
                endif
              endif
              if(ie(nie-10,m).gt.0) then
                do j = 1,ie(nie-10,m)
                  nad      = nad + 1
                  idl(nad) = gneq + j
                end do ! j
              endif
            endif
          end do ! m

c         Compute column heights

          do i = 1,nad
            ii = idl(i)
            jp(ii) = max(jp(ii),ii-mm)
          end do ! i

        end if

      end do ! n

c     Check for joints

      do n = 1,numjts

c       Add joint equations to list

        mm  = 0
        nad = 0

c       Ball and Socket

        if(jnt(1,n).eq.1) then
          do i = 2,3
            ii = jnt(i,n)
            if(ii.gt.0 .and. rbpart) then
              do j = 1,nrbdof
                if(irb(j,ii) .gt. 0) then
                  nad            = nad+1
                  idl(nad)       = irb(j,ii)
                  if(mm.eq.0) mm = idl(nad)
                  mm             = min(mm,idl(nad))
                end if
              end do ! j
            elseif(ii.lt.0) then
              do j = 1,3
                if( ndfp(j).eq.npart .and. id(j,-ii).gt.0 ) then
                  nad            = nad+1
                  idl(nad)       = id(j,-ii)
                  if(mm.eq.0) mm = idl(nad)
                  mm             = min(mm,idl(nad))
                end if
              end do ! j
            end if
          end do ! i

c         Set number of Lagrange multiplier constraints

          ii = 3

c       Revolute

        elseif(jnt(1,n).eq.2) then

          do i = 2,3
            ii = jnt(i,n)
            if(ii.gt.0 .and. rbpart) then
              do j = 4,nrbdof
                if(irb(j,ii) .gt. 0) then
                  nad            = nad+1
                  idl(nad)       = irb(j,ii)
                  if(mm.eq.0) mm = idl(nad)
                  mm             = min(mm,idl(nad))
                end if
              end do ! j
            elseif(ii.lt.0) then
              do j = 4,ndf
                if( ndfp(j).eq.npart .and. id(j,-ii).gt.0 ) then
                  nad            = nad+1
                  idl(nad)       = id(j,-ii)
                  if(mm.eq.0) mm = idl(nad)
                  mm             = min(mm,idl(nad))
                end if
              end do ! j
            end if
          end do ! i

c         Set number of Lagrange multiplier constraints

          ii = 2

c       Basic Constraint Type 2

        elseif(jnt(1,n).eq.8) then

          do i = 2,3
            ii = jnt(i,n)
            if(ii.gt.0 .and. rbpart) then
              do j = 4,nrbdof
                if(irb(j,ii) .gt. 0) then
                  nad            = nad+1
                  idl(nad)       = irb(j,ii)
                  if(mm.eq.0) mm = idl(nad)
                  mm             = min(mm,idl(nad))
                end if
              end do ! j
            elseif(ii.lt.0) then
              do j = 4,ndf
                if( ndfp(j).eq.npart .and. id(j,-ii).gt.0 ) then
                  nad            = nad+1
                  idl(nad)       = id(j,-ii)
                  if(mm.eq.0) mm = idl(nad)
                  mm             = min(mm,idl(nad))
                end if
              end do ! j
            end if
          end do ! i

c         Set number of Lagrange multiplier constraints

          ii = 1

c       Slider or Plane or Translator

        elseif(jnt(1,n).ge.3 .or. jnt(1,n).le.5) then

          do i = 2,3
            ii = jnt(i,n)
            if(ii.gt.0 .and. rbpart) then
              do j = 1,nrbdof
                if(irb(j,ii) .gt. 0) then
                  nad            = nad+1
                  idl(nad)       = irb(j,ii)
                  if(mm.eq.0) mm = idl(nad)
                  mm             = min(mm,idl(nad))
                end if
              end do ! j
            elseif(ii.lt.0) then
              do j = 1,ndf
                if( ndfp(j).eq.npart .and. id(j,-ii).gt.0 ) then
                  nad            = nad+1
                  idl(nad)       = id(j,-ii)
                  if(mm.eq.0) mm = idl(nad)
                  mm             = min(mm,idl(nad))
                end if
              end do ! j
            end if
          end do ! i

c         Set number of Lagrange multiplier constraints

          if(jnt(1,n).eq.3) then
            ii = 4   ! Slider
          elseif(jnt(1,n).eq.4) then
            ii = 1   ! Plane
          else
            ii = 5   ! Translator
          endif

c       Angle Control

        elseif(jnt(1,n).eq.6) then

          do i = 2,3
            ii = jnt(i,n)
            if(ii.gt.0 .and. rbpart) then
              do j = 4,nrbdof
                if(irb(j,ii) .gt. 0) then
                  nad            = nad+1
                  idl(nad)       = irb(j,ii)
                  if(mm.eq.0) mm = idl(nad)
                  mm             = min(mm,idl(nad))
                end if
              end do ! j
            elseif(ii.lt.0) then
              do j = 4,ndf
                if( ndfp(j).eq.npart .and. id(j,-ii).gt.0 ) then
                  nad            = nad+1
                  idl(nad)       = id(j,-ii)
                  if(mm.eq.0) mm = idl(nad)
                  mm             = min(mm,idl(nad))
                end if
              end do ! j
            end if
          end do ! i

c         Set number of Lagrange multiplier constraints

          if (hr(np(102)+9*n-1).eq.0.0d0) then
            ii = 1
          end if

c       Displacement Control

        elseif(jnt(1,n).eq.7) then

          do i = 2,3
            ii = jnt(i,n)
            if(ii.gt.0 .and. rbpart) then
              do j = 4,nrbdof
                if(irb(j,ii) .gt. 0) then
                  nad            = nad+1
                  idl(nad)       = irb(j,ii)
                  if(mm.eq.0) mm = idl(nad)
                  mm             = min(mm,idl(nad))
                end if
              end do ! j
            elseif(ii.lt.0) then
              do j = 4,ndf
                if( ndfp(j).eq.npart .and. id(j,-ii).gt.0 ) then
                  nad            = nad+1
                  idl(nad)       = id(j,-ii)
                  if(mm.eq.0) mm = idl(nad)
                  mm             = min(mm,idl(nad))
                end if
              end do ! j
            end if
          end do ! i

c         Set number of Lagrange multiplier constraints

          if(hr(np(102)+9*n-1).eq.0.0d0) then
            ii = 1
          end if

        end if

c       Lagrange multiplier equations

        if(hr(np(102)+9*n-1).eq.0.0d0) then

          do i = 1,ii
            nad = nad + 1
            idl(nad) = jnt(6,n) + i
          end do ! i

c         Set column height for Lagrange multipliers at joint

          do i = 1,nad
            ii = idl(i)
            jp(ii) = max(jp(ii),ii-mm)
          end do ! i

        end if

      end do ! n

      end
