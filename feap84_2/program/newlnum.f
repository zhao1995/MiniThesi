c$Id:$
      subroutine newlnum(lagre,lagrn,lagmn,id,ix,ie,ren,irb,jnt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add array lagmn to determine multipler equation  02/07/2007
c          numbering.  Revise algorithm.
c       2. Add array 'ie' and number multiplers based on    20/07/2007
c          active partition and region
c       3. Check for nodal Lagrange multipliers and use to  28/09/2007
c          set final numbers
c-----[--.----+----.----+----.-----------------------------------------]
c      Compute equation structure for lagrange multiplier unknowns

c      Inputs:
c        id(ndf,*)     - Solid element equation numbers
c        ix(nen1,*)    - Element connection array
c        ie(nie,*)     - Element control data
c        ren(*)        - Nodal reorder list
c        irb(nrbdof,*) - Rigid body equation numbers
c        jnt(6,numjts) - Joint equation numbers

c      Outputs
c        lagre(numel)    - First element equation number
c        lagrn(numel)    - Number equations after node number
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat1.h'
      include  'part0.h'
      include  'part1.h'
      include  'rigid1.h'
      include  'rjoint.h'
      include  'sdata.h'

      logical   eflag,pflag
      integer   eqad,neqad,i,m,ma,mm,e,n,nn
      integer   lagre(numel),lagrn(numnp),lagmn(numel)
      integer   id(ndf,numnp),ix(nen1,*),ie(nie,*)
      integer   ren(numnp),irb(nrbdof,nrbody),jnt(6,numjts)

      save

c     Zero nodal and element array

      do nn = 1,numnp
        lagrn(nn) = 0
      end do ! nn

c     Set maximum node number on element

      do e = 1,numel
        if(lagre(e).gt.0) then

c         Find maximum reordered node number on element

          nn = 0
          do i = 1,nen
            if(ix(i,e).gt.0) then
              nn = max(nn,ren(ix(i,e)))
            endif
          end do ! i
          lagmn(e) = nn
        else
          lagmn(e) = numnp + 1
        endif
      end do ! e

c     Search elements

      neqad = 0
      eqad  = 0
      eflag = .true.
      do while (eflag)

        eflag = .false.
        nn    = numnp + 1
        do mm = 1,numel
          if(lagmn(mm).lt.nn) then
            e  = mm               ! Next element to process
            nn = lagmn(mm)        ! Element maximum node number
          endif
        end do ! mm

        if(nn.le.numnp) then
          lagmn(e) = numnp + 1
          eflag    = .true.
        endif

c       Adjust all other equations

        if(eflag) then
          pflag = .true.
          ma    = ix(nen1,e)
          do m = 1,nummat
            if(ie(nie-2,m).eq.ma .and. ie(nie-8,m).gt.0 .and.
     &        (ie(nie-9,m).eq.    0     .or.
     &         ie(nie-9,m).eq.npart)  ) then
              pflag = .false.
              do n = 1,numnp
                mm = ren(n)
                if(mm.gt.nn) then
                  do i = 1,ndf
                    if(id(i,mm).gt.0 .and. ndfp(i).eq.npart .and.
     &                                     ndfl(i).eq.  0 ) then
                      id(i,mm) = id(i,mm) + lagre(e)
                    endif
                  end do ! i
                else
                  do i = 1,ndf
                    if(id(i,mm).gt.0 .and. ndfp(i).eq.npart .and.
     &                                     ndfl(i).eq.  0 ) then
                      eqad     = max(eqad,id(i,mm))
                    endif
                  end do ! i
                endif
              end do ! n

              do n = 1,numnp
                mm = ren(n)
                do i = 1,ndf
                  if(id(i,mm).gt.eqad
     &                   .and.ndfp(i).eq.npart.and.ndfl(i).ne.0) then
                    id(i,mm) = id(i,mm) + lagre(e)
                  endif
                end do ! i
              end do ! n

c             Check rigid body equation numbers

              if(rbody .and. nrbprt.eq.npart) then
                do n = 1,numjts
                  if(jnt(6,n).gt.eqad) then
                    jnt(6,n) = jnt(6,n) + lagre(e)
                  end if
                end do ! n
                do n = 1,nrbody
                  do i = 1,nrbdof
                    if(irb(i,n) .gt. eqad) then
                      irb(i,n) = irb(i,n) + lagre(e)
                      neqr     = max(neqr,irb(i,n))
                    endif
                  end do ! i
                end do ! n
              endif

c             Set numbers for element

              mm        = lagre(e)
              if(lagrn(nn).eq.0) then
                lagre(e)  = eqad  + 1
              else
                lagre(e)  = lagrn(nn) + 1
              endif
              neqad     = neqad + mm
              eqad      = eqad  + mm
              lagrn(nn) = eqad
            endif
          end do ! m
          if(pflag) then
            lagmn(e) = numnp + 1
          endif
        endif ! eflag
      end do ! while

c     Set final active equation number

      neq = neq + neqad
      if(.not.rbody .or. nrbprt.ne.npart) then
        neqr = neq
      endif

      end
