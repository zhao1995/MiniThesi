c$Id:$
      subroutine nodlnum(lagre,id,ibc,ix,ren,irb,jnt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove set of lagre(2:3,i) -- i = 0              23/08/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Compute equation structure for nodal lagrange multiplier unknowns

c      Inputs:
c        id(ndf,*)     - Nodal equation numbers
c        ibc(ndf,*)    - Nodal boundary conditions
c        ix(nen1,*)    - Element connection array
c        ren(*)        - Nodal reorder list
c        irb(nrbdof,*) - Rigid body equation numbers
c        jnt(6,numjts) - Joint equation numbers

c      Outputs
c        lagre(3,numnp)- First nodal equation number
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'iofile.h'
      include  'part0.h'
      include  'part1.h'
      include  'rigid1.h'
      include  'rjoint.h'
      include  'sdata.h'

      logical   adjust, lagmul
      integer   eqad,eqng,i,m,mm,n,nd,nn,ll
      integer   lagre(3,numnp),id(ndf,numnp),ibc(ndf,numnp),ix(nen1,*)
      integer   ren(numnp),irb(nrbdof,nrbody),jnt(6,numjts)

      save

c     Test for Lagrange multipliers

      adjust = .true.
      do i = 1,ndf
        if(ndfp(i).eq.npart .and. ndfl(i).gt.0) then
          adjust = .false.
        endif
      end do ! i
      if(adjust) return

c     Zero nodal array

      do nn = 1,numnp
        lagre(1,nn) = 0
        lagre(2,nn) = 0
        lagre(3,nn) = 0
      end do ! nn

c     Search elements

      do n = 1,numel

c       Check for active multipliers on element

        lagmul = .false.
        do i = 1,nen
          if(ix(i,n).gt.0) then
            do m = 1,ndf
              if(ndfp(m).eq.npart .and. ndfl(m).gt.0 .and.
     &                           ibc(m,ix(i,n)).eq.0) then
                lagmul = .true.
                go to 100
              endif
            end do ! m
          endif
        end do ! i

c       Find maximum reordered node number on element

100     if(lagmul) then
          nn = 0
          do i = 1,nen
            if(ix(i,n).gt.0) then
              do m = 1,ndf
                if(ndfp(m).eq.npart .and. ndfl(m).eq.0 .and.
     &                             ibc(m,ix(i,n)).eq.0) then
                  nn = max(nn,ren(ix(i,n)))
                  exit
                endif
              end do ! m
            endif
          end do ! i
          do i = 1,nen
            if(ix(i,n).gt.0) then
              do m = 1,ndf
                if(ndfp(m).eq.npart .and. ndfl(m).gt.0 .and.
     &                             ibc(m,ix(i,n)).eq.0) then
                  lagre(1,ren(ix(i,n))) = max(lagre(1,ren(ix(i,n))),nn)
                endif
              end do ! m
            endif
          end do ! i
        endif

      end do ! n

c     Set lagre(1,n) greater than number of nodes if zero

      do nd = 1,numnp
        if(lagre(1,nd).eq.0) then
          lagre(1,nd) = numnp + 1
        endif
      end do ! nd

c     Set equation numbers for non-Lagrange multipliers

      eqad  = 0
      do nd = 1,numnp
        n = numnp+1
        i = 0
        do nn = 1,numnp
          if(lagre(2,nn).eq.0 .and. lagre(1,nn).lt.n) then
            n = lagre(1,nn)
            i = nn
          endif
        end do ! nn
        if(i.gt.0 .and. i.le.numnp) then
          lagre(2, i) = 1
          lagre(3,nd) = i
        endif
        nn = ren(nd)
        do i = 1,ndf
          if(ndfp(i).eq.npart) then
            if(ndfl(i).eq.0 .and.ibc(i,nn).eq.0) then
              eqad     = eqad + 1
              id(i,nn) = eqad
            else
              id(i,nn) = 0
            endif
          endif
        end do ! i
      end do ! nd

      eqad  = 0
      do n = 1,numnp
        ll = lagre(3,n)
        if(ll.gt.0) then
          nn = lagre(1,ll)
          if(nn.gt.0 .and. nn.le.numnp) then

c           Adjust all other equations

            do m = 1,numnp
              mm = ren(m)
              if(mm.gt.nn .and. mm.le.numnp) then
                do i = 1,ndf
                  if(id(i,mm).gt. 0   .and.
     &               ndfp(i).eq.npart .and.
     &               ndfl(i).eq.  0  ) then
                    id(i,mm) = id(i,mm) + 1
                  endif
                end do ! i
              else
                do i = 1,ndf
                  if(id(i,mm).gt.0 .and. ndfp(i).eq.npart) then
                    eqad     = max(eqad,id(i,mm))
                  endif
                end do ! i
              endif
            end do ! m

c           Check rigid body equation numbers

            if(rbody .and. nrbprt.eq.npart) then
              do m = 1,numjts
                if(jnt(6,m).gt.eqad) then
                  jnt(6,m) = jnt(6,m) + 1
                end if
              end do ! m
              do m = 1,nrbody
                do i = 1,nrbdof
                  if(irb(i,m) .gt. eqad) then
                    irb(i,m) = irb(i,m) + 1
                    neqr     = max(neqr,irb(i,m))
                  endif
                end do ! i
              end do ! m
            endif

c           Set numbers for node

            do i = 1,ndf
              if(ibc(i,ll).eq.0   .and.
     &           ndfp(i).eq.npart .and.
     &           ndfl(i).gt.  0  ) then
                eqad     = eqad + 1
                id(i,ll) = id(i,ll) + eqad
              endif
            end do ! i
          endif
        endif

      end do ! n

c     Set zero equations to negative to force boundary modifications

      eqng = 0
      do n = 1,numnp
        do i = 1,ndf
          if(ndfp(i).eq.npart .and. id(i,n).eq.0) then
            eqng    = eqng - 1
            id(i,n) = eqng
          endif
        end do ! i
      end do ! n

c     Set final active equation number

      if(.not.rbody .or. nrbprt.ne.npart) then
        neqr = neq
      endif

      end
