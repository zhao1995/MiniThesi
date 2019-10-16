c$Id:$
      subroutine newnum(iad,id,nren,nrenr,irb,jnt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor            June 15, 1997            1.0
c                                           Febr 21, 2001            1.1

c      Acronym: NEW NUMbers for Lagrange multiplier equations
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'part0.h'
      include  'part1.h'
      include  'rigid1.h'
      include  'rjoint.h'
      include  'sdata.h'

      logical   check
      integer   eqad,i,j,m,mm,n,nn,nr
      integer   iad(numnp),id(ndf,numnp), nren(numnp),nrenr(numnp)
      integer   irb(nrbdof,nrbody),jnt(6,numjts)

      save

c     Adjust the computed values

      call adjleqn(iad,numnp)

      neqr = max(neq,neqr)
      do n = 1,numnp
        j  = n
        j  = nren(n)
        nn = nrenr(j)
        if(iad(j).gt.0) then
          nr = iad(j)
          check = .true.
          eqad  = neqr
          do m = 1,numnp
            mm = nrenr(m)
            if(mm.gt.nn) then
              do i = 1,ndf
                if(id(i,m).gt.0 .and. ndfp(i).eq.npart) then
                  check   = .false.
                  eqad    = min(eqad,id(i,m))
                  id(i,m) = id(i,m) + nr
                  neqr    = max(neqr,id(i,m))
                endif
              end do ! i
            endif
          end do ! m
          if(rbody .and. nrbprt.eq.npart) then
            do m = 1,numjts
              if(jnt(6,j).gt.eqad) then
                jnt(6,j) = jnt(6,j) + nr
                neq      = max(neq,jnt(6,j))
              end if
            end do ! m
            do m = 1,nrbody
              do i = 1,nrbdof
                if(irb(i,m) .gt. eqad) then
                  irb(i,m) = irb(i,m) + nr
                  neqr     = max(neqr,irb(i,m))
                  neq      = max(neq,neqr)
                endif
              end do ! i
            end do ! m
          else
            neq = neqr
          endif

c         Set numbers for node-n

          if(check) then
            neqr   = neqr + iad(j)
            iad(j) = eqad
          else
            iad(j) = eqad - 1
          endif
        endif
      end do ! n

      neq        = neqr
      nqp(npart) = neq

      end
