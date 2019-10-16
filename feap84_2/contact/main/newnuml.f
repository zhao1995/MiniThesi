c$Id:$
      subroutine newnuml(iad,id,lagre,lagrn,icx,nrenr,irb,jnt)

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
c                                           Oct. 29, 2001            1.2

c      Acronym: NEW NUMbers for Lagrange multiplier equations
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'compac.h'
      include  'part0.h'
      include  'part1.h'
      include  'rigid1.h'
      include  'rjoint.h'
      include  'sdata.h'

      integer   eqad,i,m,mm,n,nn,neqad
      integer   iad(numnp),id(ndf,numnp), nrenr(numnp)
      integer   icx(ncen1,numcels),irb(nrbdof,nrbody),jnt(6,numjts)
      integer   lagre(numel),lagrn(numnp)

      save

c     Adjust the computed values

      call adjleqn(iad,numnp)

      do n = 1,numnp
        iad(n) = - abs(iad(n))
      end do ! n

      neqad = 0
      do n = 1,numcels

c       Find maximum reorder node number on contact element
        nn = 0
        mm = 0
        do i = 1,ncen
          if(icx(i,n).gt.0) then
            if(nrenr(icx(i,n)).gt.nn) then
              nn = nrenr(icx(i,n))
              mm = n
            endif
          endif
        end do ! i

c       Adjust all other equations

        nn = mm
        if(iad(nn).lt.0) then

          eqad = 0
          do m = 1,numnp
            mm = nrenr(m)
            if(mm.gt.nn) then
              do i = 1,ndf
                if(id(i,m).gt.0 .and. ndfp(i).eq.npart) then
                  id(i,m) = id(i,m) - iad(nn)
                endif
              end do ! i
            else
              do i = 1,ndf
                if(id(i,m).gt.0 .and. ndfp(i).eq.npart) then
                  eqad = max(eqad,id(i,m))
                endif
              end do ! i
            endif
          end do ! m

c         Check rigid body equations

          if(rbody .and. nrbprt.eq.npart) then
            do m = 1,numjts
              if(jnt(6,n).gt.eqad) then
                jnt(6,n) = jnt(6,n) - iad(nn)
              end if
            end do ! m
            do m = 1,nrbody
              do i = 1,nrbdof
                if(irb(i,m) .gt. eqad) then
                  irb(i,m) = irb(i,m) - iad(nn)
                  neqr     = max(neqr,irb(i,m))
                endif
              end do ! i
            end do ! m
          endif

          do m = 1,numel*ndl
            if(lagre(m) .gt. eqad) then
              lagre(m) = lagre(m) - iad(nn)
            endif
          end do ! m

c         Set numbers for contact element

          mm      = iad(nn)
          neqad   = neqad - iad(nn)
c         iad(nn) = eqad  + lagrn(nn) + 1
          iad(nn) = eqad  + lagrn(nn)

c         Adjust other contact equations

          do i = 1,numnp
            if(iad(i).gt.iad(nn)) then
              iad(i) = iad(i) - mm
            endif
          end do ! i

        endif
      end do ! n

c     Set final active equation numbers

      neq        = neq + neqad
      nqp(npart) = neq
      if(.not.rbody .or. nrbprt.ne.npart) then
        neqr = neq
      endif

      end
