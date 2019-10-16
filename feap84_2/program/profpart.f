c$Id:$
      subroutine profpart(eq,id)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Split id(*,*,2) into eq(*,*) and id(*,*)         06/05/2009
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Compute profile for current partition

c      Inputs:
c        id(ndf,numnp) - Boundary condition codes
c        npart         - Active partition number
c        ndfp(ndf)     - Partition/dof map
c      Outputs:
c        eq(ndf,numnp) - New equation numbers
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'complx.h'
      include   'mxsiz.h'
      include   'sdata.h'
      include   'part0.h'
      include   'part1.h'
      include   'part3.h'
      include   'print.h'
      include   'rjoint.h'

      include   'pointer.h'
      include   'comblk.h'

      integer    eq(ndf,numnp), id(ndf,numnp), i,n

      save

c     Test for active partition

      if(.not.tflp(npart)) then

c       Move boundary condition codes

        do i = 1,ndf
          if(npart .eq. ndfp(i)) then
            do n = 1,numnp
              eq(i,n) = id(i,n)
            end do ! n
          endif
        end do ! i

c       Set current profile

        mxprop(npart) = 0
        mxneqp(npart) = 0
        call profil(mr(np(20+npart)),mr(np(34)),eq,mr(np(33)),1,prt)
        call profil(mr(np(20+npart)),mr(np(34)),eq,mr(np(33)),2,prt)
        nqp(npart)    = neq
        nqr(npart)    = neqr
        mxprop(npart) = max(mxprop(npart),(mr(np(20+npart)+neq-1))*ipc)
        mxneqp(npart) = max(mxneqp(npart),neq*ipc)
        mxpro         = max(mxpro,mxprop(npart))
        mxneq         = max(mxneq,mxneqp(npart))
      else
        write(*,*) ' *WARNING* Inactive partition'
      endif

      end
