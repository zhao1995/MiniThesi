c$Id:$
      subroutine pnewpr()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Recompute new profile after changes

c      Inputs:
c         none      - Transfers all through common blocks

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'cdata.h'
      include   'compas.h'
      include   'complx.h'
      include   'fdata.h'
      include   'idptr.h'
      include   'part0.h'
      include   'part1.h'
      include   'part3.h'
      include   'pointer.h'
      include   'rjoint.h'
      include   'comblk.h'

      logical    palloc, setvar
      integer    nparto, k1,k2

      save

c     Compute new equation numbers and profile

      nparto =  npart
      do k1 = 1,4
        if(.not.tflp(k1)) then
          call partpt(k1,tflp(k1),.false.)

c         Set new current profile

          call profil(mr(np(20+k1)),mr(np(34)),mr(id31),
     &                mr(np(33)),1,pfr)
          call profil(mr(np(20+k1)),mr(np(34)),mr(id31),
     &                mr(np(33)),2,pfr)
          nqp(k1)    =  neq
          nqr(k1)    =  neqr
          mxprop(k1) = (mr(np(20+k1)+neq-1))*ipc
          mxneqp(k1) =  neq*ipc

c         Compute new fill for sparse storage

          if(ittyp.gt.-3) then
            if(ittyp.eq.-1 .and. np(92).ne.0) then
              setvar = palloc(92,'OINB',0,1)
            endif
            k2 = 0
            call iters(k2,1)
          endif
        endif
      end do ! k1
      npart = nparto
      call partpt(npart,tflp(npart),.true.)

      end
