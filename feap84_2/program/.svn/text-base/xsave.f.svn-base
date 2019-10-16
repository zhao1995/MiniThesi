c$Id:$
      subroutine xsave(ir, jc, au, al, auf, alf, jp,
     &                 iunau, iunal, neq, num, iblk, nb, cfr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Store triangular parts (au/al) on disk in block form.
c               Use rfile/wfile to read/write on disk.

c      Inputs:
c         ir(*)     - row number in column for compacted form.
c         jc(*)     - pointers for columns of compacted au/al.
c         au(*)     - au in a compacted form.
c         al(*)     - al in a compacted form.
c         jp(*)     - pointers for column of profile form.
c         neq       - number of active equations.
c         num       - maximum number of terms in each block.
c         auf       - disk file for au blocks
c         alf       - disk file for al blocks
c         iunal     - logical unit associated with auf.
c         iunau     - logical unit associated with alf.
c         cfr       - unsymmetric if true

c      Outputs:
c         nb        - number of blocks au(al) take on disk.
c         iblk(3,*) - (1) first, (2) last equation in block,
c                     (3) first block for LDU decomposition.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'debugs.h'
      include  'iofile.h'

      logical   cfr
      character auf*(*), alf*(*)
      integer   iunau, iunal, neq, num, nb, i,j, ns,ne,nn
      integer   iblk(3, *), jc(neq), ir(*), jp(neq), jblk(100)
      real*8    au(*), al(*)

      save

c     Partition au/al to equal size blocks

      iblk(3,1) = 1
      nb        = 1
      ns        = 2
      nn        = 1

      do ne = 2,neq

        if(jp(ne)-jp(nn).ge.num) then
          iblk(1,nb) = ns
          iblk(2,nb) = ne - 1
          nb         = nb + 1
          ns         = ne
          nn         = ne - 1
          iblk(3,nb) = ns
        endif

c       Store smallest equation number for block

        iblk(3,nb) = min(iblk(3,nb),ne - jp(ne) + jp(ne-1))

      end do ! ne

c     Fill last block if necessary

      if(ns.le.neq) then
        iblk(1,nb) = ns
        iblk(2,nb) = neq
      else
        nb = nb - 1
      endif

      do ne = 1,nb
        jblk(ne) = jp(iblk(2,ne)) - jp(iblk(1,ne)-1)
      end do ! ne

c     Convert iblk(3,ns) to first block needed to reduce ns block

      do ns = 2,nb
        ne = iblk(3,ns)
        do nn = 1,ns-1
          if(ne.le.iblk(2,nn)) then
            iblk(3,ns) = nn
            go to 100
          endif
        end do ! ne
        iblk(3,ns) = ns - 1
100     continue
      end do ! ns

      if(debug) then
        call iprint(iblk,3,nb,3,'Block Data')
        call iprint(jblk,1,nb,1,'Block Size')
      endif

c     Write al to disk use space to fill au blocks.

      if(cfr) call wfile(iunal, alf, 1, jc(neq), al)

      do nn = 1, nb
        i = iblk(1,nn)
        j = iblk(2,nn)
        call pacau(num, neq, i, j, au, ir, jc, jp, al)
        call wfile(iunau, auf, nn, jp(j)-jp(i-1), al)
      end do ! nn

c     For unsymmetric equations read al from disk and block to disk

      if(cfr) then

        call rfile(iunal, alf, 1, jc(neq), al)

        do nn = 1, nb
          i = iblk(1,nn)
          j = iblk(2,nn)
          call pacau(num, neq, i, j, al, ir, jc, jp, au)
          call wfile(iunal, alf, nn, jp(j)-jp(i-1), au)
        end do ! nn

      endif

      end
