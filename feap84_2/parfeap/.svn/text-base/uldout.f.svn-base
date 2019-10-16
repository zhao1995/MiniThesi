c$Id$
      subroutine uldout(ldtab,ldnod,ldval, nrv,nf,partn,revx)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    05/04/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: output forces and displacements for load tables

c      Inputs:
c        ldtab(4,2,*)   - Pointer, length prop no. table
c        ldnod(*)       - List of nodes
c        ldval(ndf,*)   - List of values

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'iodata.h'
      include   'pload1.h'
      include   'prld1.h'
      include   'sdata.h'

      include   'pointer.h'
      include   'comblk.h'

      logical    lread, pflg, domnd, loadfl
      integer    i,j,n,nn,nq, nod, nd,nl,npl, nrv, nf, ng
      integer    ldtab(4,2,*),ldnod(*), partn(*), revx(2,*)
      real*8     ldval(ndf,*)

c     Loop over groups

      integer    iz
      data       iz / 0 /

      do n = 1,ldtot

c       Do forces (1) and displacements (2)

        loadfl = .true.
        do i = 1,2
          nl = ldtab(2,i,n)
          if(nl.gt.0) then

c           Set proportional load value for group

            npl    = ldtab(3,i,n)

c           Loop over nodes in group

            lread = .true.
            nd = ldtab(1,i,n)
            do nn = 1,nl
              nq  = ldnod(nd+nn)
              domnd = .false.
              do ng = 1,nf
                nod = revx(1,nrv+ng)
                if(nod.eq.nq) then
                  domnd = .true.
                  exit
                endif
              end do ! ng

c             For case where node belongs to 'nq'

              if(domnd) then
                if(partn(nod).gt.0) then
                  pflg = .false.
                  do j = 1,ndf
                    if(ldval(j,nd+nn).ne.0.0d0) then
                      pflg = .true.
                    endif
                  enddo ! j
                  if(pflg) then

c                   Write start of load table

                    if(loadfl) then
                      if(npl.le.0) then
                        write(ios,2001) 'LOAD'
                      else
                        write(ios,2001) 'LOAD PROP',npl
                      endif
                      loadfl = .false.
                    endif

c                   Write type of loading condition

                    if(lread) then
                      if(i.eq.1) then
                        write(ios,2001) '  FORCE conditions'
                      else
                        write(ios,2001) '  DISPlacement conditions'
                      endif
                      lread = .false.
                    endif
                    write(ios,2002) revx(2,nrv+ng),iz,
     &                             (ldval(j,nd+nn),j=1,ndf)
                  endif
                endif
              endif
            end do ! nn

          endif

        end do ! i

c       Write end of load table

        if(.not.loadfl) write(ios,2001) 'LOAD END'
      end do ! n

c     Formats

2001  format(/a:,i10)
2002  format(i9,i2,1p,14e14.6)

      end
