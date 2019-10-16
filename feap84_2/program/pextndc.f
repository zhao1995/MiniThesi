c$Id:$
      subroutine pextndc(ix,ie,ip, ic, ib)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change dum to dum(1)                             09/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute external nodes on mesh

c      Inputs:
c         ix(nen1,*)  - Element connection list
c         ip(*)       - Pointer for list of elements
c         ic(*)       - List of elements

c      Outputs:
c         ib(*)       - List external nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'eldata.h'
      include   'sdata.h'
      include   'corner.h'
      include   'iofile.h'
      include   'pointer.h'
      include   'comblk.h'

      integer    ix(nen1,numel), ip(numnp), ic(*), ib(numnp), jcon(100)
      integer    ie(nie,*)
      integer    ii, j,jmin,jj,jc,jcr, nn, ne, ncr
      real*8     dum(1)

      save

c     Look at each node to see if it has balances

      do nn = 1, numnp
        ib(nn) = 0
      end do ! nn
      jmin = 1
      do nn = 1,numnp
        jcr = 0
        do j = jmin,ip(nn)
          ne = ic(j)
          if(ne.gt.0) then
            ncorner = 0
            ma      = ix(nen1,ne)
            iel     = ie(nie-1,ma)
            nel     = 1
            do ii = nen,1,-1
              if(ix(ii,ne).ne.0) then
                nel = ii
                exit
              endif
            end do ! ii
            ii = nen1*(ne-1)               ! Element address
c           call elmlib(dum,dum,dum,mr(np(33)+ii),dum,
            call elmlib(dum,dum,dum,mr(npix+ii),dum,
     &                  hr(np(36)),hr(np(35)),ndf,ndm,nst,iel,26)
            do ncr = 1,ncorner
              ii = ix(icorner(ncr,1),ne)
              if(ii.eq.nn) then
                jj = ix(icorner(ncr,2),ne)
                if(jcr.eq.0) then
                  jcon(1) = jj
                  jcr     = 1
                else
                  do jc = 1,jcr
                    if(jcon(jc).eq.jj) then
                      jcon(jc) = 0
                      go to 100
                    endif
                  end do ! jc
                  jcr       = jcr + 1
                  jcon(jcr) = jj
 100              continue
                endif
              endif
            end do ! ncr
          endif ! ne > 0
        end do ! j

        do jc = 1,jcr
          if(jcon(jc).ne.0) then
            ib(nn)   = 1
            jcon(jc) = 0
          endif
        end do ! jc
        jmin = ip(nn) + 1
      end do ! nn

c     Check 2-d problems for midside nodes on boundary

      if(ndm.eq.2) then
        do nn = 1,numnp
          ip(nn) = 0
        end do ! nn

        do nn = 1,numel
          ne = 0
          do j = 1,nen
            if(ix(j,nn).gt.0) then
              ne = j
            endif
          end do ! j

c         Limit search to omit interior nodes, triangles and bars

          if(ne.le.3) then
            ne = 0
          elseif(ne.eq.9) then
            ne = 8
          elseif(ne.eq.16) then
            ne = 12
          endif
          do j = 1,ne
            jcr = ix(j,nn)
            if(jcr.gt.0) then
              ip(jcr) = ip(jcr) + 1
            endif
          end do ! j
        end do ! nn

        do nn = 1,numnp
          if(ip(nn).eq.1) then
            ib(nn) = 1
          endif
        end do ! nn
      endif

      end
