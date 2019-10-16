c$Id:$
      subroutine pnumne(ix,nen1,nen,numnp,numel, rev,ren,ip,
     &                  iq,numnp1,numel1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: New node/element numbers for mesh output
c      Inputs:
c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'compac.h'

      integer   nen1,nen, numnp,numel, numnp1,numel1, i,j,nn
      integer   ix(nen1,*), rev(*),ren(*),ip(*),iq(*)

      save

      rev(1) = 0

      if(optmsh) then

        do nn = 1,numnp
          ip(nn)         = 0
          rev(ren(nn)+1) = nn
        end do ! nn

c       Compute active number of nodes and elements

        numel1 = 0
        numnp1 = 0
        do nn = 1,numel
          if(ix(nen1-1,nn).ge.0) then
            numel1 = numel1 + 1
            do i = 1,nen
              if(ix(i,nn).gt.0) then
                j      = rev(ix(i,nn)+1)
                ip(j)  = 1
                numnp1 = max(numnp1,j)
              endif
            end do ! i
          endif
        end do ! nn

        if(opthoit) then              ! Wilson/Hoit optimizer
          do nn = 1,numnp
            ip(nn) = ren(nn)
          end do ! nn
        else                          ! Sloan optimizer
          numnp1 = 0
          do nn = 1,numnp
            iq(nn) = 1
          end do ! nn
          do nn = 1,numnp
            if(ip(nn).ne.0) then
              numnp1     = numnp1 + 1
              ip(numnp1) = ren(nn)
              iq(nn)     = 0
            endif
          end do ! nn
          do nn = numnp1+1,numnp
            ip(nn) = 0
          end do ! nn
          do nn = 2,numnp
            iq(nn) = iq(nn) + iq(nn-1)
          end do ! nn
          do nn = 1,numnp
            rev(nn+1) = rev(nn+1) - iq(rev(nn+1))
          end do ! nn
        endif

c     Mesh not optimized: Remove unused elements, nodes

      else

c       Mark nodes used in mesh

        do nn = 1,numnp
          ip(nn) = 0
        end do ! nn

        numel1 = 0
        do nn = 1,numel
          if(ix(nen1-1,nn).ge.0) then
            numel1 = numel1 + 1
            do i = 1,nen
              if(ix(i,nn).gt.0) then
                ip(ix(i,nn))  = 1
              endif
            end do ! i
          endif
        end do ! nn

c       Set new node numbers

        numnp1 = 0
        do nn = 1,numnp
          if(ip(nn).gt.0) then
            numnp1     = numnp1 + 1
            rev(nn+1)  = numnp1
            ip(numnp1) = nn
          end if
        end do ! nn

      endif

      end
