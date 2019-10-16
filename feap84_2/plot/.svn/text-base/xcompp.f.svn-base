c$Id:$
      subroutine xcompp(ie,ix,id,nie,nen1,nen,numnp,numel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change pstyp.gt.0 to pstyp.ne.0                  31/08/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Determine number of elements connected to each node
c               for use in mesh plot/outline routine xpline.

c      Inputs:
c         ie(nie,*)   - Material set assembly data
c         ix(nen1,*)  - Element nodal connection list
c         nie         - Dimension of ie array
c         nen1        - Dimension of ix array
c         nen         - Number of nodes connected to element
c         numnp       - Number of nodes in mesh
c         numel       - Number of elements in mesh

c      Outputs:
c         id(*)       - Connection data counter data
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nie,nen1,nen,numnp,numel, i,iu,i1,i2,j,ma,n, nel,pstyp
      integer   ie(nie,*),ix(nen1,*),id(*),iplt(50)

      save

      do n = 1,numnp+1
        id(n) = 0
      end do ! n

      do n = 1,numel
        if(ix(nen1,n).gt.0) then
          pstyp = ie(1,ix(nen1,n))
          if(pstyp.ne.0) then
            ma = ie(nie-1,ix(nen1,n))
            do i = nen,1,-1
              if(ix(i,n).gt.0) then
                nel = i
                exit
              endif
            end do ! i
            call plftyp(pstyp,nel,ma)
            call pltord(ix(1,n),ma,iu,iplt)
            do i = 1,iu-1
              i1 = iplt(i)
              if( i1.gt.0 .and. i1.le.nen ) then
                do j = i+1,iu
                  i2 = iplt(j)
                  if( i2.gt.0 .and. i2.le.nen ) then
                   i1     = min(ix(i1,n),ix(i2,n))
                   id(i1) = id(i1) + 1
                   go to 100
                  end if
                end do ! j
              end if
100           continue
            end do ! i
          endif
        endif
      end do ! n

      j     = 1
      i     = id(1)
      id(1) = 1
      do n = 2,numnp+1
        j     = j + i
        i     = id(n)
        id(n) = j
      end do ! n

      end
