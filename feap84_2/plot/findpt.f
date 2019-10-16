c$Id:$
      subroutine findpt(va,sne,ix,nd,nen, sv,npv)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Find intersection points with line for dplot/splot

c      Inputs:
c         va(nd,*) - Value of plot variable at nodes
c         sne(2,*) - ??
c         ix(*)    - List of element nodes
c         nd       - Number of nodes in mesh
c         nen      - Number of nodes on element

c      Outputs:
c         sv(2,*)  - Values of plot variable and location
c         npv      - Number of values in sv
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   nd,nen,npv
      integer   i,n,ne,nel,n1,ne1
      real*8    ratio,svtry

      integer   ix(nen)
      real*8    sne(2,nen),va(nd,*),sv(2,*)

      save

      do n = 1,nen
        if(ix(n).ne.0) nel = n
      end do ! n

      nel = min(nel,4)

c     Find intersections

      n1 = nel
      ne1= abs(ix(n1))

      do n = 1,nel
        ne = abs(ix(n))
        if(ne.gt.0) then
          if(sne(2,n1)*sne(2,n).lt.0.0d0) then
            ratio   = sne(2,n)/(sne(2,n) - sne(2,n1))
            svtry   = sne(1,n) + ratio*(sne(1,n1) - sne(1,n))
            if(npv.gt.0) then
              do i = 1,npv
                if(abs(svtry - sv(1,i)).lt.1.d-4) go to 110
              end do ! i
            endif
            npv = npv + 1
            sv(1,npv) = svtry
            sv(2,npv) = va(1,ne) + ratio*(va(1,ne1) - va(1,ne))
          endif
110       n1 = n
          ne1= ne
        endif
      end do ! n

      end
