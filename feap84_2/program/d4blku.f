c$Id:$
      subroutine d4blku( al , au , jp , jh, ie, j, je, neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Factor unsymmetric array A pairs of row/columns
c               Programmed for super-scalar processors

c      Inputs:
c        al(*)  - Unreduced pair of rows
c        au(*)  - Unreduced pair of columns
c        jp(*)  - Pointer array to end of columns
c        jh(2)  - height of two columns
c        ie     - Start equation for reduction
c        j      - Equation number of second column

c      Outputs:
c        al(*)  - Reduced pair of rows
c        au(*)  - Reduced pair of columns
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, imid, ie, ic, j, jmid, je, jmi, n, neq
      integer   jp(*), jh(2), ih(2), iu(4), il(4)
      real*8    dot, s00, s01, s10, s11, al(*), au(*)

      save

c     Reduce row and column in 2x2 blocks

      do i = ie, j-1, 2
        jmi  = j - i
        ic   = max(1,i-1)
        ih(1)  = jp(i  ) - jp(ic)
        ih(2)  = jp(i+1) - jp(i ) - 1

c       If all row/columns have non-zero lengths in computation
c       perform triangular decomposition step.

        imid = min( ih(1), ih(2), jh(1)-jmi, jh(2)-jmi)
        imid = max(0,imid)
        jmid  = imid + jmi

c       Locate pointers for 2 columns

        iu(1) = jp(j)+1 - jmi
        iu(2) = iu(1)   + 1
        iu(3) = jp(j+1) - jmi
        iu(4) = iu(3)   + 1

        s00   = 0.0d0
        s10   = 0.0d0
        s01   = 0.0d0
        s11   = 0.0d0

c       Reduce top parts of any columns using conventional dots

        if(i.lt.neq) then

          il(1) = jp(i)+1
          il(2) = jp(i+1)

          if( jh(1).gt.jmid ) then
            if( ih(1).gt.imid ) then
              ic    = min( ih(1), jh(1) - jmi )
              s00   = s00   + dot(au(il(1)-ic),au(iu(1)-ic),ic-imid )
            endif
            if( ih(2).gt.imid ) then
              ic    = min( ih(2), jh(1) - jmi )
              s10   = s10   + dot(au(il(2)-ic),au(iu(1)-ic),ic-imid )
            endif
          endif

          if( jh(2).gt.jmid ) then
            if( ih(1).gt.imid ) then
              ic    = min( ih(1), jh(2) - jmi )
              s01   = s01   + dot(au(il(1)-ic),au(iu(3)-ic),ic-imid )
            endif
            if( ih(2).gt.imid ) then
              ic    = min( ih(2), jh(2) - jmi )
              s11   = s11   + dot(au(il(2)-ic),au(iu(3)-ic),ic-imid )
            endif
          endif

c         Reduce remaining parts using dot product on 2x2 blocks

          if(imid.gt.0) then
            do n = 0,imid-1
              s00   = s00   + au(il(1)-imid+n)*au(iu(1)-imid+n)
              s10   = s10   + au(il(2)-imid+n)*au(iu(1)-imid+n)
              s01   = s01   + au(il(1)-imid+n)*au(iu(3)-imid+n)
              s11   = s11   + au(il(2)-imid+n)*au(iu(3)-imid+n)
            end do ! n
          endif

c         Restore reduced terms

          au(iu(1)) = au(iu(1)) - s00
          au(iu(3)) = au(iu(3)) - s01

c         Clean up last off-diagonal term in second i-row

          if(ih(2).ge.0) then
            au(iu(2)) = au(iu(2)) - (au(il(2))*au(iu(1)) + s10  )
            au(iu(4)) = au(iu(4)) - (au(il(2))*au(iu(3)) + s11  )
          else
            au(iu(2)) = au(iu(2)) - s10
            au(iu(4)) = au(iu(4)) - s11
          endif

c         Reduce top parts of any columns using conventional dots

        else
          il(1) = jp(i)+1 - je
          il(2) = jp(i+1) - je

          if( jh(1).gt.jmid ) then
            if( ih(1).gt.imid ) then
              ic    = min( ih(1), jh(1) - jmi )
              s00   = s00   + dot(al(il(1)-ic),au(iu(1)-ic),ic-imid )
            endif
            if( ih(2).gt.imid ) then
              ic    = min( ih(2), jh(1) - jmi )
              s10   = s10   + dot(al(il(2)-ic),au(iu(1)-ic),ic-imid )
            endif
          endif

          if( jh(2).gt.jmid ) then
            if( ih(1).gt.imid ) then
              ic    = min( ih(1), jh(2) - jmi )
              s01   = s01   + dot(al(il(1)-ic),au(iu(3)-ic),ic-imid )
            endif
            if( ih(2).gt.imid ) then
              ic    = min( ih(2), jh(2) - jmi )
              s11   = s11   + dot(al(il(2)-ic),au(iu(3)-ic),ic-imid )
            endif
          endif

c         Reduce remaining parts using dot product on 2x2 blocks

          if(imid.gt.0) then
            do n = 0,imid-1
              s00   = s00   + al(il(1)-imid+n)*au(iu(1)-imid+n)
              s10   = s10   + al(il(2)-imid+n)*au(iu(1)-imid+n)
              s01   = s01   + al(il(1)-imid+n)*au(iu(3)-imid+n)
              s11   = s11   + al(il(2)-imid+n)*au(iu(3)-imid+n)
            end do ! n
          endif

c         Restore reduced terms

          au(iu(1)) = au(iu(1)) - s00
          au(iu(3)) = au(iu(3)) - s01

c         Clean up last off-diagonal term in second i-row

          if(ih(2).ge.0) then
            au(iu(2)) = au(iu(2)) - (al(il(2))*au(iu(1)) + s10  )
            au(iu(4)) = au(iu(4)) - (al(il(2))*au(iu(3)) + s11  )
          else
            au(iu(2)) = au(iu(2)) - s10
            au(iu(4)) = au(iu(4)) - s11
          endif
        endif

        s00   = 0.0d0
        s10   = 0.0d0
        s01   = 0.0d0
        s11   = 0.0d0
        il(1) = iu(1) - je
        il(2) = iu(2) - je
        il(3) = iu(3) - je
        il(4) = iu(4) - je

c       Reduce top parts of any columns using conventional dots

        if( jh(1).gt.jmid ) then
          if( ih(1).gt.imid ) then
            ic    = min( ih(1), jh(1) - jmi )
            s00   = s00   + dot(au(jp(i)+1-ic),al(il(1)-ic),ic-imid )
          endif
          if( ih(2).gt.imid ) then
            ic    = min( ih(2), jh(1) - jmi )
            s10   = s10   + dot(au(jp(i+1)-ic),al(il(1)-ic),ic-imid )
          endif
        endif

        if( jh(2).gt.jmid ) then
          if( ih(1).gt.imid ) then
            ic    = min( ih(1), jh(2) - jmi )
            s01   = s01   + dot(au(jp(i)+1-ic),al(il(3)-ic),ic-imid )
          endif
          if( ih(2).gt.imid ) then
            ic    = min( ih(2), jh(2) - jmi )
            s11   = s11   + dot(au(jp(i+1)-ic),al(il(3)-ic),ic-imid )
          endif
        endif

c       Reduce remaining parts using dot product on 2x2 blocks

        if(imid.gt.0) then
          do n = 0,imid-1
            s00   = s00   + au(jp(i)+1-imid+n)*al(il(1)-imid+n)
            s10   = s10   + au(jp(i+1)-imid+n)*al(il(1)-imid+n)
            s01   = s01   + au(jp(i)+1-imid+n)*al(il(3)-imid+n)
            s11   = s11   + au(jp(i+1)-imid+n)*al(il(3)-imid+n)
          end do ! n
        endif

c       Restore reduced terms

        al(il(1)) = al(il(1)) - s00
        al(il(3)) = al(il(3)) - s01

c       Clean up last off-diagonal term in second i-row

        if(ih(2).ge.0) then
          al(il(2)) = al(il(2)) - (au(jp(i+1))*al(il(1)) + s10  )
          al(il(4)) = al(il(4)) - (au(jp(i+1))*al(il(3)) + s11  )
        else
          al(il(2)) = al(il(2)) - s10
          al(il(4)) = al(il(4)) - s11
        endif

      end do ! i

      end
