c$Id:$
      subroutine bcor3d(ixl,xl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Revise block numbering to be same as element     10/11/2008
c          numbers for 27-node Lagrange type.
c       2. Restore option for old numbering on flag oldfl   02/02/2009
c       3. Correct numbering error in emid and assignment   03/02/2009
c          of nonzero locations in ixl array
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute missing coordinate values of 27-node element

c      Inputs:
c         ixl(*)    - Nodal connection list
c         xl(3,*)   - Unadjusted coordinate array

c      Outputs:
c         xl(3,*)   - Adjusted coordinate array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'corset.h'

      integer    ixl(27), imid(12),amid(12),bmid(12),cmid(6)
      integer    dmi1(6),dmi2(6),dmi3(6),dmi4(6)
      integer    dmi5(6),dmi6(6),dmi7(6),dmi8(6)
      integer    fmi5(6),fmi6(6),fmi7(6),fmi8(6), emid(6),jmid(12)
      real*8     xl(3,27)

      integer    i,j

      save

c     New

      data       imid/ 9,10,11,12, 13,14,15,16, 17,18,19,20/
      data       amid/ 1, 2, 3, 4,  5, 6, 7, 8,  1, 2, 3, 4/
      data       bmid/ 2, 3, 4, 1,  6, 7, 8, 5,  5, 6, 7, 8/
      data       cmid/21,22,23,24,25,26/
      data       dmi1/ 1, 2, 3, 1, 1, 5/
      data       dmi2/ 4, 3, 4, 2, 2, 6/
      data       dmi3/ 8, 7, 8, 6, 3, 7/
      data       dmi4/ 5, 6, 7, 5, 4, 8/
      data       dmi5/12,10,11, 9, 9,13/
      data       dmi6/16,14,15,13,10,14/
      data       dmi7/17,18,19,17,11,15/
      data       dmi8/20,19,20,18,12,16/

c     Old

      data       jmid/13,14,15,16, 18,19,20,21,  9,10,11,12/
      data       emid/26,24,25,23,17,22/
      data       fmi5/16,14,15,13,13,18/
      data       fmi6/21,19,20,18,14,19/
      data       fmi7/ 9,10,11, 9,15,20/
      data       fmi8/12,11,12,10,16,21/

c     Numbering in 'old' form

      if(oldfl) then

c       Mid edge coordinates

        do i = 1,12
          if(ixl(jmid(i)).eq.0) then
            do j = 1,3
              xl(j,jmid(i)) = 0.5d0*(xl(j,amid(i)) + xl(j,bmid(i)))
            end do ! j
            ixl(jmid(i)) = jmid(i)
          endif
        end do ! i

c       Center face nodes

        do i = 1,6
          if(ixl(emid(i)).eq.0) then
            do j = 1,3
              xl(j,emid(i)) = 0.50d0*(xl(j,fmi5(i)) + xl(j,fmi6(i))
     &                              + xl(j,fmi7(i)) + xl(j,fmi8(i)))
     &                      - 0.25d0*(xl(j,dmi1(i)) + xl(j,dmi2(i))
     &                              + xl(j,dmi3(i)) + xl(j,dmi4(i)))
            end do ! j
            ixl(emid(i)) = emid(i)
          endif
        end do ! i

c       Center node

        if(ixl(27).eq.0) then
          do j = 1,3
            xl(j,27) = 0.125d0*(xl(j, 1) +xl(j, 2) +xl(j, 3) +xl(j, 4)
     &                        + xl(j, 5) +xl(j, 6) +xl(j, 7) +xl(j, 8))
     &               - 0.250d0*(xl(j, 9) +xl(j,10) +xl(j,11) +xl(j,12)
     &                        + xl(j,13) +xl(j,14) +xl(j,15) +xl(j,16)
     &                        + xl(j,18) +xl(j,19) +xl(j,20) +xl(j,21))
     &               + 0.500d0*(xl(j,17) +xl(j,22) +xl(j,23) +xl(j,24)
     &                        + xl(j,25) +xl(j,26))
          end do ! j
          ixl(27) = 27
        endif

c     Current numbering order on elements

      else

c       Mid edge coordinates

        do i = 1,12
          if(ixl(imid(i)).eq.0) then
            do j = 1,3
              xl(j,imid(i)) = 0.5d0*(xl(j,amid(i)) + xl(j,bmid(i)))
            end do ! j
            ixl(imid(i)) = imid(i)
          endif
        end do ! i

c       Center face nodes

        do i = 1,6
          if(ixl(cmid(i)).eq.0) then
            do j = 1,3
              xl(j,cmid(i)) = 0.50d0*(xl(j,dmi5(i)) + xl(j,dmi6(i))
     &                              + xl(j,dmi7(i)) + xl(j,dmi8(i)))
     &                      - 0.25d0*(xl(j,dmi1(i)) + xl(j,dmi2(i))
     &                              + xl(j,dmi3(i)) + xl(j,dmi4(i)))
            end do ! j
            ixl(cmid(i)) = cmid(i)
          endif
        end do ! i

c       Center node

        if(ixl(27).eq.0) then
          do j = 1,3
            xl(j,27) = 0.125d0*(xl(j, 1) +xl(j, 2) +xl(j, 3) +xl(j, 4)
     &                        + xl(j, 5) +xl(j, 6) +xl(j, 7) +xl(j, 8))
     &               - 0.250d0*(xl(j, 9) +xl(j,10) +xl(j,11) +xl(j,12)
     &                        + xl(j,13) +xl(j,14) +xl(j,15) +xl(j,16)
     &                        + xl(j,17) +xl(j,18) +xl(j,19) +xl(j,20))
     &               + 0.500d0*(xl(j,21) +xl(j,22) +xl(j,23) +xl(j,24)
     &                        + xl(j,25) +xl(j,26))
          end do ! j
          ixl(27) = 27
        endif

      endif

      end
