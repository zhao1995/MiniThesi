c$Id:$
      subroutine sethis(ie,ix,rben,nie,nen,nen1,
     &                  numel,nummat,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Correct counting of nhf for case of more than    15/11/2008
c          one material attached to an element.
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set up history addresses in ix array

c      Inputs:
c         ie(nie,*) - Material set assembly information
c         rben(*)   - List of rigid body to nodes
c         nie       - Dimension of ie array
c         nen       - Number of nodes/element
c         nen1      - Dimension of ix array
c         numel     - Number of elements in mesh
c         nummat    - Number of material sets in mesh
c         prt       - Flag, output results if true

c      Outputs:
c         ix(nen1,*)- History data pointers added to positions nen+1,
c                     nen+2 and nen+3
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'hdatam.h'
      include  'iofile.h'

      logical   flag,prt,setvar,palloc
      integer   i,n,nh0,nhf,nie,nen,nen1,numel,nummat,ma
      integer   ie(nie,*),ix(nen1,*),rben(*)

      save

c     Output amount of memory used in each user element

      if(prt) then
        flag = .true.
        do n = 1,nummat

c         Outputs

          if(ie(nie-1,n).gt.0) then
            if(flag) then
              flag = .false.
              if(ior.lt.0) write(*,2000)
              write(iow,2000)
            endif
            if(ior.lt.0) then
              write(*,2001) n,ie(nie-2,n),ie(nie-1,n),ie(nie,n),
     &                      ie(nie-5,n)
            endif
            write(iow,2001) n,ie(nie-2,n),ie(nie-1,n),ie(nie,n),
     &                      ie(nie-5,n)
          endif
        end do ! n
      endif

c     Compute maximum length necessary to store history variables

      nhmax  = 0
      nh3max = 0
      do n = 1,nummat
        nh0 = 0
        nhf = 0
        do i = 1,nummat
          if(ie(nie-2,i).eq.n) then
            ie(nie-3,i) = nh0
            ie(nie-4,i) = nhf
            nh0         = nh0 + ie(nie,  i)
            nhf         = nhf + ie(nie-5,i)
          end if
        end do ! i
        nhmax  = max(nhmax, ie(nie,n))
        nh3max = max(nh3max,ie(nie-5,n))
      end do ! n

c     Set pointers for history variables into ix-array

      nh0 = 0
      do n = 1,numel

c       Set number of history terms for rigid elements to zero

        if(rben(n).le.0) then

          ma  = ix(nen1,n)

c         Variable storage history

          nhf = 0
          do i = 1,nummat
            if(ie(nie-2,i).eq.ma) nhf = nhf + ie(nie,i)
          end do ! i
          if(nhf.gt.0) then
            ix(nen+1,n) = nh0
            nh0 = nh0 + nhf
            ix(nen+2,n) = nh0
            nh0 = nh0 + nhf
            nhf = 0
          endif

c         Fixed storage history

          do i = 1,nummat
            if(ie(nie-2,i).eq.ma) nhf = nhf + ie(nie-5,i)
          end do ! i
          if(nhf.gt.0) then
            ix(nen+3,n) = nh0
            nh0 = nh0 + nhf
          endif
        endif
      end do ! n
      nhf = nh0
      if (nhf.gt.0) then
        setvar = palloc( 49,'H    ',nhf,2)
      endif

c     Formats

2000  format(/10x,'Material    Element Tag   Element Type',
     &            '  History Terms  Element Terms')

2001  format(5i15)

      end
