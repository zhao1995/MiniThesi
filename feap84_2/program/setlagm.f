c$Id:$
      subroutine setlagm(ilagm,ix, ie)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Modify for partition setting                     13/07/2007
c       2. Revise number of call arguments on 'setlagf'     20/07/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set element lagrange multiplier equations

c      Inputs:
c        ix(nen1,*) - Element connection array
c        ie(nie,*)  - Element control data

c      Outputs:
c        ilagm(*)   - Lagrange multiplier equation numbers
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'part0.h'
      include   'sdata.h'

      integer    ilagm(*),ix(nen1,*), ie(nie,*), m,ma, n

      save

      do n = 1,numel
        ilagm(n) = 0
      end do ! n

      ndl = 0
      do n = 1,numel
        ma = ix(nen1,n)
        if(ma.gt.0) then
          do m = 1,nummat
            if(ie(nie-2,m).eq.ma) then
              if((ie(nie-9,m).eq.0 .or. ie(nie-9,m).eq.npart) .and.
     &            ie(nie-8,m).gt.0) then
                ix(nen+4,n) = n
                ix(nen+5,n) = ilagm(n)
                ilagm(n)    = ilagm(n) + ie(nie-8,m)
                ndl         = max(ndl,ilagm(n))
              endif
            endif
          end do ! m
        endif
      end do ! n

      end

      logical function setlagf(ix,ie)

      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'part0.h'
      include   'sdata.h'

      integer    n, m,ma
      integer    ix(nen1,*),ie(nie,*)

      save

c     Test for lagrange multipliers in elements

      setlagf = .false.
      do n = 1,numel
        ma = ix(nen1,n)
        if(ma.gt.0) then
          do m = 1,nummat
            if(ie(nie-2,m).eq.ma) then

c             Check for active partition and number of multiplers

              if((ie(nie-9,m).eq.0 .or. ie(nie-9,m).eq.npart) .and.
     &            ie(nie-8,m).gt.0) then
                setlagf = .true.
                return
              endif
            endif
          end do ! m
        endif
      end do ! n

      end
