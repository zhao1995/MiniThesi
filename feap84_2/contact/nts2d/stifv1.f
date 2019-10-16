c$Id:$
      subroutine stifv1 (s21,c21,csi,csi2,v,vnam)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: STIFness Vector for driver # 1

c      Purpose: form vectors N0, NC,NS,T0,TC TS

c      Inputs :
c         s21     - Sinus    of unit vector tangent to master segm 12
c         c21     - Cosinus  of unit vector tangent to master segm 12
c         csi     - normalized projection of slavenode on segment
c         csi2    - Normalized distance of nodes S - 2
c         vnam    - Requested vector

c      Outputs:
c         v(*)    - vector
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character vnam*(*)
      real*8    s21,c21,csi,csi2,v(6)

      save

      call cdebug0 ('        stifv1',-1)

      if (vnam.eq.'N0') then
        v(1) =  0.d0
        v(2) =  0.d0
        v(3) = -s21
        v(4) =  c21
        v(5) =  s21
        v(6) = -c21

      elseif (vnam.eq.'NC') then
        v(1) =  s21
        v(2) = -c21
        v(3) = -s21 * csi2
        v(4) =  c21 * csi2
        v(5) = -s21 * csi
        v(6) =  c21 * csi

      elseif (vnam.eq.'NS') then
        v(1) =  s21
        v(2) = -c21
        v(3) = -s21 * csi2
        v(4) =  c21 * csi2
        v(5) = -s21 * csi
        v(6) =  c21 * csi

      elseif (vnam.eq.'T0') then
        v(1) =  0.d0
        v(2) =  0.d0
        v(3) = -c21
        v(4) = -s21
        v(5) =  c21
        v(6) =  s21

      elseif (vnam.eq.'TC') then
        v(1) =  c21
        v(2) =  s21
        v(3) = -c21 * csi2
        v(4) = -s21 * csi2
        v(5) = -c21 * csi
        v(6) = -s21 * csi

      elseif (vnam.eq.'TS') then
        v(1) =  c21
        v(2) =  s21
        v(3) = -c21 * csi2
        v(4) = -s21 * csi2
        v(5) = -c21 * csi
        v(6) = -s21 * csi
      endif

      end
