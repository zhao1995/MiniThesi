c$Id:$
      subroutine uplagm(du,ulagr,lagre,ie,ix)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Check active partion for update of multiplers.   20/07/2007
c          Multiply increment update by gtan(1)
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Accumlate Lagrange multiplier unknowns

c      Inputs:
c         du(*)      - Increment to solution
c         lagre(*)   - Lagrange multiplier equation numbers
c         ie(nie,*)  - Element group control data
c         ix(nen1,*) - Element connection data

c      Outputs:
c         ulagr(*)   - Lagrange multiplier values
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'gltran.h'
      include   'part0.h'
      include   'sdata.h'

      integer    i, m,ma,mm, n,nn
      integer    lagre(*),ie(nie,*),ix(nen1,*)
      real*8     ulagr(*),du(*)

      save

      nn = 0
      do n = 1,numel
        ma = ix(nen1,n)
        if(ma.gt.0) then
          mm = lagre(n) - 1
          do m = 1,nummat
            if(ie(nie-2,m).eq.ma) then
              if(ie(nie-9,m).eq.0 .or. ie(nie-9,m).eq.npart) then
                do i = 1,ie(nie-8,m)
                  ulagr(nn+i) = ulagr(nn+i) + du(mm+i)*gtan(1)
                end do ! i
              endif
            endif
          end do ! m
        endif
        nn = nn + ndl
      end do ! n

      end
