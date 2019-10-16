c$Id:$
      subroutine comprol(ie,ix,lagre, neqj)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Check for Lagrange multiplier & other element equation

c      Inputs:
c         ie(nie,*)  -  Element identifier terms
c         ix(nen1,*) -  Element connection list
c         lagre(*)   -  Multiplier equation numbers

c      Outputs:
c         neqj       -  Maximum number of equations
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'sdata.h'

      integer    neqj, n, ma
      integer    ie(nie,*),ix(nen1,*),lagre(*)

      save

      do n = 1,numel
        ma   = ix(nen1,n)
        neqj = max(neqj,lagre(n) + ie(nie-8,ma) - 1)
      end do ! n

      end
