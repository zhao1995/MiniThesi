c$Id:$
      subroutine cccontab (lcc0)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: Contact Commands CONtrol TABle

c      Purpose: Set up control tables for each command

c      Inputs :

c      Outputs:
c         lcco    - Length of storage vector CC0
c         nr0     - # of rows for all tables
c         nc0(*)  - # of columns for each command
c         of0(*)  - offset within CC0 for command
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'c_dict.h'
      include  'c_comnd.h'

      integer   lcc0, k,km1

      save

      call cdebug0 ('  cccontab',-1)

c     Set dimension of each command matrix (features + user columns)

      do k=1, c_ncc
        nc0(k) = nfe(k)+nuc(k)
      end do

c     Global offsets for command array C0

      of0(1) = 1
      do k = 2 ,c_ncc
        km1 = k-1
        of0(k) = of0(km1) + nr0*(n0c(km1)+1+nc0(km1))*cck(km1)
      end do

c     Total length of commands array

      lcc0 = of0(c_ncc) + nr0*(n0c(c_ncc)+1+nc0(c_ncc))*cck(c_ncc)

c     Invert sign of # of system columns for control table dimensioning

      do k = 1 ,c_ncc
        n0c(k) = -n0c(k)
      end do

c     Provisional transfer for Fortran 77 compiler restrictions

      nc01  = nc0(1)
      nc02  = nc0(2)
      nc03  = nc0(3)

      nc09  = nc0(8)
      nc010 = nc0(9)
      nc011 = nc0(10)

      n0c1  = n0c(1)
      n0c2  = n0c(2)
      n0c3  = n0c(3)

      n0c9  = n0c(8)
      n0c10 = n0c(9)
      n0c11 = n0c(10)

      end
