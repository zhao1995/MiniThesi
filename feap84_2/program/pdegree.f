c$Id:$
      subroutine pdegree(angle, sind,cosd)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Use 'pi' from pconstant.h                        09/01/2008
c       2. Use 'pi180' from pconstant.h                     22/02/2009
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute sin and cos in terms of degree angle.

c      Input:
c         angle  - Angle in degrees

c      Outputs:
c         sind   - Sine   of angle
c         cosd   - Cosine of angle
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'pconstant.h'  ! N.B. pi180 = pi/180

      real*8     angle, sind,cosd

      sind = sin(pi180*angle)
      cosd = cos(pi180*angle)

      end
