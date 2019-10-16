c$Id:$
      subroutine uexport(r,h,neqsl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]

c     Purpose: User export interface for reduced tangent/residual from
c              slaved nodes to another program.

c     Inputs:
c          r(*)       - Residual for slaved degree of freedoms
c          h(*,*)     - Tangent for slaved degree of freedoms
c          neqsl      - Number degree of freedoms on slaved nodes

c     Outputs:
c          None

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   neqsl
      real*8    r(neqsl),h(neqsl,neqsl)

c     Output export arrays to FEAP output file.

      call mprint(r,    1,neqsl,    1,'R-export')
      call mprint(h,neqsl,neqsl,neqsl,'H-export')

      end
