c$Id:$
      subroutine ckfixed(ndtyp,id)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    06/12/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Check number of fixed degree-of-freedoms on nodes

c      Inputs:
c         ndtyp(*)    - Command option
c         id(ndf,*)   - Command parameters

c      Outputs:
c         Values reported to screen/output file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'iofile.h'
      include   'sdata.h'

      integer    i,n, nfix, ntot
      integer    ndtyp(*), id(ndf,numnp)

      if(ior.lt.0) then
        write(*,2000)
      endif
      write(iow,2000)

      do i = 1,ndf
        nfix = 0
        ntot = 0
        do n = 1,numnp
          if(ndtyp(n).ge.0) then
            if(id(i,n).ne.0) then
              nfix = nfix + 1
            endif
            ntot = ntot + 1
          endif
        end do ! n

c       Output result for DOF i

        if(ior.lt.0) then
          write(*,2001) i,nfix
        endif
        write(iow,2001) i,nfix
      end do ! i

      if(ior.lt.0) then
        write(*,2002) ntot
      endif
      write(iow,2002) ntot

c     Formats

2000  format( 5x,'Restrained Degree-of-Freedom Check'/
     &       15x,'DOF   Fixed'/10x,16('-') )

2001  format(10x,2i8)

2002  format(10x,16('-')/16x,'On',i8,' Nodes'/1x)

      end
