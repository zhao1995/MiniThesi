c$Id:$
      subroutine psregn(ix,rben,nen,nen1,ne,nf,nreg,nrig,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Add region numbers & output list of generated elements

c      Inputs:
c        ix(nen1,*)   - Element nodal connection list
c        nen          - Maximum number of nodes on element
c        nen1         - Dimension of ix array
c        ne           - Initial element number generated
c        nf           - Final   element number generated
c        nreg         - Region number
c        nrig         - Rigid  number
c        prt          - Print list if true
c        prth         - Print header if true

c      Outputs:
c        ix(nen1,*)   - Element nodal connection list with region number
c        rben(*)      - Rigid body indicator on element
c        lblend(20,*) - Global blend array
c-----[--.----+----.----+----.-----------------------------------------]

      implicit   none

      include   'iofile.h'

      logical    prt,prth
      integer    nen,nen1,ne,nf,nreg,nrig,i,in,j,ma,ix(nen1,*),rben(*)

      save

c     Set rigid body and region indicators

      do i = ne,nf
        ix(nen1-1,i) = nreg
        if(nrig.ge.1) then
          rben(i) = nrig
        else
          rben(i) = -1
        endif
      end do ! i

c     Output lists if wanted

      if(prt.and.ne.gt.0) then
        do in = ne,nf,50
          call prtitl(prth)
          write(iow,2003) (i,i=1,nen)
          if(ior.lt.0) then
            write(  *,2003) (i,i=1,nen)
          endif
          j = min(nf,in+49)
          do i = in,j
            ma = ix(nen1,i)
            write(iow,2004) i,ma,nreg,(ix(j,i),j=1,nen)
            if(ior.lt.0) then
              write(  *,2004) i,ma,nreg,(ix(j,i),j=1,nen)
            endif
          end do ! i
        end do ! in
      endif

2003  format('   E l e m e n t   C o n n e c t i o n s'//
     &   '   Elmt Mat Reg',8(i3,' node'):/(15x,8(i3,' node')))

2004  format(i7,2i4,8i8:/(15x,8i8))

      end
