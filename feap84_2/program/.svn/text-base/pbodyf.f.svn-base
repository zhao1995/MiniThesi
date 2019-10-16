c$Id:$
      subroutine pbodyf(ix,ndf,nen1,numel,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Increase allowable dof to 30                     05/07/2007
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute body forces from elements and add to load
c               vector f

c      Inputs:
c         ix(nen1,*) - Element nodal connection list
c         ndf        - Number dof/node
c         nen1       - Dimension of ix array
c         numel      - Number of elements in mesh
c         prt        - Print results if true
c         prth       - Print title/header information if true

c      Outputs:
c         f(*)       - Force vector with body forces added (via pointer)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'elbody.h'
      include  'iofile.h'
      include  'pointer.h'

      logical   allreg
      character text*15
      logical   prt,prth, pcomp, tinput, errck
      integer   i,j,n,ndf,nen1,numel
      integer   ix(nen1,*)
      real*8    td(31)

      save

c     Read data from file

1     j     = min(15,ndf+1)
      errck = tinput(text,1,td,ndf+1)
      if(ndf.gt.14) then
        j     = ndf - 14
        errck = tinput(text,0,td(16),j)
      endif

c     Set body force intensities

      do i = 1,ndf
        bodyf(i) = td(i+1)
      end do ! i

c     Check for process type

      if    (pcomp(text,'mate',4)) then
        allreg = .false.
        j      = 0
      elseif(pcomp(text,'regi',4)) then
        allreg = .false.
        j      = 1
      elseif(pcomp(text,'all' ,3)) then
        allreg = .true.
        j      = 1
      elseif(pcomp(text,'elem' ,4)) then
        j      = 2
      elseif(pcomp(text,'    ',4)) then
        return
      endif

      if(prt) then
        call prtitl(prth)
        write(iow,2001) text,int(td(1)),(i,bodyf(i),i=1,ndf)
        if(ior.lt.0) then
          write(*,2001) text,int(td(1)),(i,bodyf(i),i=1,ndf)
        endif
      endif

c     Compute and assemble body loadings

      if( j.le.1 ) then

c       Do regions or materials

        do n = 1,numel
          if(allreg .or. ix(nen1-j,n) .eq. int(td(1)) ) then
            call formfe(np(40),np(27),np(27),np(27),
     &                 .false.,.true.,.false.,.true.,15,n,n,1)
          endif
        end do ! n

c     Single element case

      elseif(j.eq.2) then

        n = nint(td(1))
        call formfe(np(40),np(27),np(27),np(27),
     &                 .false.,.true.,.false.,.true.,15,n,n,1)

      endif

      go to 1

c     Format

2001  format(/'   B o d y    F o r c e s'/
     & /'      Type:',a8,' Number',i3/(5x,3(i4,'-Value',1p,e12.4:)))
      end
