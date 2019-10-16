c$Id:$
      subroutine cslid2d(slid,cs0,ics, nd,nope,dnope,neps,ofssurf)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Auto surface descriptions in 2-d
c               Store slideline table and connection table for surfaces

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_comnd.h'
      include   'iofile.h'
      include   'print.h'

      integer    i,n, nd,nope,dnope,neps,ofssurf
      integer    slid(2,*), ics(dnope,*)
      real*8     cs0(nr0,n0c1:nc01)

      call cdebug0 ('      cslid2d',-1)

c     Store table values

      cs0(1,-1) = nd
      cs0(2,-1) = ofssurf
      cs0(3,-1) = neps
      cs0(4,-1) = dnope
      cs0(1, 0) = 1       ! Line
      cs0(2, 0) = nope

c     Store connection table for surface

      do n = 1,neps
        ics(1,n) = slid(1,n  )
        ics(2,n) = slid(1,n+1)
      end do ! n
      call csurflr(ics,dnope,neps)

c     Output table values

      if(prt) then
        write(iow,2000) nd,(n,(ics(i,n),i=1,4),n=1,neps)
      endif

c     Format

2000  format (/5x,'C o n t a c t   S u r f a c e   D a t a'//
     &         5x,'Data Set for Surface Number       ',i4//
     &         5x,'Contact facet type: 2-node line'/
     &        /5x,'S u r f a c e    F a c e t    C o n n e c t i o n s'
     &       //5x,'Facet  1 Node   2 Node   1 Side   2 Side'/(i10,4i8))

      end
