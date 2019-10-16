c$Id:$
      subroutine pblendm(isd,blend,ndm,nen1,prt,prth,eflag,nflag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change iside to iside(1) line 94                 09/01/2012
c-----[--.----+----.----+----.-----------------------------------------]

c     Purpose:  Construct interpolation using blending functions

c     Inputs:
c        isd       - Dimension for sides array
c        blend     - Dimension for blending array
c        ndm       - Spatial dimension of mesh
c        nen1      - Dimension of ix array
c        prt       - Print control
c        prth      - Print header control
c        eflag     - Element generation flag
c        nflag     - Nodal generation flag

c     Outputs stored by pointer for:
c        x(ndm,*)  - Nodal coordinates for blended patch
c        ix(nen1,*)- Element connections
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cblend.h'
      include  'cdata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'region.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   prt,prth,eflag,nflag, setvar,palloc, capfl
      integer   n,n1,isd,blend,ndm,nen1
      integer   iside(4),tblend(20)

      save

      do n = 1,numbd

c       Pointer to transformtion and surface type

        fp(1) = np(164) + blend*(n-1) - 1

c       Surface generations

        do n1 = 1,blend
          tblend(n1) = mr(fp(1)+n1)
        end do ! n1

c       Surface generations

        if(tblend(19).eq.1 .or. tblend(19).eq.4) then

          capfl = tblend(19).eq.4

          if(np(162).eq.0) then
            setvar = palloc( 162,'BSIDE',2,1)
          endif
          call pblend2a(tblend,iside,isd)
          fp(1) = np(166) + mxilr*(n-1)
          fp(2) = np(163) +    12*(n-1)
          call pblend2b(n,hr(np(161)),mr(np(162)),hr(fp(2)),tblend,
     &                  mr(fp(1)),hr(np(43)),mr(np(33)),mr(np(181)),
     &                  iside,isd,ndm,nen1,prt,prth,eflag,nflag,capfl)

c       Solid generations

        elseif(tblend(19).eq.2) then

          fp(1) = np(166) + mxilr*(n-1)
          fp(2) = np(163) +    12*(n-1)
          call pblend3(n,hr(fp(2)),tblend,mr(fp(1)),isd,ndm,nen1,
     &                 prt,prth,eflag,nflag)

c       Line generations

        elseif(tblend(19).eq.3) then

          if(np(162).eq.0) then
            setvar = palloc( 162,'BSIDE',2,1)
          endif
          call pblend1a(mr(np(162)),tblend,iside,isd)
          fp(1) = np(166) + mxilr*(n-1)
          fp(2) = np(163) +    12*(n-1)
          call pblend1b(hr(np(161)),mr(np(162)),hr(fp(2)),tblend,
     &                  mr(fp(1)),hr(np(43)),mr(np(33)),mr(np(181)),
     &                  iside(1),isd,ndm,nen1,prt,prth,eflag,nflag)

        endif

      end do ! n

      end
