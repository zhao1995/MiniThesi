c$Id:$
      subroutine pblend3(n,tb,iblend,ilr,isd,ndm,nen1,
     &                   prt,prth,eflag,nflag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add 'prt' and 'prth' to 'vblke' call             06/09/2007
c       2. Correct format spacings                          21/10/2007
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Construct 3 dimensional interpolation using blending

c     Inputs:
c        n         - Block number
c        tb(3,4)   - Transformation array
c        iblend(*) - Blending functions parameters/sides
c        ilr(*)    - Blending material numbers
c        isd       - Dimension for sides array
c        ndm       - Spatial dimension of mesh
c        nen1      - Dimension of ix array

c     Outputs through pointers to:
c        x(ndm,*)  - Nodal coordinates for blended patch
c        ix(nen1,*)- Element connections
c        rben(*)   - Rigid body element indicators
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'region.h'
      include  'trdata.h'
      include  'comblk.h'

      include  'p_int.h'
      include  'p_point.h'

      logical   prt,prth,eflag,nflag, setvar, palloc
      integer   isd,ndm,nr,ns,nt,nen1,ma,nrig,dlayer
      integer   i,ii,in, j, n,ne,ni,nf,ntyp, iblend(*),ilr(*)
      integer   lblend(20,6),mblend(20),iside(4)
      real*8    trdeto, t3(3,4),tb(3,4)

      save

c     Save current transformation determinant

      trdeto  = trdet
      do i = 1,4
        do j = 1,3
          t3(j,i) = tb(j,i)
        end do ! j
      end do ! i

c     Save iblend in case pointer changes

      do in = 1,20
        mblend(in) = iblend(in)
        do i = 1,6
          lblend(in,i) = 0
        end do ! i
      end do ! in

      do in = 9,1,-1
        if(mblend(in+10).ne.0) go to 100
      end do ! in

      write(iow,3000)
      write(ilg,3000)
      call plstop()

c     Determine the sides for the faces

100   call mkface(mblend,lblend)

c     Compute coordinates for each face

      ii = min(mblend(1),mblend(2),mblend(3))
      in = max(mblend(1),mblend(2),mblend(3))
      ii = mblend(1) + mblend(2) + mblend(3) - ii - in

      in = (ii+1)*(in+1)*3
      if(np(162).eq.0) then
        setvar = palloc(162,'BSIDE', 2, 1)
      endif
      setvar = palloc(117,'TEMP7', 6*in, 2)

      do i = 1,6
        call pblend2a(lblend(1,i),iside,isd)
        fp(i) = np(117) + in*(i-1)
        call pblend2b(n,hr(np(161)),mr(np(162)),t3,lblend(1,i),
     &                ilr,hr(fp(i)),mr(np(33)),mr(np(181)),
     &                iside,isd,3,nen1,.false.,.false.,.false.,nflag,
     &                .false.)
      end do ! i

c     Correct pointer if pblend2 added sides

      do i = 1,6
        fp(i) =np(117) + in*(i-1)
      end do ! i

c     Form coordinates for the mesh

      call pblend3x(hr(fp(1)),lblend(1,1),lblend(2,1),
     &              hr(fp(2)),lblend(1,2),lblend(2,2),
     &              hr(fp(3)),lblend(1,3),lblend(2,3),
     &              hr(fp(4)),lblend(1,4),lblend(2,4),
     &              hr(fp(5)),lblend(1,5),lblend(2,5),
     &              hr(fp(6)),lblend(1,6),lblend(2,6),
     &              mr(np(190)),hr(np(43)),ndm,mblend,nf)
      setvar = palloc(117,'TEMP7',    0, 2)

c     Output coordinates in blended block

      if(prt) then
        call prtitl(prth)
        write(iow,2001) (i,i=1,ndm)
        do i = mblend(4),nf
          point = np(43) + ndm*(i-1) -1
          write(iow,2002) i,(hr(j+point),j=1,ndm)
        end do ! i
      endif

      if(eflag) then

c       Compute current transformation determinant

        trdet   = t3(1,1)*(t3(2,2)*t3(3,3) - t3(2,3)*t3(3,2))
     &          + t3(1,2)*(t3(2,3)*t3(3,1) - t3(2,1)*t3(3,3))
     &          + t3(1,3)*(t3(2,1)*t3(3,2) - t3(2,2)*t3(3,1))

        nr     = mblend(3) + 1
        ns     = mblend(1) + 1
        nt     = mblend(2) + 1
        ni     = mblend(4)
        ne     = mblend(5)
        ma     = mblend(6)
        ntyp   = mblend(7)
        nreg   = mblend(10)
        nrig   = mblend(20)

        if(ma.lt.0) then
          dlayer = -ma
        else
          dlayer =  0
        end if
        call vblke(nr,ns,nt,hr(np(43)),mr(np(33)),ni,ne,nf,
     &             ndm,nen1,ma,ntyp,dlayer,ilr,prt,prth)

c       Set region rigid body number

        call psregn(mr(np(33)),mr(np(181)),nen,nen1,ne,nf,nreg,nrig,
     &              prt,prth)

      endif

c     Restore current transformation determinant

      trdet  = trdeto

c     Formats

2001  format('   B l e n d e d   C o o r d i n a t e s'//
     &       '     Node',4(i4,'-Coordinate':))

2002  format(i8,1x,1p,4e15.5)

3000  format(' *ERROR* PBLEND3: Incorrect number super nodes specified')

      end
