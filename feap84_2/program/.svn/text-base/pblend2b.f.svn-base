c$Id:$
      subroutine pblend2b(n,xs,is,trb,iblend,ilr,x,ix,rben,
     &                    iside,isd,ndm,nen1,prt,prth,eflag,nflag,capfl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Increase to ctype*15 for compatibility           06/09/2010
c       2. Move set of 'nf' to output correct numbers       08/02/2011
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Construct two dimensional interpolation using blending

c     Inputs:
c        n         - Block number
c        xs(3,*)   - Blending supernode connections
c        is(isd,*) - Blending side supernode lists
c        trb       - Transformation for blending coordinates
c        iblend(*) - Blending functions parameters/sides
c        ilr(*)    - Material quantities for blends
c        iside(*)  - Side connections for face
c        isd       - Dimension for sides array
c        ndm       - Spatial dimension of mesh
c        nen1      - Dimension of ix array

c     Outputs:
c        x(ndm,*)  - Nodal coordinates for blended patch
c        ix(nen1,*)- Element connections
c        rben(*)   - Rigid body markers
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'iofile.h'
      include   'pointer.h'
      include   'region.h'
      include   'trdata.h'
      include   'comblk.h'

      logical    prt,prth,eflag,nflag, setvar, palloc, capfl
      character  ctype*15
      integer    i,ii,in, j,jj, k, ma,m1,m2
      integer    n,ne,nf,ni,nm,nn,nr,ns,nodinc,ntyp, styp, dlayer
      integer    nsn(4), iside(4), isd,ndm,nen1,nrig
      integer    is(isd,*),iblend(*), ix(nen1,*),rben(*), ilr(*)
      real*8     trdeto, xs(3,*),trb(3,4),x(ndm,*)

      save

c     Set up face

      do j = 1,4
        iblend(j+10) = iside(j)
      end do ! j

c     Set for 4 sides

      in = 4

c     Check signs on sides for blend

      call mkside(n,iblend(11),is,isd)

c     Check for matching end nodes for blending functions

      do i = 1,in
        ii = iblend(i+10)
        if(ii.gt.0) then
          if(is(1,abs(ii)).eq.2) then
            do k = 3,isd
              if(is(k,abs(ii)).ne.0) m1 = k
            end do ! k
          else
            m1 = 3
          endif
        elseif(ii.lt.0) then
          m1 = 2
        else
          write(  *,3001)
          write(ilg,3001)
          write(iow,3001)
        endif

        j  = mod(i,in) + 1
        jj = iblend(j+10)
        if(jj.gt.0) then
          m2 = 2
        elseif(jj.lt.0) then
          if(is(1,abs(jj)).eq.2) then
            do k = 3,isd
              if(is(k,abs(jj)).ne.0) m2 = k
            end do ! k
          else
            m2 = 3
          endif
        else
          write(  *,3001)
          write(ilg,3001)
        endif

        if(is(m1,abs(ii)).ne.is(m2,abs(jj))) then
          if(ior.lt.0) write(*,3000) n,i,j,ii,jj
          write(ilg,3000) n,i,j,ii,jj
          write(iow,3000) n,i,j,ii,jj
        endif

      end do ! i

c     Quadrilateral region blends

      if(in.eq.4) then

c       Determine the number of nodes for each side
        do i = 1,in
          ii = abs(iblend(i+10))
          do j = isd,2,-1
            if( is(j,ii).ne.0 ) go to 110
          end do ! j
          write(  *,3002) i
          write(ilg,3002) i
          call plstop()
110       nsn(i) = j-1
        end do ! i

c       Get edge interpolations

        nr   = iblend(1)
        ns   = iblend(2)
        ntyp = iblend(6)

        setvar = palloc ( 111, 'TEMP1',(nr+1)*ndm  ,2)
        setvar = palloc ( 112, 'TEMP2',(ns+1)*ndm  ,2)
        setvar = palloc ( 113, 'TEMP3',(nr+1)*ndm  ,2)
        setvar = palloc ( 114, 'TEMP4',(ns+1)*ndm  ,2)
        setvar = palloc ( 115, 'TEMP5',max(nr,ns)+1,2)
        setvar = palloc ( 116, 'TEMP6',(nr+1)*3    ,2)

        nreg = iblend(10)
        nrig = iblend(20)
        ii   = iblend(11)
        jj   = abs(ii)
        styp = is(1,jj)

        call pside1(nr, xs, trb, ii,is(2,jj), nsn(1),ndm,
     &              hr(np(115)),hr(np(116)), hr(np(111)),styp)

        ii   = iblend(12)
        jj   = abs(ii)
        styp = is(1,jj)

        call pside1(ns, xs, trb, ii,is(2,jj), nsn(2),ndm,
     &              hr(np(115)),hr(np(116)), hr(np(112)),styp)

        ii   =-iblend(13)
        jj   = abs(ii)
        styp = is(1,jj)

        call pside1(nr, xs, trb, ii,is(2,jj), nsn(3),ndm,
     &              hr(np(115)),hr(np(116)), hr(np(113)),styp)

        ii   =-iblend(14)
        jj   = abs(ii)
        styp = is(1,jj)
        call pside1(ns, xs, trb, ii,is(2,jj), nsn(4),ndm,
     &              hr(np(115)),hr(np(116)), hr(np(114)),styp)

        ni = iblend(3)
        call pblendx(nn,nr,ns,ni,ntyp,ndm, hr(np(111)),hr(np(112)),
     &               hr(np(113)),hr(np(114)),mr(np(190)),x,
     &               nflag,prt,prth, capfl)

        setvar = palloc ( 116, 'TEMP6',0 ,2)
        setvar = palloc ( 115, 'TEMP5',0 ,2)
        setvar = palloc ( 114, 'TEMP4',0 ,2)
        setvar = palloc ( 113, 'TEMP3',0 ,2)
        setvar = palloc ( 112, 'TEMP2',0 ,2)
        setvar = palloc ( 111, 'TEMP1',0 ,2)

        if(eflag) then
          ne     = iblend(4)
          ma     = iblend(5)
          nm     = 4
          nodinc = 0
          ctype  = 'blen'

          nr     = nr + 1
          ns     = ns + 1

          trdeto = trdet
          trdet  = trb(1,1)*(trb(2,2)*trb(3,3) - trb(2,3)*trb(3,2))
     &           + trb(1,2)*(trb(2,3)*trb(3,1) - trb(2,1)*trb(3,3))
     &           + trb(1,3)*(trb(2,1)*trb(3,2) - trb(2,2)*trb(3,1))

          if(ma.lt.0) then
            dlayer = -ma
          else
            dlayer =  0
          endif

          nn = nn + 1
          call sblke(nr,ns,x,ix,ni,ne,nn,ndm,nen1,nodinc,ntyp,nm,ma,
     &               dlayer,ilr,ctype)
          nf    = nn
          nn    = nn - 1

          trdet = trdeto

        endif

      endif

c     Set region and rigid body numbers

      if(eflag) then
        do nn = ne,nf
          ix(nen1-1,nn) = nreg
          if(nrig.ge.1) then
            rben(nn) = nrig
          else
            rben(nn) = -1
          endif
        enddo ! nn

c       Print lists if wanted

        if(prt.and.ne.gt.0) then
          do nn = ne,nf,50
            call prtitl(prth)
            write(iow,2000) (i,i=1,nen)
            if(ior.lt.0) then
              write(  *,2000) (i,i=1,nen)
            endif
            j = min(nf,nn+49)
            do i = nn,j
              ma = ix(nen1,i)
              write(iow,2001) i,ma,nreg,(ix(k,i),k=1,nen)
              if(ior.lt.0) then
                write(  *,2001) i,ma,nreg,(ix(k,i),k=1,nen)
              endif
            end do ! i
          end do ! nn
        endif
      endif

c     Formats

2000  format('   E l e m e n t   C o n n e c t i o n s'//
     &   '   Elmt Mat Reg',8(i3,' node'):/(15x,8(i3,' node')))

2001  format(i7,2i4,8i8:/(15x,8i8))

3000  format(' *ERROR* PBLEND2: Blending function',i3/
     &       '         End node 2 for side',i2,' not same as'/
     &       '         End node 1 for side',i2/
     &       '         Node-2 =',i8,' Node-1 =',i8)

3001  format(' *ERROR* PBLEND2: Zero side number specified')

3002  format(' *ERROR* PBLEND2: No side nodes found for side',i4)

      end
