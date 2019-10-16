c$Id:$
      subroutine pblend1b(xs,is,trb,iblend,ilr,x,ix,rben,
     &                    iside,isd,ndm,nen1,prt,prth,eflag,nflag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change ctype to character 15                     03/11/2011
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Construct one dimensional interpolation using blending

c     Inputs:
c        xs(3,*)   - Blending supernode connections
c        is(isd,*) - Blending side supernode lists
c        trb       - Transformation for blending coordinates
c        iblend(*) - Blending functions parameters/sides
c        ilr(*)    - Material quantities for blends
c        isd       - Dimension for sides array
c        ndm       - Spatial dimension of mesh
c        nen1      - Dimension of ix array

c     Outputs:
c        x(ndm,*)  - Nodal coordinates for blended patch
c        ix(nen1,*)- Element connections
c        rben(*)   - Rigid body markers
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cblend.h'
      include   'cdata.h'
      include   'iofile.h'
      include   'pointer.h'
      include   'region.h'
      include   'trdata.h'
      include   'comblk.h'

      logical    prt,prth,eflag,nflag, setvar, palloc
      character  ctype*15
      integer    isd,ndm,nen1,nrig, nsn, iside
      integer    i, j,jj, k, ma
      integer    ne,nf,ni,nm,nn,nr,ns,nodinc,ntyp, styp, dlayer
      integer    is(isd,*),iblend(*), ix(nen1,*),rben(*), ilr(*)
      real*8     xs(3,*),trb(3,4),x(ndm,*), trdeto

      save

c     Get edge interpolations

      nr   = iblend(1)

      setvar = palloc ( 111, 'TEMP1',(nr+1)*ndm  ,2)
      setvar = palloc ( 112, 'TEMP2',(nr+1)      ,2)
      setvar = palloc ( 113, 'TEMP3',(nr+1)*3    ,2)

      nreg = iblend(10)
      nrig = iblend(20)
      jj   = abs(iside)
      styp = is(1,jj)
      do j = isd,2,-1
        if(is(j,jj).ne.0) go to 110
      end do ! j
      write(  *,3000) j
      write(ilg,3000) j
      call plstop()
110   nsn = j - 1

      call pside1(nr, xs, trb, iside,is(2,jj), nsn,ndm,
     &            hr(np(112)),hr(np(113)), hr(np(111)),styp)

      ni  = iblend(3)

      call pblend1x(nn,nr,ni,ndm, hr(np(111)),mr(np(190)),x,
     &              nflag,prt,prth)

      setvar = palloc ( 113, 'TEMP3',0 ,2)
      setvar = palloc ( 112, 'TEMP2',0 ,2)
      setvar = palloc ( 111, 'TEMP1',0 ,2)

      if(eflag) then
        ne     = iblend(4)
        ma     = iblend(5)
        ntyp   = iblend(6)
        nodinc = 0
        ctype  = 'blen'
        dlayer = 0
        ns     = max(1,iblend(2))  ! Sets line order
        nm     = 2                 ! Force line element generation
        nr     = nr + 1

        trdeto = trdet
        trdet  = trb(1,1)*(trb(2,2)*trb(3,3) - trb(2,3)*trb(3,2))
     &         + trb(1,2)*(trb(2,3)*trb(3,1) - trb(2,1)*trb(3,3))
     &         + trb(1,3)*(trb(2,1)*trb(3,2) - trb(2,2)*trb(3,1))

        call sblke(nr,ns,x,ix,ni,ne,nn,ndm,nen1,nodinc,ntyp,nm,ma,
     &             dlayer,ilr,ctype)
        trdet  = trdeto
        nf     = nn
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
        end do ! nn

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

3000  format(' *ERROR* PBLEND1: No side nodes found for side',i4)

      end
