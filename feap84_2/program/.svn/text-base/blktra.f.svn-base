c$Id:$
      subroutine blktra(ndm,nen1,x,ix,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Generate transition block of nodes in 2-d only

c      Inputs:
c         ndm        - Dimension for x array
c         nen1       - Dimension for ix array
c         prt        - Output generated data if true
c         prth       - Output title/header if true

c      Outputs:
c         x(ndm,*)   - Block of nodal coordinates
c         ix(nen1,*) - Block of elements
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cblktr.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'region.h'
      include  'comblk.h'

      logical   prt, prth, errck, pinput
      character xh*6
      integer   i,j,k,l,n,nm,nn,nr,ns,nt,nf,ng,ni,ntyp
      integer   ndm,nen1,ne,ma, me, n1,n2
      integer   ix(nen1,1),ixl(27)
      real*8    dx1,dx2, x(ndm,1),xl(3,27),xx(2),td(10)

      save

      data      xh/' coord'/

c     Block mesh generation routine

100   if(ior.lt.0) then
        write(*,5000)
        call pprint('   >')
      endif
      errck = pinput(td,7)
      if(errck) go to 100
      nn     = nint(td(1))
      nr     = nint(td(2))
      ns     = nint(td(3))
      ntyp   = nint(td(7))

c     2-d generations

      if(ntyp.lt.10) then
        nt     = 1
        ni     = nint(td(4))
        ne     = nint(td(5))
        ma     = nint(td(6))
c     3-d generations

      elseif(ntyp.lt.20) then
        nt     = nint(td(4))
        ni     = nint(td(5))
        ne     = nint(td(6))
        ma     = nint(td(7))

c     User generations

      else
        nt     = nint(td(4))
        ni     = nint(td(5))
        ne     = nint(td(6))
        ma     = nint(td(7))
      endif

c     Reset to default values if necessary

      if(ni.eq.0) ni = nio + 1
      if(ne.eq.0) ne = neo + 1
      if(ma.eq.0) ma = mao

      nr     = max(nr,1)
      ns     = max(ns,1)
      nt     = max(nt,1)
      ni     = max(ni,1)
      ma     = max(ma,1)
      if(prt) then
        call prtitl(prth)
        write(iow,2000) nr,ns,nt,ni,ne,ma,ntyp
        if(ne.le.0) write(iow,2010)
        write(iow,2002) (i,xh,i=1,ndm)
        if(ior.lt.0) then
          write(*,2000) nr,ns,nt,ni,ne,ma,ntyp
          if(ne.le.0) write(*,2010)
          write(*,2002) (i,xh,i=1,ndm)
        endif
      endif
      do n = 1,27
        do j = 1,ndm
          xl(j,n) = 0.0
          ixl(n) = 0
        end do ! j
      end do ! n
      nm = 0
      do n = 1,nn
21      if(ior.lt.0) then
          write(*,5001)
          call pprint('   >')
        endif
        errck = pinput(td,4)
        if(errck) go to 21
        l = nint(td(1))
        if(l.eq.0) l = n
        nm = max(nm,l)
        ixl(l)  = l
        xl(1,l) = td(2)
        xl(2,l) = td(3)
        xl(3,l) = td(4)
        if(prt) write(iow,2001) l,(xl(i,l),i=1,ndm)
        if(prt.and.ior.lt.0) write(*,2001) l,(xl(i,l),i=1,ndm)
      end do ! n

c     Determine last element number to be generated

      if (ntyp.eq.0) then
        nf = ne + 3*nr - 1
        if(nf.gt.numel.and.ne.gt.0) go to 401

c       Determine last node number to be generated

        ng = 3*(nr + nr/2) + ni + 1
        if(ng.gt.numnp) go to 400

c         Form the nodes

          dx1 = (xl(1,2) - xl(1,1))/nr
          dx2 = (xl(2,2) - xl(2,1))/nr
          do i = 0,nr
            nn               = ni + i
            mr(np(190)+nn-1) = 0
            x(1,nn)          = xl(1,1) + i*dx1
            x(2,nn)          = xl(2,1) + i*dx2
            if(prt) then
              write(iow,2007) nn,x(1,nn),x(2,nn)
              if(ior.lt.0) then
                write(*,2007) nn,x(1,nn),x(2,nn)
              endif
            endif
          end do ! i
          xx(1) = 0.5d0*(xl(1,1) + xl(1,4))
          xx(2) = 0.5d0*(xl(2,1) + xl(2,4))
          dx1   = (xl(1,2) + xl(1,3) - xl(1,1) - xl(1,4))/(4*nr)
          dx2   = (xl(2,2) + xl(2,3) - xl(2,1) - xl(2,4))/(4*nr)
          nn    = ni + nr + 1
          do i = 1,nr,2
            mr(np(190)+nn-1) = 0
            x(1,nn  )        = xx(1) + dx1
            x(2,nn  )        = xx(2) + dx2

            mr(np(190)+nn  ) = 0
            x(1,nn+1)        = xx(1) + dx1*2.d0
            x(2,nn+1)        = xx(2) + dx2*2.d0

            mr(np(190)+nn+1) = 0
            x(1,nn+2)        = xx(1) + dx1*3.d0
            x(2,nn+2)        = xx(2) + dx2*3.d0

            xx(1)            = xx(1) + dx1*4.d0
            xx(2)            = xx(2) + dx2*4.d0
            if(prt) then
              write(iow,2007) nn,x(1,nn),x(2,nn)
              write(iow,2007) nn+1,x(1,nn+1),x(2,nn+1)
              write(iow,2007) nn+2,x(1,nn+2),x(2,nn+2)
              if(ior.lt.0) then
                write(*,2007) nn,x(1,nn),x(2,nn)
                write(*,2007) nn+1,x(1,nn+1),x(2,nn+1)
                write(*,2007) nn+2,x(1,nn+2),x(2,nn+2)
              endif
            endif
            nn        = nn + 3
          end do ! i
          dx1 = 0.5d0*(xl(1,3) - xl(1,4))/nr
          dx2 = 0.5d0*(xl(2,3) - xl(2,4))/nr
          do i = 0,2*nr
            mr(np(190)+nn+i-1) = 0
            x(1,nn+i)          = xl(1,4) + i*dx1
            x(2,nn+i)          = xl(2,4) + i*dx2
            if(prt) then
              write(iow,2007) nn+i,x(1,nn+i),x(2,nn+i)
              if(ior.lt.0) then
                write(*,2007) nn+i,x(1,nn+i),x(2,nn+i)
              endif
            endif
          end do ! i

c       Form block of elements

        me  = ne

        n1 = nr + 1
        n2 = n1 + 3*nr/2
        do i = 1,nr,2

          nn = ni + i - 1

c         el 1

          ix(1,me  ) = nn
          ix(2,me  ) = nn + 1
          ix(3,me  ) = nn + n1 + 1
          ix(4,me  ) = nn + n1
          ix(nen1,me  ) = ma

c         el 2

          ix(1,me+1) = nn
          ix(2,me+1) = nn + n1
          ix(3,me+1) = nn + n2 + 1
          ix(4,me+1) = nn + n2
          ix(nen1,me+1) = ma

c         el 3

          ix(1,me+2) = nn + n1
          ix(2,me+2) = nn + n1 + 1
          ix(3,me+2) = nn + n2 + 2
          ix(4,me+2) = nn + n2 + 1
          ix(nen1,me+2) = ma

c         el 4

          ix(1,me+3) = nn + 1
          ix(2,me+3) = nn + 2
          ix(3,me+3) = nn + n1 + 2
          ix(4,me+3) = nn + n1 + 1
          ix(nen1,me+3) = ma

c         el 5

          ix(1,me+4) = nn + 2
          ix(2,me+4) = nn + n2 + 4
          ix(3,me+4) = nn + n2 + 3
          ix(4,me+4) = nn + n1 + 2
          ix(nen1,me+4) = ma

c         el 6

          ix(1,me+5) = nn + n1 + 1
          ix(2,me+5) = nn + n1 + 2
          ix(3,me+5) = nn + n2 + 3
          ix(4,me+5) = nn + n2 + 2
          ix(nen1,me+5) = ma

          me = me + 6
          n1 = n1 + 1
          n2 = n2 + 2

        end do ! i

      endif

c     Set old numbers

      if(ne.gt.0) neo = nf
      nio = ng
      mao = ma

c     Set region number

      do n = ne,nf
        ix(nen1-1,n) = nreg
      end do ! n

c     Print lists if wanted

      if(prt.and.ne.gt.0) then
        do n = ne,nf,50
          call prtitl(prth)
          write(iow,2005) (i,i=1,nen)
          if(ior.lt.0) write(*,2005) (i,i=1,nen)
          j = min(nf,n+49)
          do i = n,j
            write(iow,2006) i,ma,nreg,(ix(k,i),k=1,nen)
            if(ior.lt.0) write(  *,2006) i,ma,nreg,(ix(k,i),k=1,nen)
          end do ! i
        end do ! n
      endif

      return

c     Error messages

400   write(ilg,3000) ng,numnp
      write(iow,3000) ng,numnp
      if(ior.lt.0) write(  *,3000) ng,numnp
      return
401   write(ilg,3001) nf,numel
      write(iow,3001) nf,numel
      if(ior.lt.0) write(  *,3001) nf,numel

2000  format('   N o d e   T r a n s i t i o n s'//
     &   10x,'number of r-increments    ',i5/
     &   10x,'number of s-increments    ',i5/
     &   10x,'number of t-increments    ',i5/
     &   10x,'first node number         ',i5/
     &   10x,'first element number      ',i5/
     &   10x,'element material number   ',i5/
     &   10x,'node line increment       ',i5/
     &   10x,'block type (0-10)         ',i5/1x)
2001  format(i9,1p3e12.3)
2002  format(5x,'node',3(i6,a6))
2005  format('   E l e m e n t   C o n n e c t i o n s'//
     &   '   Elmt Mat Reg',8(i3,' Node')/(15x,8(i3,' Node')))
2006  format(i7,2i4,8i8/(15x,8i8))
2007  format(i8,1p,3e13.4)
2010  format(5x,'*WARNING* No elements are generated ')
3000  format(5x,'*ERROR* BLKTRA: Insufficient storage for nodes'/
     &   10x,'final node =',i5,5x,'numnp =',i5)
3001  format(5x,'*ERROR* BLKTRA: Insufficient storage for elements'/
     &   10x,'final element =',i5,5x,'numel =',i5)
5000  format(' Input: nn,nr,ns,ni,ne,ma,nodinc,ntyp')
5001  format(' Input: node, x-1, x-2, x-3')

      end
