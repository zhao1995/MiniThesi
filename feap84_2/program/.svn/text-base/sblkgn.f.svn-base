c$Id:$
      subroutine sblkgn(nsblk,ndm,nel,nen1,x,ix,
     &                  ixs,ixn,nns,xs,ts,nms,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Surface mesh description by blocks:
c               Inputs: surface
c                 nn,nr,ns,nt,ni,ne,ma,nodinc
c                 node,x_node,y_node,z_node,thick_node  (node = 1,nn)
c      Inputs:
c         nsblk     - Number of surface block
c         ndm       - Spatial dimension of mesh
c         nel       - Number of nodes on surface block
c         nen1      - Dimension of ix array
c         ixs(9)    - Block nodal connections
c         ixn       - Normal node connections
c         nns(5,*)  - Block dimensions data
c         xs(3,*)   - Block nodal coordinates
c         ts(9,*)   - Surface block nodal thickness
c         nms(3,9,*)-
c         prt       - Flag, output results if true
c         prth      - Flag, output title/header if true

c      Outputs:
c         x(ndm,*)  - Nodal coordinates for block
c         ix(nen1,*)- Element nodal connections for block
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'iofile.h'

      logical   prt,prth, errck, tinput, pinput, pcomp, zflg
      character xh*6,tx*15
      integer   i,j,k,l,n,nm,nn,nr,ns,nt,nf,ng,ni,ntyp,nodinc
      integer   ndm,nel,nen1,ne,ma, neo,nio,mao,nsblk,nb
      integer   ix(nen1,1),nns(5,nsblk),ixs(9,nsblk),ixn(9,nsblk)
      real*8    x(ndm,1),xs(ndm,9,nsblk),ts(9,nsblk),nms(3,9,*)
      real*8    td(10), xmx(2)

      save

      data      xh/' coord'/, nio,neo,mao/0,0,0/

c     Block mesh generation routine

      nb     = 0
      xmx(1) = 0.0d0
      xmx(2) = 0.0d0
      zflg   = .true.

c     Check for start of block generation

100   if(ior.lt.0) then
        write(*,5000)
        call pprint('   >')
      endif
      errck = tinput(tx,1,td,0)
      if(errck) go to 100
      if(pcomp(tx,'surf',4)) then
101     if(ior.lt.0) then
          write(*,5001)
          call pprint('   >')
        endif
        errck = pinput(td,8)
        if(errck) go to 101
        nb = nb + 1

c       3-d generations

        nn     = nint(td(1))
        nr     = nint(td(2))
        ns     = nint(td(3))
        nt     = nint(td(4))
        ni     = nint(td(5))
        ne     = nint(td(6))
        ma     = nint(td(7))
        nodinc = nint(td(8))
        ntyp   = 10

c       Save dimensioning information

        nns(2,nb) = nr
        nns(3,nb) = ns
        nns(4,nb) = nt
        nns(5,nb) = ni

c       Reset to default values if necessary

        if(nio.ne.0 .and. ni.eq.0) ni = nio + 1
        if(neo.ne.0 .and. ne.eq.0) ne = neo + 1
        if(mao.ne.0 .and. ma.eq.0) ma = mao
        nodinc = max(nodinc,0)
        nr     = max(nr,1)
        ns     = max(ns,1)
        nt     = max(nt,1)
        ni     = max(ni,1)
        ma     = max(ma,1)
        if(prt) then
          call prtitl(prth)
          write(iow,2000) nr,ns,nt,ni,ne,ma,nodinc,ntyp
          if(ne.le.0) write(iow,2010)
          write(iow,2002) (i,xh,i=1,ndm)
          if(ior.lt.0) then
            write(*,2000) nr,ns,nt,ni,ne,ma,nodinc,ntyp
            if(ne.le.0) write(*,2010)
            write(*,2002) (i,xh,i=1,ndm)
          endif
        endif
        do n = 1,9
          do j = 1,ndm
            xs(j,n,nb) = 0.0
            ixs(n,nb)  = 0
          end do ! j
        end do ! n
        nm = 0
        do n = 1,nn
21        if(ior.lt.0) then
            write(*,5002)
            call pprint('   >')
          endif
          errck = pinput(td,5)
          if(errck) go to 21
          l = nint(td(1))
          if(l.eq.0) l = n
          nm = max(nm,l)
          ixs(l,nb)   = l
          xs(1,l,nb)  = td(2)
          xs(2,l,nb)  = td(3)
          xs(3,l,nb)  = td(4)
          if(zflg) then
            xmx(2) = td(2)
            zflg   = .false.
          end if
          xmx(1) = max(xmx(1),td(2),td(3),td(4))
          xmx(2) = min(xmx(2),td(2),td(3),td(4))
          ts(l,nb)    = td(5)
          if(prt) then
            write(iow,2001) l,ts(l,nb),(xs(i,l,nb),i=1,ndm)
            if(ior.lt.0) write(*,2001) l,ts(l,nb),(xs(i,l,nb),i=1,ndm)
          endif
        end do ! n
        nns(1,nb) = nm

c       3-d generations

        nf = ne + nr*ns*nt - 1
        if(nf.gt.numel.and.ne.gt.0) go to 401
        nr = nr + 1
        ns = ns + 1
        nt = nt + 1
        ng = nr*ns*nt + ni -1
        if(ng.gt.numnp) go to 400

c       Generate elements

        call sublk(nr,ns,nt,ix,ni,ne,nen1,ma)

c       Generate normals

        call snblk(nm,ndm,ixs(1,nb),xs(1,1,nb),nms(1,1,nb))

c       Set old numbers

        if(ne.gt.0) neo = nf
        nio = ng
        mao = ma

c       Print lists if wanted

        if(prt.and.ne.gt.0) then

c         Print element lists

          do n = ne,nf,50
            call prtitl(prth)
            write(iow,2005) nb,(i,i=1,nel)
            if(ior.lt.0) write(  *,2005) nb,(i,i=1,nel)
            j = min(nf,n+49)
            do i = n,j
              write(iow,2006) i,ma,(ix(k,i),k=1,nel)
              if(ior.lt.0) write(  *,2006) i,ma,(ix(k,i),k=1,nel)
            end do ! i
          end do ! n
        endif
      endif
      tx = ' '
      if(nb.lt.nsblk) go to 100

c     Generate common normals

      call stblk(nb,ndm, nns,ixn,ixs,xs, nms, xmx(1) - xmx(2))

c     Generate coordinates for nodes

      call sxblk(nb,ndm,nns,ixs,xs,ts,nms,x,prt,prth)

      return

c     Error messages

400   write(iow,2030) ng,numnp
      if(ior.lt.0) write(  *,2030) ng,numnp
      return

401   write(iow,2031) nf,numel
      if(ior.lt.0) write(  *,2031) nf,numel

c     Formats

2000  format(/'   N o d e   G e n e r a t i o n s'//
     &   10x,'Number of r-increments    ',i5/
     &   10x,'Number of s-increments    ',i5/
     &   10x,'Number of t-increments    ',i5/
     &   10x,'First node number         ',i5/
     &   10x,'First element number      ',i5/
     &   10x,'Element material number   ',i5/
     &   10x,'Node line increment       ',i5/
     &   10x,'Block type (0-10)         ',i5/1x)

2001  format(i9,1p,4e12.3)

2002  format(5x,'node   thickness',3(i6,a6))

2005  format(/
     &   '   B l o c k',i4,'   E l e m e n t   C o n n e c t i o n s'//
     &   '    Elmt  Matl',8(i3,' Node'):/(16x,8(i3,' Node')))

2006  format(i8,i6,8i8:/(14x,8i8))

2010  format(' *WARNING* No elements are generated ')

2030  format(' *ERROR* Insufficient storage for nodes'/
     &   10x,'final node =',i5,5x,'numnp =',i5)

2031  format(' *ERROR* Insufficient storage for elements'/
     &   10x,'final element =',i5,5x,'numel =',i5)

5000  format(' Input: surface')

5001  format(' Input: nn,nr,ns,ni,ne,ma,nodinc,ntyp')

5002  format(' Input: node, thick, x-1, x-2, x-3')

      end
