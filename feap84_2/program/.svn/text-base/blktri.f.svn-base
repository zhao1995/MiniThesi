c$Id:$
      subroutine blktri(ndm,nen1,x,ix,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate triangular block of triangular elements in 2-d

c      Inputs:
c         ndm        - Dimension of x array
c         nen1       - Dimension of ix array
c         prt        - Output generated data if true
c         prth       - Output title/header if true

c      Outputs:
c         x(ndm,*)   - Block of nodal coordinates
c         ix(nen1,*) - Block of triangular elements

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cblktr.h'
      include  'cdata.h'
      include  'iofile.h'
      include  'region.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   prt,prth, errck, pinput, tinput, pcomp
      character xh*6, ctype*15
      integer   ndm,nen1,i,j,l,k,nn,nr,nod1,nuel1,ma,n,nm,nf,ng
      integer   ntyp
      integer   ix(nen1,*),ixl(6)
      real*8    x(ndm,*),xl(3,6),td(5)

      save

      data      xh/' coord'/

c     Mesh generation routine using triangular block

100   if(ior.lt.0) then
        write(*,5000)
        call pprint('   >')
      endif
      errck = tinput(ctype,1,td,5)
      if(errck) go to 100
      nr     = nint(td(1))
      nod1   = nint(td(2))
      nuel1  = nint(td(3))
      ma     = nint(td(4))
      ntyp   = nint(td(5))
      if(ntyp.eq.0) then
        ntyp = 1
      endif

c     Check for form of coordinate system

      if(.not.pcomp(ctype,'pola',4) .and.
     &   .not.pcomp(ctype,'cyli',4) .and.
     &   .not.pcomp(ctype,'sphe',4)) then
        ctype = 'cartesian'
      endif

c     Reset to default values if necessary

      if(nod1 .eq.0) nod1  = nio + 1
      if(nuel1.eq.0) nuel1 = neo + 1
      if(ma   .eq.0) ma    = mao

      ma     = max(ma,1)
      nr     = max(nr,1)
      if(prt) then
        call prtitl(prth)
        write(iow,2000) nr,nod1,nuel1,ma,ntyp
        if(nuel1.lt.0) write(iow,2005)
        write(iow,2002) ctype,(i,xh,i=1,ndm)
        if(ior.lt.0) then
          write(*,2000) nr,nod1,nuel1,ma,ntyp
          if(nuel1.eq.0) write(*,2005)
          write(*,2002) ctype,(i,xh,i=1,ndm)
        endif
      endif

c     Initialize arrays

      do n = 1,6
        ixl(n) = 0
        do j = 1,3
          xl(j,n) = 0.0
        end do ! j
      end do ! n

c     Input block coordinates

21    if(ior.lt.0) then
        write(*,5001)
        call pprint('   >')
      endif
        errck = pinput(td,4)
        if(errck) go to 21
        l = nint(td(1))
        if(l.eq.0 ) go to 22
        if(l.gt.6) go to 402
        nm = max(nm,l)
        nn = nn + 1
        ixl(l)  = l
        xl(1,l) = td(2)
        xl(2,l) = td(3)
        xl(3,l) = td(4)
      go to 21

22    if(prt) then
        do l = 1,6
          if(ixl(l).gt.0) then
            write(iow,2001) l,(xl(i,l),i=1,ndm)
            if(ior.lt.0) then
              write(*,2001) l,(xl(i,l),i=1,ndm)
            endif
          endif
        end do ! l
      endif

c     Determine last element to be generated

      if(ntyp.eq.1) then
        nf = nuel1 + nr*nr - 1
      else
        nf = nuel1 + nr*nr/4 - 1
      endif
      if(nf.gt.numel.and.nuel1.gt.0) go to 401

c     Determine last node number to be generated

      ng = nod1 + (nr+1)*(nr+2)/2 - 1
      if(ng.gt.numnp) go to 400

c     Form block of elements

      call trblk(nr,xl,ixl,x,ix,ndm,nod1,nuel1,nen1,ma,ntyp,
     &           ctype,prt,prth)

c     Set old numbers

      if(nuel1.gt.0) neo = nf
      nio = ng
      mao = ma

c     Set node type to active

      do n = nod1,ng
        mr(np(190)+n-1) = 0
      end do ! n

c     Set region number

      do n = nuel1,nf
        ix(nen1-1,n) = nreg
      end do ! n

c     Print lists if wanted

      if(prt.and.nf.ge.nuel1) then
        do n = nuel1,nf,50
          call prtitl(prth)
          write(iow,2003) (i,i=1,nen)
          if(ior.lt.0) write(  *,2003) (i,i=1,nen)
          j = min(nf,n+49)
          do i = n,j
            write(iow,2004) i,ma,nreg,(ix(k,i),k=1,nen)
            if(ior.lt.0) write(  *,2004) i,ma,nreg,(ix(k,i),k=1,nen)
          end do ! i
        end do ! n
      endif
      return

c     Error messages

400   write(ilg,4000) ng,numnp
      write(iow,4000) ng,numnp
      if(ior.lt.0) write(  *,4000) ng,numnp
      return
401   write(ilg,4001) nf,numel
      write(iow,4001) nf,numel
      if(ior.lt.0) write(  *,4001) nf,numel

402   write(ilg,4002) l
      write(iow,4002) l
      if(ior.lt.0) write(*,4002) l

c     Formats

2000  format('   N o d e   G e n e r a t i o n s'//
     1   10x,'Number of increments      ',i5/
     2   10x,'First node number         ',i5/
     3   10x,'First element number      ',i5/
     4   10x,'Element material number   ',i5/
     4   10x,'Element type              ',i5/
     &   15x,'Type 1 = 3-node triangle'/
     &   15x,'Type 2 = 6-node triangle'/1x)

2001  format(i9,1p,3e12.3)

2002  format(/5x,'Block coordinates specified as: ',a,//,
     &       5x,'node',3(i6,a6))

2003  format('   E l e m e n t   C o n n e c t i o n s'//
     1   '   Elmt Mat Reg',8(i3,' node')/(15x,8(i3,' node')))

2004  format(i7,2i4,8i8/(15x,8i8))

2005  format(' *WARNING* No elements are generated ')

4000  format(' *ERROR* BLKTRI: Insufficient storage for nodes'/
     1   10x,'final node =',i5,5x,'numnp =',i5)

4001  format(' *ERROR* BLKTRI: Insufficient storage for elements'/
     1   10x,'final element =',i5,5x,'numel =',i5)

4002  format(' *ERROR* BLKTRI: Block node has number > 6.  No. =',i8)

5000  format(' Input: nn,nr,nod1,nuel1,ma')

5001  format(' Input: node, x-1, x-2, x-3')

      end
