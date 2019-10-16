      subroutine acheck(x,y,n0,nl,nlc)
c----------------------------------------------------------------------
c      Purpose:   Parse a string to find fields separated by commas

c      Inputs:
c         x(*) -  Character string of data to parse
c         n0   -  Field width for parsed data
c         nl   -  Total number of characters in x-array
c         nlc  -  Total number of characters in y-array

c      Outputs:
c         y(*) -  Parsed data in field widths of n0
c----------------------------------------------------------------------
      character*1 x(*),y(*),macr*4

      do 100 ii = nl,1,-1
       if(x(ii).ne.' ') go to 110
100   continue

110   do 150 i = 1,nlc
       y(i) = ' '
150   continue

      k = 0
      il= 0
      do 200 i = 1,ii
        if(x(i).eq.',') then
          k  = k + n0
          if(k.gt.nlc-n0) go to 210
          il = k - i
        else
          y(i+il) = x(i)
        end if
200   continue
      k  = k + n0
cww210   call just(y,k,n0)
210   continue
c.... justify alphanumeric data in a string:- numbers right - alphanumeric remain left
c>>w
c.... do not in case of macros REST,1234 and SAVE,1234  (1234 is then a file extension)
      do im=1,4
        macr(im:im)=y(im)
      end do
      if(macr.eq.'save'.or.macr.eq.'rest'.or.
     +   macr.eq.'SAVE'.or.macr.eq.'REST') return
cww<<
      call just(y,k,n0)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine aload(x,f,ix,id,ndm,ndf,nen1,prt)
c-----------------------------------------------------------------------
c
c      Purpose:   set constant loads for defined areas
c                 (x/y[z]_min<center of element<x/y[z]_max)
c
c      Inputs:
c         x(*)        - Nodal coordinates
c         f(ndf,numnp)- Nodal force values
c         ix(nen1,*)  - Element nodal connections of mesh
c         id(ndf,*)   - Equation numbers for each active dof
c         ndm         - Spatial dimension of mesh
c         ndf         - Number dof/node
c         nen1        - Dimension for ix array
c         prt         - print flag
c
c      Outputs:
c         f(ndf)      - Nodal force values
c
c     data to read:
c     1  x_min,xmax,ymin,ymax,[zmin,zmax]
c     2  mdf:    dof to load,               NODES
c        ietyp:  4  4 node                | 1, 2, 3, 4
c               -8  8 node brick(bottom)  | 1, 2, 3, 4
c                9  9 node                | 1, 2, 3, 4, 5, 6, 7, 8, 9
c              -27 27 node brick(bottom)  | 1, 2, 3, 4,13,14,15,16,17
c                8  8 node brick(top),    | 5, 6, 7, 8
c               27 27 node brick(top),    | 5, 6, 7, 8,18,19,20,21,22
c----------------------------------------------------------------------
c     open
c     approximation: element is loaded if center is in area
c-----------------------------------------------------------------------
      USE cdata
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      dimension x(ndm,numnp),f(ndf,numnp),ix(nen1,*),id(ndf,*),td(16)
      dimension xl(ndm,9),shp(3,9),ixl(9),pl(9),ipxl(9),
     + sg(16),tg(16),wg(16),q(5)
      pi = datan(1.0d0)*4
c
      if(prt) write(iow,2000)
c.... read input
      if(ior.lt.0) write(*,3001)
3001  format('Input: x_min,x_max,y_min,y_max,[z_min,z_max] ) >',$)
c.... input coordinates
100   call dinput(td,6)
      xmin = td(1)
      xmax = td(2)
      ymin = td(3)
      ymax = td(4)
      zmin = td(5)
      zmax = td(6)
      val  = ddot(6,td,1,td,1)
      if(val.eq.0.d0) goto 200
c
c.... input other values
      if(ior.lt.0) write(*,3002)
3002  format('Input: idof,ietyp >',$)
      call dinput(td,2)
      mdf  = td(1)
      ietyp= td(2)
      if(ietyp.eq.0) ietyp=4
c.... set local node numbers for load with respect to ietyp
      call pzeroi(ipxl,9)
      if (ietyp.eq.4.or.ietyp.eq.-8.or.ietyp.eq.9.or.ietyp.eq.-27) then
        nel=4
        do i = 1,nel !      1, 2, 3, 4
          ipxl(i) = i
        end do
        if (ietyp.eq.9) then
          nel=9
          do i = 1,5 !      (1, 2, 3, 4,)  5, 6, 7, 8, 9
            ipxl(i+4) = i+4
          end do
        else if (ietyp.eq.-27) then
          nel=9
          do i = 1,5 !      (1, 2, 3, 4,) 13,14,15,16,17
            ipxl(i+4) = i+12
          end do
        end if
      else if(ietyp.eq.8.or.ietyp.eq.27) then
        nel=4
        do i = 1,nel !      5, 6, 7, 8
          ipxl(i) = i+4
        end do
        if(ietyp.eq.27) then
          nel=9
          do i = 1,5 !     (5, 6, 7, 8,) 18,19,20,21,22
            ipxl(i+4) = i+17
          end do
        end if
      else
      stop 'ALOA: ietyp not valid!' 
      end if
c.... set local node numbers for center of gravity with respect to ietyp
      npc = 8
      if (ietyp.eq.4.or.ietyp.eq.9) npc=4

c.... input load values
      call dinput(td,6)
      iqtyp = td(1)
      do i = 1,5
        q(i) = td(i+1)
      end do

c.... write data
      if(prt) then
        write(iow,2001) xmin,xmax,ymin,ymax,zmin,zmax,mdf,ietyp
        if(iqtyp.eq.1) write(iow,2002) q(1)
        if(iqtyp.eq.2) write(iow,2003) q(1),q(2),q(3)
      end if

c.... gauss points
      if(ietyp.eq.4.or.ietyp.eq.8.or.ietyp.eq.-8)   l = 2  ! linear    elements
      if(ietyp.eq.9.or.ietyp.eq.27.or.ietyp.eq.-27) l = 3  ! quadratic elements
      lint = l*l
      call pgauss(l,lint,sg,tg,wg)

c.... loop over all elements
      do ielem = 1,numel

c....   values for element/side
        call pzeroi(ixl,9)
        call pzero( xl,ndm*9)
        do iel = 1,nel
c....     nodes of element/side
          ino = ix(ipxl(iel),ielem)
          ixl(iel) = ino
c....     coordinates of element/side
          do idm = 1,ndm
            xl(idm,iel) = x(idm,ino)
          end do
        end do

c....   center of gravity (only corner nodes!!)
        xm =  0.d0
        ym =  0.d0
        zm =  0.d0
        do iel = 1,npc
          ino = ix(iel,ielem)
          xm  = xm + x(1,ino)
          ym  = ym + x(2,ino)
          if(ndm.eq.3) zm  = zm + x(3,ino)
        end do
        xm =  xm / npc
        ym =  ym / npc
        zm =  zm / npc

c....   load element if center is in area
        if(ndm.eq.2) then
          if(xmin.le.xm .and. xm .le. xmax .and.
     +       ymin.le.ym .and. ym .le. ymax) go to 101
          go to 102
        else if(ndm.eq.3) then
          if(xmin.le.xm .and. xm .le. xmax .and.
     +       ymin.le.ym .and. ym .le. ymax .and.
     +       zmin.le.zm .and. zm .le. zmax) go to 101
          go to 102
        else
          go to 102
        end if
101     continue
c....   calculate element load vector
        call pzero(pl,9)
c....   compute integrals of shape functions
        do l = 1,lint
          call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ixl,.false.)
          xsj = xsj*wg(l)
          do iel = 1,nel
            if(iqtyp.eq.1) qc =  q(1)
            if(iqtyp.eq.2) then
              xn = xl(1,iel)
              yn = xl(2,iel)
              qc = q(1)*cos(pi/(2.d0*q(2))*xn)*cos(pi/(2.d0*q(3))*yn)
            end if
            w = shp(3,iel)*xsj
            pl(iel) = pl(iel) + qc*w
          end do
        end do
c
c....   store element load into load vector f

        do iel = 1,nel
          ino = ix(ipxl(iel),ielem)
c....     from input file id=0/1
          if(ior.ge.0.and.id(mdf,ino).eq.0)
     +       f(mdf,ino)=f(mdf,ino)+pl(iel)
c....     interactive  id=-i/k via macro>mesh
          if(ior.lt.0.and.id(mdf,ino).gt.0)
     +       f(mdf,ino)=f(mdf,ino)+pl(iel)
        end do
102     continue
      end do
      goto 100
c.... print load vector
200   if(prt) then
        write(iow,2004)
        do i = 1,numnp
          if (f(mdf,i).ne.0.0d0) write(iow,2005) i,f(mdf,i)
        end do
      end if
c.... formats
2000  format(/5x,'L o a d s  o n  e l e m e n t s',/,5x,
     +'     xmin    ','     xmax    ','     ymin    ','     ymax    ',
     +'     zmin    ','     zmax    ',' dof ','Element typ')
2001  format(5x,6(1x,e12.5),2x,i1,2x,4x,i2)
2002  format(5x,'Load case 1: constant load q_0 ',e12.5)
2003  format(5x,'Load case 2: cos-load: q_0= ',e12.5,
     +          '  l_x= ',e12.5,'  l_y= ',e12.5)
2004  format(5x,'  Node',3x,' load value')
2005  format(5x,i6,3x,e12.5)
      end
c
      subroutine blkgen(ndm,ndf,nel,nel1,x,ix,id,prt,iblk)
c----------------------------------------------------------------------
c
c      Purpose: Generate a block of elements and nodes for mesh
c               descriptions.
c
c               ntyp = < 10     and 16: 2d-generation block
c               ntyp = > 10 without 16: 3d-generation block

c      Inputs:
c         x(ndm,*)   - Nodal coordinates of mesh
c         ix(nel1,*) - Element nodal connections of mesh
c         id(ndf,*)  - b.c. conditions
c         ndm        - Spatial dimension of mesh
c         ndf        - Number dof/node
c         nel        - Maximum number of nodes/element
c         nel1       - Dimension for 'ix' array
c         prt        - Print generated data if true
c         ibkl       - 0= ok, 1=error
c
c      Outputs:
c         x(ndm,*)   - Block of nodal coordinates
c         ix(nel1,*) - Block of elements
c
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE errchk
      USE iofile

      implicit double precision (a-h,o-z)
      logical prt
      character*6 xh,yyy*80
      dimension x(ndm,*),ix(nel1,*),id(ndf,*),
     *          xl(3,27),ixl(27),shp(3,9),td(10)
      data xh/' coord'/
      iblk=0
c.... block mesh generation routine
100   if(ior.lt.0) write(*,5000)
      call dinput(td,8)
      if(errck) go to 100
      nn     = td(1)
      nr     = td(2)
      ns     = td(3)
      ntyp   = td(8)
      if(ntyp.lt.10.or.ntyp.eq.16) then
c....   2-d generations
        nt     = 1
        ni     = td(4)
        ne     = td(5)
        ma     = td(6)
        nodinc = td(7)
      else
c....   3-d generations
        nt     = td(4)
        ni     = td(5)
        ne     = td(6)
        ma     = td(7)
        nodinc = 0
      end if
c.... reset to default values if necessary
      nodinc = max0(nodinc,0)
      nr     = max0(nr,1)
      ns     = max0(ns,1)
      nt     = max0(nt,1)
      ni     = max0(ni,1)
      ma     = max0(ma,1)
      if(prt) then
          if(ior.gt.0) then
            write(iow,2000)  nr,ns,nt,ni,ne,ma,nodinc,ntyp
            write(iow,2002) (i,xh,i=1,ndm)
          else
            write(  *,2000)  nr,ns,nt,ni,ne,ma,nodinc,ntyp
            write(  *,2002) (i,xh,i=1,ndm)
          end if
          if(ne.eq.0)
     +     call drawmess('warning * * * no elements are generated ',1,0)
      end if
      do 10 n = 1,27
       do 10 j = 1,3
           xl(j,n) = 0.0d0
           ixl(n) = 0
10    continue
      nm = 0
      do 20 n = 1,nn
21        if(ior.lt.0) write(*,5001)
          call dinput(td,4)
          if(errck) go to 21
          l = td(1)
          if(l.eq.0) l = n
          nm = max0(nm,l)
          ixl(l)        = l
          xl(1,l) = td(2)
          xl(2,l) = td(3)
          xl(3,l) = td(4)
          if(prt.and.ior.gt.0) write(iow,2001) l,(xl(i,l),i=1,ndm)
          if(prt.and.ior.lt.0) write(*  ,2001) l,(xl(i,l),i=1,ndm)
20    continue
c.... set generation increments of natural coordinates
      dr = 2.d0/nr
      ds = 2.d0/ns
c.... determine last element number to be generated
      if(ntyp.lt.10.or.ntyp.eq.16) then ! 2D
        if (ntyp.eq.0) then
          nf = ne + nr*ns - 1
        else if (ntyp.eq.7) then ! 6–node triangles
          nf = ne + (nr*ns)/2 - 1
        else if (ntyp.eq.8.or.ntyp.eq.9) then ! 8/9–node quadrilaterals
          nf = ne + (nr*ns)/4 - 1
        else if (ntyp.eq.16) then ! 16–node quadrilaterals
          nf = ne + (nr*ns)/9 - 1
        else
          nf = ne + 2*nr*ns - 1
        end if
cww     if(nf.gt.numel.and.ne.gt.0) go to 401
        if(nf.gt.numel.and.ne.gt.0) then
          write(yyy,2031) nf,numel
          call drawmess(yyy,1,0)
          iblk=1
        end if
c....   determine last node number to be generated
        nr = nr + 1
        ns = ns + 1
        if(ndm.eq.1) ns = 1
        ng = nr*ns + ni -1
cww     if(ng.gt.numnp) go to 400
        if(ng.gt.numnp) then
          write(yyy,2030) ng,numnp
          call drawmess(yyy,1,0)
          iblk=1
        end if
c....   form block of elements
        call sblk(nr,ns,xl,ixl,shp,x,ix,dr,ds,ni,ne,ndm,
     1           nel1,nodinc,ntyp,nm,ma,prt)
      else ! 3D
        dt = 2.d0/nt
        if(ntyp.eq.10) nf = ne + nr*ns*nt - 1    !  8–node hexahedron
        if(ntyp.eq.11) nf = ne + 6*nr*ns*nt - 1  !  4–node tetrahedron
        if(ntyp.eq.12) nf = ne + nr*ns*nt/8 - 1  ! 21–node hexahedron
        if(ntyp.eq.13) nf = ne + nr*ns*nt/8 - 1  ! 27–node hexahedron
        if(ntyp.eq.14) nf = ne + nr*ns*nt/8 - 1  ! 20–node hexahedron
        if(ntyp.eq.15) nf = ne + nr*ns*nt/4 - 1  ! 18–node hexahedron
        if(ntyp.eq.19) nf = ne + nr*ns*nt/27 - 1 ! 64 node hexahedron

c        if(ntyp.eq.15.and.nt.ne.1) then
c          write(  *,*) 'Error BLOC, tinc must be 1'
c          write(iow,*) 'Error BLOC, tinc must be 1'
c         stop
c        end if

cww     if(nf.gt.numel.and.ne.gt.0) go to 401
        if(nf.gt.numel.and.ne.gt.0) then
          write(yyy,2031) nf,numel
          call drawmess(yyy,1,0)
          iblk=1
        end if
        nr = nr + 1
        ns = ns + 1
        nt = nt + 1
        ng = nr*ns*nt + ni -1
cww     if(ng.gt.numnp) go to 400
        if(ng.gt.numnp) then
          write(yyy,2030) ng,numnp
          call drawmess(yyy,1,0)
          iblk=1
        end if
        call vblk(nr,ns,nt,xl,x,ixl,ix,id,dr,ds,dt,
     1            ni,ne,ndm,ndf,nel1,ma,ntyp,prt)
      end if
c.... print lists if wanted
      if(prt.and.ne.gt.0) then
c.... print element lists
      nf=min(nf,numel)  ! ww el. only up to numel
      do 502 n = ne,nf,50
cww     if(ior.gt.0) write(iow,2005) o,head,(i,i=1,nel)
cww     if(ior.lt.0) write(  *,2005) o,head,(i,i=1,nel)
        if(ior.gt.0) write(iow,2005)        (i,i=1,nel)
        if(ior.lt.0) write(  *,2005)        (i,i=1,nel)
        j = min0(nf,n+49)
        do 501 i = n,j
          if(ior.gt.0) write(iow,2006) i,ma,(ix(k,i),k=1,nel)
          if(ior.lt.0) write(  *,2006) i,ma,(ix(k,i),k=1,nel)
501     continue
502   continue
      end if
      return
c.... error messages
cww400   if(ior.gt.0) write(iow,2030) ng,numnp
cww      if(ior.lt.0) write(  *,2030) ng,numnp
c400   write(yyy,2030) ng,numnp
c      call drawmess(yyy,1,0)
c      iblk=1
c      return
cww401   if(ior.gt.0) write(iow,2031) nf,numel
cww      if(ior.lt.0) write(  *,2031) nf,numel
c401   write(yyy,2031) nf,numel
c      call drawmess(yyy,1,0)
c      iblk=1
c      return
cww2000  format(a1,19a4,a3//'   n o d e   g e n e r a t i o n s'//
2000  format(/'   n o d e   g e n e r a t i o n s'/
     1   10x,'number of r-increments    ',i5/
     2   10x,'number of s-increments    ',i5/
     3   10x,'number of t-increments    ',i5/
     4   10x,'first node number         ',i5/
     5   10x,'first element number      ',i5/
     6   10x,'element material number   ',i5/
     7   10x,'node line increment       ',i5/
     8   10x,'block type (0-14)         ',i5)
2001  format(i9,1p3e12.3)
2002  format(5x,'node',3(i6,a6))
cww2005  format(a1,19a4,a3//'   e l e m e n t   c o n n e c t i o n s'//
2005  format(/'   e l e m e n t   c o n n e c t i o n s'/
     1   '    elmt    matl',9(i3,' node')/(16x,9(i3,' node')))
2006  format(11i8/(16x,9i8))
cww2030  format(' **error** insufficient storage for nodes'/
cww     1        10x,'final node =',i5,5x,'numnp =',i5)
cww2031  format(' **error** insufficient storage for elements'/
cww     1        10x,'final element =',i5,5x,'numel =',i5)
2030  format(' insufficient storage: final node =',i5,' max. node =',i5)
2031  format(' insufficient storage: final element =',i5,
     +        ' max. element =',i5)
5000  format(' Input: nn,nr,ns,ni,ne,ma,nodinc,ntyp'/3x,'>',$)
5001  format(' Input: node, x-1, x-2, x-3'/3x,'>',$)
      end
c
      subroutine blkgenco(ndm,nel,nel1,x,ix,prt,iblk)
c----------------------------------------------------------------------
c
c      Purpose: Generate a block of elements and nodes for a zone of
c               shell and interface elements
c               mtyp=1: 4-node shell elmts and  8-node interface elmt
c               mtyp=2: 4-node shell elmts and 18-node interface elmt
c               mtyp=3: 4-node shell elmts and 16-node interface elmt
c
c      Inputs:
c         ndm        - Dimension of 'x' array
c         nel        - Maximum number of nodes/element
c         nel1       - Dimension for 'ix' array
c         prt        - Print generated data if true
c         ibkl       - 0= ok, 1=error
c
c      Outputs:
c         x(ndm,*)   - Block of nodal coordinates
c         ix(nel1,*) - Block of elements
c
c
c       Comments: very similar to blkgen for 4 node element
c
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      character*6 xh,yyy*80
      dimension x(ndm,*),ix(nel1,*),xl(3,9),ixl(9),shp(3,9),td(10)
      data xh/' coord'/
      iblk=0

c.... block mesh generation routine
100   if(ior.lt.0) write(*,5000)
      call dinput(td,9)
      if(errck) go to 100
      nn     = td(1)
      nr     = td(2)
      ns     = td(3)

c.... 2-d generations: 2*4-node shell elements + 8-node interface element
      nt     = 1
      ni     = td(4)
      ne     = td(5)
      ma1    = td(6)
      ma2    = td(7)
      ma3    = td(8)
      mtyp   = td(9)
      if(mtyp.eq.0) mtyp=1

c.... reset to default values if necessary
      nr     = max0(nr,1)
      ns     = max0(ns,1)
      nt     = max0(nt,1)
      ni     = max0(ni,1)

      if(prt) then
        if(ior.gt.0) then
          write(iow,2000) nr,ns,nt,ni,ne,ma1,ma2,ma3,mtyp
          write(iow,2002) (i,xh,i=1,ndm)
        else
          write(*  ,2000) nr,ns,nt,ni,ne,ma1,ma2,ma3,mtyp
          write(*  ,2002) (i,xh,i=1,ndm)
        end if
        if(ne.eq.0)
     +     call drawmess('warning * * * no elements are generated ',1,0)
      end if

      do 10 n = 1,9
       do 10 j = 1,3
           xl(j,n) = 0.0
           ixl(n) = 0
10    continue
      nm = 0

      do 20 n = 1,nn
21      if(ior.lt.0) write(*,5001)
        call dinput(td,4)
        if(errck) go to 21
        l = td(1)
        if(l.eq.0) l = n
        nm = max0(nm,l)
        ixl(l)  = l
        xl(1,l) = td(2)
        xl(2,l) = td(3)
        xl(3,l) = td(4)
        if(prt.and.ior.gt.0) write(iow,2001) l,(xl(i,l),i=1,ndm)
        if(prt.and.ior.lt.0) write(*  ,2001) l,(xl(i,l),i=1,ndm)
20    continue

c.... set generation increments of natural coordinates
      dr = 2.d0/nr
      ds = 2.d0/ns

c.... determine last element number to be generated
      nf = ne - 1 + 3*nr*ns   ! 3*times elements
      if(mtyp.eq.2) nf = ne - 1 + 2*nr*ns + nr*ns/4 ! 2*shell +1/4 interface
      if(mtyp.eq.3) nf = ne - 1 + 2*nr*ns + nr*ns/4 ! 2*shell +1/4 interface
      if(nf.gt.numel.and.ne.gt.0) then
        write(yyy,2031) nf,numel
        call drawmess(yyy,1,0)
        iblk=1
      end if

c.... determine last node number to be generated
      nr = nr + 1
      ns = ns + 1
      if(ndm.eq.1) ns = 1
      ng = 2*nr*ns + ni -1 ! 2*times nodes

      if(ng.gt.numnp) then
        write(yyy,2030) ng,numnp
        call drawmess(yyy,1,0)
        iblk=1
      end if

c.... form block of elements
      call sblkco(nr,ns,xl,ixl,shp,x,ix,dr,ds,ni,ne,ndm,
     1           nel1,nm,ma1,ma2,ma3,mtyp,prt)

c.... print lists if wanted
      if(prt.and.ne.gt.0) then
c....   print element lists
        nf=min(nf,numel)  ! ww el. only up to numel
        do 502 n = ne,nf,50
          if(ior.gt.0) write(iow,2005)        (i,i=1,nel)
          if(ior.lt.0) write(  *,2005)        (i,i=1,nel)
          j = min0(nf,n+49)
          do 501 i = n,j
            if(ior.gt.0) write(iow,2006) i,ix(nel1,i),(ix(k,i),k=1,nel)
            if(ior.lt.0) write(  *,2006) i,ix(nel1,i),(ix(k,i),k=1,nel)
501       continue
502     continue
      end if
      return
c
2000  format(/'   n o d e   g e n e r a t i o n s  for cohesive zones'/
     +   10x,'number of r-increments    ',i5/
     +   10x,'number of s-increments    ',i5/
     +   10x,'number of t-increments    ',i5/
     +   10x,'first node number         ',i5/
     +   10x,'first element number      ',i5/
     +   10x,'element material number 1 ',i5/
     +   10x,'element material number 2 ',i5/
     +   10x,'element material number 3 ',i5/
     +   10x,'mesh-typ:1=8,2=18,3=16node',i5)
2001  format(i9,1p3e12.3)
2002  format(5x,'node',3(i6,a6))
2005  format(/'   e l e m e n t   c o n n e c t i o n s cohesive zones'/
     1   '    elmt    matl',8(i3,' node')/(16x,8(i3,' node')))
2006  format(10i8/(16x,8i8))
2030  format(' insufficient storage: final node =',i5,' max. node =',i5)
2031  format(' insufficient storage: final element =',i5,
     +        ' max. element =',i5)
5000  format(' Input: nn,nr,ns,ni,ne,ma1,ma2,ma3,ityp'/3x,'>',$)
5001  format(' Input: node, x-1, x-2, x-3'/3x,'>',$)
      end
c
c
      subroutine blkgendx(ndm,nel,nel1,x,ix,prt,iblk)
c----------------------------------------------------------------------
c
c      Purpose: Generate a 2d block of elements and nodes for mesh
c                with delta x
c
c      Inputs:
c         ndm        - Dimension of 'x' array
c         nel        - Maximum number of nodes/element
c         nel1       - Dimension for 'ix' array
c         prt        - Print generated data if true
c         ibkl       - 0= ok, 1=error
c
c      Outputs:
c         x(ndm,*)   - Block of nodal coordinates
c         ix(nen1,*) - Block of elements
c
c       Comments: OPEN
c       mid side nodes, dr1,ds1 fest auf 100
c       nicht rechtecke?
c       mehr als 4 Knoten wegnehmen
c       material
c       Knotenzaehlung in pnums(feaps4)
c
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      character*6 xh,yyy*80
      dimension x(ndm,*),ix(nel1,*),xl(3,27),ixl(27),shp(3,9),td(16)
      dimension dr1(100),ds1(100),dmar(100),dmas(100)
      data xh/' coord'/
      iblk=0
c.... block mesh generation routine
100   if(ior.lt.0) write(*,5000)
      call dinput(td,8)
      if(errck) go to 100
      nn     = td(1)
      nr     = td(2)
      ns     = td(3)
      ntyp   = td(8)
c.... 2-d generations
      nt     = 1
      ni     = td(4)
      ne     = td(5)
      ma     = td(6)
      nodinc = td(7)
c.... reset to default values if necessary
      nodinc = max0(nodinc,0)
      nr     = max0(nr,1)
      ns     = max0(ns,1)
      nt     = max0(nt,1)
      ni     = max0(ni,1)
      ma     = max0(ma,1)
      if(nr.gt.100.or.ns.gt.100) stop 'only 100 increments allowed!'
      if(prt) then
          if(ior.gt.0) then
            write(iow,2000)        nr,ns,nt,ni,ne,ma,nodinc,ntyp
            write(iow,2002) (i,xh,i=1,ndm)
          else
            write(*,2000)        nr,ns,nt,ni,ne,ma,nodinc,ntyp
            write(*,2002) (i,xh,i=1,ndm)
          end if
          if(ne.eq.0)
     +     call drawmess('warning * * * no elements are generated ',1,0)
      end if
      do 10 n = 1,27
       do 10 j = 1,3
           xl(j,n) = 0.0
           ixl(n) = 0
10    continue
      nm = 0
      do 20 n = 1,nn
21        if(ior.lt.0) write(*,5001)
          call dinput(td,4)
          if(errck) go to 21
          l = td(1)
          if(l.eq.0) l = n
          nm = max0(nm,l)
          ixl(l)        = l
          xl(1,l) = td(2)
          xl(2,l) = td(3)
          xl(3,l) = td(4)
          if(prt.and.ior.gt.0) write(iow,2001) l,(xl(i,l),i=1,ndm)
          if(prt.and.ior.lt.0) write(*  ,2001) l,(xl(i,l),i=1,ndm)
20    continue
c.... read increments dx,dy for more than 16:
c     Fehler im gesamten FEAP: wenn genau 16 dann Leerzeile danach erforderlich
      il = min(nr,16)
      call dinput(dr1,il)
      if(nr.gt.16) then
        do 205 ii = 1,nr/16
          is = il+1
          il = min(is+15,nr)
203       call dinput(td,il-is+1)
          if(errck) go to 203
          do 204 k = 1,il-is+1
            dr1(k+is-1) = td(k)
204       continue
205     continue
      end if
c
      il = min(ns,16)
      call dinput(ds1,il)
      if(ns.gt.16) then
        do 215 ii = 1,ns/16
          is = il+1
          il = min(is+15,ns)
213       call dinput(td,il-is+1)
          if(errck) go to 213
          do 214 k = 1,il-is+1
            ds1(k+is-1) = td(k)
214       continue
215     continue
      end if
c.... read increments dmar,dmas for more than 16
      il = min(nr,16)
      call dinput(dmar,il)
      if(nr.gt.16) then
        do 225 ii = 1,nr/16
          is = il+1
          il = min(is+15,nr)
223       call dinput(td,il-is+1)
          if(errck) go to 223
          do 224 k = 1,il-is+1
            dmar(k+is-1) = td(k)
224       continue
225     continue
      end if
c
      il = min(ns,16)
      call dinput(dmas,il)
      if(ns.gt.16) then
        do 235 ii = 1,ns/16
          is = il+1
          il = min(is+15,ns)
233       call dinput(td,il-is+1)
          if(errck) go to 233
          do 234 k = 1,il-is+1
            dmas(k+is-1) = td(k)
234       continue
235     continue
      end if
c
      if(prt.and.ior.gt.0) write(iow,2007) (dr1(i),i=1,nr)
      if(prt.and.ior.lt.0) write(*  ,2007) (dr1(i),i=1,nr)
      if(prt.and.ior.gt.0) write(iow,2007) (ds1(j),j=1,ns)
      if(prt.and.ior.lt.0) write(*  ,2007) (ds1(j),j=1,ns)
      if(prt.and.ior.gt.0) write(iow,2008) (dmar(i),i=1,nr)
      if(prt.and.ior.lt.0) write(*  ,2008) (dmar(i),i=1,nr)
      if(prt.and.ior.gt.0) write(iow,2008) (dmas(i),i=1,ns)
      if(prt.and.ior.lt.0) write(*  ,2008) (dmas(i),i=1,ns)
2007  format('Increments:',16g10.5,/,(11x,16g10.5))
2008  format('Material  :',16f5.0,/,(11x,16f5.0))
c.... length dlx,dly
      dlx = 0.d0
      do i = 1,nr
        dlx = dlx + dr1(i)
      end do
      dly = 0.d0
      do j = 1,ns
        dly = dly + ds1(j)
      end do
c...  increments dr,ds
      do i = 1,nr
        dr1(i) = 2.d0/dlx*dr1(i)
      end do
      do j = 1,ns
        ds1(j) = 2.d0/dly*ds1(j)
      end do
c
c.... set generation increments of natural coordinates
cww   dr = 2.d0/nr
cww   ds = 2.d0/ns
c.... determine last element number to be generated
      nf = ne + nr*ns - 1
      if(nf.gt.numel.and.ne.gt.0) then
        write(yyy,2031) nf,numel
        call drawmess(yyy,1,0)
        iblk=1
      end if
c.... determine last node number to be generated
      nr = nr + 1
      ns = ns + 1
      if(ndm.eq.1) ns = 1
      ng = nr*ns + ni -1
      if(ng.gt.numnp) then
        write(yyy,2030) ng,numnp
        call drawmess(yyy,1,0)
        iblk=1
       end if
c....  form block of elements
       call sblkdx(nr,ns,xl,ixl,shp,x,ix,dr,ds,ni,ne,ndm,
     1        nel1,nodinc,ntyp,nm,ma,prt,dr1,ds1,dmar,dmas)
c.... print lists if wanted
      if(prt.and.ne.gt.0) then
c.... print element lists
      nf=min(nf,numel)  ! ww el. only up to numel
        do 502 n = ne,nf,50
          if(ior.gt.0) write(iow,2005)        (i,i=1,nel)
          if(ior.lt.0) write(  *,2005)        (i,i=1,nel)
          j = min0(nf,n+49)
          do 501 i = n,j
            if(ior.gt.0) write(iow,2006) i,ma,(ix(k,i),k=1,nel)
            if(ior.lt.0) write(  *,2006) i,ma,(ix(k,i),k=1,nel)
501       continue
502     continue
      end if
      return
2000  format(/'   n o d e   g e n e r a t i o n s'/
     1   10x,'number of r-increments    ',i5/
     2   10x,'number of s-increments    ',i5/
     3   10x,'number of t-increments    ',i5/
     4   10x,'first node number         ',i5/
     5   10x,'first element number      ',i5/
     6   10x,'element material number   ',i5/
     7   10x,'node line increment       ',i5/
     8   10x,'block type (0-14)         ',i5)
2001  format(i9,1p3e12.3)
2002  format(5x,'node',3(i6,a6))
2005  format(/'   e l e m e n t   c o n n e c t i o n s'/
     1   '    elmt    matl',8(i3,' node')/(16x,8(i3,' node')))
2006  format(10i8/(16x,8i8))
2030  format(' insufficient storage: final node =',i5,' max. node =',i5)
2031  format(' insufficient storage: final element =',i5,
     +        ' max. element =',i5)
5000  format(' Input: nn,nr,ns,ni,ne,ma,nodinc,ntyp'/3x,'>',$)
5001  format(' Input: node, x-1, x-2, x-3'/3x,'>',$)
      end
c
      subroutine blktem(ndm,t,prt,iblk)
c----------------------------------------------------------------------
c
c      Purpose: Generate a 2D/3D-block of temperatures
c
c      Inputs:
c         ndm        - Dimension of 'x' array
c         prt        - Print generated data if true
c         ibkl       - 0= ok, 1=error
c
c      Outputs:
c         t(*)       - Nodal temperature values
c
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      character yyy*80
      dimension t(*),tl(9),ixl(9),td(7)
      iblk = 0
100   if(ior.lt.0) write(*,5000)
      call dinput(td,7)
      if(errck) go to 100
      nn     = td(1)
      nr     = td(2)
      ns     = td(3)
      nt     = td(4)
      ni     = td(5)
      nodinc = td(6)
c.... 2-d generations
      if(ndm.lt.3) then
          nt     = 1
      end if
c.... reset to default values if necessary
      nodinc = max0(nodinc,0)
      nr     = max0(nr,1)
      ns     = max0(ns,1)
      nt     = max0(nt,1)
      ni     = max0(ni,1)
      if(prt) then
          if(ior.gt.0) then
cww         write(iow,2000) o,head,nr,ns,nt,ni,nodinc
            write(iow,2000)        nr,ns,nt,ni,nodinc
            write(iow,2002)
          else
cww         write(*,2000) o,head,nr,ns,nt,ni,nodinc
            write(*,2000)        nr,ns,nt,ni,nodinc
            write(*,2002)
          end if
      end if
      do 10 n = 1,9
          tl(n) = 0.0
          ixl(n) = 0
10    continue
      nm = 0
      do 20 n = 1,nn
21        if(ior.lt.0) write(*,5001)
          call dinput(td,2)
          if(errck) go to 21
          l = td(1)
          if(l.eq.0) l = n
          nm = max0(nm,l)
          ixl(l)        = l
          tl(l) = td(2)
          if(prt.and.ior.gt.0) write(iow,2001) l,tl(l)
          if(prt.and.ior.lt.0) write(*  ,2001) l,tl(l)
20    continue
c.... set generation increments of natural coordinates
      dr = 2.d0/nr
      ds = 2.d0/ns
      dt = 2.d0/nt
      nr = nr + 1
      ns = ns + 1
c.... determine last node number to be generated
      if(ndm.lt.3) then
          if(ndm.eq.1) ns = 1
          ng = nr*ns + ni -1
      else
          nt = nt + 1
          ng = nr*ns*nt + ni -1
      end if
cww   if(ng.gt.numnp) go to 400
      if(ng.gt.numnp) then
        write(yyy,2030) ng,numnp
        call drawmess(yyy,1,0)
        iblk=1
      end if
c.... form block of temperatures
      call tblk(nr,ns,nt,tl,t,ixl,dr,ds,dt,ni,ndm,nodinc,nm,prt)
      return
c.... error messages
cww400   if(ior.gt.0) write(iow,2030) ng,numnp
cww      if(ior.lt.0) write(  *,2030) ng,numnp
      return
cww2000  format(a1,19a4,a3//'   t e m p   g e n e r a t i o n s'//
2000  format(/'   t e m p   g e n e r a t i o n s'/
     1   10x,'number of r-increments    ',i5/
     2   10x,'number of s-increments    ',i5/
     3   10x,'number of t-increments    ',i5/
     4   10x,'first node number         ',i5/
     5   10x,'node line increment       ',i5/1x)
2001  format(i9,1p3e12.3)
2002  format(5x,'node',6x,' Temp.')
cww2030  format(' **error** insufficient storage for temperatures'/
cww     1        10x,'final node =',i5,5x,'numnp =',i5)
2030  format(' insuff. storage temp: final node =',i5,' max. node =',i5)
5000  format(' Input: nn,nr,ns,nt,ni,nodinc'/3x,'>',$)
5001  format(' Input: node, Temp'/3x,'>',$)
      end
c
      subroutine tblk(nr,ns,nt,tl,t,ixl,dr,ds,dt,ni,ndm,nodinc,nm,prt)
c----------------------------------------------------------------------
c      Purpose:  Generate a block of temperatures

c      Inputs:
c         nr        - Number elements in 1-local coordinate dir.
c         ns        - Number elements in 2-local coordinate dir.
c         nt        - Number elements in 3-local coordinate dir.
c         tl(*)     - Block nodal temperature array
c         ixl(*)    - Block nodal connection list
c         dr        - 1-local coordinate increment
c         ds        - 2-local coordinate increment
c         dt        - 3-local coordinate increment
c         ni        - Initial node number for block
c         ndm       - Spatial dimension of mesh
c         nodinc    - Increment array for block
c         nm        - Number master nodes on block
c         prt       - Output generated data if true
c
c      Outputs:
c         t(*)      - Nodal temperatures for block
c
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt,phd
      dimension ss(3),xl(2,9),tl(*),t(*),ixl(9),shp2(3,9),shp3(4,8)
      data xl/-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0,1.d0,
     1         0.d0,-1.d0,1.d0, 0.d0,0.d0,1.d0,-1.d0,0.d0,0.d0,0.d0/
c.... check that all corners of hexahedron are defined
      if(ndm.eq.3) then
          do 10 k = 1,nm
            if(ixl(k).ne.k) go to 900
10        continue
      end if
      n = ni
      mct = 0
      ss(3) = -1.0
      do 300 k = 1,nt
          ss(2) = -1.0
        do 200 j = 1,ns
            ss(1) = -1.0
            do 100 i = 1,nr
c....       compute shape functions and coordinates for each point
              t(n) = 0.0
              if(ndm.lt.3) then
                call shape(ss(1),ss(2),xl,shp2,xsj,2,nm,ixl,.true.)
                do 40 l = 1,nm
                      ll   = ixl(l)
                      t(n) = t(n) + shp2(3,ll)*tl(ll)
40              continue
              else if(ndm.eq.3) then
                call shp3d(ss,xsj,shp3,tl,1)
                do 50 l = 1,8
                      t(n) = t(n) + shp3(4,l)*tl(l)
50              continue
              end if
c....       output the point
              if(prt) then
                mct = mct + 1
                phd = mod(mct,50).eq.1
cww             if(phd) write(iow,2000) o,head,l
                if(phd) write(iow,2000) l
                write(iow,2001) n,t(n)
                if(ior.lt.0) then
cww                   if(phd) write(*,2000) o,head,l
                      if(phd) write(*,2000) l
                      write(*,2001) n,t(n)
                end if
              end if
              n = n + 1
            if(n.gt.numnp) goto 301
            ss(1) = ss(1) + dr
100         continue
            n = n + nodinc
          if(n.gt.numnp) goto 301
            ss(2) = ss(2) + ds
200       continue
          ss(3) = ss(3) + dt
300     continue
301   continue
      return
c.... error
900   write(iow,3000) k
      if(ior.lt.0) then
        write(*,3000) k
        return
      end if
cww   stop
      return
c.... formats
cww2000  format(a1,19a4,a3//'  n o d a l   t e m p e r a t u r e s'//
2000  format(/'  n o d a l   t e m p e r a t u r e s'/
     1    6x,'node',i7,' Temp.')
2001  format(i10,1p1e13.4)
3000  format(' **ERROR** Block node',i3,' is undefined')
      end
c
      subroutine blktri(ndm,nel,nel1,x,ix,prt,iblk)
c----------------------------------------------------------------------
c
c      Purpose:  Generate triangular block of triangular elements in 2-d

c      Inputs:
c         ndm        - Dimension of 'x' array
c         nel        - Maximum number of nodes/element
c         nel1       - Dimension for 'ix' array
c         prt        - Print generated data if true
c         ibkl       - 0= ok, 1=error
c
c      Outputs:
c         x(ndm,*)   - Block of nodal coordinates
c         ix(nel1,*) - Block of triangular elements
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE errchk
      USE iofile
      logical prt
      character*6 xh,yyy*80
      integer ndm,nel,nel1,ix,i,j,l,k,
     1        nn,nr,nod1,nuel1,ma,
     2        ixl,
     3        n,nm,nf,ng
      real*8 x,xl,td
      dimension x(ndm,*),ix(nel1,*),xl(3,6),ixl(6),td(5)
      data xh/' coord'/
      iblk=0
100   if(ior.lt.0) write(*,5000)
      call dinput(td,5)
      if(errck) go to 100
      nn     = td(1)
      nr     = td(2)
      nod1   = td(3)
      nuel1  = td(4)
      ma     = td(5)
c
c.... reset to default values if necessary
      nn     = max(nn,1)
      nr     = max(nr,1)
      nod1   = max(nod1,1)
      nuel1  = max(nuel1,1)
      ma     = max(ma,1)
      if(prt) then
          if(ior.gt.0) then
cww         write(iow,2000) o,head,nr,nod1,nuel1,ma
            write(iow,2000)        nr,nod1,nuel1,ma
            write(iow,2002) (i,xh,i=1,ndm)
          else
cww         write(*,2000) o,head,nr,nod1,nuel1,ma
            write(*,2000)        nr,nod1,nuel1,ma
            write(*,2002) (i,xh,i=1,ndm)
          end if
          if(nuel1.eq.0)
     +    call drawmess('warning * * * no elements are generated ',1,0)
      end if
c.... initialize arrays
      call pzero(xl,18)
      call pzeroi(ixl,6)
c
      nm = 0
      do 20 n = 1,nn
21      if(ior.lt.0) write(*,5001)
        call dinput(td,4)
        if(errck) go to 21
        l       = td(1)
        nm      = max(nm,l)
c
        ixl(l)  = l
        xl(1,l) = td(2)
        xl(2,l) = td(3)
        xl(3,l) = td(4)
        if(prt.and.ior.gt.0) write(iow,2001) l,(xl(i,l),i=1,ndm)
        if(prt.and.ior.lt.0) write(*,2001) l,(xl(i,l),i=1,ndm)
20    continue
c.... determine last element to be generated
        nf = nuel1 + nr*nr - 1
cww     if(nf.gt.numel.and.nuel1.gt.0) go to 401
        if(nf.gt.numel.and.nuel1.gt.0) then
        write(yyy,2031) nf,numel
        call drawmess(yyy,1,0)
        iblk=1
      end if
c
c.... determine last node number to be generated
        ng = nod1 + (nr+1)*(nr+2)/2 - 1
cww     if(ng.gt.numnp) go to 400
        if(ng.gt.numnp) then
         write(yyy,2030) ng,numnp
         call drawmess(yyy,1,0)
         iblk=1
       end if
c.... form block of elements
        call trblk(nr,xl,ixl,x,ix,ndm,nod1,nuel1,
     1            nel1,nm,ma,prt)
c
c.... print lists if wanted
      if(prt.and.nf.ge.nuel1) then
c.... print element lists
        do 502 n = nuel1,nf,50
cww       if(ior.gt.0) write(iow,2005) o,head,(i,i=1,nel)
cww       if(ior.lt.0) write(  *,2005) o,head,(i,i=1,nel)
          if(ior.gt.0) write(iow,2005)        (i,i=1,nel)
          if(ior.lt.0) write(  *,2005)        (i,i=1,nel)
          j = min0(nf,n+49)
          do 501 i = n,j
            if(ior.gt.0) write(iow,2006) i,ma,(ix(k,i),k=1,nel)
            if(ior.lt.0) write(  *,2006) i,ma,(ix(k,i),k=1,nel)
501       continue
502     continue
      end if
      return
c.... error messages
cww400   if(ior.gt.0) write(iow,2030) ng,numnp
cww      if(ior.lt.0) write(  *,2030) ng,numnp
cww      return
cww401   if(ior.gt.0) write(iow,2031) nf,numel
cww      if(ior.lt.0) write(  *,2031) nf,numel
cww      return
cww2000  format(a1,19a4,a3//'   n o d e   g e n e r a t i o n s'//
2000  format(/'   n o d e   g e n e r a t i o n s'/
     1   10x,'number of increments      ',i5/
     2   10x,'first node number         ',i5/
     3   10x,'first element number      ',i5/
     4   10x,'element material number   ',i5)
2001  format(i9,1p3e12.3)
2002  format(5x,'node',3(i6,a6))
cww2005  format(a1,19a4,a3//'   e l e m e n t   c o n n e c t i o n s'//
2005  format(/'   e l e m e n t   c o n n e c t i o n s'/
     1   '    elmt    matl',8(i3,' node')/(16x,8(i3,' node')))
2006  format(10i8/(16x,8i8))
cww2030  format(' **error** insufficient storage for nodes'/
cww     1        10x,'final node =',i5,5x,'numnp =',i5)
cww2031  format(' **error** insufficient storage for elements'/
cww     1        10x,'final element =',i5,5x,'numel =',i5)
2030  format(' insufficient storage: final node =',i5,' max. node =',i5)
2031  format(' insufficient storage: final element =',i5,
     +        ' max. element =',i5)
5000  format(' Input: nn,nr,nod1,nuel1,ma'/3x,'>',$)
5001  format(' Input: node, x-1, x-2, x-3'/3x,'>',$)
      end
c
      subroutine cfunc(shp,xl,ixl,ndm,x)
c--------------------------------------------------------------------
c      Purpose: Compute coordinates for point defined by local arrays,
c               node generation for blockgen
c
c      Inputs:
c         shp(3,*)  - Shape function array
c         xl(ndm,*) - Array of element coordinates
c         ixl(*)    - Element node numbers
c         ndm       - Spatial dimension of mesh
c
c      Outputs:
c         x(ndm)    - Coordinates of point

c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension shp(3,*),xl(3,*),ixl(*),x(*)
      do 55 l = 1,ndm
        x(l) = 0.0
        do 50 k = 1,9
          m = ixl(k)
          if(m.gt.0) x(l) = x(l) + shp(3,m)*xl(l,m)
50      continue
55    continue
      return
      end
c
      subroutine chkblk(yy,n0,nt)
c----------------------------------------------------------------------
c
c      Purpose: Add zero entry to character array which has blank
c               field.

c      Inputs:
c         y(*) - array to check
c         n0   - Field width of data
c         nt   - Size of array to check

c      Outputs:
c         y(*) - Blank fields have zero (0) added to field 'n0'
c
c----------------------------------------------------------------------
      character*1 y(75)
      character*75 yy
      do 200 i=1,75
200   y(i) = yy(i:i)
c.... patch for the plot device type if necessary
      if(y(1).eq.'t'.and.y(2).eq.'e'.and.y(3).eq.'k'.and.
     1   y(4).eq.'t') then
        if(y(16).eq.' ') y(16) = y(27)
        if(y(17).eq.' ') y(17) = y(28)
        if(y(18).eq.' ') y(18) = y(29)
        if(y(19).eq.' ') y(19) = y(30)
      end if
c.... add a character if y(nt) is blank
      do 100 n = n0,nt,n0
        if(y(n).eq.' ') y(n) = '0'
100   continue
      return
      end
c
      subroutine ckbrk8 ( n, ix, xl, ndm, nel, shp )
c----------------------------------------------------------------------
c      Purpose: Check 8-node hexahedron for bad data.
c               Write message to file on errors located.

c      Inputs:
c         n         - Number of element being checked
c         ix(*)     - List of nodes for element
c         xl(ndm,*) - Coordinate array
c         ndm       - Spatial dimension of mesh
c         nel       - Number of nodes on element

c      Outputs:
c         None

c      Scratch:
c         shp(*)    - Array to store shape functions
c----------------------------------------------------------------------
c....  Declare variable types
      USE iofile
      integer ndm, nel,  i, l, n, ineg
      real*8  detj
c....  Declare array types
      integer ix(*), ic(16)
      real*8  rst(3,8), xl(ndm,*), shp(4,8)
c..... Intrinsics
      intrinsic abs
c
      data rst/-1.d0,-1.d0,-1.d0,   1.d0,-1.d0,-1.d0,
     1          1.d0, 1.d0,-1.d0,  -1.d0, 1.d0,-1.d0,
     2         -1.d0,-1.d0, 1.d0,   1.d0,-1.d0, 1.d0,
     3          1.d0, 1.d0, 1.d0,  -1.d0, 1.d0, 1.d0/
c....  Check the element for input errors
      ineg = 0
      do 100 l = 1,nel
        if(xl(1,l).eq. -999.0d0 .and. ix(l).ne.0) then
          ic(ineg+1) = l
          ic(ineg+2) = abs(ix(l))
          ineg = ineg + 2
        end if
100   continue
c....  Node numbering errors
      if(ineg.gt.0) then
        write(iow,2000) n,(ic(i),i=1,ineg)
        if(ior.lt.0) write(*,2000) n,(ic(i),i=1,ineg)
      else
c....  Compute the jacobian at each corner of the element
        do 110 l = 1,nel
          call bjac3d ( rst(1,l) , xl, ndm, shp, detj )
          if(detj.le.0.0d0) then
            ic(ineg+1) = l
            ic(ineg+2) = abs(ix(l))
            ineg = ineg + 2
          end if
  110   continue
        if(ineg.gt.0) then
          write(iow,2001) n,(ic(i),i=1,ineg)
          if(ior.lt.0) write(*,2001) n,(ic(i),i=1,ineg)
        end if
      end if
c
2000  format(' >Element',i4,' coordinates not input for nodes:'/
     1      ('                Local =',i3,' Global =',i4))
2001  format(' >Element',i4,' has negative jacobian at nodes:'/
     1      ('                Local =',i3,' Global =',i4))
      end
c
      subroutine bjac3d ( rst , xl, ndm, shp, detj )
c----------------------------------------------------------------------
c
c      Purpose: Compute jacobian determinant and shape functions
c               with natural coord. derivatives for an 8-node hexahedron.

c      Inputs:
c         rst(3)    - Natural coordinate location
c         xl(ndm,*) - Array of element coordinates
c         ndm       - Space dimension of mesh

c      Outputs:
c         shp(4,8)  - Shape functions and derivatives w/r natural coords.
c         detj      - Determinant of jacobian determinant
c
c----------------------------------------------------------------------
c....  Declare variable types
      integer ndm, i, j, k
      real*8  detj, xii, eti, zti
c....  Declare array types
      real*8  rst(3), xl(ndm,*), shp(4,8), xs(3,3)
      real*8  xi(8), eta(8), zta(8)
c
      data xi /-0.5d0, 0.5d0, 0.5d0,-0.5d0,-0.5d0, 0.5d0, 0.5d0,-0.5d0/
      data eta/-0.5d0,-0.5d0, 0.5d0, 0.5d0,-0.5d0,-0.5d0, 0.5d0, 0.5d0/
      data zta/-0.5d0,-0.5d0,-0.5d0,-0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0/
c
c....  Compute the shape functions and their derivatives
      do 100 i = 1,8
        xii = 0.5 +  xi(i)*rst(1)
        eti = 0.5 + eta(i)*rst(2)
        zti = 0.5 + zta(i)*rst(3)
        shp(1,i) =  xi(i)*eti*zti
        shp(2,i) = eta(i)*xii*zti
        shp(3,i) = zta(i)*xii*eti
        shp(4,i) =  xii*eti*zti
100   continue
c
c....  Compute the jacobian matrix
      call pzero(xs,9)
      do 210 i = 1,3
      do 210 j = 1,3
        do 200 k = 1,8
          xs(i,j) = xs(i,j) + xl(i,k)*shp(j,k)
200     continue
210   continue
c
c....  Compute the jacobian determinant
      detj = xs(1,1)*(xs(2,2)*xs(3,3) - xs(2,3)*xs(3,2))
     1     + xs(1,2)*(xs(2,3)*xs(3,1) - xs(2,1)*xs(3,3))
     1     + xs(1,3)*(xs(2,1)*xs(3,2) - xs(2,2)*xs(3,1))
      end
c
      subroutine ckisop(ix,xl,shp,ndm)
c----------------------------------------------------------------------
c
c      Purpose: Check isoparametric elements for data input errors
c
c      Inputs:
c         ix(*)     - List of nodes connected to element
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c
c      Outputs:
c         None

c      Scratch:
c         shp(*)    - Storage for shape functions
c
c----------------------------------------------------------------------
      USE eldata
      USE iofile
      implicit real*8 (a-h,o-z)
c.... check isoparametric elements
      integer xn(9),yn(9),ic(18),ix(*)
      real*8  shp(3,*),xl(ndm,*)
      data xn/-1,1,1,-1,0,1,0,1,0/, yn/-1,-1,1,1,-1,0,1,0,0/
c.... check the element for input errors
      ineg = 0
      do 100 l = 1,nel
        if(xl(1,l).eq. -999.0d0 .and. ix(l).ne.0) then
          ic(ineg+1) = l
          ic(ineg+2) = abs(ix(l))
          ineg = ineg + 2
        end if
100   continue
      if(ineg.gt.0) then
        write(iow,2000) n,(ic(i),i=1,ineg)
        if(ior.lt.0) write(*,2000) n,(ic(i),i=1,ineg)
      else
        do 110 l = 1,nel
          ss = xn(l)
          tt = yn(l)
          call  shape (ss,tt,xl,shp,xsj,ndm,nel,ix,.false.)
          if(xsj.le.0.0d0) then
            ic(ineg+1) = l
            ic(ineg+2) = abs(ix(l))
            ineg = ineg + 2
          end if
  110   continue
        if(ineg.gt.0) then
          write(iow,2001) n,(ic(i),i=1,ineg)
          if(ior.lt.0) write(*,2001) n,(ic(i),i=1,ineg)
        end if
      end if
      return
2000  format(' >Element',i4,' coordinates not input for nodes:'/
     1      ('                Local =',i3,' Global =',i4))
2001  format(' >Element',i4,' has negative jacobian at nodes:'/
     1      ('                Local =',i3,' Global =',i4))
      end
c
      subroutine copyhis(h1,h2,h1e,h2e,nh,isw)
c----------------------------------------------------------------------
c
c      Purpose:
c
c      Inputs/Outputs:
c
c      isw = 1 copy history variables from global arrays h1e,h2e
c              to element working arrays h1,h2
c
c      isw = 2 copy history variables from global arrays h3e
c              to element working arrays h3  (in h1e,h1)
c
c      isw = 3 copy history variables from element working arrays h1,h2
c              to  global arrays h1e,h2e
c
c      isw = 4 copy history variables from element working arrays h3
c              to  global arrays h3e (in h1,h1e
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension h1e(nh),h2e(nh),h1(nh),h2(nh)
      if(isw.eq.1) then
c....   copy element values into working array of actual element
        do n = 1,nh
          h1(n)=h1e(n)
          h2(n)=h2e(n)
        end do
      elseif(isw.eq.2) then
c....   copy element values into working array of actual element
        do n = 1,nh
          h1(n)=h1e(n)
        end do
      elseif(isw.eq.3) then
c....   copy working array of actual element into element array
c       reset element values into working array of actual element
c       (= save in case need to recompute tangent)
        do n = 1,nh
          temp   = h1e(n)
          h1e(n) = h1(n)
          h1(n)  = temp

          temp   = h2e(n)
          h2e(n) = h2(n)
          h2(n)  = temp
        end do
      elseif(isw.eq.4) then
c....   copy working array of actual element into element array
c       reset element values into working array of actual element
c       (= save in case need to recompute tangent)
        do n = 1,nh
          temp   = h1e(n)
          h1e(n) = h1(n)
          h1(n)  = temp
        end do
      end if
      return
      end
c
      subroutine dasbly(s,p,ld,jp,ns,alfl,aufl,bfl,b,al,au,ad)
c----------------------------------------------------------------------
c
c      Purpose: Assemble symmetric/unsymmetric arrays for 'DATRI'

c      Inputs:
c         s(ns,ns) - Element array to assemble
c         p(ns)    - Element vector to assemble
c         ld(ns)   - Local to Global equation numbers for assemble
c         jp(*)    - Pointer array for upper/lower parts of A array.
c         ns       - Size of element arrays
c         aufl     - If true, assemble A array
c         alfl     - If true, array A is unsym
c         bfl      - If true, assemble B vector

c      Outputs:
c         b(*)     - Assembled right hand side B vector
c         al(*)    - Assembled lower part of A array
c         au(*)    - Assembled upper part of A array
c         ad(*)    - Assembled diagonal part of A array
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical alfl,aufl,bfl
      dimension al(*),au(*),ad(*),b(*),s(ns,ns),p(ns),ld(ns),jp(*)
c.... loop through the rows to perform the assembly
      do 200 i = 1,ns
        ii = ld(i)
        if(ii.gt.0) then
          if(aufl) then
c....       loop through the columns to perform the assembly
            do 100 j = 1,ns
              if(ld(j).eq.ii) then
                !$OMP ATOMIC
                ad(ii) = ad(ii) + s(i,j)
              else if(ld(j).gt.ii) then
                jc = ld(j)
                jj = ii + jp(jc) - jc + 1
                !$OMP ATOMIC
                au(jj) = au(jj) + s(i,j)
                if(alfl) then
                  !$OMP ATOMIC
                  al(jj) = al(jj) + s(j,i)
                end if
              end if
100         continue
          end if
          if(bfl) then
            !$OMP ATOMIC
            b(ii)  = b(ii)  + p(i)
          end if
        end if
200   continue
      return
      end
c
      subroutine dasol (al,au,ad,dr,jdiag,neq,aengy)
c-----------------------------------------------------------------------
c      Purpose: Solution of the problem Ax=b
c
c      Inputs:
c         al(*)    - Factored lower triangular terms
c         au(*)    - Factored upper triangular terms
c         ad(*)    - Factored diagonal terms
c         dr(*)    - right hand side vector b
c         jdiag(*) - Pointer to row/column ends in 'al' and 'au'.
c         neq      - Number of equations
c
c      Outputs:
c         dr(*)    - Solution vector x
c         aengy    - Energy residual
c
c      Comments:
c      istyp = 0  => standard solver FEAP
c            = 1  => SM       solver (CSR sparse matrix without minimum degree)
c            = 2  => SM       solver (CSR sparse matrix with    minimum degree)
c            = 3  => SuperLU  solver CSR
c            = 4  => Pardiso  solver CSR
c            = 5  => PBCG     solver CSR
c            = 6  => PGMRES   solver CSR
c            = 7  => PGMRES2  solver CSR
c            = 8  => Pardiso  solver CSR iterative
c
c-----------------------------------------------------------------------

      USE iofile
      USE soltyp
      implicit double precision (a-h,o-z)
      dimension dr(*),ad(*),al(*),au(*),jdiag(*)

      if ( istyp.eq.0 ) then
c-----------------------------------------------------------------------
c...    standard profil solver
          call dasol1 ( al,au,ad,dr,jdiag,neq,aengy )

      else if ( istyp.eq.1 .or. istyp.eq.2 ) then
c-----------------------------------------------------------------------
c...    SM solver
        call dasol2(dr,jdiag,neq,aengy)

      else if ( istyp.eq.3 ) then
c-----------------------------------------------------------------------
c...    SuperLU solver
        call dasol3 ( ad,dr,jdiag,neq,aengy )

      else if ( istyp.eq.4.or.istyp.eq.8 ) then
c-----------------------------------------------------------------------
c...    Pardiso solver
        call dasol4 ( ad,dr,jdiag,neq,aengy )

      else if ( istyp.eq.5 ) then
c-----------------------------------------------------------------------
c...    PBCG solver
        call dasol5 ( ad,dr,jdiag,neq,aengy )

      else if ( istyp.eq.6 ) then
c-----------------------------------------------------------------------
c...    PGMRES solver
        call dasol6 ( ad,dr,jdiag,neq,aengy )

      else if ( istyp.eq.7 ) then
c-----------------------------------------------------------------------
c...    PGMRES2 solver
        call dasol7 ( ad,dr,jdiag,neq,aengy )

      end if

      return
      end
c
      subroutine dasol1 (al,au,ad,b,jp,neq, energy)
c----------------------------------------------------------------------
c
c      Purpose: Solution of the problem Ax=b  for the standard solver
c
c      Inputs:
c         al(*)    - Factored lower triangular terms
c         au(*)    - Factored upper triangular terms
c         ad(*)    - Factored diagonal terms
c         b(*)     - right hand side vector b
c         jp(*)    - Pointer to row/column ends in 'al' and 'au'.
c         neq      - Number of equations
c
c      Outputs:
c         b(*)     - Solution vector x
c         energy   - Energy residual
c
c       Comment:
c         Coefficient matrix must be decomposed into its triangular
c         factors using datri before using dasol.
c
c----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension al(*),au(*),ad(*),b(*),jp(*)

c.... find the first non-zero entry in the right hand side
      do 100 is = 1,neq
          if(b(is).ne.0.0d0) go to 200
100   continue
      call drawmess(
     +      ' ***DASOL warning*** zero right-hand-side vector',1,0)
      return
c.... reduce the right hand side
200   continue
      call perform(0,neq,1,26)
      do 300 j = is+1,neq
        jr = jp(j-1)
        jh = jp(j) - jr
        if(jh.gt.0) then
          b(j) = b(j) - ddot(jh,al(jr+1),1,b(j-jh),1)
        end if
        call perform(j,neq,2,26)
300   continue
      call perform(j,neq,3,26)
c.... multiply by inverse of diagonal elements
      energy = 0.0d0
      do 400 j = is,neq
          bd = b(j)
          b(j) = b(j)*ad(j)
          energy = energy + bd*b(j)
400   continue
c.... backsubstitution
      call perform(0,neq,1,27)
      do 500 j = neq,2,-1
          jr = jp(j-1)
          jh = jp(j) - jr
          if(jh.gt.0) then
c           call colred(au(jr+1),b(j),jh, b(j-jh))
            call daxpty(jh,-b(j),au(jr+1),b(j-jh))
          end if
      call perform(neq-j,neq,2,27)
500   continue
      call perform(neq,neq,3,27)
      return
      end
c

      subroutine datest(au,jh,daval)
c----------------------------------------------------------------------
c
c      Purpose: Check if equations are singular when zero diagonal
c               exists

c      Inputs:
c         au(*) - Column of A array
c         jh    - Height of column

c      Outputs:
c         daval - Sum of absolute values of column.
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension au(jh)
      daval = 0.0d0
      do 100 j = 1,jh
         daval = daval + dabs(au(j))
100   continue
      return
      end
c
      subroutine datri ( al,au,ad,jd,neq,cfr )
c-----------------------------------------------------------------------
c
c      Purpose: Triangular decomposition of matrix stored in profile or
c               sparse form.
c
c      Inputs:
c         al(*)  - Unfactored lower triangular terms
c         au(*)  - Unfactored upper triangular terms
c         ad(*)  - Unfactored diagonal terms
c         jd(*)  - Pointer to row/column ends in 'al' and 'au'.
c         neq    - Number of symmetric equations
c         cfr    - .false:=sym .true.=unsym
c
c      Outputs:
c         al(*)  - Factored lower triangular terms
c         au(*)  - Factored upper triangular terms
c         ad(*)  - Factored diagonal terms
c
c      Comments:
c      istyp = 0  => standard solver FEAP
c            = 1  => SM       solver (CSR sparse matrix without minimum degree)
c            = 2  => SM       solver (CSR sparse matrix with    minimum degree)
c            = 3  => SuperLU  solver CSR
c            = 4  => Pardiso  solver CSR
c            = 5  => PBCG     solver CSR
c            = 6  => PGMRES   solver CSR
c            = 7  => PGMRES2  solver CSR
c            = 8  => Pardiso  solver CSR iterative
c            0 9  => simpex algorithm
c
c-----------------------------------------------------------------------
      USE conv
      USE dii
      USE iofile
      USE soltyp
      implicit double precision (a-h,o-z)
      logical cfr
      dimension ad(*),al(*),au(*),jd(*)

c---------------------------------
c.... reset values for plotting neg.Dii
      call pzeroi(ndii,150)
      call pzeroi(ii,3)
c--------------------------------

      if ( istyp.eq.0 ) then
c-----------------------------------------------------------------------
c....   standard profil solver
        call datri1 (al,au,ad,jd,neq,cfr)

      else if ( istyp.eq.1 .or. istyp.eq.2 ) then
c-----------------------------------------------------------------------
c....   SM solver
        call datri2(ad,jd)

      else if ( istyp.eq.3 ) then
c-----------------------------------------------------------------------
c....   SuperLU solver
        call datri3(ad,jd,neq,nneg)

      else if ( istyp.eq.4.or.istyp.eq.8 ) then
c-----------------------------------------------------------------------
c....   Pardiso solver
        call datri4(ad,jd,neq,nneg)

      else if ( istyp.eq.5 ) then
c-----------------------------------------------------------------------
c....   PBCG solver
c       no triangular decomposition !!

      else if ( istyp.eq.6 ) then
c-----------------------------------------------------------------------
c....   PGMRES solver
c       no triangular decomposition !!

      else if ( istyp.eq.7 ) then
c-----------------------------------------------------------------------
c....   PGMRES2 solver
c       no triangular decomposition !!

      end if

c.... print number of negative diagonals (if there are any)
      if(nneg.gt.0) then
                     write(iow,1000) nneg
        if(ior.lt.0) write(*  ,1000) nneg
      end if
      return
cww1000  format(5x,'factorization message : number of neg. diag.=',i5)
1000  format(5x,'###### factorization warning:',i5,' neg. diagonal el.')
      end
c
      subroutine datri1(al,au,ad,jp,neq,flg)
c----------------------------------------------------------------------
c
c      Purpose: Triangular decomposition of matrix stored in profile
c               form - standard feap solver.
c
c      Inputs:
c         al(*)  - Unfactored lower triangular terms
c         au(*)  - Unfactored upper triangular terms
c         ad(*)  - Unfactored diagonal terms
c         jp(*)  - Pointer to row/column ends in 'al' and 'au'.
c         neq    - Number of symmetric equations
c         flg    - .false:=sym .true.=unsym
c
c      Outputs:
c         al(*)  - Factored lower triangular terms
c         au(*)  - Factored upper triangular terms
c         ad(*)  - Factored diagonal terms, reciprocals
c         nneg   - number of negative diagonal entries
c
c----------------------------------------------------------------------
      USE arcext
      USE conv
      USE fdata
      USE iofile
      implicit double precision (a-h,o-z)
      logical flg
      dimension al(*),au(*),ad(*),jp(*)
c.... n.b.  tol should be set to approximate half-word precision.
      data zero,one/0.0d0,1.0d0/, tol/0.5d-07/
c.... set initial values for conditioning check
      dimx = zero
      dimn = zero
      nneg = 0
      do 50 j = 1,neq
        dimn = max(dimn,abs(ad(j)))
50    continue
      dfig = zero
c.... loop through the columns to perform the triangular decomposition
cc---------------------------------
cwwc.... reset values for plotting neg.Dii
cww      call pzeroi(ndii,150)
cww      call pzeroi(ii,3)
cc--------------------------------
C.... find max. difference of off-diagonal elements (only for utang)
      if(flg) then
       diffmax = 0.d0
       do j=1,jp(neq)
         diff= abs(al(j)-au(j))
         diffmax = max(diff,diffmax)
       end do
cww    if (ior .le. 0) write(*,  1002) diffmax
cww                    write(iow,1002) diffmax
cww1002 format(' max. difference of off-diagonals :', 2x,e18.10)
      end if
C
      jd = 1
c.... show state
      call perform(0,neq,1,25)
      do 200 j = 1,neq
        jr = jd + 1
        jd = jp(j)
        jh = jd - jr
        if(jh.gt.0) then
          is = j - jh
          ie = j - 1
c....     if diagonal is zero compute a norm for singularity test
          if(ad(j).eq.zero) call datest(au(jr),jh,daval)
          do 100 i = is,ie
            jr = jr + 1
            id = jp(i)
            ih = min0(id-jp(i-1),i-is+1)
            if(ih.gt.0) then
              jrh = jr - ih
              idh = id - ih + 1
              au(jr) = au(jr) - ddot(ih,au(jrh),1,al(idh),1)
              if(flg) al(jr) = al(jr) - ddot(ih,al(jrh),1,au(idh),1)
            end if
100       continue
        end if
c....   reduce the diagonal
        dd = ad(j)
        if(jh.ge.0) then
          jr = jd - jh
          jrh = j - jh - 1
          call dredu(al(jr),au(jr),ad(jrh),jh+1,flg  ,ad(j))
c....     extended system - add exeps*e(k)*e(k)
          if(kflg) then
            if(j.eq.kex) then
              if(exeps.eq.zero) exeps = dimn
              ad(j) = ad(j) + exeps
              write(iow,3000) kex,exeps
              if(ior.lt.0) write(*,3000) kex,exeps
            end if
          end if
c....     count number of negative diagonals
          if(ad(j).lt.zero) nneg=nneg+1
c....     check for possible errors and print warnings
          if(ior.gt.0) then
            if(dabs(ad(j)).lt.tol*dabs(dd))     call prtdii1(j,3)
            if(dd.lt.zero.and.ad(j).gt.zero)    call prtdii1(j,1)
            if(dd.gt.zero.and.ad(j).lt.zero)    call prtdii1(j,1)
            if(ad(j) .eq.  zero)                call prtdii1(j,2)
c....       complete rank test for a zero diagonal case
            if(dd.eq.zero.and.jh.gt.0) then
              if(dabs(ad(j)).lt.tol*daval)   write(iow,2003) j
            end if
          else
            if(dabs(ad(j)).lt.tol*dabs(dd))     call prtdii1(j,3)
            if(dd.lt.zero.and.ad(j).gt.zero)    call prtdii1(j,1)
            if(dd.gt.zero.and.ad(j).lt.zero)    call prtdii1(j,1)
            if(ad(j) .eq.  zero)                call prtdii1(j,2)
c....       complete rank test for a zero diagonal case
            if(dd.eq.zero.and.jh.gt.0) then
              if(dabs(ad(j)).lt.tol*daval)   write(    *,2003) j
            end if
          end if
        end if
c....   store reciprocal of diagonal
        if(abs(ad(j)).gt.1.d-300) then
          dimx = max(dimx,abs(ad(j)))
          dimn = min(dimn,abs(ad(j)))
          dfig = max(dfig,dd/abs(ad(j)))
          ad(j) = one/ad(j)
        end if
        call perform(j,neq,2,25)
200   continue
      call perform(neq,neq,3,25)
c.... print conditioning information
      dd = zero
      if(abs(dimn).gt.1.d-300) dd = dimx/dimn
cww>  cases added for mixed elements
      if(dfig.gt.0.d0) then
        ifig = dlog10(dfig) + 0.6d0
      else
        ifig = 999
      end if
cww<
      if(pfr)              write(iow,2004) dimx,dimn,dd,ifig
      if(pfr.and.ior.lt.0) write(  *,2004) dimx,dimn,dd,ifig

      return
c.... formats
c2000  format('**datri warning 1** loss of at least 7 digits in',
c     1 ' reducing diagonal of eq.',i5)
c2001  format('**datri warning 2** sign of diagonal changed when',
c     1 ' reducing eq.',i5)
c2002  format(' **datri warning 3** reduced diagonal is zero for',
c     1 ' eq.',i5)
2003  format('**datri warning 4** rank failure for zero unreduced',
     1 ' diagonal in eq.',i5)
2004  format(' Condition check: D-max',e11.4,'; D-min',e11.4,
     1 '; Ratio',e11.4/' Maximum no. diagonal digits lost:',i3)
cww2005  format('   FEAP solve equations: '/)
cww2006  format('  current eq. ',i6,'  max. no. of eqs. ',i6)
3000  format(' Singularity in factoring matrix by datri.'/,
     1       '   kex =',i5,', eps =',g12.5)
      end
c
      subroutine dcheck(x,vd,nt,error)
c----------------------------------------------------------------------
c      Purpose: Internal input of values from string data

c      Inputs:
c         x(*)  - Character array from inputs
c         nt    - Length of character string

c      Outputs:
c         vd(*) - Numerical values from string inputs
c         error - True of error occurs during inputs
c----------------------------------------------------------------------
      USE conval
      USE iofile
      implicit double precision (a-h,o-z)
      logical error
      character*1 x(nt)
      character*255 y
      character*1 yy(255)
      dimension ivd(2,16),vd(*)
c
c.... expression substitution using predefined constants
c
      n0 = 15
      nn = 1
      do 50 i = 1,nt
        if(x(i).eq.',') nn = nn + 1
50    continue
      call acheck(x,yy,n0,nt,255)
      call pcharr(yy,ivd,n0,nn*n0)
      error = .false.
      do 12 i=1,255
12    y(i:i) = yy(i)
      read(y,'(17f15.0)',err=200) (vd(i),i=1,nn)
      do 60 i = 1,nn
        if(ivd(1,i).gt.0) then
          vd(i) = vvv(ivd(1,i),ivd(2,i))
        else if(ivd(1,i).lt.0) then
          vd(i) = www(-ivd(1,i))
        end if
 60   continue

      return
c.... error
 200    write(iow,2000) nn,y(1:nn*15)
      if(ior.lt.0) then
        write(  *,2000) nn,y(1:nn*15)
      end if
      call errclr('DCHECK')
      error = .true.
      return
2000  format(/' *ERROR* attempting to input',i3,' value(s).',
     &        '  Input is:'/a)
      end
c
      subroutine detkt(ad,neq,ifdet)
c----------------------------------------------------------------------
c
c      Purpose: calculate determinant of stiffness matrix
c
c      Inputs:
c         ad(neq)  - Factored diagonal terms of A-array
c         neq(*)   - Length of A array
c         ifdet    - 1: set (det K) to (det K_0) see ARCLEN

c      Outputs:
c         detc     - Scaled determinant of K(_fact)
c         nneg     - Number of negative diagonals  in K(_fact)
c
c      Comments:
c      Only correct, if ad are terms of factored(!) A-array!
c
c     W.Wagner IBNM/UH 6/92
c----------------------------------------------------------------------
      USE arcl
      USE iofile
      USE soltyp
      implicit double precision(a-h,o-z)
      dimension ad(*)

c.... calculate determinant and neg. dii
      det1 = 0.d0
      nneg = 0

c.... calculate value det1 = LN(det K)

      if(istyp.eq.0) then                   ! standard solver 0

        do i = 1,neq
          adi = 1.d0/ad(i)       ! Diags are stored 1/adi
          if(abs(adi).gt.1.d-30) then
            if(adi.lt.0.d0) then
              xx = -adi
              nneg = nneg + 1
            else
              xx = adi
            end if
            det1 = det1 + log(xx)  ! add on LN-basis, sign later
          end if

c          write(iow,*) 'dof', i,' D_ii',xx,' LN D_ii',log(xx)

        end do

      else if(istyp.eq.1.or.istyp.eq.2) then  ! SM
        call detkt2(neq,nneg,det1)

      else if(istyp.eq.3) then   ! SuperLU
        call detkt3(neq,nneg,det1)

      else if(istyp.eq.4.or.istyp.eq.8) then   ! Pardiso

        call detkt4(nneg,det1)

      else if(istyp.eq.5) then   ! PCG

c       write(*,*) 'determinant for solver not implemented'

      else if(istyp.eq.6) then   ! PGMRES

c       write(*,*) 'determinant for solver not implemented'

      else if(istyp.eq.7) then   ! PGMRES2

c       write(*,*) 'determinant for solver not implemented'

      end if

c.... for all solvers:

c.... initial value of LN(detK)
      if(det0.eq.0.d0 .or. ifdet.eq.1) det0=det1

c.... scaling
      det1 = det1 - det0

c.... overflow det K_T
cww      if(det1.gt.230) write(*,*)
cww     + '#### Warning: Overflow Determinant possible LN DET K =', det1

c.... back from detKq=e^LN(detKq)
      det1 = exp(det1)

c.... add sign by number of negative diagonal elements
      detc=det1*(-1)**nneg

      return
      end
c
      subroutine dinput(d,nn)
c----------------------------------------------------------------------
c
c      Purpose: Data input subprogram for real values
c
c      Inputs:
c         nn    - Number of data items to input
c
c      Outputs:
c         d(nn) - Values for nn items input
c
c      Comments:
c         possible characters in line: 200
c         max number of items:          16
c         length(print format) of item: 15
c         for input with variables
c           max. length in line:        75
c
c----------------------------------------------------------------------
      USE errchk
      USE iofile
      USE iosave
      implicit double precision(a-h,o-z)
cww      character xxx(80)*1
      character xxx(201)*1
      dimension d(nn)
c.... check on number of items
      errck = .false.
      if(nn.gt.16) then
          if(ior.lt.0) then
            write(*,2000) nn
            errck = .true.
            return
          else
            write(iow,2000) nn
cww         stop
          return
          end if
      end if
      do 50 n = 1,nn
50      d(n) = 0.0d0
51    if(ior.gt.0) then
          read (ior,1000,err=901,end=902) xxx
      else
          read (*  ,1000,err=901,end=902) xxx
      end if
c.... check length of string
      if(xxx(201).ne.' ') then
         write(*,*) 'Warning: string too long'
         write(*,*) xxx
         stop
      end if
52    if(lsave) write(lfile,1000) xxx
c.... if no characters in list return
      if (xxx(1).eq.char(0))  return
cww      do 60 nl = 80,1,-1
      do 60 nl = 201,1,-1
          if(xxx(nl).ne.' ') go to 70
60    continue
c.... if blank characters in list return
      if (xxx(1).eq.' ')  return
70    no = 1
      nv = 1
c.... skip leading blank characters
      do 80 n = 1,nl
          if(xxx(n).ne.' ') go to 90
80    continue
c.... format separators are blanks or commas
90    if(xxx(n).eq.' '.or.xxx(n).eq.',') then
          if(n.gt.no) call setval(xxx(no),n-no,d(nv))
92        n = n + 1
          if(n.lt.nl.and.xxx(n).eq.' ') go to 92
          no = n
          nv = nv + 1
      else
          n  = n + 1
      end if
      if(n.le.nl.and.nv.le.nn) go to 90
c.... fill in last value if needed
      if(n.gt.no.and.nv.le.nn) call setval(xxx(no),n-no,d(nv))
      return
c.... read error encoutered
901   call  errclr ('DINPUT')
      goto  51
c.... eof encountered
902   call  endclr ('DINPUT',xxx(1))
      goto  52
c
cww1000  format(80a1)
1000  format(201a1)
2000  format(' ** ERROR ** too many items requested, limit = 16')
      end
      subroutine dinput2(d,nn)
c----------------------------------------------------------------------
c
c      Purpose: Data input subprogram for real values=dinput with ior=32
c
c      Inputs:
c         nn    - Number of data items to input
c
c      Outputs:
c         d(nn) - Values for nn items input
c
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      character xxx(201)*1
      dimension d(nn)
      ior = 32
c.... check on number of items
      if(nn.gt.16) then
        write(*,2000) nn
        return
      end if
      do 50 n = 1,nn
50      d(n) = 0.0d0
      read (ior,1000) xxx
c.... check length of string
      if(xxx(201).ne.' ') then
         write(16,*) 'Warning: string too long'
         write(16,*) xxx
         stop
      end if
c.... if no characters in list return
      if (xxx(1).eq.char(0))  return
      do 60 nl = 201,1,-1
          if(xxx(nl).ne.' ') go to 70
60    continue
c.... if blank characters in list return
      if (xxx(1).eq.' ')  return
70    no = 1
      nv = 1
c.... skip leading blank characters
      do 80 n = 1,nl
          if(xxx(n).ne.' ') go to 90
80    continue
c.... format separators are blanks or commas
90    if(xxx(n).eq.' '.or.xxx(n).eq.',') then
        if(n.gt.no) call setval2(xxx(no),n-no,d(nv))
92      n = n + 1
        if(n.lt.nl.and.xxx(n).eq.' ') go to 92
        no = n
        nv = nv + 1
      else
        n  = n + 1
      end if
      if(n.le.nl.and.nv.le.nn) go to 90
c.... fill in last value if needed
      if(n.gt.no.and.nv.le.nn) call setval2(xxx(no),n-no,d(nv))
      return
c
1000  format(201a1)
2000  format(' ** ERROR ** too many items requested, limit = 16')
      end
