c$Id:$
      subroutine peforc(x,f,ndtyp,ndm,ndf,numnp,type,prt,prth,name)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change default gap by dividing by sqrt(numnp)    05/12/2006
c       2. Allow up to 30 dof/node on inputs                26/06/2007
c       3. Move prtitl inside prt check                     24/09/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set force or displacement values based on specified
c               edge coordinates.

c      Inputs:
c         x(ndm,*)   - Nodal coordinates of mesh
c         ndtyp(*)   - Nodal type (>=0 exists; <0 tied)
c         ndm        - Spatial dimension of mesh
c         ndf        - Number dof/node
c         numnp      - Number of nodes in mesh
c         prt        - Output generated results if true
c         prth       - Output title/header data if true
c         name       - File name for inputs

c      Outputs:
c         f(ndf,*)   - Nodal force or displacement values
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   prt,prth, errck, pinput, tinput, vinput, pcomp, gapfl
      character name*5, text*15, type*4
      integer   ndm,ndf,numnp,i,j,n,nn, ndtyp(*)
      real*8    pdiff, x0,dx, gap, x(ndm,numnp),f(ndf,numnp),td(32)

      save

c     Read input of boundary edge forces

      gap   = 1.d-3/sqrt(dble(max(1,numnp)))
      gapfl = .false.
100   if(ior.lt.0) then
        write(*,3001) ndf
        call pprint('   >')
      endif
      if(ndf.le.14) then
        nn = ndf + 1
      else
        nn = 15
      endif
      errck = tinput(text,1,td(2),nn)
      if(errck) go to 100
      if(pcomp(text,'gap',3)) then
        gap   = td(2)
        gapfl = .true.
        if(prt) then
          write(iow,2002) gap
          if(ior.lt.0) then
            write(*,2002) gap
          endif
        endif
        go to 100
      else
        errck = vinput(text,15,td(1),1)
        i  = abs(nint(td(1)))
      endif
      if(ndf.gt.14) then
        nn    = ndf - 14
        errck = pinput(td(17),nn)
      endif
      if(i.le.0.or.i.gt.ndm) go to 4
      x0 = td(2)

c     Locate nodes on which to apply loading

      if(gapfl) then
        dx = gap
      else
        dx  = pdiff(x,i,ndm,numnp)*1.d-3
      endif
      do n = 1,numnp
        if(ndtyp(n).ge.0.and.abs(x(i,n)-x0).le.dx) then
          if(pcomp(type,'add',3)) then
            do j = 1,ndf
              f(j,n) = f(j,n) + td(j+2)
            end do ! j
          else
            do j = 1,ndf
              f(j,n) = td(j+2)
            end do ! j
          endif
        endif
      end do ! n
      go to 100

4     if(prt) call prtitl(prth)
      if(prt) write(iow,2000) (i,name,i=1,ndf)
      if(prt.and.ior.lt.0) write(*,2000) (i,name,i=1,ndf)
      do n = 1,numnp
        do i = 1,ndf
          if(f(i,n).ne.0.0d0) go to 400
        end do ! i
        go to 410
400     if(prt .and. ndtyp(n).ge.0) then
          write(iow,2001) n,(f(i,n),i=1,ndf)
          if(ior.lt.0) then
            write(*,2001) n,(f(i,n),i=1,ndf)
          endif
        endif
410     continue
      end do ! n

c     Formats

2000  format('  E d g e    N o d a l    V a l u e'/
     &       /4x,'Node',6(i6,'-',a5):/(8x,6(i6,'-',a5)))

2001  format(i8,1p,6e12.4:/(8x,1p,6e12.4))

2002  format(/'     Search gap = ',1p,1e12.4/)

3001  format(' Input: ndir,x(ndir),(fl(i),i=1,',i2,')')

      end
