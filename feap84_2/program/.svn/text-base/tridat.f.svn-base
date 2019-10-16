c$Id:$
      subroutine tridat(x,ix,numnp,numel,nen,ndm,nen1,intr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Output mesh for tri2d

c      Inputs:
c         x(ndm,*)  - Nodal coordinates for mesh
c         ix(nen1,*)- Element nodal connection list
c         numnp     - Number of nodes in mesh
c         numel     - Number of elements in mesh
c         nen       - Number of nodes on element
c         ndm       - Spatial dimension of mesh
c         nen1      - Dimension of ix array
c         intr      - Flag, inquire if true.

c      Outputs:
c         none      - Outputs go to file for use by tri2d
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iodata.h'
      include  'iofile.h'
      include  'prflag.h'
      include  'prstrs.h'
      include  'psdat4.h'
      include  'comblk.h'

      character text1*80, y*1
      logical   lfil,intr
      integer   i,j,numnp,numel,ndm,nen,nen1
      integer   npoig,neleg,nfn,nbcs,nreg,nmat,ndof,ncon,ix(nen1,*)
      real*8    x(ndm,*)

      save

c     If error analysis performed, generate new tri2d mesh

      if (trifl) then

c       Inquire if mesh is to be generated

        if(intr) then
          call pprint('/ *Generate mesh for tri2d (y or n)? -> ')
          read (*,'(a)') y
          if(y.ne.'y' .and. y.ne.'Y' .and. y.ne. ' ') return
        endif

c       Open file for tri2d mesh

        i     = 1
        fname      = 'Atria'
100     fname(5:5) = char(96+i)
        i     = i + 1
        inquire(file = fname , exist = lfil )
        if(lfil) then
          if(i.gt.26) then
            write(*,*) 'Cannot create file: ',fname
            call plstop()
          else
            go to 100
          endif
        endif

        open(unit = iwd, file = fname, access = 'sequential',
     &       form = 'formatted', status = 'unknown')

c       Output title and parameters

        inquire(file = 'NAME', exist = lfil)
        if(lfil) then
          open(unit = ird, file = 'NAME', access = 'sequential',
     &         form = 'formatted', status = 'old')
          read (ird,'(a)') text1
          read (ird,   *) npoig,neleg,nfn,nbcs,nreg,nmat,ndof,ncon
          write(iwd,2005) text1,numnp,numel,nfn,nbcs,nreg,nmat,ndof,ncon
        else
          write(  *,3002)
          write(ilg,3002)
          call plstop()
        end if

c       Output nodal coordinates

        write(iwd,2000)
        do i = 1,numnp
          write(iwd,2001) i,x(1,i),x(2,i)
        end do ! i

c       Output element connections

        write(iwd,2002)
        do i = 1,numel
          write(iwd,2003) i,ix(nen1,i),(ix(j,i),j=1,nen)
        end do ! i

c       Output error data

        write(iwd,2004)
        do i = 1, numnp
          write(iwd,2001) i, hr(ner-1+i)
        end do ! i

c       Output segment and region data

        do i = 1, 5000
          read (ird,'(a)',end=200) text1
          write(iwd,'(a)') text1
        end do ! i

c       Close all files and update tri2d execution defaults

200     close(ird)
        close(iwd)

        open(unit = iwd, file = 'triname', access = 'sequential',
     &       form = 'formatted', status = 'unknown')

        text1(1:1)  = 'I'
        text1(2:10) = fname(1:9)

        write(iwd,2006) fname,text1(1:10)
        close(iwd)

        write(  *,3001) fname

      endif

c     Formats

2000  format(/'coor')
2001  format(i7,'  0',2e13.5)
2002  format(/'elem')
2003  format(i7,'  0',10i7)
2004  format(/'erro')
2005  format(a/8i7)
2006  format(2(a10,2x))

3001  format(/' *End of <TRI2D> mesh generation, File: ',a/)
3002  format(/'  *ERROR* TRIDAT: Could not find file: NAME.')

      end
