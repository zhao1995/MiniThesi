c$Id:$
      subroutine pdelfl()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add delete of mesh tie file 'IXfilename'         29/05/2009
c       2. Set length of file extension to 8                25/01/2011
c       3. Add 'body' for ftyp data to delete               08/06/2011
c          Add 'setext' to set correct file name
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Clean up old files: Erase temporary mesh input files

c      Inputs:
c         none

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'comfil.h'
      include   'compas.h'
      include   'cornum.h'
      include   'iodata.h'
      include   'part1.h'
      include   'part3.h'
      include   'part8.h'
      include   'setups.h'

      include   'pointer.h'
      include   'comblk.h'

      include   'p_int.h'

      character  fname*132,fext*8,ftyp(19)*8, flnk(4)*8
      logical    exst,isopen, flags(5)
      integer    i,m,mmx, mtest

      save

      data       fp   / 10*0 /

c     data       ftyp /'an0','bn0','co0','ds0','ep0','fr0',
c    &                 'gb0','in0','ld0','sl0','ud0','ww0',
c    &                 'yp0'/
      data       ftyp /'eang','ebou','edis','efor','epro','ebas',
     &                 'csur','cbou','cang','cdis','cfor','cpro',
     &                 'cdam','cmas','csti','cbas','curv','lfor',
     &                 'body'/
c     data       flnk /'eln','lnk','jnt','tem'/
      data       flnk /'elin','link','jnts','temp'/

c     Check if 'ios' open

      inquire(unit = ios, opened = exst)
      if(exst) then
        close(unit = ios)
      endif

c     Check for maximum number of files active

      mmx = max(nsurf,nbouf,ndisf,nforf,nangf)

c     Delete mesh generation files

      do m = 0,mmx

        mtest = m
        do i = 1,19
          fname = fsav
          call setext(ftyp(i), mtest, fext, .false.)
          call addext(fname,fext,128,8)
          inquire(file=fname,exist=exst,opened=isopen)
          if(exst) then
            if(.not.isopen) then
              open (unit = ios, file = fname, status = 'old')
            endif
            close(unit = ios, status = 'delete')
          endif
        end do ! i
      end do ! m

c     Delete mesh manipulation files

      do m = 1,4
        fname = fsav
        mtest = 0
        call setext(flnk(m), mtest, fext, .false.)
        call addext(fname,fext,128,8)
        inquire(file=fname,exist=exst,opened=isopen)
        if(exst) then
          if(.not.isopen) then
            open (unit = ios, file = fname, status = 'old')
          endif
          close(unit = ios, status = 'delete')
        endif
      end do ! m

c     Delete 'feaploop' mesh blocks

      fname = ' '   ! clear name
      fname = 'feaploop000.0'
      if(rank.gt.0) then
        if(rank.lt.10) then
          write(fname(11:11),'(i1)') rank
        elseif(rank.lt.100) then
          write(fname(10:11),'(i2)') rank
        else
          write(fname( 9:11),'(i3)') rank
        endif
      endif
      do m = 0,9
        write(fname(13:13),'(i1)') m
        inquire(file=fname,exist=exst,opened=isopen)
        if(exst) then
          if(.not.isopen) then
            open (unit = ios, file = fname, status = 'old')
          endif
          close(unit = ios, status = 'delete')
        endif
      end do ! m

c     Delete mesh tie file

      fname      = ' '  ! Clear name
      fname(1:2) = 'IX'
      i          = index(finp,' ')
      fname(3:i) = finp(2:i-1)
      inquire(file=fname,exist=exst,opened=isopen)
      if(exst) then
        if(.not.isopen) then
          open (unit = ios, file = fname, status = 'old')
        endif
        close(unit = ios, status = 'delete')
      endif

c     Delete 'Material' file

      inquire(file= fmtl ,exist=exst,opened = isopen)
      if(exst) then
        if(.not.isopen) then
          open (unit = ios, file =  fmtl , status = 'old')
        endif
        close(unit = ios, status = 'delete')
      endif

c     Delete out-of-core blocks (additional to look for contacts)

      do m = 1,min(999,maxbl+10)

c       Check upper

        fname = '      '
        if(m.lt.10) then
          write(fname,'(a6,i1)') 'Aupper',m
        elseif(m.lt.100) then
          write(fname,'(a6,i2)') 'Aupper',m
        else
          write(fname,'(a6,i3)') 'Aupper',m
        endif
        inquire(file=fname,exist=exst,opened=isopen)
        if(exst) then
          if(.not.isopen) then
            open (unit = ios, file = fname, status = 'old')
          endif
          close(unit = ios, status = 'delete')
        endif

c       Check lower

        fname = '      '
        if(m.lt.10) then
          write(fname,'(a6,i1)') 'Alower',m
        elseif(m.lt.100) then
          write(fname,'(a6,i2)') 'Alower',m
        else
          write(fname,'(a6,i3)') 'Alower',m
        endif
        inquire(file=fname,exist=exst,opened=isopen)
        if(exst) then
          if(.not.isopen) then
            open (unit = ios, file = fname, status = 'old')
          endif
          close(unit = ios, status = 'delete')
        endif

      end do ! m

c     Purge solver storage if necessary

      do i = 1,5
        if(.not.tflp(i)) then
          if(nsolver(i)) then
            fp(1) = nap (i)
            fp(2) = naup(i)
            fp(3) = nalp(i)
            fp(4) = np(20+i)
            call psolve(nittyp(i),hr(np(26)),fp,
     &                 .false.,.false.,.false.,.false.)
          else
            flags(1) = .false.
            flags(2) = .false.
            flags(3) = .false.
            flags(4) = .false.
            flags(5) = .true.
            call usolve(flags,hr(np(26)))
          endif
        endif
      end do ! i

      end
