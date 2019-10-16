c$Id:$
      subroutine umacr9(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Convert a mesh of 10-node tetrahedra into a mesh of
c                4-node tetrahedra.  Each quadratic tet is divided into
c                8 linear ones.

c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters

c      Outputs:
c         N.B.  Users are responsible for command actions.  See
c               programmers manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'iofile.h'
      include  'sdata.h'
      include  'umac1.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   pcomp
      character lct*15
      real*8    ctl(3)

      save

c     Set command word

      if(pcomp(uct,'mac9',4)) then      ! Usual    form
        uct = 'qtol'                    ! Specify 'name'
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation

c     Convert a quadratic order element mesh to a linear order one

      call pqdtolin(hr(np(27)),mr(np(31)),mr(np(33)),hr(np(43)))

      endif

      end

      subroutine pqdtolin(f,id,ix,x)

      implicit   none

      include   'bdata.h'
      include   'cdata.h'
      include   'chdata.h'
      include   'comfil.h'
      include   'iodata.h'
      include   'iofile.h'
      include   'pglob1.h'
      include   'refng.h'
      include   'sdata.h'

      logical    fflag,oflag
      integer    i,e,n, nel
      integer    id(ndf,numnp,2),ix(nen1,numel)
      integer    it(4,8)
      real*8     f(ndf,numnp,2),x(ndm,numnp)

c     Transfer 10-node tet to 8 4-node tets

      data       it / 1, 5,7,8, 2,6, 5,9, 3,7, 6,10, 4, 9,8,10,
     &                5,10,9,6, 5,7,10,6, 5,9,10, 8, 5,10,7, 8/

c     Output the new mesh

      nel = 4
      open(unit = ios, file = 'Itet4')
      write(ios,'(a)') 'NOCOunt'
      write(ios,'(20a4/6i8)') head,numnp,numel*8,nummat,ndm,ndf,nel

c     Check for active global parameters

      if(gtypfl.or.gdeffl.or.gomgfl.or.
     &   gtdofl.or.grayfl.or.groufl) then
        write(ios,'(/a)') 'GLOBal parameters'

c       Output 2-d solution type

        if(gtypfl) then
          if(g2type.eq.1) then
            write(ios,2008) '  PLANe STREss'
          elseif(g2type.eq.2) then
            write(ios,2008) '  PLANe STRAIN'
          elseif(g2type.eq.3) then
            write(ios,2008) '  AXISymmetric'
          elseif(g2type.eq.8) then
            write(ios,2008) '  AXISymmetric TORSion'
          endif
        endif

c       Output deformation type: Small or Finite

        if(gdeffl) then
          if(gdtype.eq.1) then
            write(ios,2008) '  SMALl deformation'
          elseif(gdtype.eq.-1) then
            write(ios,2008) '  FINIte deformation'
          endif
        endif

c       Output rotational velocity parameters

        if(gomgfl) then
            write(ios,2008) '  OMEGa radians    ',gomega(1)
            write(ios,2008) '  OMEGa COORdinate ',(gomex(n),n=1,ndm)
            write(ios,2008) '  OMEGa VECTor     ',(gomev(n),n=1,ndm)
        endif

c       Output thermal dof

        if(gtdofl) then
          write(iow,2008) '  TEMPerature DOF ',gtdof
        endif

c       Output rayleigh damping parameters
        if(grayfl) then
          write(iow,2008) '  RAYLeigh damping ',gray(1),gray(2)
        endif

c       Output ground/group parameters

        if(groufl) then
          write(iow,2008) '  GROUp  factors ',(gfac(n),n=1,ndf)
        endif
      endif

c     Output coordinate data

      write(ios,'(a)') 'NOPRint'
      write(ios,'(/a)') 'COORdinate ALL'
      do n = 1,numnp
        write(ios,'(i8,i3,1p,3e15.5)') n,0,(x(i,n),i=1,ndm)
      end do ! n

      write(ios,'(/a)') 'ELEMents ALL'
      nel = 0
      do n = 1,numel
        do e = 1,8
          nel = nel + 1
          write(ios,'(10i8)') nel,0,ix(nen1,n),(ix(it(i,e),n),i=1,4)
        end do ! e
      end do ! n

      oflag = .true.
      do n = 1,numnp
        fflag = .false.
        do i = 1,ndf
          if(f(i,n,1).ne.0.0d0) then
            fflag = .true.
          endif
        end do ! i
        if(oflag .and. fflag) then
          write(ios,'(/a)') 'FORCes'
          oflag = .false.
        endif
        if(fflag) then
          write(ios,'(i8,i3,1p,3e15.5)') n,0,(f(i,n,1),i=1,ndf)
        endif
      end do ! n

      oflag = .true.
      do n = 1,numnp
        fflag = .false.
        do i = 1,ndf
          if(f(i,n,2).ne.0.0d0) then
            fflag = .true.
          endif
        end do ! i
        if(oflag .and. fflag) then
          write(ios,'(/a)') 'DISPlacements'
          oflag = .false.
        endif
        if(fflag) then
          write(ios,'(i8,i3,1p,6e15.5)') n,0,(f(i,n,2),i=1,ndf)
        endif
      end do ! n

      oflag = .true.
      do n = 1,numnp
        fflag = .false.
        do i = 1,ndf
          if(id(i,n,2).ne.0) then
            fflag = .true.
          endif
        end do ! i
        if(oflag .and. fflag) then
          write(ios,'(/a)') 'BOUNdary codes'
          oflag = .false.
        endif
        if(fflag) then
          write(ios,'(i8,i3,6i8)') n,0,(id(i,n,2),i=1,ndf)
        endif
      end do ! n

c     Retrieve material data set entries from file "fmtl"
c     Output material/parameter lists

      open(unit=iwd,file=fmtl,status='old')
      fflag = .true.
      write(ios,'(a)') ' '
      do while(fflag)
        read(iwd,'(a)',end=200) xxx
        do e = 255,1,-1
          if(xxx(e:e).ne.' ') go to 100
        end do ! e
        e = 1
100     if(xxx(e:e).eq.char(13)) xxx(e:e) = ' '
        write(ios,'(a)') xxx(1:e)
      end do ! while
200   close(iwd)

      write(ios,'(/a)') 'END mesh'
      write(ios,'(/a)') 'INTEractive'
      write(ios,'(/a)') 'STOP'

      close(ios,status = 'keep')

c     Format

2008  format( a,1p,3e14.6)

      end
