c$Id:$
      subroutine umacr8(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Convert mesh of 10-node tetrahedra into mesh of 4-node
c                tetrahedra

c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters

c      Outputs:
c               Feap mesh
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'iofile.h'
      include  'sdata.h'
      include  'umac1.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   pcomp, setvar, palloc
      character lct*15
      real*8    ctl(3)

      save

c     Set command word

      if(pcomp(uct,'mac8',4)) then      ! Usual    form
        uct = 'two2'                    ! Specify 'name'
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation

c     Convert a quadratic order element mesh to a linear order one

      setvar = palloc(111,'TEMP1',numnp+1,1) ! used to store linear nodes

      call ptwo2one(hr(np(27)),mr(np(31)),mr(np(33)),hr(np(43)),
     &              mr(np(111)))

      endif

      end

      subroutine ptwo2one(f,id,ix,x,ip)

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
      integer    i,e,n, nnp,nel,nel1
      integer    id(ndf,numnp,2),ix(nen1,numel), ip(0:numnp)
      real*8     f(ndf,numnp,2),x(ndm,numnp)

c     Determine lower order nodes

      do n = 0,numnp
        ip(n) = 0
      end do ! n

c     Loop through elements

      nel  = 0
      nel1 = 0
      do n = 1,numel
        do i = 1,nen
          if(ix(i,n).gt.0) then
            nel = i
          endif
        end do ! i

c       Check element type

        if(ndm.eq.2) then
          if(nel.eq.6) then                      ! quadratic triangle
            nel = 3
          elseif(nel.eq.8 .or. nel.eq.9) then    ! quadratic quadrilateral
            nel = 4
          endif
        elseif(ndm.eq.3) then
          if(nel.eq.10) then                     ! quadratic tetrahedron
            nel = 4
          elseif(nel.ge.20 .and. nel.le.27) then ! quadratic brick
            nel = 8
          endif
        endif
        nel1 = max(nel1,nel)

c       Mark ip array for linear nodes

        do i = 1,nel
          if(ix(i,n).gt.0) then
            ip(ix(i,n)) = 1
          endif
        end do ! i
      end do ! n

c     Convert ip to new node numbers

      nnp = 0
      do n = 1,numnp
        if(ip(n).gt.0) then
          nnp = nnp + 1
          ip(n) = nnp
        endif
      end do ! n

c     Output the new mesh

      open(unit = ios, file = 'Ilinear')
      write(ios,'(a)') 'NOCOunt'
      write(ios,'(20a4/6i8)') head,nnp,numel,nummat,ndm,ndf,nel1

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
        if(ip(n).gt.0) then
          write(ios,'(i8,i3,1p,3e15.5)') ip(n),0,(x(i,n),i=1,ndm)
        endif
      end do ! n

      write(ios,'(/a)') 'ELEMents ALL'
      do n = 1,numel
        write(ios,'(10i8)') n,0,ix(nen1,n),(ip(ix(i,n)),i=1,nel1)
      end do ! n

      oflag = .true.
      do n = 1,numnp
        if(ip(n).gt.0) then
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
            write(ios,'(i8,i3,1p,3e15.5)') ip(n),0,(f(i,n,1),i=1,ndf)
          endif
        endif
      end do ! n

      oflag = .true.
      do n = 1,numnp
        if(ip(n).gt.0) then
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
            write(ios,'(i8,i3,1p,6e15.5)') ip(n),0,(f(i,n,2),i=1,ndf)
          endif
        endif
      end do ! n

      oflag = .true.
      do n = 1,numnp
        if(ip(n).gt.0) then
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
            write(ios,'(i8,i3,6i8)') ip(n),0,(id(i,n,2),i=1,ndf)
          endif
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
