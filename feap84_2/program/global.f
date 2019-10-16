c$Id:$
      subroutine global()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c     2. Add 'augm,<on,off>' option                         14/03/2007
c     3  Add global 'equa'tion option                       27/03/2008
c     4. Add warning if number global equations is zero     02/05/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set global solution parameter

c      Inputs:
c         none

c      Outputs:
c         Global parameters output through common blocks
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'iofile.h'
      include  'modcon.h'
      include  'pconstant.h'
      include  'pglob1.h'
      include  'sdata.h'
      include  'refng.h'

      logical   pcomp, tinput, errck, palloc, setval
      character type(2)*15, gqd(2)*5
      integer   i
      real*8    td(14), gvec

      save

      data      gqd / 'Gauss', 'Nodal' /

c     Input a record

      write(iow,2000)
      if(ior.lt.0) then
        write(*,2000)
      end if
100   if(ior.lt.0) then
        write(*,3000)
        call pprint('          >')
      endif
      errck = tinput(type,2,td,14)

c     Plane stress/strain option

      if(pcomp(type(1),'plan',4)) then

        gtypfl = .true.
        if(ndm.le.2) then

c         Plane stress option

          if(pcomp(type(2),'stre',4)) then

            g2type = 1
            write(iow,2001)
              if(ior.lt.0) then
              write(*,2001)
            end if

c         Plane strain option (default for 'plane')

          else

            g2type = 2
            write(iow,2002)
            if(ior.lt.0) then
              write(*,2002)
            end if

          endif

c       Error

        else

          write(iow,4000)
          if(ior.lt.0) then
            write(*,4000)
          end if

        endif

c     Axisymmetric option

      elseif(pcomp(type(1),'axis',4)) then

        gtypfl = .true.
        if(ndm.le.2) then

          if(pcomp(type(2),'tors',4)) then
            g2type = 8
            write(iow,2014)
            if(ior.lt.0) then
              write(*,2014)
            end if
          else
            g2type = 3
            write(iow,2003)
            if(ior.lt.0) then
              write(*,2003)
            end if
          endif

c       Error

        else

          write(iow,4000)
          if(ior.lt.0) then
            write(*,4000)
          end if

        endif

c     Kinematics: Small Deformation

      elseif(pcomp(type(1),'smal',4)) then

        gdeffl = .true.
        gdtype = 1
        write(iow,2004)
        if(ior.lt.0) then
          write(*,2004)
        end if

c     Kinematics: Finite Deformation

      elseif(pcomp(type(1),'fini',4)) then

        gdeffl = .true.
        gdtype = -1
        write(iow,2005)
        if(ior.lt.0) then
          write(*,2005)
        end if

c     Thermal-Mechanical Coupling: Temperature DOF

      elseif(pcomp(type(1),'temp',4)) then
        if(pcomp(type(2) ,'dof',3)) then

          gtdofl = .true.
          gtdof  = nint(td(1))
          if(gtdof.gt.0 .and. gtdof.le.ndf) then
            write(iow,2006) gtdof
            if(ior.lt.0) then
              write(*,2006) gtdof
            end if
          else
            if(ior.lt.0) then
              write(*,4001) gtdof
              return
            end if
            write(iow,4001) gtdof
            call plstop()
          end if
        end if

c     Define reference node/vector for element use

      elseif(pcomp(type(1),'refe',4)) then

        if(pcomp(type(2),'node',4)) then

          gref = 1
          do i = 1,ndm
            grefx(i) = td(i)
          end do ! i

          write(iow,2007) (grefx(i),i=1,ndm)
          if(ior.lt.0) then
            write(*,2007) (grefx(i),i=1,ndm)
          end if

        elseif(pcomp(type(2),'vect',4)) then

          gref = 2
          do i = 1,ndm
            gtref(i) = td(i)
          end do ! i

          write(iow,2008) (gtref(i),i=1,ndm)
          if(ior.lt.0) then
            write(*,2008) (gtref(i),i=1,ndm)
          end if

        elseif(pcomp(type(2),'pola',4)) then

          gref = 3
          do i = 1,ndm
            gtref(i) = td(i)
          end do ! i

          write(iow,2009) (gtref(i),i=1,ndm)
          if(ior.lt.0) then
            write(*,2009) (gtref(i),i=1,ndm)
          end if

        elseif(pcomp(type(2),'axia',4)) then

          gref = 4
          do i = 1,ndm
            gtref(i) = td(i)
          end do ! i

          write(iow,2010)
          if(ior.lt.0) then
            write(*,2010)
          end if

        else

          gref = 0
          write(iow,2011)
          if(ior.lt.0) then
            write(*,2011)
          end if

        endif

c     Global dof factors

      elseif(pcomp(type(1),'grou',4)) then

        groufl = .true.
        do i = 1,ndf
          gfac(i) = td(i)
        end do ! i

        write(iow,2012) (i,gfac(i),i=1,ndf)
        if(ior.lt.0) then
          write(*,2012) (i,gfac(i),i=1,ndf)
        end if

c     Global proportional load factors

      elseif(pcomp(type(1),'prop',4)) then

        groupl = .true.
        do i = 1,ndf
          gprop(i) = td(i)
        end do ! i

        write(iow,2021) (i,gprop(i),i=1,ndf)
        if(ior.lt.0) then
          write(*,2021) (i,gprop(i),i=1,ndf)
        end if

c     Rayleigh Damping factors

      elseif(pcomp(type(1),'rayl',4)) then

        grayfl  = .true.
        gray(1) = td(1)
        gray(2) = td(2)
        write(iow,2013) gray(1),gray(2)
        if(ior.lt.0) then
          write(*,2013) gray(1),gray(2)
        endif
        rayla0 = gray(1)
        rayla1 = gray(2)

c     Angular velocity

      elseif(pcomp(type(1),'omeg',4)) then

        gomgfl  = .true.

c       Coordinate option

        if(pcomp(type(2),'coor',4) .or. pcomp(type(2),'node',4)) then

          do i = 1,ndm
            gomex(i) = td(i)
          end do ! i

          write(iow,2016) (gomex(i),i=1,ndm)
          if(ior.lt.0) then
            write(*,2016) (gomex(i),i=1,ndm)
          end if

        elseif(pcomp(type(2),'vect',4)) then

          gvec = 0.0d0
          do i = 1,ndm
            gomev(i) = td(i)
            gvec     = gvec + gomev(i)**2
          end do ! i
          if(gvec.gt.0.0d0) then
            gvec = 1.d0/sqrt(gvec)
            do i = 1,ndm
              gomev(i) = gomev(i)*gvec
            end do ! i
          else
            write(iow,4002) (i,gomev(i),i=1,ndm)
            if(ior.lt.0) then
              write(*,4002) (i,gomev(i),i=1,ndm)
            endif
            call plstop()
          endif

          write(iow,2017) (gomev(i),i=1,ndm)
          if(ior.lt.0) then
            write(*,2017) (gomev(i),i=1,ndm)
          end if

        else
          if(pcomp(type(2),'cycl',4)) then
            gomega(1) = td(1)*2.d0*pi
            gomega(2) = td(2)*2.d0*pi
            gomega(3) = td(3)*2.d0*pi
          else
            gomega(1) = td(1)
            gomega(2) = td(2)
            gomega(3) = td(3)
          endif
          write(iow,2015) gomega(1)
          if(ior.lt.0) then
            write(*,2015) gomega(1)
          endif
        endif

c     Quadrature options

      elseif(pcomp(type(1),'quad',4)) then
        if(pcomp(type(2),'node',4) .or. pcomp(type(2),'noda',4)) then
          gquadn = 1.d0
        elseif(pcomp(type(2),'gaus',4)) then
          gquadn = 0.d0
        endif
        write(iow,2018) gqd(nint(gquadn)+1)
        if(ior.lt.0) then
          write(*,2018) gqd(nint(gquadn)+1)
        endif

c     Augmenting options

      elseif(pcomp(type(1),'augm',4)) then
        if(pcomp(type(2),'off',3)) then
          gaugm = 0.0d0
          write(iow,2019) 'off'
          if(ior.lt.0) then
            write(*,2019) 'off'
          endif
        else
          gaugm = 1.0d0
          write(iow,2019) 'on'
          if(ior.lt.0) then
            write(*,2019) 'on'
          endif
        endif

c     Global equation numbers

      elseif(pcomp(type(1),'equa',4)) then

        geqnum = nint(td(1))
        gpart  = nint(td(2))
        write(iow,2020) geqnum,gpart
        if(ior.lt.0) then
          write(*,2020) geqnum,gpart
        endif
        if(geqnum.gt.0) then
          setval = palloc( 258,'GUVAL', geqnum, 2)
        else
          write(iow,*) ' *WARNING* Number global equations is zero'
        endif

        if(geqnum.gt.nst-nen*ndf) then
          write(iow,4003) geqnum,nst-nen*ndf
          call plstop()
        endif

c     User: Global parameters

      elseif(.not.pcomp(type(1),' ',1)) then

        call uglobl(type(1),td)

      elseif(pcomp(type(1),' ',1)) then

        return

      end if

      go to 100

c     Formats

2000  format(/5x,'G l o b a l   P a r a m e t e r s'/1x)

2001  format(10x,'Plane Stress Analysis'/1x)

2002  format(10x,'Plane Strain Analysis'/1x)

2003  format(10x,'Axisymmetric Analysis'/1x)

2004  format(10x,'Kinematics: Small Deformation'/1x)

2005  format(10x,'Kinematics: Finite Deformation'/1x)

2006  format(10x,'Thermo-mechanical Coupling: Temperature DOF =',i3/1x)

2007  format(10x,'Reference node coordinates'/15x,3(1p,e14.5:))

2008  format(10x,'Reference vector components'/15x,3(1p,e14.5:))

2009  format(10x,'Reference polar components'/15x,3(1p,e14.5:))

2010  format(10x,'Reference axial')

2011  format(10x,'No reference state set')

2012  format(10x,'Ground acceleration components'/(15x,i10,1p,e14.5:))

2013  format(10x,'Rayleigh Damping Ratios'/
     &       15x,'Mass  value: a0',1p,1e14.5/
     &       15x,'Stiff value: a1',1p,1e14.5)

2014  format(10x,'Axisymmetric Analysis with Torsion'/1x)

2015  format(10x,'Angular Velocity (radians/time)',1p,1e16.5)

2016  format(10x,'Angular Velocity Reference Coordinate'/
     &       15x,3(1p,e14.5:))

2017  format(10x,'Angular Velocity Vector'/15x,3(1p,e14.5:))

2018  format(10x,'Quadrature Type: ',a)

2019  format(10x,'Augmenting ',a)

2020  format(10x,'Global number of equations =',i5/
     &       10x,'Global partion number      =',i5)

2021  format(10x,'Ground proportional components'/(15x,i10,1p,e14.5:))

3000  format(/5x,'Input Global Parameter')

4000  format(10x,'*WARNING* Can not set plane/axisymmetric option in ',
     &           'this mode.')

4001  format(10x,'*ERROR* GLOBAL: Temperature degree-of-freedom input',
     &           ' as ',i3/1x)

4002  format(10x,'*ERROR* GLOBAL: Zero length angular velocity vector'/
     &      (18x,'X(',i1,') =',1p,1e12.4:/))

4003  format(10x,'*ERROR* GLOBAL: Number of global equations too large'/
     &           '          Global  =',i5/
     &           '          Element =',i5/
     &           '        Set NAD (field 7) on control record')

      end
