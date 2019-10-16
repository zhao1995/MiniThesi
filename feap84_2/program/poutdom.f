c$Id:$
      subroutine poutdom (ld,id,x,ix,ang,f,t,ip,iq,dn, numnp1,numel1)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set search from character 256                    21/12/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output a mesh with links/tie effects removed

c      Inputs:
c         d(*)      - Material set parameters
c         id(*)     - Equation numbers for active dof
c         x(*)      - Nodal coordinates for mesh
c         u(*)      - Nodal coordinates for mesh
c         ix(*)     - Element nodal connection list
c         ang(*)    - Nodal angles
c         f(*)      - Interpolated nodal force/displacement
c         t(*)      - Nodal temperature values

c      Scratch:
c         ld(*)     - Element local/global equation numbers
c         ip(*)     - Active node numbers
c         iq(*)     - Active element numbers
c         dn(*)     - Active region (domain) numbers

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'
      include  'bdata.h'
      include  'cdata.h'
      include  'chdata.h'
      include  'comfil.h'
      include  'compac.h'
      include  'ddata.h'
      include  'endata.h'
      include  'fdata.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'mdata.h'
      include  'mxsiz.h'
      include  'ndata.h'
      include  'part0.h'
      include  'pdata2.h'
      include  'pdata6.h'
      include  'plflag.h'
      include  'pointer.h'
      include  'prlod.h'
      include  'rdata.h'
      include  'rigid1.h'
      include  'rjoint.h'
      include  'sdata.h'
      include  'tdata.h'
      include  'comblk.h'

      logical   lread, pflg
      character fnamr*134,fext*8, eformat*26, size*2
      integer   i,j, iz, ne,nn, numnp1, numel1, maxdom, numnpd,numeld
      integer   ld(*),id(ndf,*),ix(nen1,*),ip(*),iq(*),dn(*)
      real*8    x(ndm,*),ang(*),f(ndf,numnp,2),t(*)

      save

      data      iz / 0 /
c     Determine number of regions (domains)

      maxdom = 0
      do j = 1,numel
        maxdom = max(maxdom,ix(nen1-1,j))
      end do ! j

      do j = 1,maxdom
        do nn = 1,numnp1
          dn(nn) = 0
        end do ! nn

c       Set active node numbers for domain

        numeld = 0
        do ne = 1,numel1
          if(ix(nen1-1,ne).eq.j) then
            numeld = numeld + 1
            do  nn = 1,nen
              if(iq(ix(nn,ne)+1).gt.0) then
                dn(iq(ix(nn,ne)+1)) = 1
              endif
            end do ! nn
          endif
        end do ! ne

        numnpd = 0
        do nn = 1,numnp1
          if(dn(nn).eq.1) then
            numnpd = numnpd + 1
            dn(nn) = numnpd
          endif
        end do ! nn

c       Name and open output file

        fnamr = finp
        if(j.lt.10) then
          write(fext,'(a3,i1)') '000',j
        elseif(j.lt.100) then
          write(fext,'(a2,i2)') '00',j
        elseif(j.lt.1000) then
          write(fext,'(a1,i3)') '0',j
        elseif(j.lt.10000) then
          write(fext,'(i4)') j
        endif

        call addext(fnamr,fext,128,8)
        write(*,*) 'Text mode of mesh output: filename = ',fnamr
        call opnfil(fext,fnamr,-1,ios,lread)
        rewind ios

c       Output problem parameters

        write(ios,2000) head,numnpd,numeld,nummat,ndm,ndf,nen

c       Coordinates - nodal list

        write(ios,2001) 'COORdinates ALL'
        do nn = 1,numnp1
          if(dn(nn).gt.0) then
            write(ios,2002) dn(nn),iz,(x(i,ip(nn)),i=1,ndm)
          endif
        end do ! nn

c       Build output format to fit data

c       Format:    1...5...10...15...20...25...30
        eformat = '(i10,i10,i6,13i10/(16i10))'

        nn      = nint(log10(dble(max(numnp1,numel)))) + 2
        write(size,'(i2)') nn
        eformat( 3: 4) = size
        eformat( 7: 8) = size
        eformat(16:17) = size
        eformat(23:24) = size

c       Element connection list

        write(ios,2001) 'ELEMents ALL'
        numeld = 0
        do nn = 1,numel
          if(ix(nen1-1,nn).eq.j) then
            numeld = numeld + 1
            write(ios,eformat) numeld,iz,ix(nen1,nn),
     &                     (dn(iq(ix(i,nn)+1)),i=1,nen)
          endif
        end do ! nn

c       Angle conditions lists

        lread = .true.
        do nn = 1,numnp1
          if(dn(nn).gt.0 .and. ang(ip(nn)).ne.0.0d0) then
            if(lread) then
              write(ios,2001) 'ANGLe conditions'
              lread = .false.
            endif
            write(ios,2002) dn(nn),iz,ang(ip(nn))
          end if
        end do ! nn

c       Boundary conditions lists

        lread = .true.
        do nn = 1,numnp1
          if(dn(nn).gt.0) then
            pflg = .false.
            do i = 1,ndf
              if(id(i,ip(nn)).ne.0) then
                pflg = .true.
                ld(i) = 1
              else
                ld(i) = 0
              end if
            end do ! i

            if(pflg) then
              if(lread) then
                write(ios,2001) 'BOUNdary conditions'
                lread = .false.
              endif
              write(ios,2003) dn(nn),iz,(ld(i),i=1,ndf)
            endif
          end if
        end do ! nn

c       Forced conditions lists

        lread = .true.
        do nn = 1,numnp1
          if(dn(nn).gt.0) then
            pflg = .false.
            do i = 1,ndf
              if(f(i,ip(nn),1).ne.0.0d0) then
                pflg = .true.
              end if
            end do ! i
            if(pflg) then
              if(lread) then
                write(ios,2001) 'FORCe conditions'
                lread = .false.
              endif
              write(ios,2002) dn(nn),iz,(f(i,ip(nn),1),i=1,ndf)
            endif
          end if
        end do ! nn

c       Displacement conditions lists

        lread = .true.
        do nn = 1,numnp1
          if(dn(nn).gt.0) then
            pflg = .false.
            do i = 1,ndf
              if(f(i,ip(nn),2).ne.0.0d0) then
                pflg = .true.
              end if
            end do ! i
            if(pflg) then
              if(lread) then
                write(ios,2001) 'DISPlacement conditions'
                lread = .false.
              endif
              write(ios,2002) dn(nn),iz,(f(i,ip(nn),2),i=1,ndf)
            end if
          end if
        end do ! nn

c       Temperature array list

        lread = .true.
        do nn = 1,numnp1
          if(dn(nn).gt.0) then
            pflg = .false.
            if(t(ip(nn)).ne.0.0d0) then
              pflg = .true.
            end if
            if(pflg) then
              if(lread) then
                write(ios,2001) 'TEMPerature conditions'
                lread = .false.
              endif
              write(ios,2002) dn(nn),iz,t(ip(nn))
            end if
          end if
        end do ! nn

c       Material/Params - list

        open(unit=iwd,file=fmtl,status='old')
        pflg = .true.
        write(ios,'(a)') ' '
        do while(pflg)
          read(iwd,'(a)',end=200) xxx
          do i = 256,1,-1
            if(xxx(i:i).ne.' ') go to 100
          end do ! i
          i = 1
100       if(xxx(i:i).eq.char(13)) xxx(i:i) = ' '
          write(ios,'(a)') xxx(1:i)
        end do ! while
200     close(iwd)

c       Contact list

c       call contact(315)

c       Closing list

        write(ios,2004)
        close(ios)
      end do ! j

c     Formats

2000  format('NOCOunt'/20a4/6i10/' '/'NOPArse'/'NOPRint'/' ')

2001  format(/a)

!2002  format(i9,i2,1p,10e23.15)
2002  format(i9,i2,1p,10e15.7)

2003  format(i9,i3,14i3/(16i3))

2004  format(/'END'//'INTEractive'//'STOP')

      end
