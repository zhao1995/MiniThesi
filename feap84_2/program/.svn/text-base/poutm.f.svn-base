c$Id:$
      subroutine poutm (ld,ie,d,id,x,u,ix,ang,f,t,fpro,ndam,nmas,nsti,
     &                  ip,iq, lct,ct, j)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set search from character 256                    21/12/2008
c       2. Add output of load groups to flat files          05/04/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output a mesh with links/tie effects removed

c      Inputs:
c         ie(*)     - Assembley information for material sets
c         d(*)      - Material set parameters
c         id(*)     - Equation numbers for active dof
c         x(*)      - Nodal coordinates for mesh
c         u(*)      - Nodal coordinates for mesh
c         ix(*)     - Element nodal connection list
c         ang(*)    - Nodal angles
c         f(*)      - Interpolated nodal force/displacement
c         t(*)      - Nodal temperature values
c         fpro(*)   - Nodal proportional load numbers
c         ndam(*)   - Nodal dampers
c         nmas(*)   - Nodal masses
c         nsti(*)   - Nodal stiffness
c         lct       - Command option
c         ct(*)     - Command data
c         j         - Command number in this routine

c      Scratch:
c         ld(*)     - Element local/global equation numbers
c         ip(*)     - Active node numbers
c         iq(*)     - Active element numbers

c      Outputs:
c         Depends on command number j
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'
      include  'bdata.h'
      include  'cdata.h'
      include  'cdat1.h'
      include  'chdata.h'
      include  'comfil.h'
      include  'compac.h'
      include  'ddata.h'
      include  'endata.h'
      include  'fdata.h'
      include  'idptr.h'
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
      include  'prflag.h'
      include  'prlod.h'
      include  'rdata.h'
      include  'rigid1.h'
      include  'rjoint.h'
      include  'sdata.h'
      include  'tdata.h'
      include  'comblk.h'

      logical   lread,pcomp, pflg
      character lct*15, fnamr*134,fext*8, eformat*26, size*2
      integer   i, iz, j, k, nn, numnp1, numel1
      integer   ld(*),ie(nie,*),id(ndf,*),ix(nen1,*),fpro(ndf,*)
      integer   ip(*),iq(*)
      real*8    d(ndd,*),x(ndm,*),u(ndf,*),ang(*),f(ndf,numnp,2),t(*)
      real*8    ndam(ndf,*),nmas(ndf,*),nsti(ndf,*), ct(3),x0(3)

      save

c     Determine number of active nodes/elements

      call pnumne(ix,nen1,nen,numnp,numel,iq,mr(np(89)),ip,
     &            iq(numnp+2),numnp1,numel1)

c     [outm]    Output a Renumbered Input File

      if(j.eq.1) then

c     1. Output Mesh -

c       Control data

        iz     = 0
        fnamr  =  finp
        if(pcomp(lct,'bina',4)) then

c         Binary write of mesh data

          fext   =  'binary'
          call addext(fnamr,fext,128,8)
          write(*,*) 'Binary mode of mesh output: filename = ',fnamr
          call opnfil(fext,fnamr,-3,ios,lread)
          rewind ios
          write(ios) head,numnp1,numel,nummat,ndm,ndf,nen,ndd,nud
          write(ios) ((x(i,ip(nn)),i=1,ndm),nn=1,numnp1)
          write(ios) ((iq(ix(i,nn)+1),i=1,nen),ix(nen1,nn),nn=1,numel)
          write(ios) ((id(i,ip(nn)),i=1,ndf),nn=1,numnp1)
          write(ios) (( f(i,ip(nn),1),i=1,ndf),nn=1,numnp1)
          write(ios) (( f(i,ip(nn),2),i=1,ndf),nn=1,numnp1)
          write(ios) (  t(ip(nn)) ,nn=1,numnp1)
          write(ios) (  hr(np(45)+ip(nn)-1) ,nn=1,numnp1)
          write(ios) ((ie(i,nn),i=1,nie),nn=1,nummat)
          write(ios) (( d(i,nn),i=1,ndd),nn=1,nummat)
          write(ios) ia,inord,ipord

c       Text (ASCII) write of domain decomposition mesh

        elseif(pcomp(lct,'doma',4)) then

          call poutdom(ld,id,x,ix,ang,f,t,ip,iq,iq(numnp+2),
     &                 numnp1,numel1)

c       Text (ASCII) write of mesh data for element "nn"

        elseif(pcomp(lct,'elem',4)) then

          nn = max(1,min(nint(ct(1)),numel))
          write(*,2001) '  Output for element ',nn
          i  = index(fnamr,' ')
          if(nn.lt.10) then
            write(fnamr(i:i+1),'(a1,i1)') '_',nn
          elseif(nn.lt.100) then
            write(fnamr(i:i+2),'(a1,i2)') '_',nn
          elseif(nn.lt.1000) then
            write(fnamr(i:i+3),'(a1,i3)') '_',nn
          elseif(nn.lt.10000) then
            write(fnamr(i:i+4),'(a1,i4)') '_',nn
          elseif(nn.lt.100000) then
            write(fnamr(i:i+5),'(a1,i5)') '_',nn
          elseif(nn.lt.1000000) then
            write(fnamr(i:i+6),'(a1,i6)') '_',nn
          elseif(nn.lt.10000000) then
            write(fnamr(i:i+7),'(a1,i7)') '_',nn
          elseif(nn.lt.100000000) then
            write(fnamr(i:i+8),'(a1,i8)') '_',nn
          endif
          fext = 'elm'
          call opnfil(fnamr,fnamr,-1,ios,lread)
          rewind ios
          write(ios,2000) head,nen,1,nummat,ndm,ndf,nen

c         Coordinates - nodal list

          do i = 1,ndm
            x0(i) = 0.0d0
            iz    = 0
            do k = 1,nen
              if(ix(k,nn).ne.0) then
                x0(i) = x0(i) + x(i,ix(k,nn))
                iz    = iz + 1
              endif
            end do ! k
            x0(i) = x0(i)/dble(iz)
          end do ! i
          iz = 0
          write(ios,2001) 'COORdinates ALL ',nen
          do k = 1,nen
            if(ix(k,nn).ne.0) then
              write(ios,2002) k,iz,((x(i,ix(k,nn))-x0(i)),i=1,ndm)
            endif
          end do ! k

c         Element connection list

          k  = 1
          write(ios,2001) 'ELEMent'
          write(ios,2003) k,iz,ix(nen1,nn),(i,i=1,nen)

c         Material/Params - list

          open(unit=iwd,file=fmtl,status='old')
          pflg = .true.
          write(ios,'(a)') ' '
          do while(pflg)
            read(iwd,'(a)',end=200) xxx
            do i = 256,1,-1
              if(xxx(i:i).ne.' ') go to 100
            end do ! i
            i = 1
100         if(xxx(i:i).eq.char(13)) xxx(i:i) = ' '
            write(ios,'(a)') xxx(1:i)
          end do ! while
200       close(iwd)

c         Closing list

          write(ios,2004)

c       Text (ASCII) write of mesh data

        else

          if(optmsh) then
            fext   =  'opt'
          else
            fext   =  'rev'
          endif
          call addext(fnamr,fext,128,6)
          i = index(fnamr,' ') - 1
          write(*,*) ' Text mode of mesh output: filename = ',fnamr(1:i)
          call opnfil(fext,fnamr,-1,ios,lread)
          rewind ios

c         Output problem parameters

          write(ios,2000) head,numnp1,numel1,nummat,ndm,ndf,nen

c         Coordinates - nodal list

          write(ios,2001) 'COORdinates ALL'
          if(.not.pcomp(lct,'defo',4)) then
            do nn = 1,numnp1
              write(ios,2002) nn,iz,(x(i,ip(nn)),i=1,ndm)
            end do ! nn
          else
            do nn = 1,numnp1
              write(ios,2002) nn,iz,
     &                       (x(i,ip(nn))+u(i,ip(nn)),i=1,min(ndm,ndf)),
     &                       (x(i,ip(nn)),i=min(ndm,ndf)+1,ndm)
            end do ! nn
          endif

c         Build output format to fit data

c         Format:    1...5...10...15...20...25...30
          eformat = '(i10,i10,i6,13i10/(16i10))'

          if(.not.pcomp(lct,'fixe',4)) then
            nn      = nint(log10(dble(max(numnp1,numel)))) + 2
            write(size,'(i2)') nn
            eformat( 3: 4) = size
            eformat( 7: 8) = size
            eformat(16:17) = size
            eformat(23:24) = size
          endif

c         Element connection list

          write(ios,2001) 'ELEMents ALL'
          numel1 = 0
          do nn = 1,numel
            if(ix(nen1-1,nn).ge.0) then
              numel1 = numel1 + 1
              do i = 1,nen
                ld(i) = iq(ix(i,nn)+1)
              end do ! i

c             Remove repeated node numbers

              if(pcomp(lct,'repe',4)) then
                do i = 1,nen-1
                  do j = i+1,nen
                    if(ld(j).eq.ld(i)) then
                      ld(j) = -1
                    endif
                  end do ! j
                end do ! i
                j = 1
                do i = 2,nen
                  if(ld(i).gt.0) then
                    j = j + 1
                    ld(j) = ld(i)
                  endif
                end do ! i
                do i = j+1,nen
                  ld(i) = 0
                end do ! i
              endif

c             Output element list

              write(ios,eformat) numel1,iz,ix(nen1,nn),(ld(i),i=1,nen)
            endif
          end do ! nn

c         Region list

          lread  = .true.
          numel1 = 0
          do nn = 1,numel
            if(ix(nen1-1,nn).ge.0) then
              numel1 = numel1 + 1
              if(ix(nen1-1,nn).gt.0) then
                if(lread) then
                  write(ios,2001) 'EREGion data'
                  lread = .false.
                endif
                write(ios,2007) numel1,iz,ix(nen1-1,nn)
              endif
            endif
          end do ! nn

c         Angle conditions lists

          lread = .true.
          do nn = 1,numnp1
            if(ang(ip(nn)).ne.0.0d0) then
              if(lread) then
                write(ios,2001) 'ANGLe conditions'
                lread = .false.
              endif
              write(ios,2002) nn,iz,ang(ip(nn))
            end if
          end do ! nn

c         Boundary conditions lists

          lread = .true.
          do nn = 1,numnp1
            pflg = .false.
            do i = 1,ndf
c             if(id(i,ip(nn)).le.0) then
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
              write(ios,2003) nn,iz,(ld(i),i=1,ndf)
            endif
          end do ! nn

c         Forced conditions lists

          lread = .true.
          do nn = 1,numnp1
            pflg = .false.
            do i = 1,ndf
              if(f(i,ip(nn),1).ne.0.0d0) then
                pflg = .true.
              end if
            end do ! i
            if(pflg) then
              if(lread) then
                lread = .false.
                write(ios,2001) 'FORCe conditions'
              endif
              write(ios,2002) nn,iz,(f(i,ip(nn),1),i=1,ndf)
            endif
          end do ! nn

c         Displacement conditions lists

          lread = .true.
          do nn = 1,numnp1
            pflg = .false.
            do i = 1,ndf
              if(f(i,ip(nn),2).ne.0.0d0) then
                pflg = .true.
              end if
            end do ! i
            if(pflg) then
              if(lread) then
                lread = .false.
                write(ios,2001) 'DISPlacement conditions'
              endif
              write(ios,2002) nn,iz,(f(i,ip(nn),2),i=1,ndf)
            endif
          end do ! nn

c         Proportional load lists

          lread = .true.
          do nn = 1,numnp1
            pflg = .false.
            do i = 1,ndf
              if(fpro(i,ip(nn)).ne.0) then
                pflg = .true.
              end if
            end do ! i
            if(pflg) then
              if(lread) then
                lread = .false.
                write(ios,2001) 'FPROportional load'
              endif
              write(ios,2003) nn,iz,(fpro(i,ip(nn)),i=1,ndf)
            endif
          end do ! nn

c         Add loads from ldtab

          if(np(265).ne.0) then
            call pldout(mr(np(265)),mr(np(266)),hr(np(267)), iq)
          endif

c         Nodal dampers, mass and stiffness

          if(nmfl) then

c           Nodal dampers

            lread = .true.
            do nn = 1,numnp1
              pflg = .false.
              do i = 1,ndf
                if(ndam(i,ip(nn)).ne.0.0d0) then
                  pflg = .true.
                end if
              end do ! i
              if(pflg) then
                if(lread) then
                  lread = .false.
                  write(ios,2001) 'DAMPERS nodal'
                endif
                write(ios,2002) nn,iz,(ndam(i,ip(nn)),i=1,ndf)
              endif
            end do ! nn

c           Nodal masses

            lread = .true.
            do nn = 1,numnp1
              pflg = .false.
              do i = 1,ndf
                if(nmas(i,ip(nn)).ne.0.0d0) then
                  pflg = .true.
                end if
              end do ! i
              if(pflg) then
                if(lread) then
                  lread = .false.
                  write(ios,2001) 'MASS nodal'
                endif
                write(ios,2002) nn,iz,(nmas(i,ip(nn)),i=1,ndf)
              endif
            end do ! nn

c           Nodal stiffness values

            lread = .true.
            do nn = 1,numnp1
              pflg = .false.
              do i = 1,ndf
                if(nsti(i,ip(nn)).ne.0.0d0) then
                  pflg = .true.
                end if
              end do ! i
              if(pflg) then
                if(lread) then
                  lread = .false.
                  write(ios,2001) 'STIFfness nodal'
                endif
                write(ios,2002) nn,iz,(nsti(i,ip(nn)),i=1,ndf)
              endif
            end do ! nn
          endif ! nmfl = .true.

c         Material/Params - list

          open(unit=iwd,file=fmtl,status='old')
          pflg = .true.
          write(ios,'(a)') ' '
          do while(pflg)
            read(iwd,'(a)',end=210) xxx
            do i = 256,1,-1
              if(xxx(i:i).ne.' ') go to 110
            end do ! i
            i = 1
110         if(xxx(i:i).eq.char(13)) xxx(i:i) = ' '
            write(ios,'(a)') xxx(1:i)
          end do ! while
210       close(iwd)

c         Contact list

          call contact(315)

c         Closing list

          write(ios,2004)
        end if
        close(ios)
      end if

c     [renu] Output a renumbered list

      if(j.eq.2) then
        write(iow,2005)
        do nn = 1,numnp1
          write(iow,2006) nn,ip(nn),(x(i,ip(nn)),i=1,ndm)
        end do ! nn
      end if

c     Formats

2000  format(20a4/6i10/' ')

2001  format(/a:,i10)

2002  format(i9,i2,1p,14e14.6)
!2002  format(i9,i2,1p,10e24.15)

2003  format(i10,i4,14i10/(16i10))

2004  format(/'END'//'INTEractive'//'STOP')

2005  format('  New node','  Old node','   X_i Coordinates')

!2006  format(2i10,2x,1p,3e24.15)
2006  format(2i10,2x,1p,3e14.6)

2007  format(i10,i4,i10)

      end
