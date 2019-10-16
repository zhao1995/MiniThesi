c$Id:$
      subroutine pelink(id,x,ndm,ndf,numnp,neq,prt, elnk,ip)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Allow up to 30 dof/node on inputs                26/06/2007
c       2. Correct setting for multiple settings; add 'ip'  13/12/2007
c       3. Set value of 'elnk' for periodic boundary use    25/12/2007
c       4. Add read of 'ldnum,ldprp,ldflg' and 'pload1.h'   25/01/2009
c          Change m1 to ndf + 2 after statement 10 and later
c          set m1 = ndf - 13 (instead of 14).
c       5. Add read of 'spnum,spflg'                        09/03/2009
c       6. Increase idl to 44 and td to 48 (same as plink)  04/05/2010
c       7. Set length of file extension to 8                25/01/2011
c       8. Change file extender to 'elin'                   17/03/2011
c       9. Remove extra test on values for reducing the id  19/05/2011
c      10. Call 'setext' to set 'fext' value                09/06/2011
c      11. Do loop m1 from 1 to numnp, not numnp-1          07/01/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform link of degree of freedom based on input
c               of edge coordinate values and direction.

c      Inputs:
c         x(ndm,*)    - Nodal coordinates of mesh
c         ndm         - Spatial dimension of mesh
c         ndf         - Number dof/node
c         numnp       - Number of nodes in mesh
c         prt         - Output results if true

c      Scratch::
c         ip(*)       - Used to flag nodes that are renumbered

c      Outputs:
c         id(ndf,*)   - Modified equation numbers from links
c         neq         - Number of equations after links
c         elnk(ndm,*) - Edge link direction
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdat2.h'
      include  'comfil.h'
      include  'conval.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'ioincl.h'
      include  'part0.h'
      include  'pload1.h'
      include  'trdata.h'

      logical   prt,lsave,errck,pinput,tinput,vinput,pcomp
      logical   oprt,prth, equalfl, dolink
      character text*15, fnamr*132,fext*8, type*4
      integer   ndm,ndf,numnp,neq, iosfil, i,ii,i1,i2, j, m1,m2,nmax
      integer   id(ndf,*),idl(45), elnk(ndm,numnp), ip(numnp)
      real*8    gap, x1, x2, td(48),x(ndm,*)

      save

c     Routine to link degrees of freedom together

      gap   =  1.d-04
      fnamr =  fsav
      call setext('elin', 0, fext, .false.)
      call addext(fnamr,fext,128,8)
      call opnfil(fext,fnamr,-2,ios,lsave)
      if(lsave) then
        iosfil = ior
        ior    = ios

        oprt = prt
        do i = 0,36
          do j = 1,26
            vvsave(j,i) = vvv(j,i)
          end do ! j
        end do ! i
        do i = 1,3
          x0sav(i) = x0(i)
        end do ! i
        read(ior,1000) type,fincld(isf),irecrd(isf),prt,prth
        read(ior,1001) vvv
        read(ior,1001) tr,xr,trdet,x0
        read(ior,1002) ldnum,ldprp,spnum,ldflg,spflg

10      m1 = min(ndf+2,15)
        errck =  tinput(text,1,td(2),m1)
        if(pcomp(text,'gap',3)) then
          gap   = td(2)
          if(prt) then
            write(iow,2001) gap
            if(iosfil.lt.0) then
              write(*,2001) gap
            endif
          endif
          go to 10
        else
          errck = vinput(text,15,td(1),1)
          i1  = nint(td(1))
        endif
        if(ndf.gt.13) then
          if(ndf.le.29) then
            errck = pinput(td(17),ndf-13)
          else
            errck = pinput(td(17),16)
            errck = pinput(td(33),min(16,ndf-29))
            if(ndf.gt.45) then  ! Dimension limitation
              write(iow,3001)
              call plstop()
            endif
          endif
        endif
        i1  = min(i1,ndm)
        x1  = td(2)
        x2  = td(3)
        do i = 1,ndf
          idl(i) = nint(td(i+3))
        end do ! i

c       End of inputs

        if(i1.eq.0) then

c         Compute elnk array

          call psetlnk(elnk,id,ndm,ndf,numnp)

c         Compute number of active equations

          neq = 0
          do ii = 1,numnp
            do i = 1,ndf
              neq = max(neq,id(i,ii))
            end do ! i
          end do ! ii

          close(ior)
          ior = iosfil

          prt = oprt
          do i = 0,36
            do j = 1,26
              vvv(j,i) = vvsave(j,i)
            end do ! j
          end do ! i
          do i = 1,3
            x0(i) = x0sav(i)
          end do ! i
          return
        endif

        do m1 = 1,numnp
          ip(m1) = 0
        end do ! m1

c       Look for matching coordinates in i1-direction

        do m1 = 1,numnp
          if(abs(x(i1,m1)-x1) .lt. gap ) then
            ip(m1) = 1
            do m2 = 1,numnp
              if( abs(x(i1,m2)-x2) .lt. gap .and. ip(m2).eq.0) then

c               Check that all other values are same

                if(abs(x1-x2) .gt. gap) then
                  equalfl = .false.
                  do i2 = 1,ndm
                    if(i2.ne.i1) then
                      if(abs(x(i2,m1)-x(i2,m2)) .gt. gap ) go to 20
                    endif
                  end do ! i2
                else
                  equalfl = .true.
                endif

c               Match found print result

                ip(m2) = 1
                if(prt) then
                  write(iow,2000) m1,m2,(idl(i),i=1,ndf)
                  if(iosfil.lt.0) then
                    write(*,2000) m1,m2,(idl(i),i=1,ndf)
                  endif
                endif

c               Modify equations and set dolink

                dolink = .false.
                do j = 1,ndf
                  if(ndfp(j).eq.npart .and. idl(j).eq.0) then

c                   Both nodes active

                    if(id(j,m1).gt.0  .and. id(j,m2).gt.0) then

                      dolink = .true.

                      if(id(j,m1).ne.id(j,m2)) then

c                       Select node to renumber dof

                        if(id(j,m1).lt.id(j,m2)) then
                          nmax     = id(j,m2)
                          id(j,m2) = id(j,m1)
                        else
                          nmax     = id(j,m1)
                          id(j,m1) = id(j,m2)
                        endif
                        do ii = 1,numnp
                          if(id(j,ii).eq.nmax) then
                            id(j,ii) = id(j,m1)
                            end if
                        end do ! ii

c                       Loop through all nodes to reduce eq numbers

                        do i = 1,ndf
                          if(ndfp(i).eq.npart) then
                            do ii = 1,numnp
                              if(id(i,ii).gt.nmax) then
c                             if(id(i,ii).gt.nmax     .and.
c    &                           id(i,ii).ne.id(i,m1)) then
                                id(i,ii) = id(i,ii) - 1
                              endif
                            end do ! ii
                          endif
                        end do ! i

                      endif ! id's not equal
                    endif ! id check
                  end if
                end do ! j

c               Set e-link

                if(dolink) then
c                 if(x(i1,m2).gt.x(i1,m1)) then
c                   elnk(i1,m2) =  m1
c                   elnk(i1,m1) = -m2
c                 else
c                   elnk(i1,m1) =  m2
c                   elnk(i1,m2) = -m1
c                 endif
                  dxlnk(i1) = abs(x(i1,m2) - x(i1,m1))
                endif

              endif
20            continue
            end do ! m2
            if(equalfl) go to 10
          endif
        end do ! m1
        go to 10
      else
        write(iow,3000)
        write(ilg,3000)
        call plstop()
      endif

c     Formats

1000  format(a4,2x,a12,i8,2l5)
1001  format(1p,4e20.12)
1002  format(3i8,2l3)

2000  format(5x,'Link node',i8,' to',i8,' DOF = 0 to link:',12i4)
2001  format(/5x,'Search gap =',1p,1e12.4/)

3000  format(5x,'*ERROR* PELINK: Link file does not exist')
3001  format(5x,'*ERROR* PELINK: Program cannot link more than 45 DOF')

      end

      subroutine psetlnk(elnk,id,ndm,ndf,numnp)

      implicit   none

      integer    ndm,ndf,numnp, i,m,n,nn, ii
      integer    elnk(ndm,numnp), id(ndf,numnp)

      do i = 1,ndm
        do n = 1,numnp
          ii = min(i,ndf)
          if(id(ii,n).gt.0) then
            nn = id(ii,n)
            do m = n+1,numnp
              if(id(ii,m).eq.nn .and. elnk(i,m).eq.0) then
                elnk(i,m) = n
              endif
            end do ! m
          endif
        end do ! n
      end do ! i

      end
