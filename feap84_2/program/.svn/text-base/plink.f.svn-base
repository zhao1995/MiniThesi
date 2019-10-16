c$Id:$
      subroutine plink(id,ntyp,ndf,numnp,neq,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Move 'go to 40' from line 103 to line 101        21/02/2007
c          Modify line with 'sign' to check for error combinations.
c       2  Allow up to 48 dof/node for links                26/06/2007
c       3. Remove the go to 40 statement after endif.       29/06/2007
c          Change 'm1.lt.n1' and 'm2.lt.n2' to 'le' in
c          the 'go to 20' expression.
c       4. Return checks in 3. to 'lt'; add 'exit' so that  08/07/2007
c          last non-zero input record is processed.
c       5. Add read of ldnum,ldprp,ldflg                    06/03/2009
c       6. Add read of spnum,spflg                          09/03/2009
c       7. Increase dimension of idl and jdl to 44          04/05/2010
c       8. Set length of file extension to 8                25/01/2011
c          Call setext to set 'fext'
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Link degrees of freedom to have same solution value

c      Inputs:
c         id(ndf,*) - Equation numbers before link
c         ntyp(*)   - Node type: >= 0 exist; < 0 not active
c         ndf       - Number dof/node
c         numnp     - Number of nodes in mesh
c         prt       - Output links performed if true

c      Outputs:
c         id(ndf,*) - Equation numbers for each dof after link
c         neq       - Number of equations active after link
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

      logical   prt,lsave,errck, pinput, oprt,prth, exit
      character fnamr*132,fext*8, type*4
      integer   ndf, numnp, neq, iosfil
      integer   i, ii, i1,i2, j, j1,j2
      integer   m1,m2, n1,n2, nmax, id(ndf,*),idl(44),jdl(44),ntyp(*)
      real*8    td(48)

      save

c     Routine to link degrees of freedom together

      fnamr =  fsav
      call setext('link',0,fext,.false.)
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

c       Output header information

        if(prt) then
          write(iow,2000) (i,i=1,ndf)
          if(iosfil.lt.0) then
            write(*,2000) (i,i=1,ndf)
          endif
        endif

c       Start link search

        m1 = 0
10      if(ndf.le.12) then
          errck = pinput(td,4+ndf)
        else
          errck = pinput(td,16)
          if(ndf.le.28) then
            errck = pinput(td(17),ndf-12)
          else
            errck = pinput(td(17),16)
            errck = pinput(td(33),min(16,ndf-28))
            if(ndf.gt.44) then   ! Dimension limitation
              write(iow,3005)
              call plstop()
            endif
          endif
        endif
        n1 = nint(td(1))
        n2 = nint(td(2))
        i1 = nint(td(3))
        i2 = nint(td(4))
        do i = 1,ndf
          idl(i) = nint(td(i+4))
        end do ! i
        if(n1.eq.0 .or. n1.gt.numnp .or.
     &     n2.eq.0 .or. n2.gt.numnp) then
          close(ior)
          ior = iosfil
          if(n1.gt.numnp. or. n2.gt.numnp) then
            if(prt) then
              write(iow,3001) n1,n2
              if(ior.lt.0) write(*,3001) n1,n2
            endif
            n1 = 0
            n2 = 0
          endif
          exit = .true.
        else
          exit = .false.
        endif
        if(m1.gt.0) then
20        if(ntyp(m1).lt.0 .or. ntyp(m2).lt.0) then
            if(prt) then
              write(iow,3001) m1,m2
              if(iosfil.lt.0) write(*,3001) m1,m2
            endif
          elseif(m1.eq.m2) then
            if(prt) then
              write(iow,3002) m1,m2
              if(iosfil.lt.0) write(*,3002) m1,m2
            endif
          else
            if(prt) then
              write(iow,2001) m1,m2,(jdl(i),i=1,ndf)
              if(iosfil.lt.0) then
                write(*,2001) m1,m2,(jdl(i),i=1,ndf)
              endif
            endif

c           Check that node pair has not already linked d.o.f.

            do j = 1,ndf
              if(ndfp(j).eq.npart .and. jdl(j).eq.0) then
                if(id(j,m1).gt.0 .and. id(j,m2).gt.0) then

c                 Select node to renumber dof

                  if(id(j,m1).eq.id(j,m2)) then
                    if(prt) then
                      write(iow,3004) m1,m2,j
                      if(iosfil.lt.0) then
                        write(*,3004) m1,m2,j
                      endif
                    endif
                    go to 30
                  elseif(id(j,m1).lt.id(j,m2)) then
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

c                 Loop through all nodes to reduce equation numbers

                  errck = .false.
                  do i = 1,ndf
                    if(ndfp(i).eq.npart) then
                      do ii = 1,numnp
                        if(id(i,ii).gt.nmax) then
                          id(i,ii) = id(i,ii) - 1
                          errck    = .true.
                        endif
                      end do ! ii
                    endif
                  end do ! i
                  if(errck) neq = neq - 1
                else

c                 Error

                  write(ilg,3003) m1,m2,j
                  write(iow,3003) m1,m2,j
                  if(ior.lt.0) then
                    write(*,3003) m1,m2,j
                    go to 40
                  endif
                  call plstop()

                endif
              endif
30            continue
            end do ! j
            if(exit) go to 40
          endif
          m1 = m1 + j1
          m2 = m2 + j2
          if( (j1.gt.0 .and. m1.lt.n1)  .or.
     &        (j2.gt.0 .and. m2.lt.n2) ) go to 20

c         Check for error combinations

          if(n1.gt.0) then
            j1 =  1
          else
            j1 = -1
          endif
          if(n2.gt.0) then
            j2 =  1
          else
            j2 = -1
          endif
          if(j1*j2.le.0) go to 40
        endif
        m1 = n1
        m2 = n2
        j1 = i1
        j2 = i2
        do i = 1,ndf
          jdl(i) = idl(i)
        end do ! i

        go to 10
      else
        write(ilg,3000)
        write(iow,3000)
        call plstop()
      endif

c     Reset parameter values

40    do i = 0,36
        do j = 1,26
          vvv(j,i) = vvsave(j,i)
        end do ! j
      end do ! i
      do i = 1,3
        x0(i) = x0sav(i)
      end do ! i
      prt = oprt

c     Check that the number of equations is correct

      neq = 0
      do i = 1,ndf
        if(ndfp(i).eq.npart) then
          do ii = 1,numnp
            neq = max(neq,id(i,ii))
          end do ! ii
        endif
      end do ! i

c     Formats

1000  format(a4,2x,a12,i8,2l5)
1001  format(1p,4e20.12)
1002  format(3i8,2l3)

2000  format('    N o d a l    L i n k    O u t p u t s'//
     &    '     Linked  Pairs     DOF Link Pattern (0=link; 1=no link)'/
     &    '    1-Node   2-Node',8(i3,'-dof')/(18x,8(i3,'-dof')))

2001  format(2i9,8i7:/(18x,8i7))

3000  format(5x,'*ERROR* PLINK: Link file does not exist')

3001  format(5x,'*WARNING* Nodes',i8,' and',i8,' are not active.')

3002  format(5x,'*WARNING* Nodes',i8,' and',i8,' are same.')

3003  format(5x,'*ERROR* PLINK: Attempt to link restrained DOF',
     &          ' to active DOF'/
     &       5x,'        Nodes are',i8,' and',i8,'; DOF =',i4)

3004  format(5x,'*WARNING* Nodes',i8,' and',i8,
     &          ' already linked for DOF =',i4)

3005  format(5x,'*ERROR* PLINK: Program cannot link more than 44 DOF')

      end
