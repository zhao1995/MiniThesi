c$Id:$
      subroutine pnorml(ie,ix,x,norm,ip,nie,ndm,nen,nen1,numnp,numel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Include check on surface type for 3-d problems   22/01/2011
c       2. Force computation of normals for all 2-d node    15/02/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute normals to locate external nodes

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'pdata5.h'
      include 'pdata6.h'
      include 'pointer.h'
      include 'comblk.h'

      include 'p_int.h'

      logical  lclip,flag
      integer  nie,ndm,nen,nen1,numnp,numel
      integer  iel,iiel,iu,i,j,k,l,m,n, i1,i2,i3,i4, ma, nel,pstyp
      integer  ie(nie,*),ix(nen1,*),ip(*), iq(4,6), it(3,4), iplt(30)
      real*8   leng,xmin,xmax,tol
      real*8   dx1(3), dx2(3), dn(3), x(ndm,*), norm(3,*)

      save

c     8-node brick faces

      data iq/3,2,1,4, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 5,6,7,8/

c     4-node tet faces

      data it/1,2,4, 2,3,4, 3,1,4, 1,3,2/


c     Zero normal vector

      do n = 1,numnp
        do i = 1,3
          norm(i,n) = 0.0d0
        end do ! i
      end do ! n

c     Loop through nodes to compute maximum and minimum

      flag = .true.
      fp(1)  = np(190) - 1
      do n = 1,numnp
        if(mr(fp(1)+n).ge.0) then
          do i = 1,ndm
            if(flag) then
              xmin = x(i,n)
              xmax = x(i,n)
              flag = .false.
            endif
            xmin = min(xmin,x(i,n))
            xmax = max(xmax,x(i,n))
          end do ! i
        endif
      end do ! n

      tol = 1.d-8*(xmax-xmin)

c     Loop through elements

      do n = 1,numel
        ma = ix(nen1,n)
        pstyp = ie(1,ma)
        if(ix(nen1-1,n).ge.0 .and. ma.ge.0 .and. pstyp.gt.0) then
          iel = ie(nie-1,ma)
          do i = nen,1,-1
            if(ix(i,n).gt.0) then
              nel = i
              exit
            endif
          end do ! i
          call plftyp(pstyp,nel,iel)
          if(iel.gt.0) then
            iiel = inord(iel)
          else
            iiel = exord(-iel)
          endif

c         2-D problems

          if(ndm.eq.2) then

            call pltord(ix(1,n),iel, iu,iplt)
            do i = 1,iu-1
              i1 = ix(iplt(i  ),n)
              i2 = ix(iplt(i+1),n)

              dn(1) = x(2,i2) - x(2,i1)
              dn(2) = x(1,i1) - x(1,i2)

              leng  = sqrt(dn(1)**2 + dn(2)**2)
              if(leng.gt.tol) then
                dn(1) = dn(1)/leng
                dn(2) = dn(2)/leng
                norm(1,i1) = norm(1,i1) + dn(1)
                norm(2,i1) = norm(2,i1) + dn(2)
                norm(1,i2) = norm(1,i2) + dn(1)
                norm(2,i2) = norm(2,i2) + dn(2)
              endif
            end do ! i

c         3-D problems

          elseif(ndm.eq.3) then

            if (iiel .eq. 9 .and. pstyp.eq.3) then
              if( lclip(ix(1,n),4,x,ndm) ) then
                do m = 1,4
                  i = 1
                  do j = 2,3
                    if(ix(it(j,m),n).lt.ix(it(i,m),n)) i = j
                  end do ! j
                  j = mod(i,3) + 1
                  k = mod(j,3) + 1

                  i1 = ix(it(i,m),n)
                  i2 = ix(it(j,m),n)
                  i3 = ix(it(k,m),n)

                  do j = 1,3
                    dx1(j) = x(j,i2) - x(j,i1)
                    dx2(j) = x(j,i3) - x(j,i1)
                  end do ! j

                  dn(1) = dx1(2)*dx2(3) - dx1(3)*dx2(2)
                  dn(2) = dx1(3)*dx2(1) - dx1(1)*dx2(3)
                  dn(3) = dx1(1)*dx2(2) - dx1(2)*dx2(1)

                  leng  = sqrt(dn(1)**2 + dn(2)**2 + dn(3)**2)

                  if(leng.gt.tol) then

                    leng = 1.d0/leng
                    do j = 1,3
                      dn(j) = dn(j)*leng
                    end do ! j

                    do j = 1,3
                      norm(j,i1) = norm(j,i1) + dn(j)
                      norm(j,i2) = norm(j,i2) + dn(j)
                      norm(j,i3) = norm(j,i3) + dn(j)
                    end do ! j

                  endif
                end do ! m
              end if

            elseif (iiel .gt. 10 .and. pstyp.eq.3) then

              if( lclip(ix(1,n),8,x,ndm) ) then
                do m = 1,6
                  i = 1
                  do j = 2,4
                    if(ix(iq(j,m),n).lt.ix(iq(i,m),n)) i = j
                  end do ! j
                  j = mod(i,4) + 1
                  k = mod(j,4) + 1
                  l = mod(k,4) + 1
                  i1 = ix(iq(i,m),n)
                  i2 = ix(iq(j,m),n)
                  i3 = ix(iq(k,m),n)
                  i4 = ix(iq(l,m),n)

                  do j = 1,3
                    dx1(j) = x(j,i3) - x(j,i1)
                    dx2(j) = x(j,i4) - x(j,i2)
                  end do ! j

                  dn(1) = dx1(2)*dx2(3) - dx1(3)*dx2(2)
                  dn(2) = dx1(3)*dx2(1) - dx1(1)*dx2(3)
                  dn(3) = dx1(1)*dx2(2) - dx1(2)*dx2(1)

                  leng  = sqrt(dn(1)**2 + dn(2)**2 + dn(3)**2)
                  if(leng.ne.0.0d0) then
                    leng = 1.d0/leng
                    do j = 1,3
                      dn(j) = dn(j)*leng
                    end do ! j
                    do j = 1,3
                      norm(j,i1) = norm(j,i1) + dn(j)
                      norm(j,i2) = norm(j,i2) + dn(j)
                      norm(j,i3) = norm(j,i3) + dn(j)
                      norm(j,i4) = norm(j,i4) + dn(j)
                    end do ! j
                  endif
                end do ! m
              end if

            else  ! Use surface quadrilateral

              call pltord(ix(1,n),iel, iu,iplt)
              if(iu.gt.3) then

                iu = min(nel,4)

                do j = 1,ndm
                  dx1(j) = x(j,ix( 3,n)) - x(j,ix(1,n))
                  dx2(j) = x(j,ix(iu,n)) - x(j,ix(2,n))
                end do ! j

                dn(1) = dx1(2)*dx2(3) - dx1(3)*dx2(2)
                dn(2) = dx1(3)*dx2(1) - dx1(1)*dx2(3)
                dn(3) = dx1(1)*dx2(2) - dx1(2)*dx2(1)

                leng  = sqrt(dn(1)**2 + dn(2)**2 + dn(3)**2)
                if(leng.ne.0.0d0) then
                  leng = 1.d0/leng
                  do j = 1,3
                    dn(j) = dn(j)*leng
                  end do ! j
c                 do m = 1,iu-1
                  do m = 1,nel
                    i1 = ix(m,n)
                    if(i1.gt.0) then
                      do j = 1,3
                        norm(j,i1) = norm(j,i1) + dn(j)
                      end do ! j
                    endif
                  end do ! m
                endif
              endif ! iu > 3 (not rod)

            endif
          endif
        endif

      end do ! n

c     Set indicator for nodes

      do n = 1,numnp
        ip(n) = 0
        if(max(abs(norm(1,n)),abs(norm(2,n)),abs(norm(3,n)))
     &                                           .gt. 1.d-05) then
          leng = 1.d0/sqrt(norm(1,n)**2 + norm(2,n)**2 + norm(3,n)**2)
          do j = 1,ndm
            norm(j,n) = norm(j,n)*leng
          end do ! j
          ip(n) = 1
        endif
      end do ! n

      end
