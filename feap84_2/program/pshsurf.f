c$Id:$
      subroutine pshsurf(x, xlm, norm, ix, ie, mo)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change use of rotation update for non-smooth to  06/02/2009
c          type -1 (This needs to be revisited)
c       2. Recode for non-smooth surfaces                   01/11/2010
c       3. Don't allow non-smooth shell if ndf = 5          21/06/2013
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Compute directors for shells

c      Inputs:
c      Outputs:
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'erotas.h'
      include   'iofile.h'
      include   'sdata.h'

      integer    i,j, ii,jj,kk, ma, n,nel
      real*8     x1,x2,x3, vn, v1(3),v2(3),v3(3)
      integer    ix(nen1,*), ie(nie,*), mo(*)
      real*8     x(ndm,*), xlm(9,6,*), norm(3,*)

      save

c     Zero norm

      do n = 1,numnp
        do i = 1,3
          norm(i,n) = 0.0d0
        end do ! i
      end do ! n

c     Compute director for shell elements

      do n = 1,numel

c       Check if finite deformation shell element

        ma = ix(nen1,n)
        if(ie(nie-1,ma).eq.-5) then
          nel = 0
          do i = nen,1,-1
            if(ix(i,n).gt.0) then
              nel = i
              exit
            endif
          end do ! i

          kk = ix(nel,n)
          do i = 1,nel
            ii = ix(i,n)
            if(ii.gt.0) then
              if(mo(ii).eq.-1) then
                x1 = xlm(7,1,ii)
                x2 = xlm(8,1,ii)
                x3 = xlm(9,1,ii)
                if(max(abs(x1),abs(x2),abs(x3)).eq.0.0d0) then
                  jj = ix(mod(i,nel)+1,n)
                  do j = 1,3
                    v1(j) = x(j,jj) - x(j,ii)
                    v2(j) = x(j,kk) - x(j,ii)
                  end do ! j
                  call vecp(v1,v2,v3)
                  vn = sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)
                  if(vn.gt.0.0d0) then
                    do j = 1,3
                      v3(j) = v3(j)/vn
                    end do ! j
                    vn = sqrt(norm(1,ii)**2+norm(2,ii)**2+norm(3,ii)**2)
                    if(vn.gt.0.0d0) then
                      do j = 1,3
                        v1(j) = norm(j,ii)/vn
                      end do ! j

                      vn = v1(1)*v3(1) + v1(2)*v3(2) + v1(3)*v3(3)

c                     Non-smooth shell use 6 dof/node

                      if(abs(vn).lt.0.95d0 .and. ndf.ge.6) then
                        mo(ii) = -5

c                     Smooth shell average the normal

                      else
                        do j = 1,3
                          norm(j,ii) = norm(j,ii) + v3(j)
                        end do ! j
                      endif

c                   First occurrance

                    else
                      do j = 1,3
                        norm(j,ii) = v3(j)
                      end do ! j
                    endif

c                 Error

                  else
                    write(iow,4000) ii,n
                    write(ilg,4000) ii,n
                    call plstop()
                  endif

                endif
              endif ! mo = -1
            endif
            kk = ii
          end do ! i

        endif ! is shell

      end do ! numel

c     Compute averaged normals

      do n = 1,numnp
        vn = sqrt(norm(1,n)**2+norm(2,n)**2+norm(3,n)**2)*2.d0
        if(vn.gt.0.0d0) then
          do j = 1,3
            xlm(j+6,1,n) = x(j,n) + norm(j,n)/vn
          end do ! j
        endif

      end do ! numnp

c     Formats

4000  format(' *ERROR* Shell normal computation. Node',i8,
     &       ' on element',i8,' has zero area.')

      end
