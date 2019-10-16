c$Id:$
      subroutine aslid2a(con,ib,ip,iq,ix, norm,
     &                   ma,nen,nen1,numnp,numel)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Auto surface descriptions in 2-d
c               Compute boundary patches and surface connection array

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_contac.h'
      include   'iofile.h'

      integer    ma,nen,nen1,numnp,numel
      integer    i,ii, j,jj, k, n,nel
      real*8     dotn

      integer    con(2,numnp)
      integer    ib(numnp),ip(numnp),iq(*),ix(nen1,numel)
      real*8     norm(3,numnp)

      call cdebug0 ('    aslid2a',-1)

c     1.) Build boundary element patches

      do n = 1,ip(numnp)
        iq(n) = 0
      end do ! n

      do n = 1,numel
        if(ma.eq.0 .or. ma.eq.ix(nen1,n)) then
          do i = 1,nen
            ii = ix(i,n)
            if(ii.gt.0) then
              if(ib(ii).eq.1) then
                if(ii.gt.1) then
                  jj = ip(ii-1) + 1
                else
                  jj = 1
                end if
                do j = jj,ip(ii)
                  if(iq(jj).eq.0 .or. iq(jj).eq.n) then
                    iq(jj) = n
                    go to 300
                  else
                    jj = jj + 1
                  end if
                end do ! j
                write(*,*) 'ERROR - AUTOCN: N,I,J',n,i,j
  300           continue
              end if
            end if
          end do ! i
        endif
      end do ! n

      if(ifdb) then
        call iprint(  iq,1,ip(numnp),1,'E-PATCH')
      endif

c     2.) Loop over patches to compute boundary connection array

      call pzeroi(con,2*numnp)
      jj  = 0
      do n = 1,numnp
        if(ib(n).gt.0) then

c         Form boundary connections for patch

          do j = jj+1,ip(n)

c           Get contribution from each element

            if(iq(j).gt.0) then

c             Determine number of nodes on element

              nel = 0
              do i = 1,nen
                if(ix(i,iq(j)) .gt.0) then
                  nel = i
                endif
              end do ! i

c             Find boundary segment

              if(nel.gt.2) then
                do i = 1,nen
                  ii = ix(i,iq(j))
                  if(ii.eq.n) then
                    k = ix(mod(i,nel)+1,iq(j))

c                 Check for active boundary node and normal possible

                    if(ib(k).gt.0) then
                      dotn = norm(1,n)*norm(1,k) + norm(2,n)*norm(2,k)
                      if(dotn .gt. -0.1d0) then
                        con(1,n) = k
                        con(2,n) = ix(nen1,iq(j))
                      endif
                    endif

                  endif
                end do ! i
              endif
            endif
          end do ! j
        end if

        jj = ip(n)

      end do ! n

      if(ifdb) then
        call iprint(con,2,numnp,2,'CONS?')
      endif

      end
