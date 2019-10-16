c$Id:$
      subroutine aslid3db(iq,con,ib,ip,ix,norm,
     &                    ma,nen,nen1,numnp,numel,ns)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Auto surface descriptions in 2-d
c               Set element face node numbers

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_contac.h'
      include   'iofile.h'
      include   'pointer.h'
      include   'comblk.h'

      integer             slidn
      common     /aslids/ slidn

      logical    flagng,flagch, setvar,palloc
      integer    ma,nen,nen1,numnp,numel,ns,nslid
      integer    i,ii,is, j,jj, k,kk, n,ne,nf, nfaces
      real*8     norm13, norm24

      integer    con(2,*),iq(numnp),ib(numnp),ip(numnp),ix(nen1,numel)
      integer    faces(4,6),ic(4)
      real*8     norm(3,numnp)

      save

      call cdebug0 ('    aslid3db',-1)

      data       faces / 1,4,3,2, 1,2,6,5, 2,3,7,6,
     &                   3,4,8,7, 4,1,5,8, 5,6,7,8 /

c     Set face/element numbers into nodes

      call pzeroi(con,2*ip(numnp))

      do n = 1,numel
        if(ma.eq.0 .or. ma.eq.ix(nen1,n)) then
          jj = 0
          do j = 1,min(8,nen)
            if(ix(j,n).gt.0) jj = jj + 1
          end do ! j

c         8-node brick elements

          if(jj.eq.8) then
            nfaces = 6
            do i = 1,nfaces
              jj = 0
              do j = 1,4
                ic(j) = ix(faces(j,i),n)
                if(ib(ic(j)).gt.0) jj    = jj + ib(ic(j))
              end do ! j
              if(jj.eq.4) then

c               Check normal to see if facet can be an interior surface!

                norm13 = norm(1,ic(1))*norm(1,ic(3))
     &                 + norm(2,ic(1))*norm(2,ic(3))
     &                 + norm(3,ic(1))*norm(3,ic(3))

                norm24 = norm(1,ic(2))*norm(1,ic(4))
     &                 + norm(2,ic(2))*norm(2,ic(4))
     &                 + norm(3,ic(2))*norm(3,ic(4))

                if(norm13 .gt. -0.25d0 .and. norm24 .gt. -0.25d0 ) then
                  do j = 1,4
                    if(ic(j).eq.1) then
                      ns = 1
                    else
                      ns = ip(ic(j)-1) + 1
                    endif
                    do k = ns,ip(ic(j))
                      if(con(1,k).eq.0) then
                        con(1,k) = n
                        con(2,k) = i
                        go to 100
                      end if
                    end do ! k
100                 continue
                  end do ! j
                endif
              endif
            end do ! i
          endif
        endif
      end do ! n

      if(ifdb) then
        call iprint(con,2,ip(numnp),2,'CONS?')
      endif

c     Number slide lines

      do n = 1, numnp
        iq(n) = - abs(ib(n))
      end do ! n

      ns = 0
      do n = 1, numnp

        if(iq(n).lt.0) then
          ns    = ns + 1
          iq(n) = ns
          do kk = 1,numnp
            flagch = .true.
            do ii = n,numnp
              if(iq(ii).eq.ns) then
                if(ii.eq.1) then
                  is = 1
                else
                  is = ip(ii-1) + 1
                endif
                flagng = .false.
                do k = is,ip(ii)
                  ne = con(1,k)
                  nf = con(2,k)
                  do j = 1,4
                    jj = ix(faces(j,nf),ne)
                    if(iq(jj).lt.0) then
                      flagng = .true.
                    end if
                  end do ! j
                end do ! k
                if(flagng) then
                  flagch = .false.
                  do k = is,ip(ii)
                    ne = con(1,k)
                    nf = con(2,k)
                    do j = 1,4
                      jj = ix(faces(j,nf),ne)
                      iq(jj) = ns
                    end do ! j
                  end do ! k
                endif
              endif
            end do ! ii
            if(flagch) go to 200
          end do ! kk
200       continue
        end if

      end do ! n

      if(ifdb) then
        call iprint(iq,1,numnp,1,'SLIDE-LINE NO.')
      endif

c     Set up facets

      nslid = 0
      do n = 1,numel
        jj = 0
        do j = 1,min(8,nen)
          if(ix(j,n).gt.0) jj = jj + 1
        end do ! j
        if(jj.eq.8) then

          nfaces = 6
          do i = 1,nfaces
            do j = 1,4
              ic(j) = ix(faces(j,i),n)
            end do ! j
            jj = min(iq(ic(1)),iq(ic(2)),iq(ic(3)),iq(ic(4)))
            kk = max(iq(ic(1)),iq(ic(2)),iq(ic(3)),iq(ic(4)))
            if(jj.eq.kk .and. jj.ne.0) then
              nslid  = nslid + 1
              setvar = palloc( 221,'ASLD2', 5*nslid, 1)
              call autoslid(mr(np(221)),ic,jj,nslid)
            endif
          end do ! i
        endif
      end do ! n

      if(ifdb) then
        write(iow,*) ' Number of slidelines/facets =',ns,nslid
        if(ior.lt.0) then
          write(*,*) ' Number of slidelines/facets =',ns,nslid
        endif
        call iprint(mr(np(221)),5,nslid,5,'SLIDE-LINES')
      endif

      slidn = nslid

      end
