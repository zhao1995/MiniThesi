c$Id:$
      subroutine xpline(x,ie,ix,id,ic,ip,numnp,numel,ndm,
     &                 nen1,nen,nie,ct,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change pstyp.gt.0 to pstyp.ne.0                  31/08/2008
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Determine plot line sequence to display 2-d mesh or
c               outline of 2-d mesh

c      Inputs:
c         x(ndm,*)  - Nodal coordinates for mesh
c         ie(nie,*) - Material set assembly data
c         ix(nen1,*)- Element nodal connection list
c         ip(8,*)   - Symmetry sorts for element sequences to plot
c         numnp     - Number of nodes in mesh
c         numel     - Number of elements in mesh
c         ndm       - Dimension of x array
c         nen1      - Dimension of ix array
c         nen       - Number of nodes/element
c         nie       - Dimension of ie array
c         ct        - Plot by material numbers if negative
c         isw       - Flag, plot mesh if true, otherwise do outline

c      Scratch:
c         id(*)     - Number of elements connected to nodes
c         ic(*)     - Element numbers connected to each node

c      Outputs:
c         none      - Plot to screen/file
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      include  'pbody.h'
      include  'pdatas.h'
      include  'pdatay.h'
      include  'pdatxt.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   ifl,iend,isw
      integer   numnp,numel,ndm,nen1,nen,nie,nume, pstyp,nel
      integer   i, j, k, ii, jj, ij, iju, n, n1, n2, ni, nn, nsy
      real*8    ct, x2, x3, tol, rtol

      integer   ie(nie,*),ix(nen1,*),ic(*),ip(8,numel),jplt(50),id(*)
      real*8    x(ndm,*)

      save

      data      rtol / 1.d-06/

c     Maximum connections to any node

      do nsy = 1,nsym
        lsym = isym(nsy)
        nume = nfac(lsym)
        call pltsym(x,ndm,numnp,lsym)

c       Initialize connection array

        do i = 1,id(numnp+1)
          ic(i) = 0
        end do ! i

c       Loop through elements to set up list

        do nn = 1,nume
          n  = ip(lsym,nn)
          if(ix(nen1-1,n).ge.nreg1 .and. ix(nen1-1,n).le.nreg2) then
            ii = ix(nen1,n)
            pstyp = ie(1,ii)

c           Plot material number: maplt (0 = all); ii > 0 active matl

            jj = maplt
            if(pstyp.ne.0 .and. (jj.eq.0 .or. ii.eq.jj)) then
              if(ii.eq.jj .and. ct.lt.0.0d0) call pppcol(jj,1)
              do i = nen,1,-1
                if(ix(i,n).gt.0) then
                  nel = i
                  exit
                endif
              end do ! i
              call plftyp(pstyp,nel,ie(nie-1,ii))
              call pltord(ix(1,n),ie(nie-1,ii), iju,jplt)

c             Look up element nodes

              ii = abs(ix(jplt(1),n))
              do ij = 2,iju
                j = jplt(ij)
                if((j.le.nen).and.(j.gt.0).and.(ix(j,n).ne.0)) then
                  jj = abs(ix(j,n))
                  if(jj.ne.ii) then
                    n1 = min(ii,jj)
                    n2 = max(ii,jj)
                    do k = id(n1),id(n1+1)-1
                      if(ic(k).eq.0) then
                        ic(k) =  n2
                        go to 100
                      elseif(abs(ic(k)).eq.n2) then
                        ic(k) = -abs(n2)
                        go to 100
                      endif
                    end do ! k
100                 ii = jj
                  endif
                endif
              end do ! ij
            endif
          endif
        end do ! nn

c       Change signs to permit mesh plot

        if(isw) then
          do n = 1,numnp
            do i = id(n),id(n+1)-1
              ic(i) = abs(ic(i))
            end do ! i
          end do ! n

c       Check for removal of symmetry conditions

        else
          fp(1) = npty - 1
          do ii = 1,ndm
            if(isymm(ii,1).gt.0) then
              x2 = x(ii,1)
              x3 = x(ii,1)
              do n = 1,numnp
                if(mr(fp(1)+n).ge.0 ) then
                  x2 = min(x2,x(ii,n))
                  x3 = max(x3,x(ii,n))
                endif
              end do ! n
              tol = rtol*(x3-x2)
              do n = 1,numnp
                if(abs(x(ii,n)-xsyc(ii)).lt.tol) then
                  do i = id(n),id(n+1)-1
                    if(ic(i).gt.0) then
                      if(abs(x(ii,ic(i))-xsyc(ii)).lt.tol) then
                        ic(i) = - abs(ic(i))
                      endif
                    endif
                  end do ! i
                endif
              end do ! n
            endif
          end do ! ii
        endif

c       Plot outline of part with continuous lines

        x3 = 0.0d0
        do ni = 1,numnp
          iend = .true.
          do n = 1,numnp
            ifl = .true.
            n1  =  n
101         do i = id(n1),id(n1+1)-1
              if(ic(i).gt.0) then
                go to 102
              elseif(ic(i).eq.0) then
                go to 103
              endif
            end do ! i
            go to 103
102         iend = .false.
            if(ifl) then
              if(ndm.ge.3) x3 = x(3,n1)
              call plotl(x(1,n1),x(2,n1),x3,3)
              ifl = .false.
            endif
            n2    =  ic(i)
            ic(i) = -n2
            if(ndm.ge.3) x3 = x(3,n2)
            call plotl(x(1,n2),x(2,n2),x3,2)
            n1 = n2
            go to 101
103         continue
          end do ! n
          if(iend) go to 104
        end do ! ni

104     call pltsym(x,ndm,numnp,lsym)

      end do ! nsy

      end
