c$Id:$
      subroutine pout3d(iface,ip,elface,x,xl,norm,econ,
     &                  n1,cosa,nen,ndm,nface,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Plot edge definitions for 3-d objects

c      Inputs:
c         iface(7,*)- Face nodal connections, materials, regions
c         ip(8,*)   - Sort information for symmetry plots
c         elface(*) - Pointer for element faces
c         x(ndm,*)  - Nodal coordinates
c         xl(ndm,*) - storage for face coordinates
c         norm(3,*) - Normal vectors to facets
c         econ(*)   - Faces connected to each node
c         n1        - Color for face plot
c         cosa      - Value of angle cosine for node to face normal
c         nen       - Dimension of xl array
c         ndm       - Dimension of x and xl arrays
c         nface     - Number of faces
c         numnp     - Number of nodes in mesh

c      Outputs:
c         none      - Output is plot to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pbody.h'
      include  'pdata2.h'
      include  'pdatay.h'

      integer   n1,nen,ndm,nface,numnp
      integer   ma,nsy,nume,nn,ii,jj,icol, i,j,k,n,nel
      real*8    ang, dot, cosa

      integer   iface(7,nface),ip(8,*),elface(*),econ(*)
      real*8    x(ndm,numnp),xl(ndm,nen),norm(3,nface)

      save

      if(cosa.le.0.0d0) cosa = 0.96d0

      do nsy = 1,nsym
        lsym = isym(nsy)
        call pltsym(x,ndm,numnp,lsym)
        nume = nfac(lsym)
        do nn = 1,nume

          n  = ip(lsym,nn)

c         Check active material set

          if(iface(6,n).ge.nreg1 .and. iface(6,n).le.nreg2) then

            ma = iface(7,n)
            if(maplt.eq.0 .or. ma.eq.maplt) then

c             Set panel color

              if(n1.gt.0) then
                icolr = n1
              elseif(n1.eq.0) then
                icolr = ma + 1
              else
                icolr = -1
              endif
              call pppcol(icolr,0)

c             Set coordinates for face

              do i = 1,4
                ii = iface(i,n)
                if(ii.gt.0) then
                  nel = i
                  do j = 1,ndm
                    xl(j,i) = x(j,ii)
                  end do ! j
                else
                  do j = 1,ndm
                    xl(j,i) = 0.0d0
                  end do ! j
                endif
              end do ! i

c             Check for line elements

              if(iface(1,n).eq.iface(4,n) .and.
     &           iface(2,n).eq.iface(3,n)) nel = 2

              call plot9(1,iface(1,n),xl,ndm,nel,-1)

c             Plot boundaries in appropriate color

              icol = 5
              call pppcol(icol,0)

              do i = 1,4

                j = mod(i,4) + 1
                do ii = 1,ndm
                  xl(ii,1) = x(ii,iface(i,n))
                  xl(ii,2) = x(ii,iface(j,n))
                end do ! ii

c               Find edges

                ii = iface(i,n)
                if(ii.gt.0) then

                  do jj = elface(ii)+1,elface(ii+1)
                    if(econ(jj).gt.0 .and. econ(jj).ne.n) then

                      do k = 1,4
                        if(iface(k,econ(jj)).eq.iface(j,n)) then

                          ang = abs(dot(norm(1,n),norm(1,econ(jj)),3))
                          if(ang.lt.cosa) then

c                           Draw Line

                            call plotl(xl(1,1),xl(2,1),xl(3,1),3)
                            call plotl(xl(1,2),xl(2,2),xl(3,2),2)
                            go to 200

                          end if
                        end if
                      end do ! k

                    end if
                  end do ! jj
                end if
200             continue
              end do ! i

            end if
          end if

        end do ! nn

        call pltsym(x,ndm,numnp,lsym)

      end do ! nsy

      end
