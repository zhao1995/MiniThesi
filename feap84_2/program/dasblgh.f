c$Id:$
      subroutine dasblgh(ix,id,ixl,igl,ige,p,s, r,hh,g)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Assemble all but block for passing to other program

c     Inputs:
c         ix(*)     - Element connection array
c         id(*,*)   - Active equation numbers
c         ixl(*)    - Local working array
c         igl(*)    - Local working array
c         ige(*)    - Local working array
c         p(*)      - Element right hand side
c         s(*,*)    - Element matrix

c     Outputs:
c         r(*)      - Reduced right hand side
c         hh(*,*)   - Diagonal block of reduced matrix
c         g(*,*,*)  - Cloupling block of reduced matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'eqslv.h'
      include  'eqsym.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'comblk.h'

      integer   i,ii,j,jj,jg,m1,mm,n,nn,eqm,eqn,numlg
      integer   ix(*),id(ndf,*),ixl(*),igl(*),ige(*)
      real*8    p(*),s(nst,*),r(*),hh(neqg,*),g(neq,neqg,2)

      save

c     Mark element export nodes and slave nodes on element

      do mm = 1,nen
        ige(mm) = 0
        igl(mm) = 0
        ixl(mm) = 0
      end do ! mm

      numlg = 0
      eqm  = 0
      do nn = 0,numg-1
        do mm = 1,nen
          if(ix(mm).eq.mr(np(205)+nn)) then
            numlg      = numlg + 1
            ige(numlg) = eqm
            igl(numlg) = mr(np(205)+nn)
            ixl(mm)    = 1
          endif
        end do ! mm
        do i = 1,ndf
          if(id(i,mr(np(205)+nn)).le.0) eqm = eqm + 1
        end do ! i
      end do ! nn

c     Loop through element nodes

      do mm = 1,nen

c       Check for slave node

        ii   = (mm-1)*ndf
        do m1 = 1,numlg

c         Assemble R

          if(igl(m1).eq.ix(mm) .and. ix(mm).gt.0) then

            eqm = ige(m1)
            do i = 1,ndf
              if(id(i,igl(m1)).le.0) then
                eqm    = eqm + 1
                r(eqm) = r(eqm) + p(ii+i)

                do nn = 1,nen
                  jj   = (nn-1)*ndf

c                 Current node is slaved

                  if(ixl(nn).ne.0) then

                    do n = 1,numlg

c                     Assemble H contributions from elements

                      if(igl(n).eq.ix(nn) .and. ix(nn).gt.0) then

                        eqn = ige(n)
                        do j = 1,ndf
                          if(id(j,igl(n)).le.0) then
                            eqn         = eqn + 1
                            hh(eqm,eqn) = hh(eqm,eqn) + s(ii+i,jj+j)
                          endif
                        end do ! j

                      endif
                    end do ! n

                  endif
                end do ! nn

              endif
            end do ! i

          endif

        end do ! m1

      end do ! mm

c     Loop through element nodes and assemble each G-vector

      do mm = 1,nen

c       Check for slave node on element

        ii   = (mm-1)*ndf
        do m1 = 1,numlg

          if(igl(m1).eq.ix(mm) .and. ix(mm).gt.0) then

            eqm = ige(m1)

c           Do individual RHS for G's

            do i = 1,ndf
              if(id(i,igl(m1)).le.0) then
                eqm    = eqm + 1

                do nn = 1,nen
                  jj   = (nn-1)*ndf

c                 Assemble G

                  if(ix(nn).gt.0) then
                    do j = 1,ndf
                      jg = id(j,ix(nn))
                      if(jg.gt.0) then
                        g(jg,eqm,1) = g(jg,eqm,1) + s(jj+j,ii+i)
                        g(jg,eqm,2) = g(jg,eqm,2) + s(ii+i,jj+j)
                      endif
                    end do ! j
                  endif

                end do ! nn

              endif ! Equation eqn
            end do ! i
          endif
        end do ! m1
      end do ! mm

      end
