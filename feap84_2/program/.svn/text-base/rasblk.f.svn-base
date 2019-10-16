c$Id:$
      subroutine rasblk(s,p,ld,ix,rixt,rlink,x,aufl,bfl,jp,b,al,au,ad)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Modify and assemble rigid dof blocks

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'eldata.h'
      include  'eqsym.h'
      include  'sdata.h'

      logical   bfl,aufl, mmfl,nnfl
      integer   i,j, mm,nn, i1,j1, slavem,slaven, masterm,mastern
      integer   ld(*),ix(*),rixt(*),rlink(ndf,*),jp(*)
      real*8    ti(3,3),tj(3,3), s11(3,3),s12(3,3),s21(3,3),s11tj(3,3)
      real*8    s(nst,*), p(*), x(ndm,*),b(*),al(*),au(*),ad(*)

      save

c     3-D case with 6 degrees of freedom

      if(ndm.eq.3 .and. ndf.ge.6) then

        i1 = 0
        do mm = 1,nel

          do i = 1,3
            ti(i,1) = 0.0d0
            ti(i,2) = 0.0d0
            ti(i,3) = 0.0d0
          end do ! i

          if(ix(mm).gt.0) then
            slavem  = ix(mm)
            masterm = rixt(slavem)
            if(masterm.gt.0) then
              ti(1,2) = -(x(3,masterm) - x(3,slavem))
              ti(1,3) =  (x(2,masterm) - x(2,slavem))
              ti(2,3) = -(x(1,masterm) - x(1,slavem))
              ti(2,1) = -ti(1,2)
              ti(3,1) = -ti(1,3)
              ti(3,2) = -ti(2,3)
              do j = 1,3
                if(rlink(j+ndm,slavem).ne.0) then
                  ti(1,j) = 0.0d0
                  ti(2,j) = 0.0d0
                  ti(3,j) = 0.0d0
                end if
              end do ! j

c             Transform residual coefficients

              if(bfl) then
                do i = 1,3
                  p(i+i1+ndm) = p(i+i1+ndm) + ti(1,i)*p(1+i1)
     &                                      + ti(2,i)*p(2+i1)
     &                                      + ti(3,i)*p(3+i1)
                end do ! i

              endif ! bfl
              mmfl = .true.
            else
              mmfl = .false.
            endif

c           Transform stiffness coefficients

            if(aufl) then
              j1 = 0
              do nn = 1,nel

                do i = 1,3
                  tj(i,1) = 0.0d0
                  tj(i,2) = 0.0d0
                  tj(i,3) = 0.0d0
                end do ! i

                if(ix(nn).gt.0) then
                  slaven  = ix(nn)
                  mastern = rixt(slaven)
                  if(mastern.gt.0) then
                    tj(1,2) = -(x(3,mastern) - x(3,slaven))
                    tj(1,3) =  (x(2,mastern) - x(2,slaven))
                    tj(2,3) = -(x(1,mastern) - x(1,slaven))
                    tj(2,1) = -tj(1,2)
                    tj(3,1) = -tj(1,3)
                    tj(3,2) = -tj(2,3)
                    do j = 1,3
                      if(rlink(j+ndm,slaven).ne.0) then
                        tj(1,j) = 0.0d0
                        tj(2,j) = 0.0d0
                        tj(3,j) = 0.0d0
                      end if
                    end do ! j
                    nnfl = .true.
                  else
                    nnfl = .false.
                  endif

c                 Extract stiffness coefficients

                  do j = 1,3
                    do i = 1,3
                      s11(i,j) = s(i+i1    ,j+j1    )
                      s12(i,j) = s(i+i1    ,j+j1+ndm)
                      s21(i,j) = s(i+i1+ndm,j+j1    )
                    end do ! i
                  end do ! j

                  do j = 1,3
                    do i = 1,3
                      s11tj(i,j) = s11(i,1)*tj(1,j)
     &                           + s11(i,2)*tj(2,j)
     &                           + s11(i,3)*tj(3,j)
                    end do ! i
                  end do ! j

                  if(mmfl) then
                    do j = 1,3
                      do i = 1,3
                        s(i+i1+ndm,j+j1    ) = s(i+i1+ndm,j+j1    )
     &                                       + ti(1,i)*s11(1,j)
     &                                       + ti(2,i)*s11(2,j)
     &                                       + ti(3,i)*s11(3,j)
                        s(i+i1+ndm,j+j1+ndm) = s(i+i1+ndm,j+j1+ndm)
     &                                       + ti(1,i)*s12(1,j)
     &                                       + ti(2,i)*s12(2,j)
     &                                       + ti(3,i)*s12(3,j)
                      end do ! i
                    end do ! j
                  endif
                  if(nnfl) then
                    do j = 1,3
                      do i = 1,3
                        s(i+i1    ,j+j1+ndm) = s(i+i1    ,j+j1+ndm)
     &                                       + s11tj(i,j)
                        s(i+i1+ndm,j+j1+ndm) = s(i+i1+ndm,j+j1+ndm)
     &                                       + s21(i,1)*tj(1,j)
     &                                       + s21(i,2)*tj(2,j)
     &                                       + s21(i,3)*tj(3,j)
                      end do ! i
                    end do ! j
                    if(mmfl) then
                      do j = 1,3
                        do i = 1,3
                          s(i+i1+ndm,j+j1+ndm) = s(i+i1+ndm,j+j1+ndm)
     &                                         + ti(1,i)*s11tj(1,j)
     &                                         + ti(2,i)*s11tj(2,j)
     &                                         + ti(3,i)*s11tj(3,j)
                        end do ! i
                      end do ! j
                    endif
                  endif
                endif
                j1 = j1 + ndf
              end do ! nn
            endif ! aufl
          endif
          i1 = i1 + ndf
        end do ! mm

      else

        write(*,*) ' *WARNING* No master-slave modify performed'

      end if ! ndm & ndf test

c     Perform assembly into global arrays

      call dasble(s,p,ld,jp,nst,neqs,aufl,bfl,b,al,au,ad)

      end
