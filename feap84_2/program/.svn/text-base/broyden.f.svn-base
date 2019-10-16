c$Id:$
      subroutine broyden(vtil,delx,v,numvec)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Separate 'id' & 'eq' on call to pload            27/04/2009
c       2. Add convergence check to exit                    03/04/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Broyden update for unsymmetric arrays.
c      Call order: > utan,,1
c                  > broy,"iter",vecs

c      Inputs:
c         vtil(*)   - Broyden vector set 1
c         delx(*)   - Broyden vector set 2
c         v(*)      - Broyden working vector
c         numvec    - Number vector iterations to perform

c      Outputs:
c                   - Solution through pointers
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'compas.h'
      include   'fdata.h'
      include   'hdatam.h'
      include   'idptr.h'
      include   'iofile.h'
      include   'ndata.h'
      include   'part0.h'
      include   'pfeapb.h'
      include   'pointer.h'
      include   'print.h'
      include   'prlod.h'
      include   'rdata.h'
      include   'rdat1.h'
      include   'sdata.h'
      include   'comblk.h'

      include   'p_int.h'

      logical    fa,tr
      integer    numvec,iter,k,j, toteq
      real*4     tary(2), etime, tt
      real*8     w,z,dot
      real*8     vtil(neq,numvec),delx(neq,numvec),v(neq)

      save

c     Set parameters

      if(pfeap_on) then
        toteq = numteq
      else
        toteq = neq
      endif

      fa = .false.
      tr = .true.

c     Iterate for updates

      do iter = 2,numvec

c       Allow for history updates

        hflgu  = tr
        h3flgu = tr

c       Get external load vector

        call pload(mr(np(31)+nneq),mr(id31),hr(np(40)),hr(np(30)),
     &             hr(np(26)),prop,tr,fa)

c       Get total negative residual in hr(np(26))

        call formfe(np(40),np(26),np(26),np(26),fa,tr,fa,fa,6,
     &              1,numel,1)

        rnorm = sqrt(dot(hr(np(26)),hr(np(26)),neq))
        if(ior.lt.0 .and. prt) then
          tt = etime(tary)
          write(*,2000) iter,rnorm,rnorm/rel0,tary
        endif

c       Check convergence

        if(abs(rnorm)/dble(toteq) .lt.                   ! Force meas.
     &     abs(rnorm1/rnormn)*sqrt(tol)*1.d-3) then
          write(*,*) ' --> Solution Converged'
          return
        endif


c       Set flags for u-symm solve & solve w result found in hr(np(26))

        fp(1) = na
        fp(2) = nau
        fp(3) = nal
        fp(4) = np(20+npart)
        call psolve(ittyp,hr(np(26)),fp,fa,tr,tr,prnt)

        do k=1,neq
          v(k) = -hr(np(26)+k-1)
        end do ! k

c       Broyden update computation first one

        if(iter.eq.2) then

          w = dot(v,delx(1,iter-1),neq)
          z = dot(delx(1,iter-1),delx(1,iter-1),neq)
          z = -1.d0/(w+z)
          do k = 1,neq
            vtil(k,iter-1) =   v(k)*z
            delx(k,iter)   = -(v(k) + w*vtil(k,iter-1))
          end do ! k

c       Broyden update computation generic form iter > 2

        else
          do j = 1,iter-2
            w = dot(v,delx(1,j),neq)
            do k = 1,neq
              v(k) = v(k) + vtil(k,j)*w
            end do ! k
          end do ! j

          w = dot(v,delx(1,iter-1),neq)
          z = dot(delx(1,iter-1),delx(1,iter-1),neq)
          z = -1.d0/(w+z)
          do k = 1,neq
            vtil(k,iter-1) = v(k)*z
            delx(k,iter)   = -(v(k) + w*vtil(k,iter-1))
          end do ! k

        endif

c       Actual update of solution vector

        call update(mr(id31),hr(np(30)),hr(np(40)),hr(np(42)),
     &              delx(1,iter),fl(9),2)

      end do ! iter

      if(ior.lt.0) then
        write(*,*) ' --> *WARNING*: BROYDEN did not convergence'
      endif

c     Format

2000  format(4x,'Iter. =',i3,' Resid. Norm =',1p,2e12.4,2x,
     &       't=',0p,2f9.2)

      end
