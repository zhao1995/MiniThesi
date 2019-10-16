c$Id:$
      subroutine arclen(u,du,u1,u2,f,id,ndf,time)

c     * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Save aengy to prevent print error of convergence 17/10/2010
c       2. Change: r = -sign(1,r) proposed by Stanciulescu  11/09/2012
c          for complex load curves (instead of r = rr)
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:     Perform arc length solution to find limit states

c      Programmed - 08/19/86: Peter Wriggers
c      Modified   - 10/25/94: Sanjay Govindjee
c                   12/20/95: R.L. Taylor

c      Inputs:
c         u(*)    - Current solution
c         du(neq) - Solution increment from last solution
c         u1(*)   - Displacement before first iteration (N/R) ! np(84)
c         u2(*)   - Displacement for load level 1.0 (N/R)     ! np(85)
c         f(*)    - Load vector (N/R)
c         id(*)   - Equation numbers for each degree-of-freedom
c         ndf     - Number dof/node
c         time    - Time

c      Outputs:
c         u(*)   - Converged solution

c      Scratch:
c         u1(neq) - Working vector space
c         u2(neq) - Working vector space

c      Input record:
c        arclen,,kflag

c      Solution options:
c        kflag    = 0: Mod.  N/R
c                 = 1: Updated normal plane iteration Mod.  N/R
c                 = 2: Orig. N/R
c                 = 3: Updated normal plane iteration Orig. N/R
c                 = 4: Displacement control Mod.  N/R
c                 = 5: Displacement control Orig. N/R
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'arclei.h'
      include  'arcler.h'
      include  'cdata.h'
      include  'compas.h'
      include  'counts.h'
      include  'endata.h'
      include  'eqsym.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'ndata.h'
      include  'part0.h'
      include  'pointer.h'
      include  'print.h'
      include  'prlod.h'
      include  'pscal.h'
      include  'tdata.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   udis
      integer   i,j,jj,ndf
      integer   id(ndf,*)
      real*8    time, cs,ds,dcs,damp,dot,rr,rat, aengysv
      real*8    alfa,alfa1,alfold, dtold,tol, undis
      real*8    u(ndf,*),du(*),u1(*),u2(*),f(ndf,*)

      save

      data      tol /1.d-08/, dtold /0.0d0/

c     Exit on near zero proportional loading

      if( abs(prop).lt.1.0d-10) return

c     Initialization at beginning of load step

      if(time.ne.timold) then

        ite    = 0
        timold = time

c       Calculate arc length at beginning of time step

        ds = sqrt(dot(du,du,neq))

c       Accumulate nodal forces into 'u1'

        do i = 1,neq
          u1(i) = 0.0d0
        end do ! i

        do j = 1,ndf
          if(ndfp(j).eq.npart) then
            do i = 1,numnp
              jj = id(j,i)
              if (jj.gt.0) u1(jj) = u1(jj) + f(j,i)*prop
            end do ! i
          endif
        end do ! j

c       Accumulate surface loads into 'u1'

        call ploads(u,u1,prop,.false.,.false.)

c       Accumulate element loads into 'u1'

        call formfe(np(40),np(84),np(84),np(84),.false.,.true.,
     &             .false.,.false.,23,1,numel,1)

c       Check for limit points

        rr     = dot(du,u1,neq)/prop

c       Set arc length constant (in first time step only)
c       Safe energy value at beginning of loading

c       First load/time step only

        if (nstep.le.1) then
          r     = 1.0d0
          ds0   = ds
          c0    = rr
          dtold = dt
          undis = alfa0
        endif

c       Rescale arclength for changing time steps

        if(dt.ne.dtold.and.dtold.ne.0.0d0) then
          ds0   = ds0*abs(dt/dtold)
          dtold = dt
        end if

c       Current stiffness parameter of current load increment

        cs = c0/rr

c       Change sign if stiffness param passes infinity or limit point

        if (nstep.gt.2) then
          rr = sign(1.d00,rr)
          if (cs02.gt.0.0d0.and.cs.lt.0.0d0) then
            dcs = cs02 - cs01
c           if (dcs.lt.0.0d0) r = rr
            if (dcs.lt.0.0d0) r = -sign(1.d0,r)
          elseif (cs02.lt.0.0d0.and.cs.gt.0.0d0) then
            dcs = cs02 - cs01
c           if (dcs.gt.0.0d0) r = rr
            if (dcs.gt.0.0d0) r = -sign(1.d0,r)
          endif
        endif

c       Save current stiffness values

        if (nstep.eq.1) then
          cs01 = cs
        else
          cs01 = cs02
        endif
        cs02 = cs

c       Constant arc length method

        if (kflag.le.3.or.kflag.eq.6) then

c         Calculate reducing factor and change direction for limit point

          alfa0 = ds0/ds * r

c         Actual load level

          rlnew = rlnew + alfa0 - dt

c         Save displacement vector for load step

          do i = 1,neq
            du(i) = du(i)*alfa0
            u1(i) = du(i)
          end do ! i

          if(prnt) then
            write (iow,2001) time,ds,alfa0
            write (iow,2002) rlnew*prop
            write (iow,2004) cs,cs01,cs02
            if(ior.lt.0.and.pfr) then
              write (*,2001) time,ds,alfa0
              write (*,2002) rlnew*prop
              write (*,2004) cs,cs01,cs02
            endif
          endif

c       Displacement control: Calculate displacements for load level 1.0

        elseif (kflag.eq.4.or.kflag.eq.5) then

          do i = 1,neq
            u1(i) = u1(i)/prop
          end do ! i

c         Load vector zero case

          if(dot(u1,u1,neq).eq.0.0d0) then
            udis  = .false.
            alfa  = 0.0d0
            alfa1 = 0.0d0

c         Load vector non-zero case

          else
            udis    = .true.
            fp(1)   = na
            fp(2)   = nau
            fp(3)   = nal
            fp(4)   = np(20+npart)
            aengysv = aengy   ! Save energy to avoid convergence error
            call psolve(ittyp,u1,fp,.false.,.true.,.true.,prnt)
            aengy   = aengysv

c           Load factor for displacement control

            ndis  = id(nddis,nodis)
            alfa1 = (alfa0 - du(ndis))/u1(ndis)
            alfa  = alfa1
          endif

c         Update displacements

          do i = 1,neq
            du(i) = du(i) + alfa1*u1(i)
          end do ! i

c         Update load level

          rlnew = rlnew + alfa
          if(prnt) then
            write (iow,2005) alfa*prop,rlnew*prop
            write (iow,2004) cs,cs01,cs02
            if(ior.lt.0) then
              if(pfr) then
                write (*,2005) alfa*prop,rlnew*prop
                write (*,2004) cs,cs01,cs02
              else
                write (*,2002) rlnew*prop
              endif
            endif
          else
            write (iow,2002) rlnew*prop
          endif
        endif

c     Iteration starts

      else

        ite = ite + 1

c       Rescale arclength for changing time steps

        if(dt.ne.dtold.and.dtold.ne.0.0d0) then
          ds0   = ds0*abs(dt/dtold)
          dtold = dt
        end if

        if(kflag.gt.1.and.kflag.ne.4) then

c         Calculate Newton displacement for load level 1.0

          do i = 1,neq
            u2(i) = 0.0d0
          end do ! i

          do j = 1,ndf
            if(ndfp(j).eq.npart) then
              do i = 1,numnp
                jj = id(j,i)
                if (jj.gt.0) u2(jj) = u2(jj) + f(j,i)*prop
              end do ! i
            endif
          end do ! j
          call ploads(u,u2,prop,.false.,.false.)
          call formfe(np(40),np(85),np(85),np(85),.false.,.true.,
     &               .false.,.false.,23,1,numel,1)

          if( dot(u2,u2,neq).eq.0.0d0 ) then
            udis     = .false.
          else
            udis     = .true.
          endif
          fp(1)   = na
          fp(2)   = nau
          fp(3)   = nal
          fp(4)   = np(20+npart)
          aengysv = aengy   ! Save energy to avoid convergence error
          call psolve(ittyp,u2,fp,.false.,.true.,.true.,prnt)
          aengy   = aengysv
          do jj = 1,neq
            u2(jj) = u2(jj)/prop
          end do ! jj
        endif

c-----[--.----+----.----+----.-----------------------------------------]
c       Update 1. Modified and original Newton Raphson
c-----[--.----+----.----+----.-----------------------------------------]

c       Calculate reducing factor by iteration on normal plane
c         (unscaled displ. vector)

        if (kflag.eq.0) then

          alfa  =  dot(du,u1,neq)
          alfa  = -alfa*alfa0/(ds0*ds0)

c       Iteration on updated normal plane

        elseif (kflag.eq.1) then

          if (ite.eq.1) then
            do i = 1,neq
              u2(i) = u1(i)
            end do ! i
          endif
          alfa  = dot(du,u2,neq)
          alfa1 = dot(u1,u2,neq)
          alfa  = alfa/alfa1

c       Single displacement control

        elseif (kflag.eq.4) then

          ndis  =  id(nddis,nodis)
          alfa  = -du(ndis)/u1(ndis)

c       Displacement control

        elseif (kflag.eq.5) then

          ndis = id(nddis,nodis)
          if(udis) then
            alfa = -du(ndis)/u2(ndis)
          else
            alfa =  0.0d0
            do j = 1,ndf
              if(ndfp(j).eq.npart) then
                do i = 1,numnp
                  if(id(j,i).eq.ndis) then
                    alfa = u(i,j)
                  endif
                end do ! i
              endif
            end do ! j
            alfa  = - du(ndis)/alfa
          endif

c       Arc length (unscaled displacement vector)

        else

          alfa  = dot(du,u1,neq)
          alfa1 = dot(u1,u2,neq)
          if(abs(alfa1).gt.tol*abs(alfa)) then
            alfa = -alfa/alfa1    ! previously, only possibility
          else
            alfa = -alfa/(alfa1 + sign(tol,alfa1))
          endif
        endif

c       Numerical damping if sign of alfa is varying

        damp = 1.d0
        if (ndamp.eq.0) then
          if (ite.eq.1) alfold = alfa
          if(abs(alfold).gt.tol) then
            rat = alfa/alfold
          else
            rat = alfa/(alfold + tol)
          endif
          if (-1.0d0.lt.rat.and.rat.lt.0.0d0) damp = 0.5d0
          alfa   = alfa*damp
          alfold = alfa
        endif

c       Update new load level

        rlnew = rlnew + alfa
        if(prnt) then
          write (iow,2003) ite,alfa*prop,rlnew*prop,damp
          if(ior.lt.0) then
            if(pfr) then
              write (*,2003) ite,alfa*prop,rlnew*prop,damp
            else
              write (*,2002) rlnew*prop
            endif
          endif
        else
          write (iow,2002) rlnew*prop
        endif

c       Update displacements in iteration step

        if(kflag.le.1.or.kflag.eq.4) then
          do i = 1,neq
            du(i) = du(i)*damp + u1(i)*alfa
          end do ! i

c         Updated normal plane iteration

          if (kflag.eq.1) then
            do i = 1,neq
              u2(i) = u1(i) + du(i)
            end do ! i
          endif

c       Update displacements in iteration step

        else

          do i = 1,neq
            du(i) = u2(i)*alfa + du(i)*damp
          end do ! i

c         Updated normal plane iteration

          if (kflag.eq.3) then
            do i = 1,neq
              u1(i) = du(i) + u1(i)
            end do ! i

c         New tangent plane

          elseif (kflag.eq.6) then
            do i = 1,neq
              u1(i) =  du(i)
            end do ! i
          endif
        endif

      endif

c     Formats

 2001 format(/3x,'A r c   L e n g t h   M e t h o d'//3x,'Time =',f12.2,
     &           ' Arc Length = ',g12.5,' Red. factor = ',g12.5)

 2002 format( 3x,'Load level = ',g12.5)

 2003 format(/3x,'Iteration:',i5,' Reduct.factor = ',g12.5,
     &           ' Load level = ',g12.5/ 3x,'Damp. fact = ',g12.5)

 2004 format( 3x,'Cs-param.  = ',g12.5,/,3x,'Cs01       = ',g12.5,/,
     &        3x,'Cs02       = ',g12.5)

 2005 format( 3x,'Displacement Control Parameters',/,
     &        3x,'Reduction factor  = ',g12.5,' Load level = ',g12.5)

      end
