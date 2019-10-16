c$Id:$
      subroutine pmacr2(lct,ct,j)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Set length on assignment to uct                  17/04/2007
c       2. Echo solution step number to screen.             31/10/2007
c          Flush log file data to permit monitoring.
c       3. Change 'Macro' to 'Command' in format            07/01/2008
c       4. Call 'updrot' for static problems with isw = 1   18/11/2008
c       5. Set tranfl (transient) and intlfl (initial)      01/02/2009
c       6. Add funciton pflush to do flushes                25/02/2009
c       7. Separate id31 and np(31) on profpart call        29/04/2009
c       8. Reduce counter value on auto dt checks           25/11/2009
c       9. Add set of 'convfl' for convergence outputs      16/02/2010
c      10. Add reset of partition number at line            15/03/2010
c      11. Add 'tolr' set and use in residual convergence   10/08/2010
c      12. Add  common /cfreq/ and set omega1,omega2        18/08/2010
c      13. Remove 'and' dt=0 from ljump check on next       18/03/2011
c      14. Add initial acceleration sets                    24/04/2011
c      15. Change 'i' to 'n' in [newf loop on 'n'           04/05/2013
c      16. Set 'flconv' to export convergence result        06/05/2013
c      17. Add 'loop,rve' option                            22/05/2013
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Command language instruction subprogram: Part 2

c      Inputs:
c         lct(*)     - Command option
c         ct(3,*)    - Command parameters
c         j          - Command number in this routine

c      Outputs:
c         Depends on command number j
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'arclel.h'
      include  'arcler.h'
      include  'augdat.h'
      include  'auto1.h'
      include  'auto2.h'
      include  'cdata.h'
      include  'cfreq.h'
      include  'chdata.h'
      include  'comfil.h'
      include  'compas.h'
      include  'corfil.h'
      include  'counts.h'
      include  'crotas.h'
      include  'ddata.h'
      include  'debugs.h'
      include  'elauto.h'
      include  'endata.h'
      include  'evdata.h'
      include  'fdata.h'
      include  'hdatam.h'
      include  'idptr.h'
      include  'iofile.h'
      include  'ldata.h'
      include  'mxsiz.h'
      include  'ndata.h'
      include  'part0.h'
      include  'part1.h'
      include  'part7.h'
      include  'pconstant.h'
      include  'pointer.h'
      include  'print.h'
      include  'prflag.h'
      include  'prld1.h'
      include  'prlod.h'
      include  'rdata.h'
      include  'rdat0.h'
      include  'rdat1.h'
      include  'sdata.h'
      include  'tdata.h'
      include  'tdato.h'
      include  'umac1.h'
      include  'comblk.h'

      include  'setups.h'

      logical   pcomp,err, errck, tinput, lexpr,setvar,palloc,oprt
      logical   vinput,crflg
      integer   i,j,n,larept,nparta, k1, npl(2)
      real*4    etime, tary(2)
      real*8    propld, dtnew,dtsav, v1,v2,v3,v4, ct(3,*),ctl(25)
      character lct(*)*15,lctl(2)*15, tname*21, tmeth(2)*6, expr*80

      save

      data      tmeth/ 'Newton','Secant'/

c     Transfer to correct process

      go to (1,2,3,4,5,6,7,8,8,10,11,12,13,14,15,16,17,18,19,20,
     &       21,22,10,24,25,26,27,28), j

c     Set solution tolerance
c     [tol,,eval,rval]     - set relative convergence tolerance to tol
c     [tol,ener,value]     - set absolute convergence tolerance to enzer
c     [tol,emax,value]     - set comparison value of convergence to emax
c     [tol,iter,itol,atol] - set residual & absolute convergence tols

1     if(pcomp(lct(l),'ener',4)) then

c       Set energy assumed zero value

        enzer = abs(ct(1,l))
        write(iow,2008) enzer
        if(ior.lt.0 .and. prnt) write(*,2008) enzer

      elseif(pcomp(lct(l),'emax',4)) then

c       Set convergence limit

        if(ct(1,l).ne.0.0d0) then
          rnmax = ct(1,l)
          write(iow,2007) rnmax
          if(ior.lt.0 .and. prnt) write(*,2007) rnmax
        end if

      elseif(pcomp(lct(l),'iter',4)) then

c       Set residual and absolute convergence limit

        itol            = ct(1,l)
        nitolp(npart)   = itol
        if(ct(2,l).ne.0.0d0) then
          atol          = ct(2,l)
          natolp(npart) = atol
        endif
        if(ct(3,l).ne.0.0d0) then
          dtol          = ct(3,l)
        endif
        err = vinput(lzz(l),80,ctl,4)
        if(nint(ctl(4)).gt.0) then
          pmaxit = nint(ctl(4))
        endif

        write(iow,2022) itol,atol,dtol
        if(pmaxit.gt.0) write(iow,2023) pmaxit
        if(ior.lt.0 .and. prnt) then
          write(*,2022) itol,atol,dtol
          if(pmaxit.gt.0) write(*,2023) pmaxit
        endif

      else

c       Set normal tolerance value

        if(ct(1,l).ne.0.0d0) then
          tol  = ct(1,l)
          if(ct(2,l).ne.0.0d0) then
            tolr = ct(2,l)
          else
            tolr = sqrt(tol)*100.d0
          endif
          write(iow,2013) tol,tolr
          if(ior.lt.0 .and. prnt) write(*,2013) tol,tolr
        else
          write(iow,2021) tol
          if(ior.lt.0 .and. prnt) write(*,2021) tol
        endif

      endif

      return

c     Set time increment
c     [dt,,value]

2     dtold = dt
      dt    = ct(1,l)
      return

c     Set loop start indicators
c     [loop,,number]

3     lv            =  lv + 1
      lvs(lv)       =  l
      lve(lv)       =  nint(ct(2,l))
      ct(1,lve(lv)) =  1.d0
      ct(3,lve(lv)) =  1.d0
      if(lv.gt.1) then
        floop(2)      = .true.
      endif
      flncon   = .false.
      flconv   = .false.


c     Check for infinite looping

      if(pcomp(lct(l),'infi',4)) then
        if(lv.le.2) then
          if(ct(1,l).le.0.0d0) then
            ct(1,l) = 1.0
          endif
        else
          write(ilg,4003) lv-1
          write(iow,4003) lv-1
          call plstop()
        endif
      endif

c     Auto dt adjustment

      if(autofl .and. lvauto.eq.lv ) then
        larept = 0
      end if

c     Set problem type if not outer default loop.

      if(l.ne.1) then
        if(linear .and. prnt) then
          if(ior.lt.0) write(  *,2006)
          write(iow,2006)
        end if
        linear  = .false.
      end if
      return

c     Loop terminator control
c     [next]

4     n       = nint(ct(2,l))   ! Location of matching LOOP
      ct(1,l) = ct(1,l) + 1.0d0
      ct(3,l) = ct(3,l) + 1.0d0

c     Auto dt adjustment

      if(autofl .and. lvauto.eq.lv ) then
        if(ct(1,l).le.ct(1,n)) then

c         Solution is diverging, reduce and resolve (NEEDS BETTER ALGO)

          if (autr(1).gt.1000.d0*autr0 .and. niter.ge.3) then

            call autdt(dt,dtnew, .false.)
            call autbac(dtnew)
            ct(1,l) = 1.0d0
            ct(3,l) = 1.0d0
            larept  = larept + 1
            iautl   = 5

          end if

c         Rattle flag turned off when residual drops by 1.d-04

          if ( aratfl .and. autr(1) .lt. autr0*1.d-4 ) then
            aratfl  = .false.
            ct(3,l) = iaopt
          end if

          l = n

        else
          if(larept.ge.iarept .and. .not.autcnv) then
            if(ior.lt.0)
     &        write(*,*) 'BYE - there is something really wrong!'
            write(iow,*) 'BYE - there is something really wrong!'
            call plstop()

c         NO CONVERGENCE: Need to reduce time step and redo solution

          elseif(ct(3,l).gt.ct(1,n) .and. .not.autcnv) then

c           First Check for rattling residual before reducing time step;
c           if residual is rattling, then increase time step.

            if(ct(3,l) .gt. 4) then
              v1 = autr(1)-autr(3)
              v2 = autr(2)-autr(4)
              v3 = autr(1)-autr(2)+autr(3)-autr(4)
              v4 = v3*v3/(v1*v1+v2*v2)
              v4 = v1*v1+v2*v2
              if(v4.eq.0.d0) then
                 v4 = 401.d0
              else
                 v4 = v3*v3/v4
              end if
              ct(1,l) = 1.0d0
              ct(3,l) = 1.0d0
              if(ior.lt.0) write(*,*) 'Rattle Check: ',v4
              write(iow,*) 'Rattle Check: ',v4
              v4=399.d0!SW:added this line
              if(v4 .gt. 400.d0) then
                if(aratfl) then
                 if(ior.lt.0) write(*,*) 'BYE - too many rattles'
                 write(iow,*) 'BYE - too many rattles'
                 call plstop()
                end if
                aratfl = .true.
                ct(1,l) = -5.d0
                ct(3,l) = 1.0d0
              else
                aratfl = .false.
              end if
            else
              ct(1,l) = 1.0d0
              ct(3,l) = 1.0d0
              aratfl = .false.
            end if

            call autdt(dt,dtnew, .false.)
            call autbac(dtnew)
            l       = n
            larept  = larept + 1

          elseif(autcnv) then

c           CONVERGENCE: Need to reduce time step and do next solution

            aratfl = .false.
            if(int(ct(3,l))-1.gt.int(afact)) then

              call autdt(dt,dtnew, .false.)
              dt = dtnew

c           CONVERGENCE: Need to increase time step and do next solution

            elseif(int(ct(3,l))-1.lt.iaopt) then

              iautl = iautl - 2
              if(iautl.lt.0) then
                call autdt(dt,dtnew, .true.)
                dt = dtnew
                iautl = 0
              endif

            end if

c           Terminate current loop

            floop(2) = .false.
            lv       =  lv -1

          end if
        end if

      elseif(autofli .and. lvautoi.eq.lv ) then

        if(ct(1,l).le.ct(1,n)) then

           if (rmeas .gt. 2.0d0) then   ! Iterations too large, backup
!             write(*,*) 'rmeas.gt.2.d0'
            write(iow,*) 'Old dt',dt
            dtnew = dt*0.85d0/rmeas
            write(iow,*) 'New dt',dtnew
            call autbac(dtnew)
            ct(1,l) = 1.0d0          ! Reset iteration counter
            ct(3,l) = 1.0d0          ! Reset iteration counter
           end if

           l = n                      ! Goto top of loop

        else

          if(autcnv .and. rmeas.le.1.25d0) then
            write(iow,*) 'Old dt',dt
            if(rmeas.le.0.5d0) then
              dt = dt*1.5d0
            elseif(rmeas.le.0.8d0) then
              dt = dt*1.25d0
            else
              dt = dt/rmeas
            endif
            dt = min(dt,dtmax)
            write(iow,*) 'New dt',dt
            floop(2) = .false.
            lv       = lv -1
          else
            if(rmeas.gt.1.25d0) then
              dtnew = dt*0.85d0/rmeas
            else
              dtnew = dt/3.d0
            endif
            write(iow,*) 'Old dt',dt
            write(iow,*) 'New dt',dtnew
            call autbac(dtnew)
            ct(1,l) = 1.0d0          ! Reset iteration counter
            ct(3,l) = 1.0d0          ! Reset iteration counter
            l = n                    ! Goto top of loop
          endif

        endif
        rmeas = 0.d0             ! Accpt conv val reset for next step

c     Augmented loop check:
c     N.B. Do not check if augm & solution in same loop

      elseif(lv.eq.lvaug .and. lv.ne.lvsol) then

        if((auggfl .and. augmfl .and. naugm.gt.1) .or.
     &                          ct(1,l).gt.ct(1,n)) then
          lv       = lv - 1              ! terminate loop
          lvaug    = 0
          if(auggfl .and. augmfl) then
            write(iow,2019) naugm, augg
            naugm  = 0
          endif
        else
          l        = n                   ! increment loop
        endif

c     Conventional loop checks

      else

c       Check for infinite looping

        if(pcomp(lct(n),'infi',4)) then
          ct(1,l) = 1.0d0          ! Reset iteration counter
          ct(3,l) = 1.0d0          ! Reset iteration counter
        endif

c       Check for a jump

        if(ljump .and. pcomp(lct(l),cjump,4)  .and.
     &            .not.pcomp('    ',cjump,4)) then
          floop(2) = .false.
          l        = njump
          lv       = vjump
          ljump    = .false.

        elseif(nint(ct(1,l)) .le. nint(ct(1,n))) then

c         if(ljump .and. dt.eq.0.0d0) then  ! Test zero time increment
          if(ljump) then  ! Reached end of time
            floop(2) = .false.
            lv       = lv - 1
            if(rank.eq.0) flconv   = .true.
            return
          endif

          l  = n                         ! increment loop

        else

c         No tangent convergence warning

          if(floop(1) .and. floop(2) .and. .not.linear .and.
     &       lv.eq.lvcn .and. niter.gt.1) then
            write(iow,2014) ttim,rnorm
            if(ior.lt.0) write(*,2014) ttim,rnorm
            floop(1) = .false.
            lvcn     = -1
            flncon   = .true.  ! Signal no convergence
          else
            if(rank.eq.0) flconv   = .true.
          end if

c         Exit loop

          floop(2) = .false.
          lv       = lv - 1

        end if
      end if

      return

c     Input proportional load table
c     [prop,,num1]      - input 1 to num1 proportional loads
c     [prop,,num1,num2] - input proportional loads num1 to num2
c     [prop,off]        - disable proportional loading

5     if(pcomp(lct(l),'off',3)) then
        npld = 0
        prop = 1.0d0
        ttim = ct(1,l)
        dt   = ct(2,l)
        if(np(121).ne.0) setvar = palloc(121,'PROP0',0,1) ! Remove prop
        if(np(122).ne.0) setvar = palloc(122,'PROP1',0,2) ! table data
        write(iow,2001) ttim,dt
        if(ior.lt.0 .and. prnt) write(*,2001) ttim,dt
      else
        uct    = lct(l)(1:4)
        npl(1) = nint(min(ct(1,l),ct(2,l)))
        npl(2) = nint(max(ct(1,l),ct(2,l)))
        if(npl(1).eq.0) then
          npl(2) = max(1,npl(2))
          npl(1) = npl(2)
        elseif(npl(1).lt.0) then
          npld = 0
          prop = 1.0d0
          write(iow,2001) ttim,dt
          if(ior.lt.0 .and. prnt) write(*,2001) ttim,dt
          return
        endif
        npld = max(npld,min(npl(2),50))
        prop = propld (ttim,npl)
        uct  = ' '
      endif

      return

c     Data command
c     [data,tol] : Set 'tol' as data (change during execution)
c     [data,dt]  : Set 'dt'  as data (change during execution)

6     if(ior.lt.0 .and. prnt) then
        write(xxx,3000) lct(l)
        call pprint(xxx)
      endif
      errck = tinput(lctl,2,ctl,3)
      if(errck) go to 6

c     Error diagnostics

      if(.not.pcomp(lct(l),lctl(1),4)) then
        write(ilg,4001)
        write(iow,4001)
        if(ior.lt.0) then
          write(*,4001)
          go to 6
        end if
        call plstop()
      end if

c     Set appropriate data values

      if(pcomp(lctl(1),'tol ',4)) tol  = ctl(1)
      if(pcomp(lctl(1),'tolr',4)) tolr = ctl(1)
      if(pcomp(lctl(1),'dt  ',4)) dt   = ctl(1)
      return

c     [time],,<tmax>      : Advance time by 'dt', quit after time > tmax
c     [time,set,ttim_new] : Time set to 'ttim_new'
c     [time,expl,<tmax>,c : Explicit critical time: dt = c * dt_cr

7     if(pcomp(lct(l),'set ',4)) then
        ttim = ct(1,l)
        write(iow,2016) ttim
        if(ior.lt.0) then
          write(*,2016) ttim
        endif
        return
      end if

c     Do an update if necesasry

      if(niter.le.1) then
        dtsav  = dt
        dt     = dtold
        hflgu  = .true.
        h3flgu = .true.
        call formfe(np(40),np(26),np(26),np(26),
     &             .false.,.false.,.false.,.false.,6,1,numel,1)
        dt     = dtsav
      endif

c     Write solution status to log file

      if(nstep.gt.0) then
        if(niter.gt.1) then
          write(ilg,3001) nstep,niter,nform,ttim,dt,rel0,rnorm,rnmax,
     &                    aengy,etime(tary)
        else
          write(ilg,3002) nstep,niter,nform,ttim,dt,rel0,rnorm,rnmax,
     &                    etime(tary)
        endif
        call pflush(ilg)
      endif

      call ptimpl()

c     Check critical time step for explicit transient problems

      if(pcomp(lct(l),'expl',4)) then
        crflg  = .false.
        nparta = npart
        dtcr   = 0.0d0
        do  i = 1,maxpart
          npart = i
          err   = .false.
          call partpt(npart,err,.true.)

c         Update dynamic vectors for step

          if(fl(9) .and. expflg) then
            call dsetci(.true.)
            call formfe(np(40),np(26),np(26),np(26),
     &             .false.,.false.,.false.,.false.,21,1,numel,1)
            if(dtcr.gt.0.0d0) then
              crflg = .true.
            endif
          end if

        end do ! i
        npart = nparta
        call partpt(npart,err,.true.)

c       Set critical time

        if(crflg) then
          if(ct(2,l).eq.0.0d0) then
            v1 = 0.9d0
          else
            v1 = min(1.d0,ct(2,l))
          endif
          dt = dtcr*v1
          if(ior.lt.0) then
            write(*,2024) v1
          endif
        endif
      endif

c     Increment time

      ttim  = ttim + dt

c     Set increment to end at tmax

      if(ct(1,l).gt.0.0d0 .and. ttim+dt*1.d-8.ge.ct(1,l)) then
        write(*,*) ' REACHED END TIME'
        ljump         = .true.
        v1            = dt
        dt            = dt + (ct(1,l) - ttim)
        ttim          = ct(1,l)
        ct(1,lve(lv)) = ct(1,lvs(lv))
        if(dt .le. 0.999999d0 * v1) then
          write(iow,2018) ttim,v1,dt
          if(ior.lt.0) then
            write(*,2018) ttim,v1,dt
          endif
        endif
      end if

c     Initialize step value and flags

      do i = 1,4
        ofl9(i)  = flp(9,i)
        fstep(i) = .true.
        steps(i) = 1.d0
      end do ! i

c     Update Database to n+1 (hflgu is true)

      hflgu  = .true.
      h3flgu = .true.
      call formfe(np(40),np(26),np(26),np(26),
     &           .false.,.false.,.false.,.false.,12,1,numel,1)

      npl(1) = 0
      propo = prop
      if(npld.gt.0) prop = propld(ttim,npl)
      if(arcf) then
        if(abs(prop-propo).gt.1.d-8) then
          write(ilg,4002)
          write(iow,4002)
          if(ior.lt.0) then
            write(*,4002)
          endif
          call plstop()
        endif
      endif
      if(prnt) then
        if(npld.gt.1) then
          k1 = npld
        else
          k1 = 0
        endif
        write(iow,2002) ttim,prop,(i,prldv(i),i=1,k1)
        if(ior.lt.0) then
          write(*,2002) ttim,prop,(i,prldv(i),i=1,k1)
        endif
      endif
      augf  = 1.0d0
      rnmax = 0.0d0

c     Zero displacement increment for time step

      call pzero(hr(np(40)+nneq),nneq)

c     Update monolithic solution for dynamics and contact

      if(solnfl .and. .not.monofl) then
        solnfl = .false.
      endif

      if(solnfl) then

        if(fl(9)) then
          call dsetci(.true.)
          call update(mr(id31),hr(np(30)),hr(np(40)),hr(np(42)),
     &                hr(np(26)),fl(9),1)
        end if

c       Update contact arrays

        call contact (301)

c     Update partition solution for dynamics and contact

      else

c       Check if using monolithic solution

        if(monofl) then
          solnfl = .true.
          do i = 1,ndf
            ndfp(i) = ndfg(i)
          end do ! i
        else
          solnfl = .false.
        endif

        nparta = npart
        do  i = 1,maxpart
          npart = i
          err   = .false.
          call partpt(npart,err,.true.)

c         Update dynamic vectors for step

          if(fl(9)) then
            call dsetci(.true.)
            call update(mr(id31),hr(np(30)),hr(np(40)),hr(np(42)),
     &                  hr(np(26)),fl(9),1)
          end if

c         Update contact arrays

          call contact (301)

        end do ! i

c       Update rotational arrays for static problem

        if(.not.fl(9) .and. frotas) then
          call updrot(hr(np(40)+2*nneq),ndf,hr(np(82)),mr(np(81)),
     &                numnp,1)
        endif

c       Restore active partition or set monolithic algorithm

        if(monofl) then
          call partpt(5,.false.,.true.)
          npart = 1
          do i = 1,ndf
            ndfp(i) = 1
          end do ! i
          oprt = prt
          prt  = .false.
          call profpart(mr(id31),mr(np(31)+nneq))
          prt  = oprt
        else
          npart = nparta
          call partpt(npart,err,.true.)
        endif
      endif

      fl( 8) = .false.
      fl(10) = .true.

c     Move interpolated force vector

      call pmove (hr(np(30)       ),hr(np(30)+  nneq), nneq)
      call pmove (hr(np(30)+2*nneq),hr(np(30)+3*nneq), nneq)

c     Reset history variables

      call reshis(mr(np(33)+nen),nen1,numel,2, 1)

c     Parallel time update

      call parsend()

c     Set iteration counter to indicate begining of time step

      nstep  = nstep + 1
      titer  = titer + niter
      niter  = 0
      taugm  = taugm + naugm
      naugm  = 0
      iaugm  = 0
      tform  = tform + nform
      nform  = 0
      iform  = 0
      intvc  = intvc + 1
      dtold  = dt

      if(echo .and. ior.gt.0) then
        write(*,3003) nstep,ttim
      endif

      return

c     Set print flag for print outputs
c     [prin/nopr,data]      - Permits data from inputs to be on or off
c     [prin/nopr,comm]and   - Controls print/noprint commands to screen
c     [prin/nopr,on/off]    - Old form for 'comm'and option
c     [prin/nopr,less]      - Prints shorter prompts to screen
c     [prin/nopr]           - Old form for 'less' option
c     [prin,mass/cmas/geom] - Print consistent mass array diagonals
c     [prin,iden/lmas]      - Print identity or lump mass array diags
c     [prin,tang/utan]      - Print tangent stiffness array diagonals
c     [prin,resi]           - Print residual array (diagonal)

8     if    (pcomp(lct(l),'data',4)) then
        prt  = j.eq.8
      elseif(pcomp(lct(l),'comm',4)) then
        prnt = j.eq.8
      elseif(pcomp(lct(l),'less',4)) then
        pfr  = j.eq.8
      elseif(pcomp(lct(l),'cont',4)) then
        conprt = j.eq.8
      elseif(pcomp(lct(l),'off', 3)) then
        prnt = .false.
      elseif(pcomp(lct(l),'on',  2)) then
        prnt = .true.
      elseif(pcomp(lct(l),'lmas',4).or.pcomp(lct(l),'iden',4)) then
        call mprint(hr(nl),1,neq,1,'Lmas/Iden ' )
      elseif(pcomp(lct(l),'cmas',4).or.pcomp(lct(l),'mass',4)
     &                             .or.pcomp(lct(l),'geom',4)) then
        call mprint(hr(nm),1,neq,1,'Cmas/Geom ' )
      elseif(pcomp(lct(l),'tang',4).or.pcomp(lct(l),'utan',4)) then
        call mprint(hr(na),1,neq,1,'Tang-diag ' )
      elseif(pcomp(lct(l),'resi',4)) then
        call mprint(hr(np(26)),1,neq,1,'Residuals ' )
      elseif(pcomp(lct(l),'    ',4)) then
        pfr  = j.eq.8
      end if
      return

c     Input integration parameters and initialize vectors
c     [beta,xxxx,beta,gamma,alpha]: see dparam for 'xxxx' options
c     [tran,xxxx,beta,gamma,alpha]: see dparam for 'xxxx' options

10    ofl9(npart)  =  flp(9,npart)
      fl(9)        = .true.
      flp(9,npart) = .true.
      call dparam(ct(1,l),lct(l))
      setvar = palloc( 42,'VEL  ',nneq*nrt ,2)
      tranfl = .true.
      return

c     Input initial conditions for dynamic integration
c     [init,disp] - set initial displacements
c     [init,rate] - set initial rates
c     [init,velo] - set initial velocity (same as rate)
c     [init,acce] - set initial accelerations
c     [init,spin,omg1,omg2,omg3] - set initial spins
c     [init,mate,v1,v2,v3] - Set velocity for material to a constant
c     [init,regi,v1,v2,v3] - Set velocity for region to a constant

11    call pinitl(lct(l),ct(1,l),err)
      intlfl = .true.
      return

c     Define an identity vector for stiffness eigen computation
c     [iden,,n1,n2]  - set dof n1 to n2 to unity

12    write(tname,'(4hLMAS,i1)') npart
      setvar = palloc(npart+12,tname,neq,2)
      nl    = np(npart+12)
      imtyp = 1
      idenf = .true.
      fl(1) = .false.
      fl(2) = .true.
      fl(5) = .false.
      nx    = npart+12
      call pzero(hr(nl),neq)
      n = nint(ct(1,l))
      n = max(1,(n-1)*ndf+1)
      i = nint(ct(2,l))*ndf
      if(i.eq.0) i = neq
      call piden(hr(nl),n,i)
      return

c     Update current force vector f0
c     [newf]
c     [newf,zero]

13    if(pcomp(lct(l),'zero',4)) then
        do i = 0,ndf-1
          if(ndfp(i+1).eq.npart) then
            do n = i,nneq-1,ndf
              hr(np(28)+n+2*nneq) = 0.0d0  ! f0(n,1) = fixed forces
              hr(np(28)+n+3*nneq) = 0.0d0  ! f0(n,2) = fixed displ.
            end do ! n
          end if
        end do ! i
        write(iow,2009) npart
        if(ior.lt.0 .and. prnt) write(*,2009) npart
      else
        call pload0(hr(np(27)),hr(np(28)+2*nneq),hr(np(40)),
     &              nneq,prop*rlnew )
      end if
      rlnew = 0.0d0
      return

c     Backup a time step
c     [back,,dt] - back-up to beginning of time step reset dt.

14    dtnew = max(0.0d0,ct(1,l))
      call autbac(dtnew)
      return

c     Debug flag on/off
c     [debug,on,level] or [debug,off] or [debug,,level]

15    if(pcomp(lct(l),'    ',4)) debug = .true.
      if(pcomp(lct(l),  'on',2)) debug = .true.
      if(pcomp(lct(l), 'off',3)) debug = .false.
      ndebug = nint(ct(1,l))
      if(debug) then
        if(ior.lt.0) write(  *,2003) ndebug
        if(ior.gt.0) write(iow,2003) ndebug
      else
        if(ior.lt.0) write(  *,2004)
        if(ior.gt.0) write(iow,2004)
      end if
      return

c     Linear problem - no test on convergence
c     [line]ar

16    if(.not.linear) then
        if(prnt) then
          if(ior.lt.0) write(  *,2005)
          write(iow,2005)
        end if
      end if
      linear = .true.
      return

c     Non-linear problem - test on convergence
c     [nonl]inear

17    if(linear) then
        if(prnt) then
          if(ior.lt.0) write(  *,2006)
          write(iow,2006)
        end if
      end if
      linear = .false.
      return

c     Set for auto time stepping control
c     [auto,time,iopt,fact,repeat]
c     [auto,dt,dtmin,dtmax]
c     [auto,mate] (old [auto,intv])
c     [auto,off]

18    if( pcomp('off',lct(l),3) ) then
        autofl  = .false.
        autofli = .false.
      elseif( pcomp('time',lct(l),4) ) then
        aratfl = .false.
        autofl = .true.
        lvauto = lv + 1
        iaopt  = nint(ct(1,l))
        afact  = ct(2,l)
        iarept = nint(ct(3,l))
        if(dtmax.eq.0.0d0) dtmax  = dt
        if(prnt) then
          if(ior.lt.0) then
            write(*,2010) iaopt,int(afact),iarept
          end if
          write(iow,2010) iaopt,int(afact),iarept
        end if
      elseif( pcomp('mate',lct(l),4) .or. pcomp('intv',lct(l),4) ) then
        autofli = .true.
        lvautoi = lv + 1
        rmeas   = 0.d0
        if(ct(1,l).eq.0.0d0) then
          rvalu(1) = 1.d0
        else
          rvalu(1)= ct(1,l)
        end if
        rvalu(2)= ct(2,l)
        rvalu(3)= ct(3,l)

        if(dtmax.eq.0.0d0) dtmax  = dt
        if(prnt) then
          if(ior.lt.0) then
            write(*,2020) rvalu
          end if
          write(iow,2020) rvalu
        end if
      elseif( pcomp('dt',lct(l),2) ) then
        dtmin = ct(1,l)
        dtmax = ct(2,l)
        if(prnt) then
          if(ior.lt.0) then
            write(*,2017) dtmin,dtmax
          endif
          write(iow,2017) dtmin,dtmax
        endif
      end if
      return

c     [method,type] - Set solution method for linear equations

19    nmeth = max(1,min(2,abs(int(ct(1,l)))))
      if(noi.ne.0) then
        if(ior.lt.0 .and. prnt) write(*,2012)
        write(iow,2012)
        nmeth = 1
      endif
      nmethp(npart) = nmeth
      if(ior.lt.0 .and. prnt) write(*,2011) tmeth(nmeth)
      write(iow,2011) tmeth(nmeth)
      return

c     [if,expression]   - begin if/then/else
c     [else,expression] - begin if/then/else

20    li      = li + 1
      lie(li) = int(ct(3,l))
      lexpr   = .false.
21    if(lexpr) then
        l = lie(li)
      else
        expr = lct(l)
        call evalex(expr,ctl,v1,80,errck)
        if(v1.lt.0.0d0) then
          l     = int(ct(2,l))
        else
          lexpr = .true.
        endif
      endif
      return

c     [endif] - endif if/then/else

22    li = li - 1
      return

c     [step],load,ds

24    if(steps(npart).eq.1.d0 .and. fstep(npart)) then
        steps(npart) =  0.0d0
        fstep(npart) = .false.
      endif

      if(pcomp(lct(l),'load',4) .or. pcomp(lct(l),'    ',4)) then
        if(ct(1,l).eq.0.0d0) ct(1,l) = 0.2d0
        steps(npart) = min(1.d0,steps(npart)+abs(ct(1,l)))
        write(iow,2015) steps(npart)
        if(ior.lt.0) then
          write(*,2015) steps(npart)
        endif
      endif
      return

c     [jump <label>] - Jump label point (does nothing)

25    continue
      return

c     [echo] -- Echo commands to screen

26    if(pcomp(lct(l),'off',3)) then
        echo = .false.
        write(*,*) ' -> ECHO off'
      else
        echo = .true.
        write(*,*) ' -> ECHO on'
      endif
      return

c     [conv] -- Output convergence details to 'log' file

27    if(pcomp(lct(l),'off',3)) then
        convfl = .false.
        if(ior.lt.0) write(*,*) ' -> CONVERGENCE outputs off'
      else
        convfl = .true.
        if(ior.lt.0) write(*,*) ' -> CONVERGENCE outputs on'
      endif
      return

c     [omeg]a,<hz>   - Frequency set in Hertz (convert to rad/sec)

28    if(pcomp(lct(l),'hz',2)) then
        omega1 = 2.d0*pi*ct(1,l)
        write(iow,2025) ct(1,l),omega1
        if(ior.lt.0) write(*,2025) ct(1,l),omega1
      else
        omega1 = ct(1,l)
        write(iow,2025) ct(1,l)/(2.d0*pi),omega1
        if(ior.lt.0) write(*,2025) ct(1,l)/(2.d0*pi),omega1
      endif
      omega2 = omega1**2
      return
c     formats

2001  format(' Number of proportional loads set to zero & prop = 1.0'/
     &       ' Time =',1p,1e12.4,' Dt =',1p,1e12.4)

2002  format(/,'   Computing solution at time ',1p,1e11.4,
     &         ': Total proportional load ',1p,1e11.4:/
     &         '   Individual factors: '/(3x,4(i4,' =',1p,1e12.4)))

2003  format(/'   Debug flag is set to .true. - Printing is on'/
     &        '   Debug level =',i4)

2004  format(/'   Debug flag is set to .false. - Printing is off'/)

2005  format(/'   Linear Problem -- no test on convergence with TOL.'/)

2006  format(/'   Non-linear Problem - test on convergence with TOL.'/)

2007  format(/'   Non-linear Problem - RNMAX =',1p,1e12.5/)

2008  format(/'   Energy assumed zero when less than ',1p,1e12.5/)

2009  format(/'   F0 vector set to zero in partition',i3/)

2010  format(/'   Auto time stepping :  Optimal iteration window'/
     &        '      Min. iterations :', i4/
     &        '      Max. iterations :', i4/
     &        '      Max. No. repeats:', i4/)

2011  format(/'   Solution of incremental equations by ',a,' method.'/)

2012  format(/'   Only Newton method allowed for transient solutions.')

2013  format(/'   Solution energy   tolerance = ',1p,1e12.5/
     &        '   Solution residual tolerance = ',1p,1e12.5/)

2014  format(/'   *WARNING* NO CONVERGENCE: Time =',1p,1e13.5,
     &        ': Residual =',1p,1e13.5/)

2015  format(/'   Load step factor = ',1p,1e12.4)

2016  format(/'   Time reset to = ',1p,1e12.4)

2017  format(/'   Auto time stepping:  Time increment limits'/
     &        '      Min. dt = ',1p,1e12.4/
     &        '      Max. dt = ',1p,1e12.4/)

2018  format( '   *WARNING* Maximum time',1p,1e12.4,' reached:',
     &        ' Time increment reset:'/
     &        '      dt_old = ',1p,1e12.4/
     &        '      dt_new = ',1p,1e12.4)

2019  format(/'   Augmented iteration converged in loop =',i4,
     &        ': Maximum gap =',1p,1e12.4/)

2020  format(/'   Auto time stepping :  Material R factor control'/
     &        '      Parameter 1     :', 1p,1e12.4/
     &        '      Parameter 2     :', 1p,1e12.4/
     &        '      Parameter 3     :', 1p,1e12.4/)

2021  format(/' *WARNING* Zero tolerance input - tol not reset:'/
     &        '   Solution tolerance = ',1p,1e12.5/)

2022  format('   Iteration: Residual tolerance = ',1p,1e13.6/
     &       '              Absolute tolerance = ',1p,1e13.6/
     &       '              Diverge  tolerance = ',1p,1e13.6)

2023  format('              Iterations         = ',i10)

2024  format('   Critical time factor times CFL:',1p,1e13.6)

2025  format('   Frequency:',1p,1e13.6,' Hz.',1p,1e13.6,' rad/sec')

3000  format(' Input ',a4,' Command >')

3001  format(2i6,i5,1p,1e11.3,1p,5e10.2,    0p,1f10.2)

3002  format(2i6,i5,1p,1e11.3,1p,4e10.2,10x,0p,1f10.2)

3003  format(' --> Step',i6,' Solution: Time =',1p,1e12.5)

4001  format(/' *ERROR* PMACR2: Command label mismatch on data',
     &        ' command')

4002  format(/' *ERROR* PMACR2: Variable proportional loading not',
     &        ' allowed with ARCLength.')

4003  format(/' *ERROR* LOOP: Infinite loop allowed in first loop only'/
     &        '         Current level = ',i4)
      end
