c$Id:$
      subroutine rinput (isw,flg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Control program for FEAP problem input and solution.

c      Inputs:
c        isw    - Switch parameter for operation
c         1       [rigid] - compute the rigid properties
c         2       [joint] - input joint definitions
c         3       [rload] - input rigid loads
c         4       [rboun] - specify boundary restraints on rigid body
c         5       [rdisp] - specify joint displacement values

c      Outputs:
c        flg    - Flag for return of action/error
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'allotd.h'
      include  'bdata.h'
      include  'cblend.h'
      include  'cdata.h'
      include  'cdat1.h'
      include  'chdata.h'
      include  'codat.h'
      include  'corset.h'
      include  'cornum.h'
      include  'comfil.h'
      include  'compac.h'
      include  'complx.h'
      include  'conval.h'
      include  'crotas.h'
      include  'edgdat.h'
      include  'errchk.h'
      include  'hlpdat.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'ioincl.h'
      include  'iosave.h'
      include  'linka.h'
      include  'mdata.h'
      include  'modreg.h'
      include  'mxsiz.h'
      include  'part3.h'
      include  'pdata2.h'
      include  'pdata5.h'
      include  'pdata6.h'
      include  'pointer.h'
      include  'prflag.h'
      include  'print.h'
      include  'region.h'
      include  'rigid1.h'
      include  'rigid2.h'
      include  'rjoint.h'
      include  'sdata.h'
      include  'trdata.h'
      include  'vdata.h'
      include  'comblk.h'

      include  'p_int.h'

      character titl*4,dnam*15, fext*8, type*4
      logical   flg, errs, pcomp, pinput,tinput, setvar,palloc
      logical   oprt,oprth
      integer   isw, i, j,jj,jt,jx, l1, neqrb,nrevs
      real*8    td(11)

      save

c     Input rigid body data and generate inertial properties

c     [rigid,nrbdof,npart]  - rigid body specification,dof/rb,partition

      if    (isw.eq.1) then

        rbody = .true.
        call acheck(xxx,yyy,15,80,80)
        read(yyy,1000,err=900,end=900) titl,nrbdof,nrbprt,neqrb
        if(nrbdof.le.0) then
          if(ndm.eq.1) then
            nrbdof = 1
          elseif(ndm.eq.2) then
            nrbdof = 3
          elseif(ndm.eq.3) then
            nrbdof = 6
          endif
        endif
        if(nrbprt.le.0) then
          nrbprt = 1
        endif
        if(neqrb.eq.0) then
          neqrb = -3
        endif
        write(iow,2001) nrbdof,nrbprt,neqrb,(jj,jj=1,ndm)
        if(ior.lt.0) then
          write(*,2001) nrbdof,nrbprt,neqrb,(jj,jj=1,ndm)
        end if
        do j = 1,nrbody
          if(rbcen(j).eq.0) then
            write(iow,2002) j,rbtype(j)
            if(ior.lt.0) then
              write(*,2002) j,rbtype(j)
            end if
          else
            write(iow,2003) j,rbtype(j),(rbx0(jj,j),jj=1,ndm)
            if(ior.lt.0) then
              write(*,2003) j,rbtype(j),(rbx0(jj,j),jj=1,ndm)
            end if
          endif
        end do ! j

        setvar = palloc( 96,'REQRB',nrbody*2, 1) ! Equation update type

c       Set location of update matrices

        do j = 1,nrbody
          mr(np(96) + nrbody + j - 1) = j
        end do ! j

c       Output modal body desciptions

        if(nmbody.gt.0) then
          do j = 1, nmbody
            write(iow,2010) modbod(nmbody)
            if(ior.lt.0) then
              write(*,2010) modbod(nmbody)
            endif
          end do ! j
        endif

        call rignod(mr(np(100)),mr(np(32)),mr(np(33)),mr(np(181)),
     &              nie,nen,nen1,numnp,numel)

c       Allocate inertia / R arrays

        setvar = palloc( 95,'RCG  ', nrbody*33    ,   2) ! CG locations
        setvar = palloc( 98,'RINER', nrbody*9     ,   2) ! Intertia body
        setvar = palloc( 99,'RIRB ', nrbody*nrbdof,   1) ! Eqn numbers
        setvar = palloc(104,'RLAMB', nrbody*54    ,   2) ! Orientation
        setvar = palloc(107,'RMASS', nrbody*1     ,   2) ! Mass body
        setvar = palloc(108,'RUROT', nrbody*6     ,   2) ! Rotn params
        setvar = palloc(182,'RBOU ', nrbody*nrbdof,   1) ! Boundary cond
        if(neqrb.eq.-5 .or.neqrb.eq.-6) then
          setvar = palloc(109,'REXMS', ndf*numnp  ,   2) ! Explicit Mass
          setvar = palloc(110,'REXIN', nrbody*36  ,   2) ! Expl Inertia
        endif

c       Compute location for center of mass in initial configuration

        call formrb(hr(np(41)),hr(np(44)),hr(np(39)),hr(np(35)),
     &              hr(np(36)),mr(np(32)),hr(np(25)),hr(np(43)),
     &              mr(np(33)),mr(np(181)),hr(np(104)),mr(np(96)),
     &              nrbody,neqrb)

c       Create list of joints
c       [joint]  - joint specifications
c       [type,node1,node2,body1,body2] - type data
c         type = 'ball'&socket - needs node1,node2
c         type = 'revo'lute    - needs node1,node2,body1,body2

      elseif(isw.eq.2) then

        write(iow,2004)

c       Input all joint data to determine storage allocation

        call plinka('jnts','set','   ')

c       Allocate storage for joint data

        setvar = palloc(101,'RJNT ',6*iclink,  1) ! Joint data
        setvar = palloc(102,'RJNX ',9*iclink,  2) ! Jt coord/angle/disp

        fext  = 'jnts'
        call pinpfl('RINPUT',fext, type, 1)
        oprt  = prt
        oprth = prth

c       Input joint data from temporary file

        jt     =  0
        jx     =  0
        nrevs  =  0
        numjts = -1

        dnam   = 'start'
        do while(.not.pcomp(dnam,'    ',4))

          errs = tinput(dnam,1,td,11)

c         Set coordinate array

          do i = 0,8
            hr(np(102)+jx+i) = td(i+3)
          end do ! i

c         Type 1: Ball and Socket Joint

          if(pcomp(dnam,'ball',4) .or.pcomp(dnam,'sphe',4)) then
            mr(np(101)+jt  ) = 1
            mr(np(101)+jt+1) = nint(td(1))
            mr(np(101)+jt+2) = nint(td(2))
            write(iow,2005) numjts+2,'Ball & Socket',
     &           (mr(np(101)+jt+i),i=1,2),(hr(np(102)+jx+i),i=0,2)

c         Type 2: Revolute Joint

          elseif(pcomp(dnam,'revo',4)) then
            nrevs            = nrevs + 1
            mr(np(101)+jt  ) = 2
            mr(np(101)+jt+1) = nint(td(1))
            mr(np(101)+jt+2) = nint(td(2))
            write(iow,2005) numjts+2,'Revolute',
     &           (mr(np(101)+jt+i),i=1,2),(hr(np(102)+jx+i),i=0,5)

c         Type 8: Basic Constraint Type 2

          elseif(pcomp(dnam,'bct2',4)) then
            nrevs            = nrevs + 1
            mr(np(101)+jt  ) = 8
            mr(np(101)+jt+1) = nint(td(1))
            mr(np(101)+jt+2) = nint(td(2))
            write(iow,2005) numjts+2,'Revolute',
     &           (mr(np(101)+jt+i),i=1,2),(hr(np(102)+jx+i),i=0,8)

c         Type 3: Slider Joint

          elseif(pcomp(dnam,'slid',4)) then
            nrevs            = nrevs + 1
            mr(np(101)+jt  ) = 3
            mr(np(101)+jt+1) = nint(td(1))
            mr(np(101)+jt+2) = nint(td(2))
            write(iow,2005) numjts+2,'Slider',
     &           (mr(np(101)+jt+i),i=1,2),(hr(np(102)+jx+i),i=0,5)

c         Type 4: Plane Joint

          elseif(pcomp(dnam,'plan',4)) then
            nrevs            = nrevs + 1
            mr(np(101)+jt  ) = 4
            mr(np(101)+jt+1) = nint(td(1))
            mr(np(101)+jt+2) = nint(td(2))
            write(iow,2005) numjts+2,'Plane',
     &           (mr(np(101)+jt+i),i=1,2),(hr(np(102)+jx+i),i=0,5)

c         Type 5: Translation Joint

          elseif(pcomp(dnam,'tran',4)) then
            nrevs            = nrevs + 1
            mr(np(101)+jt  ) = 5
            mr(np(101)+jt+1) = nint(td(1))
            mr(np(101)+jt+2) = nint(td(2))
            write(iow,2005) numjts+2,'Translation',
     &           (mr(np(101)+jt+i),i=1,2),(hr(np(102)+jx+i),i=0,5)

c         Type 6: Angular Control

          elseif(pcomp(dnam,'angl',4)) then
            nrevs            = nrevs + 1
            mr(np(101)+jt  ) = 6
            mr(np(101)+jt+1) = nint(td(1))
            mr(np(101)+jt+2) = nint(td(2))
            write(iow,2005) numjts+2,'Angle Control',
     &           (mr(np(101)+jt+i),i=1,2),(hr(np(102)+jx+i),i=0,8)

c         Type 7: Displacement Control

          elseif(pcomp(dnam,'disp',4)) then
            nrevs            = nrevs + 1
            mr(np(101)+jt  ) = 7
            mr(np(101)+jt+1) = nint(td(1))
            mr(np(101)+jt+2) = nint(td(2))
            write(iow,2005) numjts+2,'Displacement Control',
     &           (mr(np(101)+jt+i),i=1,2),(hr(np(102)+jx+i),i=0,8)

c         Error: No match on type

          elseif(.not.pcomp(dnam,'    ',4)) then
            write(iow,3001) dnam
            flg = .true.
          endif

c         Increase counters

          jt     = jt + 6
          jx     = jx + 9
          numjts = numjts + 1

        end do ! while

c       Set final memory and basis vectors for rotational joints

        if(nrevs.gt.0) then
          setvar = palloc( 97,'REVO ',nrevs*9,   2)! Jt revolute data
          jt = 0
          jx = 0
          l1 = 0
          do i = 1,numjts
            if(mr(np(101)+jt).ne.1) then
              call rbasis(hr(np(102)+jx), hr(np(97)+9*l1))
              l1 = l1 + 1
            end if
            jt = jt + 6
            jx = jx + 9
          end do ! i
        endif
        setvar = palloc(103,'RJTU ',numjts*15,   2)  ! Jt soln params

c       Close temporary file and restore logical input unit number

        call pinpfl('RINPUT',fext, type, 2)

        prt  = oprt
        prth = oprth

c     Rigid body loads
c     [rloa]  - Rigid load specifications
c     [body,1-comp,2-comp,etc.]

      elseif(isw.eq.3) then

        write(iow,2006) (j,j=1,6)
        j     = 7*nrlds
        errs  = pinput(td,7)
        l1    = nint(td(1))
        do while(l1.gt.0)
          setvar = palloc(106,'RLOAD',7*(nrlds+1),  2) ! Rigid body load
          hr(np(106)+j  ) = dble(l1)
          hr(np(106)+j+1) = td(2)
          hr(np(106)+j+2) = td(3)
          hr(np(106)+j+3) = td(4)
          hr(np(106)+j+4) = td(5)
          hr(np(106)+j+5) = td(6)
          hr(np(106)+j+6) = td(7)
          write(iow,2007) l1,(hr(np(106)+j+i),i = 1,6)
          j     = j + 7
          nrlds = nrlds + 1
          errs  = pinput(td,7)
          l1    = nint(td(1))
        end do ! while

c     Rigid body boundary conditions
c     [rbou]  - Rigid b.c. specifications
c     [body,1-comp,2-comp,etc.]

      elseif(isw.eq.4) then

        write(iow,2008) (j,j=1,nrbdof)
        errs = pinput(td,nrbdof+1)
        j = nint(td(1))
        do while(j.gt.0)
          fp(1) = np(182) + (j-1)*nrbdof - 1
          do i = 1,nrbdof
            mr(i+fp(1)) = nint(td(i+1))
          end do ! i
          write(iow,2009) j,(mr(i+fp(1)),i=1,nrbdof)
          errs = pinput(td,nrbdof+1)
          j = nint(td(1))
        end do ! while

      elseif(isw.eq.5) then

      endif

      return

900   call  endclr ('RINPUT',record)
      call pdelfl()

c     Formats

1000  format(a4,11x,3i15)

2001  format(/'   R i g i d   B o d y    D a t a'/
     &        '      ndf       Partition   Rot. Update'/ i8,2i12//
     &        '                Rigid Body     Type'/
     &        '                  Number    (0=R,1=M)',3(i6,'-Coord':))

2002  format(i22,i12,9x,'At center of mass')
2003  format(i22,i12,3x,1p,3e12.4)

2004  format(/'   R i g i d   B o d y    J o i n t   D a t a'//
     & '     Joint      Type        1-Body  2-Body',
     &         2x,'1-Coord     2-Coord     3-Coord'/
     &        44x,'1-Coord     2-Coord     3-Coord'/
     &        44x,'1-Value     2-Value     3-Value'/)

2005  format(/i8,3x,a13,2i8,1p,3e12.4:/(40x,1p,3e12.4))

2006  format(/'   R i g i d   B o d y    F o r c e   D a t a'//
     &        '   Body',6(i7,'-Comp'))
2007  format(i7,1p,6e12.4)

2008  format(/'   R i g i d   B o d y    B.  C.    D a t a'//
     &        '   Body',6(i3,'-b.c.'))
2009  format(i7,6i8)

2010  format(/'  Modal region number =',i4/)

3001  format(/' *ERROR* Incorrect joint type:',a15)

      end
