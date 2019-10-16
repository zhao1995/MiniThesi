c$Id:$
      subroutine pmacr7 (lct,ct,j)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Increase size of GH-array to *4 at line 570      16/03/2007
c       2. Add 'glob'al and 'loca'l commands                13/04/2007
c       3. Remove argument 'prt7'                           13/12/2008
c       4. Change 'include/finclude' to 'finclude'          23/01/2009
c       5. Add output of initial conditions to domain file  01/02/2009
c          in subprogram uoutdom.
c       6. Remove common 'pfeapa' (values in 'setups')      05/02/2009
c       7. Separate eq and id on call to psetid             18/11/2009
c          Change mr(31) to id31 on call to profil, idreset
c          Change 'id' to 'eq' on call to ueqn1
c       8. Add output for turning PETSc ON for ntasks equal 24/11/2009
c          one (1) processor
c       9. Change FEAPHOME8_2 to FEAPHOME8_3                01/12/2009
c      10. Correct format 2006                              06/04/2011
c      11. In 'uoutdom' extend formats 2003 and 2006 to     02/02/2012
c          larger numbers of dofs.
c      12. Change FEAPHOME8_3 to FEAPHOME8_4                05/11/2012
c      13. Add output of load groups to domains             06/11/2012
c      14. Change outd options to aij, baij with option     05/01/2013
c          for including bc-equations in aij format
c      15. Remove double allocation                         07/01/2013
c          Change processor to ntasks
c      16. Set npnp and nper for glot                       09/01/2013
c          Change mpirun call for graph,part
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Command language instruction subprogram: Part 7

c      Inputs:
c         lct      - Command option for current command
c         ct(3)    - Command parameters for current command
c         j        - Number of command to execute

c      Outputs:
c         Depends on value of command j
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

#     include  "finclude/petsc.h"

#     include   "cdata.h"
#     include   "cdat1.h"
#     include   "comfil.h"
#     include   "compas.h"
#     include   "complx.h"
#     include   "debugs.h"
#     include   "evdata.h"
#     include   "fdata.h"
#     include   "idptr.h"
#     include   "iodata.h"
#     include   "iofile.h"
#     include   "part0.h"
#     include   "part1.h"
#     include   "part3.h"
#     include   "pdata3.h"
#     include   "pfeapb.h"
#     include   "plist.h"
#     include   "print.h"
#     include   "prstrs.h"
#     include   "rdata.h"
#     include   "region.h"
#     include   "rjoint.h"
#     include   "sdata.h"
#     include   "setups.h"
#     include   "strnum.h"
#     include   "tdata.h"
#     include   "xtout.h"

#     include   "c_0.h"
#     include   "c_contac.h"
#     include   "c_ccp.h"

#     include   "pointer.h"
#     include   "comblk.h"

      logical    pcomp, setval, palloc, gprt, exec, pmetis
      logical    prntsv, sfl, initf, opfl, cknon0

      character  gfile*89,gplf*128,fext*5,c*1

      integer    i,j,k, domains,nftot, nodglen, ixinvlen,numedg, nparto
      integer    ndlast,newnel,oldnel, ndi
      integer    imas, mad, pneq, pqmx, gnode, bserchi
      integer    ii,ista
      integer    skips,skipm,nslaves,nmaster

      character  lct*15, tname*5
      real*8     ct(3), stol
      real*4     etime,tary(2)

      save

c     Macro instruction subprogram - part 7.

c     [graph,<node>,d] d = number of domains
c     [graph,file]     Input graph partion from existing file

      if(j.eq.1) then

c       Set flag for graph prints

        if(ct(2).gt.0.0d0) then
          gprt = .true.
        else
          gprt = .false.
        endif

c       Compute nodal graph partition vector using metis

        if(pcomp(lct,'    ',4) .or. pcomp(lct,'node',4)) then

c         Set number of domains (regions in feap data sets)

          if(debug) then
            write(  *,*) ' Start IX_INV domains',etime(tary)
            write(iow,*) ' Start IX_INV domains',etime(tary)
          endif
          domains     = nint(ct(1))
          if(domains.le.0) then
            write(*,2000)
            read (*,   *) domains
            write(*,2001) domains
          endif

c         Set maximum number of domains in problem

          mxreg       = domains

c         Use Metis to build nodal partitions

          if(np(254).eq.0) then
            setval = palloc(254,'NODPT',numnp,1)
          endif

c         Case with 1 domain

          if(domains.eq.1) then

            do i = 0,numnp-1
              mr(np(254)+i) = 1
            end do ! i

c         Case with more than 1 domain

          elseif(domains.gt.1) then

c           Set up contact
            if(ifct) then
              setval = palloc(140,'CTEM5',numnp,1)
              setval = palloc(143,'CTEM8',cck(3)+1,1) ! slave pointers
              setval = palloc(144,'CTEM9',cck(3)+1,1) ! master pointers
              mr(np(143)) = 1
              mr(np(144)) = 1

c             Get the count of total number of slave and master nodes
              do i = 1,cck(3)
c
c               np(140)+(k-1) = +contact_pair or -contact_pair if node
c               k is part of the contact pair
c
                call crsurfinv(i,hr(ccp(1)),hr(ccp(3)),mr(np(140)))
                nslaves = 0
                nmaster = 0
                do k = 0,numnp-1
                  if(mr(np(140)+k).eq.-i) nslaves = nslaves + 1
                  if(mr(np(140)+k).eq. i) nmaster = nmaster + 1
                end do
c
c               pointers to start of slave/master side nodes per
c               contact pair in np(141)/np(142)
c
                mr(np(143)+i) = mr(np(143)+i-1) + nslaves
                mr(np(144)+i) = mr(np(144)+i-1) + nmaster
              end do
              setval = palloc(141,'CTEM6',mr(np(143)+cck(3))-1,1) ! slave list
              setval = palloc(142,'CTEM7',mr(np(144)+cck(3))-1,1) ! master list

c             Fill arrays
              skips = 0
              skipm = 0
              do i = 1,cck(3)
                call crsurfinv(i,hr(ccp(1)),hr(ccp(3)),mr(np(140)))
                do k = 0,numnp-1
                  if(mr(np(140)+k).eq.-i) then
                    mr(np(141)+skips) = k+1
                    skips = skips + 1
                  elseif(mr(np(140)+k).eq. i) then
                    mr(np(142)+skipm) = k+1
                    skipm = skipm + 1
                  endif
                end do
              end do
              setval = palloc(140,'CTEM5',0,1)
             if(debug) then
              call iprint(mr(np(143)),1,cck(3)+1,1,'SLAVE-POINTERS')
              call iprint(mr(np(141)),1,mr(np(143)+cck(3))-1,1,'SLAVES')
              call iprint(mr(np(144)),1,cck(3)+1,1,'MASTER-POINTERS')
              call iprint(mr(np(142)),1,mr(np(144)+cck(3))-1,1,'MASTER')
             endif
            endif

c           Build inverse ix list

            setval  = palloc(111,'TEMP1',numnp+1,1)
            call size_ixinv(mr(np(33)),mr(np(111)),ixinvlen)

            setval  = palloc(112,'TEMP2',ixinvlen,1)
            setval  = palloc(113,'TEMP3',numnp,1)
            call create_ixinv(mr(np(33)),mr(np(111)),mr(np(113)),
     &                        mr(np(112)))
            setval  = palloc(113,'TEMP3',0,1)
            if(debug) then
              write(*,*) 'IXINVLEN ',ixinvlen
              call iprint(mr(np(111)),1,numnp+1,1,'V')
              call iprint(mr(np(112)),1,ixinvlen,1,'IXINV')
              write(  *,*) ' End   IX_INV domains',etime(tary)
              write(iow,*) ' End   IX_INV domains',etime(tary)
            endif

c           Build nodal graph and partition

            if(debug) then
              write(  *,*) ' Start METIS  domains',etime(tary)
              write(iow,*) ' Start METIS  domains',etime(tary)
            endif

c           Allocate memory for graph pointer array

            setval  = palloc(252,'XADJN',numnp+1,1)
            nodglen = numnp*30     ! Estimate for nodal graph list
            oldnel  = 30
            if(debug) then
              write(  *,*) ' STARTING NODG STORAGE OLDNEL =',oldnel
              write(iow,*) ' STARTING NODG STORAGE OLDNEL =',oldnel
            endif
100         setval = palloc(253,'NODG',nodglen,1)
            call pzeroi(mr(np(253)),nodglen)
c                                IX          V           IXINV
c                                XADJN       NODG        LEN
            call make_node_graph(mr(np( 33)),mr(np(111)),mr(np(112)),
     &                           mr(np(252)),mr(np(253)),nodglen,numedg,
     &                           ndlast)

c           If list not completed need to redo with more storage

            if(ndlast.lt.numnp) then
              newnel  =(nodglen+ndlast-1)/ndlast
              if(newnel.le.oldnel) then
                newnel = oldnel + 2
              endif
              nodglen = numnp*newnel
              oldnel  = newnel
              if(debug) then
                write(  *,*) ' INCREASE NODG STORAGE OLDNEL =',oldnel
                write(iow,*) ' INCREASE NODG STORAGE OLDNEL =',oldnel
              endif
              setval  = palloc(253,'NODG',0,1)
              go to 100
            endif
            setval = palloc(112,'TEMP2',0,1)
            setval = palloc(111,'TEMP1',0,1)
            setval = palloc(253,'NODG',nodglen,1)  ! Final storage


            if(debug) then
              write(  *,*) ' MEMORY USE: IXINVLEN =',ixinvlen,
     &                                 ' NODGLEN = ',nodglen
              write(iow,*) ' MEMORY USE: IXINVLEN =',ixinvlen,
     &                                 ' NODGLEN = ',nodglen
              write(  *,*) ' NODG =',nodglen/numnp,oldnel,
     &                     ' time =',etime(tary)
              write(iow,*) ' NODG =',nodglen/numnp,oldnel,
     &                     ' time =',etime(tary)
              call iprint(mr(np(252)),1,numnp+1,1,'XADJN')
              call iprint(mr(np(253)),1,nodglen,1,'NODEG')
            endif

c           Call Metis

            call pzeroi(mr(np(254)),numnp)
            exec = pmetis(domains,numedg)
          endif

          if(debug .and. ndebug.gt.1) then
            call iprint(mr(np(254)),1,numnp,1,' NODAL REGIONS ')
          endif

          if(debug) then
            write(  *,*) ' End   METIS  domains',etime(tary)
            write(iow,*) ' End   METIS  domains',etime(tary)
          endif

c       Graph partition vector in disk file "graph.filename"

        elseif(pcomp(lct,'file',4) .or. pcomp(lct,'part',4)) then

c         Create graph file using ParMETIS

          if(pcomp(lct,'part',4)) then

c           Copy input file to Gfile

            i                =  index(finp,' ')
            gfile(1:3)       = 'cp '
            gfile(4:i+4)     = finp(1:i)
            gfile(i+3:i+6)   = '.rev'
            gfile(i+7:i+12) = ' Gfile'
            call system (gfile)

c           Set file command to run partitioning using ParMETIS

            gfile( 1:15) = 'mpirun -n   2  '
            write(gfile(11:13),'(i3)') nint(ct(1)) ! Number of processors

c           Location of file 'partition using FEAPHOME8_4 definition'

            gfile(16:79) = '$FEAPHOME8_4/parfeap/partition/partition '
            gfile(80:89) = '   2 Gfile'
            write(gfile(81:83),'(i3)') nint(ct(1)) ! Number of domains

c           Execute 'partition'

            call system (gfile)

c           Remove scratch Gfile and finp.rev files

            gfile          = 'rm -f Gfile'
            call system (gfile)
            gfile(  7:i+7) = finp(1:i)
            gfile(i+6:i+9) = '.rev'
            call system (gfile)

c           Set graph filename

            gfile = 'graph.file'

c         Extract graph filename from input file name

          else
            gfile = 'graph. '
            i     =  index(finp,' ')
            gfile(7:4+i) = finp(2:i-1)
          endif
          if(ior.lt.0) then
            write(*,*) '  Input graph from ',gfile(1:4+i)
          endif

c         Perform graph input from file 'gfile'

          inquire(file = gfile,exist = initf)
          if(initf) then
            open(unit = ios, file=gfile, status = 'old')
            setval = palloc(254,'NODPT',numnp,1)
            setval = palloc(252,'XADJN',numnp+1,1)
            call getgfile(ios,mr(np(254)),numnp)    ! Input node partitions
            call getgfile(ios,mr(np(252)),numnp+1)  ! Input graph pointers
            nodglen = mr(np(252)+numnp)-1           ! Set node graph length
            setval = palloc(253,'NODG',nodglen,1)
            call getgfile(ios,mr(np(253)),nodglen)  ! Input graph
            close(unit = ios)

c           Compute number of domains to output

            domains = 0
            do i = 0,numnp-1
              domains = max(mr(np(254)+i),domains)
            end do ! i
            mxreg = domains
            if(debug) then
              write(  *,*) ' Partition into ',domains,' domains'
              write(iow,*) ' Partition into ',domains,' domains'
            endif
          else
            write(*,*) 'FILE = ',gfile(1:4+i),' does not exist'
            write(*,*) '        Run partition to create file.'
            if(ior.gt.0) then
              call plstop()
            else
              return
            endif
          endif

c         Remove graph file

          if(pcomp(lct,'part',4)) then
            gfile = 'rm -f graph.file'
            call system(gfile)
          endif
        endif

c       Output lists if requested

        if(gprt) then
          call iprint(mr(np(254)),1,numnp,1,'NODPT')
          call iprint(mr(np(252)),1,numnp+1,1,'XADJN')
          call iprint(mr(np(253)),1,nodglen,1,'NODG')
          call ioprof(mr(np(252)),mr(np(253)),numnp,'Node Graph',2)
        endif

c     [outd,<aij,baij>,<bcin_flag>,<nsbk>] - Output of domain meshes

      elseif(j.eq.2) then

c       Set number of domains (regions in feap data sets)

        domains     = max(1,mxreg)

c       Set blocking type

        pfeap_bcin = .false.        ! Default no BC eqs in the assembly
        if(.not.pcomp('baij',lct,4) .and. .not.pcomp('aij',lct,3)) then
          pfeap_blk  = .false.
          if(debug) then
            write(  *,*) ' Equations AIJ format'
            write(iow,*) ' Equations AIJ format'
          endif
        else                         ! BAIJ or AIJ
          if(pcomp('baij',lct,4)) then
            pfeap_blk  = .true.
            pfeap_bcin = .true.       ! BC eqs are in the assemble for BAIJ form
            if(debug) then
              write(  *,*) ' Equations BAIJ format'
              write(iow,*) ' Equations BAIJ format'
            endif
          else
            pfeap_blk = .false.
            if(nint(ct(1)).ne.0) pfeap_bcin = .true.
            if(debug) then
              write(  *,*) ' Equations AIJ format'
              write(iow,*) ' Equations AIJ format'
            endif
          endif

          if(pfeap_bcin) then
            if(nint(ct(2)).eq.0) then
              nsbk = ndf
            else
              nsbk = nint(ct(2))
              if(mod(ndf,nsbk).ne.0) then
                if(mod(ndf,ndm).eq.0) then
                  i    = ndf/ndm
                  nsbk = ndf/i
                else
                  nsbk = ndf
                endif
                if(debug) then
                  write(  *,3011) nsbk
                  write(iow,3011) nsbk
                endif
              endif
            endif
          endif

        endif

c       Output the domain data

        setval = palloc(111,'TEMP1',numnp+1  ,1)  ! nodx local node array
        setval = palloc(112,'TEMP2',domains*4,1)  ! revp array

c       Determine size for renumber list

        if(debug) then
          write(  *,*) ' Start UOUTD1 domains',etime(tary)
          write(iow,*) ' Start UOUTD1 domains',etime(tary)
        endif
        call uoutd1 (mr(np(33)),mr(np(254)),mr(np(111)),mr(np(112)),
     &               domains, nftot)
        setval = palloc(113,'TEMP3',nftot*2,1)  ! revx array
        if(debug) then
          write(  *,*) ' End   UOUTD1 domains',etime(tary)
          write(iow,*) ' End   UOUTD1 domains',etime(tary)
        endif

c       Set renumber list

        if(debug) then
          write(  *,*) ' Start UOUTD2 domains',etime(tary)
          write(iow,*) ' Start UOUTD2 domains',etime(tary)
        endif
        call uoutd2 (mr(np(33)),mr(np(254)),mr(np(111)),mr(np(112)),
     &               mr(np(113)), domains)

c       Allocate additional arrays used in output of partition meshes

        setval = palloc(114,'TEMP4',numel*2,1)
        setval = palloc(115,'TEMP5',nneq*2,1)
        setval = palloc(116,'TEMP6',domains*2,1)
        setval = palloc(117,'TEMP7',numnp,1)
        setval = palloc(118,'TEMP8',domains*domains*5,1)
        if(debug) then
          write(  *,*) ' End   UOUTD2 domains',etime(tary)
          write(iow,*) ' End   UOUTD2 domains',etime(tary)
        endif

c       Output domain data files

        call uoutdom(mr(id31),mr(np(31)+nneq),mr(np(33)),mr(np( 34)),
     &               hr(np( 43)),hr(np( 45)),hr(np( 27)),mr(np( 29)),
     &               hr(np( 86)),hr(np( 87)),hr(np( 88)),mr(np(254)),
     &               mr(np(111)),mr(np(112)),mr(np(113)),mr(np(114)),
     &               mr(np(115)),mr(np(116)),mr(np(117)),mr(np(118)),
     &               domains)
        setval = palloc(118,'TEMP8',0,1)
        setval = palloc(117,'TEMP7',0,1)
        setval = palloc(116,'TEMP6',0,1)
        setval = palloc(115,'TEMP5',0,1)
        setval = palloc(114,'TEMP4',0,1)
        setval = palloc(113,'TEMP3',0,1)
        setval = palloc(112,'TEMP2',0,1)
        setval = palloc(111,'TEMP1',0,1)

c       Reset ID equation numbers

        setval = palloc(111,'TEMP1',nneq,  1)
        call psetid(mr(id31),mr(np(31)+nneq),mr(np(33)),mr(np(111)),
     &              mr(np(32)),mr(np(240)),nie,ndf,nen,nen1,
     &              numel,numnp,nummat)
        setval = palloc(111,'TEMP1',   0,  1)

c       Set region indicator so all active region elements are assembled

        nreg = -1

c       Compute new equation numbers and profile

        nparto =  npart
        do i = 1,4
          if(.not.tflp(i)) then
            call partpt(i,tflp(i),.false.)

c           Set new current profile

            call profil(mr(np(20+i)),mr(np(34)),mr(id31),
     &                  mr(np(33)),1,pfr)
            call profil(mr(np(20+i)),mr(np(34)),mr(id31),
     &                  mr(np(33)),2,pfr)
            nqp(i)    =  neq
            nqr(i)    =  neqr
            mxprop(i) = (mr(np(20+i)+neq-1))*ipc
            mxneqp(i) =  neq*ipc

c           Compute new fill for sparse storage

            if(ittyp.gt.-3) then
              if(ittyp.eq.-1 .and. np(92).ne.0) then
                setval = palloc(92,'OINB',0,1)
              endif
              nftot = 0
              call iters(nftot,1)
            endif
          endif
        end do ! i
        npart = nparto
        call partpt(npart,tflp(npart),.true.)

c     [pets]c,<on,off,view,novi> - Parallel solution using PETSc

      elseif(j.eq.3) then

        if(pcomp(lct,'view',4)) then      ! Turn on PETSc array viewer
          pfeap_view = .true.
          if(ior.lt.0) then
            write(*,2003) 'ON'
          endif
        elseif(pcomp(lct,'novi',4)) then  ! Turn off PETSc array viewer
          pfeap_view = .false.
          if(ior.lt.0) then
            write(*,2003) 'OFF'
          endif
        elseif(pcomp(lct,'on',2)) then
          solver     = .false.
          pfeap_on   = .true.
          pfeap_dstr = .false.
          if(pfr) then
            write(iow,2002)
            if(ior.lt.0) then
              write(*,2002)
            endif
          endif
          if(pfeap_bcin) then
            call pidreset(mr(id31))
          endif
        elseif(pcomp(lct,'off',3)) then   ! Turn off PETSc solution
          solver   = .true.
          pfeap_on = .false.
          write(tname,'(4hTANG,i1)') npart
          setval = palloc(npart,tname, 0,2)
        endif

c     Subspace eigencomputations (for: mass,iden,geom)
c     [psub,,k1,k2]     - subspace for k1 eigenpairs
c     [psub,prin,k1,k2] - subspace for k1 eigenpairs - print matrices
c                       - k2 used to overwrite default no. guard vects.

      elseif(j.eq.4) then

c       Check for tangent matrix

        if(fl(4)) then
          write(ilg,3002) 'PSUB'
          write(iow,3002) 'PSUB'
          if(ior.lt.0) write(*,3002) 'PSUB'
          return
        endif

c       Check for mass matrix

        if(fl(5).and.fl(6)) then
          write(ilg,3003)
          write(iow,3003)
          if(ior.lt.0) then
            write(*,3003)
            return
          endif
          call plstop()
        endif

        if(fl(1)) then         ! Use mass matrix
          imas = 1
        else                   ! Use geometric stiffness
          imas = 2
        endif

        mf  = nint(ct(1))           ! Number of eigenpairs
        mad = nint(ct(2))           ! Number of extra vectors
        stol=      ct(3)            ! Solution tolerance

c       Mask tolerance if it is too small

        if(stol.eq.0.d0) then
          stol = max(tol, 1.d-12)
        endif

        if(mad.le.0) mad = max(mf,8)
        mf = min(numteq,max(1,mf))
        mq = min(mf+mad,numteq)

c       Number of possible terms positive

        if(mq.gt.0) then
          prntsv = prnt
          if(ittyp.eq.-2) then
            prnt = .false.
          endif
          if(mq.lt.mf) then
            write(iow,2004) mq
            if(ior.lt.0) then
              write(*,2004) mq
            endif
          endif
          mf     = min(mf,mq)
          pneq   = numpn*ndf               ! Size of active domain
          vneq   = numnp*ndf               ! Size of domain vector
          setval = palloc( 76,'EVAL',mq     ,2)       ! d  - array
          setval = palloc( 77,'EVEC',mq*vneq,2)       ! v  - array

          setval = palloc(111,'TEMP1',vneq       ,2)  ! t  - temps
          setval = palloc(112,'TEMP2',mq*(mq+1)*4,2)  ! gh - array
          setval = palloc(113,'TEMP3',mq         ,2)  ! dp - array
          setval = palloc(114,'TEMP4',mq         ,2)  ! dtl- array
          setval = palloc(115,'TEMP5',mq*mq      ,2)  ! p  - array
          setval = palloc(116,'TEMP6',ntasks+1,1)  ! proceq

          sfl    = pcomp(lct,'prin',4)
          setval = pfr
          pfr    = sfl
          pqmx   = (mq*(mq+1))/2
          call psubsp(hr(np(77)) ,hr(np(111)),
     &                hr(np(112)),hr(np(76)) ,hr(np(113)),
     &                hr(np(114)),hr(np(115)),mr(np(116)),
     &                mf,mq,pqmx,pneq,imas,shift,stol,sfl,25)
          pfr    = setval

          setval = palloc(116,'TEMP6', 0,1)
          setval = palloc(115,'TEMP5', 0,2)
          setval = palloc(114,'TEMP4', 0,2)
          setval = palloc(113,'TEMP3', 0,2)
          setval = palloc(112,'TEMP2', 0,2)
          setval = palloc(111,'TEMP1', 0,2)

c         Restore print flag

          prnt = prntsv

c       There are no possible modes

        else

          if(imtyp.eq.1) then
            write(ilg,3009)
            write(iow,3009)
            if(ior.lt.0) then
              write(*,3009)
            endif
          elseif(imtyp.eq.2) then
            write(ilg,3010)
            write(iow,3010)
            if(ior.lt.0) then
              write(*,3010)
            endif
          endif
          if(ior.gt.0) call plstop()

        endif

c     [u_in]<inpu,outp> - initial solution inputs/outputs

      elseif(j.eq.5) then

        if(pcomp(lct,'outp',4)) then
          call u_init_out(hr(np(40)))
        elseif(pcomp(lct,'inpu',4)) then
          call u_init_inp(hr(np(40)))
        endif

c     [glis]t,,# -- Input list of global node numbers to output

      elseif(j.eq.6) then

        i = max(1,min(3,nint(ct(1))))
        write(iow,2005) i
        if(ior.lt.0) then
          write(*,2005) i
        endif
        call prlist(iolist(1,i),niols(i))

c       Convert list to global list

        if(pfeap_on .and. ntasks.gt.1) then
          k = 0
          do j = 1,niols(i)
            gnode = bserchi(mr(np(244)),numpn, iolist(j,i))
            if(gnode.gt.0) then
              k           = k + 1
              iolist(k,i) = gnode
            endif
          end do ! j
          niols(i) = k
          if(niols(i).gt.0) then
            write(iow,2006) (j,j=1,min(niols(i),8))
            do k = 1,niols(i),8
              write(iow,2007) (iolist(j,i),j=k,min(niols(i),k+7))
            end do ! k
            if(ior.lt.0) then
              write(*,2006) (j,j=1,min(niols(i),8))
              do k = 1,niols(i),8
                write(*,2007) (iolist(j,i),j=k,min(niols(i),k+7))
              end do ! k
            endif
          endif
        endif

c     [gplo]t,disp,#  -- Output displacement comp. # for global plot with ndata
c     [gplo]t,stre,#  -- Output nodal stress comp. # for global plot with ndata
c     [gplo]t,pstr,#  -- Output nodal stress comp. # for global plot with ndata

      elseif(j.eq.7) then

       opfl = pcomp('disp',lct,4) .or. pcomp('stre',lct,4) .or.
     &        pcomp('pstr',lct,4)
       if(opfl) then

        if(pcomp('disp',lct,4)) c = 'd'
        if(pcomp('stre',lct,4)) c = 's'
        if(pcomp('pstr',lct,4)) c = 'p'

c       Open up next G-file for writing
        inquire(unit=ios,opened = opfl)
        if(opfl) then
           close(ios)
        endif

        ii = 1
55      gplf = fplt
        gplf(1:1) = 'G'
        fext(1:5) = ' '
        if(ii.lt.10) then
          write(fext,'(a1,a3,i1)') c,'000',ii
        elseif(ii.lt.100) then
          write(fext,'(a1,a2,i2)') c,'00',ii
        elseif(ii.lt.1000) then
          write(fext,'(a1,a1,i3)') c,'0',ii
        elseif(ii.lt.10000) then
          write(fext,'(a1,i4)') c,ii
        endif
        call addext(gplf,fext,128,5)
        inquire(file=gplf, exist=opfl )
        if(opfl) then
          ii = ii + 1
          go to 55
        endif
        open(unit=ios,file=gplf,form='formatted',status='new')
        write(*,*) 'Creating global plot file: ',gplf
        write(ios,*) ntasks

c       Write data
        if(pcomp('disp',lct,4)) then
          ndi = max(1,min(nint(ct(1)),ndf)) - 1
          write(ios,*) ttim, '  DISPLACEMENT',ndi+1
          do ii = 1,numpn
           write(ios,*) mr(np(244)+ii-1),
     &                  hr(np(40)+ndi+(ii-1)*ndf)
          end do
        elseif(pcomp('stre',lct,4)) then
          if(.not.fl(11)) then
           if (plfl) then
             setval = palloc ( 58,'NDNP',numnp*npstr,2)
             setval = palloc ( 57,'NDER',numnp*8    ,2)
             setval = palloc ( 60,'NDNS',max(nen*npstr,nst*nst),2)
             setval = palloc (207,'NSCR',numel      ,2)
             plfl = .false.
           endif
           nph  = np(58)
           ner  = np(57)
           call pjstrs(.false.)
           fl(11) = .true.
          endif
c         Compute number of non-zero stress components
          ista = 0
          do ii = 1,abs(istv)
            if(cknon0(hr(np(58)+numnp*ii),numnp) ) ista = ii
          end do
          ndi = max(1,min(nint(ct(1)),ista))
          write(ios,*) ttim, ' STRESS',ndi
          do ii = 1,numpn
            write(ios,*) mr(np(244)+ii-1),
     &       hr(np(58)+numnp*ndi+(ii-1))
          end do
        elseif(pcomp('pstr',lct,4)) then
          if(.not.fl(11)) then
           if (plfl) then
             setval = palloc ( 58,'NDNP',numnp*npstr,2)
             setval = palloc ( 57,'NDER',numnp*8    ,2)
             setval = palloc ( 60,'NDNS',max(nen*npstr,nst*nst),2)
             setval = palloc (207,'NSCR',numel      ,2)
             plfl = .false.
           endif
           nph  = np(58)
           ner  = np(57)
           npnp = nph
           nper = ner
           call pjstrs(.false.)
           fl(11) = .true.
          endif
          ndi = max(1,min(nint(ct(1)),7))
          write(ios,*) ttim, ' PSTRESS',ndi
           do ii = 1,numpn
             write(ios,*) mr(np(244)+ii-1),
     &        hr(np(57)+numnp*ndi+(ii-1))
           end do
        endif

c       Close
        close(ios)
       endif

c     [global <on,off, > - Global output flag

      elseif(j.eq.8) then

        if(pcomp(lct,'off',3)) then
          pfeap_glob = .false.
        else
          pfeap_glob = .true.
        endif

c     [local             - Local output flag

      elseif(j.eq.9) then

        pfeap_glob = .false.

      endif

c     Formats

2000  format(/5x,'Specify number of domains. > ',$)

2001  format(/5x,'Use METIS to partition mesh into',i5,' domains'/)

2002  format( '  PETSC Solver'/)
2003  format( '  PETSC Viewer: ',a)

2004  format('   Number of eigenpairs reduced to',i4,' by number of',
     &       ' nonzero matrix diagonal terms')

2005  format('   O u t p u t   G l o b a l   N o d e   L i s t',i5)

2006  format('   O u t p u t   L o c a l   N o d e   L i s t'/
     &       2x,8(i5,'Node':))

2007  format(2x,8i9)
3002  format(/' *ERROR* ',a,': No stiffness matrix, use TANG or UTAN'/)

3003  format(/' *ERROR* PSUB: No mass matrix, use MASS or IDEN'/)

3009  format(/' *ERROR* PSUB: Mass has no terms. May need to specify',
     &        ' a density.'/)

3010  format(/' *ERROR* PSUB: Geometric tangent has no terms.  Compute'
     &       /'               solution before GEOM to get values.'/)

3011  format(/' *WARNING* Block size not correct.  Reset to size =',i4)

      end

      subroutine size_ixinv(ix,v,ixinvlen)

c-----[--+---------+---------+----------------------------------------]
c     Purpose: Compute size of inverse list

c     Inputs:  ix(*)    - Array of nodes on elements

c     Outputs: v(*)     - Array of pointers to start of records for ixinv
c              ixinvlen - Length of ixinv array
c-----[--+---------+---------+----------------------------------------]
      implicit   none

#     include   "cdata.h"
#     include   "sdata.h"

#     include   "c_0.h"
#     include   "c_contac.h"
#     include   "c_ccp.h"

#     include   "pointer.h"
#     include   "comblk.h"

      logical   setvar, palloc

      integer   e,j,n,ixinvlen,k
      integer   ix(nen1,numel),v(*)

      do e = 1, numel
        do j = 1,nen
          n = ix(j,e)
          if (n.ne.0) v(n+1) = v(n+1) + 1
        end do ! j
      end do ! e

      if(ifct) then
        do j = 0,mr(np(143)+cck(3))-2  ! -1 -1
          n = mr(np(141)+j)
          v(n+1) = v(n+1) + 1
        end do
        do j = 0,mr(np(144)+cck(3))-2  ! -1 -1
          n = mr(np(142)+j)
          v(n+1) = v(n+1) + 1
        end do
      endif

      ixinvlen = 0
      v(1)     = 1
      do j = 2, numnp
        ixinvlen = ixinvlen + v(j)
        v(j)     = v(j-1) + v(j)
      end do ! j
      ixinvlen   = ixinvlen + v(numnp+1)
      v(numnp+1) = ixinvlen + 1

      end

      subroutine create_ixinv(ix,v,skip,ixinv)

c-----[--+---------+---------+----------------------------------------]
c     Purpose: Create ix-inverse in compressed form from pointer list
c              in v(*)

c     Inputs:  ix(*)    - Array of nodes on elements
c              v(*)     - Array of pointers to start of records for ixinv

c     Output:  ixinv(*) - Array of inverse ix list
c-----[--+---------+---------+----------------------------------------]
      implicit   none

#     include   "cdata.h"
#     include   "sdata.h"

#     include   "c_0.h"
#     include   "c_contac.h"
#     include   "c_ccp.h"

#     include   "pointer.h"
#     include   "comblk.h"

      logical   setvar, palloc


      integer   e,j,n,k
      integer   ix(nen1,numel),v(*),skip(*),ixinv(*)

      do e = 1, numel
        do j = 1, nen
          n = ix(j,e)
          if (n.ne.0) then
            ixinv(v(n)+skip(n)) = e
            skip(n)             = skip(n) + 1
          endif
        end do ! j
      end do ! e

      if(ifct) then
        do j = 1,cck(3)
          do k = mr(np(143)+j-1),mr(np(143)+j)-1
            n = mr(np(141)+k-1)
            ixinv(v(n)+skip(n)) = numel+j
            skip(n) = skip(n) + 1
          end do
          do k = mr(np(144)+j-1),mr(np(144)+j)-1
            n = mr(np(142)+k-1)
            ixinv(v(n)+skip(n)) = numel+j
            skip(n) = skip(n) + 1
          end do
        end do
      endif

      end

      subroutine make_node_graph(ix,v,ixinv,point,nodeg,len,numedg,
     &                           ndlast)

c-----[--+---------+---------+----------------------------------------]
c     Purpose: Build nodal graph for use by Metis

c     Inputs:  ix(*)    - Node connection list
c              v(*)     - Pointer list
c              ixinv(*) - Inverse ix(*) list

c     Outputs: nodeg(*) - Node graph
c              point(*) - Pointer array for node graph
c              len      - Length of node graph array
c              numedg   - Number of edges in graph  (for Windows only)
c-----[--+---------+---------+----------------------------------------]
      implicit   none

#     include   "cdata.h"
#     include   "sdata.h"

#     include   "c_0.h"
#     include   "c_contac.h"
#     include   "c_ccp.h"

#     include   "pointer.h"
#     include   "comblk.h"

      logical   setvar, palloc, slaveflg

      logical   insert
      integer   pairn
      integer   e,n,i,cn,k, len,skip,t, numedg, ndlast, istrt, iend
      integer   ipass
      integer   ix(nen1,numel),v(*),ixinv(*),point(*),nodeg(*)

      point(1) = 1
      numedg   = 0
      ndlast   = 1


      do n = 1, numnp
        skip = 0
        do k = v(n), v(n+1)-1

          e = ixinv(k)
          if(e.le.numel) then
            do i = 1, nen
              insert = .true.
              cn     = ix(i,e)

              do t = 0, skip
                if(cn.eq.nodeg(point(n)+t) .or.
     &             cn.eq.n .or. cn.eq.0) insert=.false.
              end do ! t

              if(insert) then
                nodeg(point(n)+skip) = cn
                skip                 = skip + 1
                if(cn.gt.n) then
                  numedg = numedg + 1
                endif
              endif

            end do ! i
          else
            pairn = e - numel

c           Check if node n is on the slave side or the master side
            slaveflg = .true.
c           do i = mr(np(144)+pairn-1),mr(np(144)+pairn)-1    XX Switch to conne
c             if(mr(np(142)+i-1).eq.n) slaveflg = .false.     XX Switch to conne
c           end do                                            XX Switch to conne

            do ipass=1,2

c           Set loop limits
            if(slaveflg) then
              istrt = mr(np(144)+pairn-1)   ! When slave, loop master
              iend  = mr(np(144)+pairn)-1
            else
              istrt = mr(np(143)+pairn-1)   ! When master, loop slave
              iend  = mr(np(143)+pairn)-1
            endif

            do i = istrt,iend
              insert = .true.
              if(slaveflg) then
                cn = mr(np(142)+i-1)
              else
                cn = mr(np(141)+i-1)
              endif
              do t = 0, skip
                  if(cn.eq.nodeg(point(n)+t) .or.
     &               cn.eq.n .or. cn.eq.0) insert=.false.
              end do ! t

              if(insert) then
                nodeg(point(n)+skip) = cn
                skip                 = skip + 1
                if(cn.gt.n) then
                  numedg = numedg + 1
                endif
              endif
            end do ! i
            slaveflg = .false. ! switch from slave side to master side
            end do ! pass over slave,master sides
          endif
        end do ! k

        if(point(n)+skip.gt.len) then
          len    = point(n)
          ndlast = n-1
          return
        endif
        point(n+1) = point(n) + skip

      end do ! n

      len    = point(numnp+1) - 1
      ndlast = numnp

      end

      subroutine uoutdom(eq,id,ix,ld,x,ang,f,fpro,ndam,nmas,nsti,
     &                   partn,nodx,revp,revx,domn,domeq,domd,
     &                   igets,dact,domains)
c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove dead code associated with bflag           05/01/2013
c       2. Separate blocking from bc-eqs being in           05/01/2013
c       3. Call contact (317) to output contact info        29/01/2013
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Output individual domain input files for parallel solution

c     Inputs:
c       eq(*)      - Equation numbers
c       id(*)      - Boundary condition data
c       ix(*)      - Element connection data
c       x(*)       - Nodal coordinate data
c       ang(*)     - Angle of boundary condition slope
c       f(*)       - Nodal force/displacements
c       fpro(*)    - Proportional load values for nodes
c       ndam(*)    - Nodal damping   values
c       nmas(*)    - Nodal mass      values
c       nsti(*)    - Nodal stiffness values
c       domains    - Number of domains

c     Scratch:
c       ld(*)      - Local array for list
c       nodx(*)    - List of local nodes for domain
c       partn(*)   - Nodal partition list
c       revp(4,*)  - Number of (1) domain nodes
c                              (2) total nodes
c                              (3) Pointer for GET  data
c                              (4) Pointer for SEND data
c       revx(2,*)  - List of (1) global nodes for all domains
c                            (2) local  nodes for all domains
c       domn(2,*)  - List of (1) Elements in domain
c                            (2) Ghost elements
c       domeq(2,*) - Number of equations (1) In current domain
c                                        (2) In other domains
c       domd(2,*)  - Range of equation numbers for each domain
c       igets(*)   - List of values for get/send
c       dact(*,*)  - Number get/sends for each domain pair

c     Outputs:
c                  - Files with input data for individual domains
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

#     include   "cdata.h"
#     include   "chdata.h"
#     include   "comfil.h"
#     include   "compac.h"
#     include   "corset.h"
#     include   "ddata.h"
#     include   "debugs.h"
#     include   "iodata.h"
#     include   "iofile.h"
#     include   "pglob1.h"
#     include   "pfeapb.h"
#     include   "prflag.h"
#     include   "refng.h"
#     include   "sdata.h"
#     include   "setups.h"
#     include   "pointer.h"
#     include   "comblk.h"

#     include   "p_int.h"

      logical    dflag,eflag,gflag, lread,pflg, setval,palloc
      logical    binary
      character  fnamr*134,fext*8, eformat*26, size*2, str*141
      character  fnamt*128
      character  twd(10)*4
      integer    domains, d, e,ee, n,ne,nf,ng,nq, nfm,nstr,nstf, ze
      integer    neg,nep,nrv, dneg,dnep, dned, nnu, matd,mato, j,jj
      integer    tnsbk
      integer    eq(ndf,numnp),id(ndf,numnp),ix(nen1,*), ld(*)
      integer    fpro(ndf,*)
      integer    partn(*),nodx(0:*),revp(4,*),revx(2,*),domn(2,*)
      integer    domeq(2,*),domd(2,*), igets(*), dact(domains,domains)
      real*8     x(ndm,*), ang(*), f(ndf,numnp,*)
      real*8     ndam(ndf,*),nmas(ndf,*),nsti(ndf,*), eula(3)

      real*4     etime,tary(2)

      save

      data       ze     /  0  /
      data       binary / .false. /
      data       twd    / 'NEWM', 'BACK', 'HHT ', 'EXPL', 'CONS',
     &                    'STAT', 'GEN1', 'CENT', 'BDF2', 'EULE' /

      nrv = 0
      do n = 1,domains
        nrv = nrv + revp(2,n)
      end do ! n
      if(debug .and. ndebug.gt.1) then
        call iprint(revx,2,nrv,2,'REVX-tot')
      endif

c     Set equation numbers for domains

      neg = 0
      do d = 1,domains
        domd(1,d) = neg + 1
        if(pfeap_bcin) then             ! Equation numbers blocked
          do n = 1,numnp
            if(partn(n).eq.d) then
c             Active equation case: Number all equations at node
              do e = 1,ndf
                neg = neg + 1
                eq(e,n) = neg
              end do ! e
            endif
          end do ! n
        else                            ! Equation numbers not blocked
          do n = 1,numnp
            if(partn(n).eq.d) then
              do e = 1,ndf
                if(id(e,n).eq.0) then
                  neg = neg + 1
                  eq(e,n) = neg
                else
                  eq(e,n) = 0
                endif
              end do ! e
            endif
          end do ! n
        endif
        domd(2,d) = neg
      end do ! d

      neq = neg  ! Set total number of equations

      if(debug .and. ndebug.gt.1) then
        call iprint(domd,2,domains,2,'DOMD_eq')
      endif

c     Compute active domain table

      if(debug) then
        write(  *,*) ' Start DACT domains',etime(tary)
        write(iow,*) ' Start DACT domains',etime(tary)
      endif
      call udact(partn,revp,revx,dact, domains)
      if(debug) then
        write(  *,*) ' End   DACT domains',etime(tary)
        write(iow,*) ' End   DACT domains',etime(tary)
      endif

c     Zero element partition

      do n = 1,numel
        domn(1,n) = 0
        domn(2,n) = 0
      end do ! n

c     Loop over domains

      neg = 0
      nep = 0
      nrv = 0
      do d = 1,domains
        ndomn = d

        dneg = 0
        dnep = 0

c       Assign name for domain mesh files

        fnamr = finp
        do n = 1,128
          if(fnamr(n:n).eq.'.') then
            fnamr(n:128) = ' '
            exit
          endif
        end do ! n
        nstr = index(fnamr,' ')

c       Set string for solution file name

        str( 1: 14)     = 'INCLude solve.'
        str(15:13+nstr) =  fnamr(2:nstr)

c       Add extender to domain mesh files

        if(d.lt.10) then
          write(fext,'(a3,i1)') '000',d
        elseif(d.lt.100) then
          write(fext,'(a2,i2)') '00',d
        elseif(d.lt.1000) then
          write(fext,'(a1,i3)') '0',d
        elseif(d.lt.10000) then
          write(fext,'(i4)') d
        endif

c       Add extender and replace '.' by '_' in file name

        call addext(fnamr,fext,128,6)
        n = index(fnamr,'.')
        if(n.gt.0) then
          fnamr(n:n) = '_'
        endif

        n    = index(fnamr,' ')
        nstf = n-1
        fnamt(1:nstf)   =  fnamr(1:nstf)
        write(  *,*) ' '
        write(iow,*) ' '
        write(  *,*) ' Mesh output:  Filename = ',fnamr(1:n)
        write(iow,*) ' Mesh output:  Filename = ',fnamr(1:n)

c       Open domain mesh file

        call opnfil(fext,fnamr,-1,ios,lread)
        rewind ios

c       Set node numbers for this domain

        if(debug) then
          write(  *,*) ' Start NODE numbers',etime(tary)
          write(iow,*) ' Start NODE numbers',etime(tary)
        endif

        do n = 0,numnp
          nodx(n) = 0
        end do ! n

        ne = 0
        do n = 1,numel
          eflag   = .false.
          do e = 1,nen
            if(ix(e,n).gt.0) then
              if(partn(ix(e,n)).eq.d) then
                eflag         = .true.
              endif
            endif
          end do ! e

c         Count number of elements and set ghost nodes

          if(eflag) then
            ne = ne + 1
          endif
        end do ! n

c       Set number of nodes in partition

        nq = revp(1,d)       !  Nodes for active partition
        nf = revp(2,d)       !  Nodes for ghost values

        write(  *,*) ' Number partition nodes = ',nq,
     &               ' Number total nodes =',numnp
        write(iow,*) ' Number partition nodes = ',nq,
     &               ' Number total nodes =',numnp

        if(debug) then
          write(  *,*) ' End   NODE numbers',etime(tary)
          write(iow,*) ' End   NODE numbers',etime(tary)
        endif

c       Output control data

        write(ios,2000) 'NOCOunt ! # nodes, elements matl specified'

        write(ios,2000) 'feap * * Domain ',d, nq,
     &                   nf,ne,nummat,ndm,ndf,nen

        write(ios,2000) 'NOPRint ! Do not echo data to output file'

c       Check for active global parameters

        if(gtypfl.or.gdeffl.or.gomgfl.or.
     &     gtdofl.or.grayfl.or.groufl) then
          write(ios,2000) 'GLOBal parameters'

c         Output 2-d solution type

          if(gtypfl) then
            if(g2type.eq.1) then
              write(ios,2008) '  PLANe STREss'
            elseif(g2type.eq.2) then
              write(ios,2008) '  PLANe STRAIN'
            elseif(g2type.eq.3) then
              write(ios,2008) '  AXISymmetric'
            elseif(g2type.eq.8) then
              write(ios,2008) '  AXISymmetric TORSion'
            endif
          endif

c         Output deformation type: Small or Finite

          if(gdeffl) then
            if(gdtype.eq.1) then
              write(ios,2008) '  SMALl deformation'
            elseif(gdtype.eq.-1) then
              write(ios,2008) '  FINIte deformation'
            endif
          endif

c         Output rotational velocity parameters

          if(gomgfl) then
              write(ios,2008) '  OMEGa radians    ',gomega(1)
              write(ios,2008) '  OMEGa COORdinate ',(gomex(n),n=1,ndm)
              write(ios,2008) '  OMEGa VECTor     ',(gomev(n),n=1,ndm)
          endif

c         Output thermal dof

          if(gtdofl) then
            write(iow,2008) '  TEMPerature DOF ',gtdof
          endif

c         Output rayleigh damping parameters
          if(grayfl) then
            write(iow,2008) '  RAYLeigh damping ',gray(1),gray(2)
          endif

c         Output ground/group parameters

          if(groufl) then
            write(iow,2008) '  GROUp  factors ',(gfac(n),n=1,ndf)
          endif
        endif

c       Output node list

        if(debug) then
          write(  *,*) ' Start MESH numbers',etime(tary)
          write(iow,*) ' Start MESH numbers',etime(tary)
        endif

        write(  *,*) ' Start COORD: TIME =',etime(tary)
        write(iow,*) ' Start COORD: TIME =',etime(tary)

c       Output node list

        if(debug) then
          write(  *,*) ' Start MESH numbers',etime(tary)
          write(iow,*) ' Start MESH numbers',etime(tary)
        endif

        write(  *,*) ' Start COORD: TIME =',etime(tary)
        write(iow,*) ' Start COORD: TIME =',etime(tary)

c       Binary file with nodal data

        if(binary) then
          fnamt(1:1) =  'C'
          write(ios,2007) 'COORDinate BINAry = ',fnamt(1:nstf)
          open(unit = iop, file = fnamt(1:nstf), form = 'unformatted')
          write(iop) ((x(e,revx(1,nrv+n)),e=1,ndm),n=1,nf)
          close(iop)

c       ASCII form of nodal data

        else
          write(ios,2000) 'COORdinate ALL ! All nodes in list'
          do n = 1,nf
            ng = revx(1,nrv+n)
            write(ios,2001) revx(2,nrv+n),ze,(x(e,ng),e=1,ndm)
          end do ! n
        endif

        write(  *,*) ' End   COORD: TIME =',etime(tary)
        write(iow,*) ' End   COORD: TIME =',etime(tary)

c       Build remap for element outputs

        do n = 1,numnp
          nodx(n) = 0
        end do ! n

        do n = 1,nf
          ng       = revx(1,nrv+n)
          nodx(ng) = revx(2,nrv+n)
        end do ! n

c       Output element list

        write(  *,*) ' Start ELEMT: TIME =',etime(tary),nodx(0)
        write(iow,*) ' Start ELEMT: TIME =',etime(tary),nodx(0)

c       Binary files with element data (Needs corrections)

        if(binary) then
          fnamt(1:1) = 'E'
          write(ios,2000) 'ELEMents BINAry = ',fnamt(1:nstf)
          open(unit = iop, file = fnamt(1:nstf), form = 'unformatted')
          write(iop) ((ix(e,revx(1,nrv+n)),e=1,ndm),n=1,nf)
          close(iop)

c       ASCII data in current file

        else

c         Build output format for elements that fits data

c         Format:    1...5...10...15...20...25...30
          eformat = '(i10,i10,i6,13i10/(16i10))'

          nfm = nint(log10(dble(max(nf,ne)))) + 2
          write(size,'(i2)') nfm
          eformat( 3: 4) = size
          eformat( 7: 8) = size
          eformat(16:17) = size
          eformat(23:24) = size

c         Output element list

          write(ios,2000) 'ELEMents ALL ! All elements in list'
          ne = 0
          do n = 1,numel
            eflag = .false.
            dflag = .false.
            gflag = domn(1,n).eq.0
            domn(2,n) = 0
            dned = 0
            do e = 1,nen
              if(ix(e,n).gt.0) then
                ee = e
                if(partn(ix(e,n)).eq.d) then
                  dned  = dned + 1
                  eflag = .true.
                elseif(gflag) then
                  dflag = .true.
                endif
              endif
            end do ! e
            if(eflag) then
              if(dflag) then
                dneg = dneg + 1
              elseif(dned.eq.ee) then
                dnep = dnep + 1
              endif
              ne = ne + 1
              write(ios,eformat) ne,ze,ix(nen1,n),(nodx(ix(e,n)),e=1,ee)
     &                                ,(0,e=ee+1,nen)
              if(domn(1,n).eq.0 .and. dflag) then
                domn(1,n) = -d
                domn(2,n) = ne
              endif
            endif
          end do ! n
        endif

        write(  *,*) ' End   ELEMT: TIME =',etime(tary),nodx(0),ne
        write(iow,*) ' End   ELEMT: TIME =',etime(tary),nodx(0),ne


        write(  *,*) ' Ghost  Stress Elements = ',dneg,
     &               ' Partn Stress Elements = ',dnep
        write(iow,*) ' Ghost  Stress Elements = ',dneg,
     &               ' Partn Stress Elements = ',dnep
        neg = neg + dneg
        nep = nep + dnep

        if(debug .and. ndebug.gt.1) then
          call iprint(domn,2,numel,2,'DOMN')
        endif

c       Element region for stress output

        eflag = .true.
        do n = 1,numel
          if(domn(1,n).eq.-d) then
            if(eflag) then
              eflag = .false.
              write(ios,2000) 'EREGions  ! Region 21 ignored for STREss'
            endif
            write(ios,2002) domn(2,n),ze,21
          endif
        end do ! n

c       Angle conditions lists

        if(anglefl) then
          lread = .true.
          do n = 1,nf
            ng = revx(1,nrv+n)
            if(partn(ng).gt.0 .and. ang(ng).ne.0.0d0) then
              if(lread) then
                write(ios,2000) 'ANGLe conditions'
                lread = .false.
              endif
              write(ios,2001) revx(2,nrv+n),ze,ang(ng)
            end if
          end do ! n
        endif

c       Euler angle conditions lists

        if(eulerfl) then
          lread = .true.
          do n = 1,nf
            ng = revx(1,nrv+n)
            if(partn(ng).gt.0) then
              do e = 1,3
                eula(e) = hr(np(242)+3*(ng-1)+e)
              end do ! e
              if(eula(1).ne.0.0d0 .or.
     &           eula(2).ne.0.0d0 .or.
     &           eula(3).ne.0.0d0) then
                if(lread) then
                  write(ios,2000) 'EULEr angle conditions'
                  lread = .false.
                endif
                write(ios,2001) revx(2,nrv+n),ze,(eula(e),e=1,3)
              end if
            endif
          end do ! n
        endif

c       Boundary conditions lists

        if(debug .and. ndebug.gt.1) then
          call iprint(id,ndf,numnp,ndf,'ID_bc')
        endif
        lread = .true.
        do n = 1,nf
          ng = revx(1,nrv+n)
          if(partn(ng).gt.0) then
            pflg = .false.
            do e = 1,ndf
              if(id(e,ng).ne.0) then
                pflg = .true.
                ld(e) = 1
              else
                ld(e) = 0
              end if
            end do ! e

            if(pflg) then
              if(lread) then
                write(ios,2000) 'BOUNdary conditions'
                lread = .false.
              endif
              write(ios,2002) revx(2,nrv+n),ze,(ld(e),e=1,ndf)
            endif
          end if
        end do ! n

c       Forced conditions lists

        lread = .true.
        do n = 1,nq
          ng = revx(1,nrv+n)
          if(partn(ng).gt.0) then
            pflg = .false.
            do e = 1,ndf
              if(f(e,ng,1).ne.0.0d0) then
                pflg = .true.
              end if
            end do ! e
            if(pflg) then
              if(lread) then
                write(ios,2000) 'FORCe conditions'
                lread = .false.
              endif
              write(ios,2003) revx(2,nrv+n),ze,(f(e,ng,1),e=1,ndf)
            endif
          end if
        end do ! n

c       Displacement conditions lists

        lread = .true.
        do n = 1,nf
          ng = revx(1,nrv+n)
          if(partn(ng).gt.0) then
            pflg = .false.
            do e = 1,ndf
              if(f(e,ng,2).ne.0.0d0) then
                pflg = .true.
              end if
            end do ! e
            if(pflg) then
              if(lread) then
                write(ios,2000) 'DISPlacement conditions'
                lread = .false.
              endif
              write(ios,2003) revx(2,nrv+n),ze,(f(e,ng,2),e=1,ndf)
            end if
          end if
        end do ! n

c       Proportional load conditions lists

        lread = .true.
        do n = 1,nf
          ng = revx(1,nrv+n)
          if(partn(ng).gt.0) then
            pflg = .false.
            do e = 1,ndf
              if(fpro(e,ng).ne.0.0d0) then
                pflg = .true.
              end if
            end do ! e
            if(pflg) then
              if(lread) then
                write(ios,2000) 'FPROportional load'
                lread = .false.
              endif
              write(ios,2002) revx(2,nrv+n),ze,(fpro(e,ng),e=1,ndf)
            end if
          end if
        end do ! n

c       Base proportional load conditions lists

        if(np(125).ne.0) then
          if(debug .and. ndebug.gt.1) then
            call iprint(mr(np(125)),ndf,numnp,ndf,'Base')
          endif
          lread = .true.
          do n = 1,nf
            ng = revx(1,nrv+n)
            if(partn(ng).gt.0) then
              pflg = .false.
              fp(1) = np(125) + ndf*(ng - 1) - 1
              do e = 1,ndf
                if(mr(fp(1)+e).ne.0) then
                  pflg = .true.
                  ld(e) = 1
                else
                  ld(e) = 0
                end if
              end do ! e

              if(pflg) then
                if(lread) then
                  write(ios,2000) 'BASE conditions'
                  lread = .false.
                endif
                write(ios,2002) revx(2,nrv+n),ze,(ld(e),e=1,ndf)
              endif
            end if
          end do ! n
        end if

c       Add loads from ldtab

        if(np(265).ne.0) then
          call uldout(mr(np(265)),mr(np(266)),hr(np(267)),
     &                nrv, nf, partn, revx)
        endif

c       Check for existing nodal dampers, masses and stiffnesses

        if(nmfl) then

c         Output nodal damper values

          lread = .true.
          do n = 1,nq
            ng = revx(1,nrv+n)
            if(partn(ng).gt.0) then
              pflg = .false.
              do e = 1,ndf
                if(ndam(e,ng).ne.0.0d0) then
                  pflg = .true.
                end if
              end do ! e
              if(pflg) then
                if(lread) then
                  write(ios,2000) 'DAMPer nodal'
                  lread = .false.
                endif
                write(ios,2003) revx(2,nrv+n),ze,(ndam(e,ng),e=1,ndf)
              end if
            end if
          end do ! n

c         Output nodal mass values

          lread = .true.
          do n = 1,nq
            ng = revx(1,nrv+n)
            if(partn(ng).gt.0) then
              pflg = .false.
              do e = 1,ndf
                if(nmas(e,ng).ne.0.0d0) then
                  pflg = .true.
                end if
              end do ! e
              if(pflg) then
                if(lread) then
                  write(ios,2000) 'MASS nodal'
                  lread = .false.
                endif
                write(ios,2003) revx(2,nrv+n),ze,(nmas(e,ng),e=1,ndf)
              end if
            end if
          end do ! n

c         Output nodal stiffness values

          lread = .true.
          do n = 1,nq
            ng = revx(1,nrv+n)
            if(partn(ng).gt.0) then
              pflg = .false.
              do e = 1,ndf
                if(nsti(e,ng).ne.0.0d0) then
                  pflg = .true.
                end if
              end do ! e
              if(pflg) then
                if(lread) then
                  write(ios,2000) 'STIFness nodal'
                  lread = .false.
                endif
                write(ios,2003) revx(2,nrv+n),ze,(nsti(e,ng),e=1,ndf)
              end if
            end if
          end do ! n

        endif ! nmfl = .true.

        if(debug) then
          write(  *,*) ' End   MESH numbers',etime(tary)
          write(iow,*) ' End   MESH numbers',etime(tary)
        endif

c       Retrieve material data set entries from file "fmtl"
c       Output material/parameter lists

        open(unit=iwd,file=fmtl,status='old')
        pflg = .true.
        write(ios,'(a)') ' '
        do while(pflg)
          read(iwd,'(a)',end=200) xxx
          do e = 255,1,-1
            if(xxx(e:e).ne.' ') go to 100
          end do ! e
          e = 1
100       if(xxx(e:e).eq.char(13)) xxx(e:e) = ' '
          write(ios,'(a)') xxx(1:e)
        end do ! while
200     close(iwd)

c       Write end of mesh

        write(ios,2000) 'END MESH'

c       Output contact surfaces and pairs

        call contact (317)

c       Output domain information

        write(ios,2005) 'DOMAin number ',d,domains,nq,numnp,neq

c       Write blocked equation form and block size

        if(pfeap_blk) then
          write(ios,2000) 'FORMat BAIJ'
        else
          write(ios,2000) 'FORMat AIJ'
        endif
        if(pfeap_bcin) then
          write(ios,2000) 'BCIN BLOCked',nsbk
        endif

c       Output global/local node map

        write(ios,2000) 'LOCAl to GLOBal node numbers'
        do n = 1,nf
          ng = revx(1,nrv+n)
          write(ios,2004) n,revx(1,nrv+n),eq(ndf,ng)/ndf
        end do ! n

c       Compute get/send lists

        if(debug) then
          write(  *,*) ' Start GETS numbers',etime(tary)
          write(iow,*) ' Start GETS numbers',etime(tary)
        endif

        call ugetsend(d,nrv,nodx,revp,revx,igets,dact,domains)

        if(debug) then
          write(  *,*) ' End   GETS numbers',etime(tary)
          write(iow,*) ' End   GETS numbers',etime(tary)
        endif

c       Compute global equation numbers

        if(debug) then
          write(  *,*) ' Start EQNS numbers',etime(tary)
          write(iow,*) ' Start EQNS numbers',etime(tary)
        endif

c       Case for multiple domains

        if(domains.gt.1) then
          call ueqn1(eq,partn,domd,mr(np(252)),mr(np(253)), domeq,d,
     &               ndf,numnp,neq)

c       Case for a single domain

        else
          ubycol = .false.      ! Store by columns
          udiag  = .false.      ! Diagonal in separate part of array
          uall   = .true.       ! Assemble all of column

          call uiters(nnu,-1)
          call ueqn2(domd,mr(np(225)),mr(np(226)), domeq, d,neq)
          setval = palloc(226,'USOL2',0,1)
          setval = palloc(225,'USOL1',0,1)
        endif

c       Global matrix storage

        if(debug .and. ndebug.gt.1) then
          call iprint(domeq,2,neq,2,'DOM_eq')
        endif

c       Output matrix storage

        write(ios,2000) 'MATRix  storage'
        matd = 0
        mato = 0

c       Boundary condition equations in the assembly
c       BAIJ and AIJ

        if(pfeap_bcin) then

c         Divide current equation numbers by block size "nsbk" for BAIJ
c         else by 1 for AIJ since AIJ needs all equations

          tnsbk = nsbk
          if(.not.pfeap_blk) tnsbk = 1

          do n = tnsbk,nneq,tnsbk
            domeq(1,n) = domeq(1,n)/tnsbk
            domeq(2,n) = domeq(2,n)/tnsbk
          end do ! n

c         Set number of blocks/node

          jj = ndf/tnsbk

c         Check for nodes with active equations

          do n = 1,nq
            ng = revx(1,nrv+n)

            do j = 1,jj
              ee = eq(j*tnsbk,ng)

c             Test for boundary condition where "id(e,ng) .ne. 0"

              do e = (j-1)*tnsbk+1,j*tnsbk
                if(id(e,ng).eq.0) then
                  go to 300
                endif
              end do ! e

c             Reset domain equations for boundary block

c             domeq(1,ee) = 1
c             domeq(2,ee) = 0

c             Output equation information for node; count blocks

300           if(ee.gt.0) then
                write(ios,2006) ee/tnsbk,domeq(1,ee),domeq(2,ee)
                matd = matd + domeq(1,ee)
                mato = mato + domeq(2,ee)
              endif
            end do ! j

          end do ! n

c         Report blocking information
          if(pfeap_blk) then
            write(  *,*) ' Block Matrix Store: Diagonal =',matd,
     &                   ' Off-Diagonal =',mato
            write(iow,*) ' Block Matrix Store: Diagonal =',matd,
     &                   ' Off-Diagonal =',mato
          else
            write(  *,*) ' Matrix Store: Diagonal =',matd,
     &                   ' Off-Diagonal =',mato
            write(iow,*) ' Matrix Store: Diagonal =',matd,
     &                   ' Off-Diagonal =',mato
          endif

c       AIJ case w/o bc-eqs: Report equation and numbers of terms on
c                            current processor  "domeq(1,ee)" and on
c                            other   processors "domeq(2,ee)"
        else
          do n = 1,nq
            ng = revx(1,nrv+n)
            do e = 1,ndf
              ee = eq(e,ng)        ! Equation number
              if(ee.gt.0) then
                write(ios,2006) ee,domeq(1,ee),domeq(2,ee)
                matd = matd + domeq(1,ee)
                mato = mato + domeq(2,ee)
              endif
            end do ! e
          end do ! n

c         Report equation storage information

          write(  *,*) ' Matrix Store: Diagonal =',matd,
     &                 ' Off-Diagonal =',mato
          write(iow,*) ' Matrix Store: Diagonal =',matd,
     &                 ' Off-Diagonal =',mato

c         Global equation number map

          if(debug .and. ndebug.gt.1) then
            call iprint(eq(1,1),ndf,numnp,ndf,'ID_eq')
          endif

c         Output list of equation numbers for each node

          write(ios,2000) 'EQUAtion numbers'
          do n = 1,nf
            ng = revx(1,nrv+n)
            write(ios,2006) revx(2,nrv+n),(eq(e,ng),e=1,ndf)
          end do ! n
        endif

c       Write final data

        write(ios,2000) 'END DOMAIN'

        if(debug) then
          write(  *,*) ' End   EQNS numbers',etime(tary)
          write(iow,*) ' End   EQNS numbers',etime(tary)
        endif

c       Output lists initial conditions

        if(intlfl) then

c         Initial displacements

          if(idisfl) then
            write(ios,2000) 'BATCh  ! Initial conditions'
            write(ios,2010) 'INITial DISPlacement'
            write(ios,2008) 'END BATCH   ! Initial conditions'
            call uinitld(revx, hr(np(40)), nf,nrv, ndf, ios)
          endif

c         Initial velocity conditions

          if(ivelfl) then
            write(ios,2000) 'BATCh  ! Initial conditions'
            if(tranfl .and. noi.gt.0) then
              write(ios,2010) 'TRAN',twd(noi),(theta(n),n=1,3)
            endif
            write(ios,2010) 'INITial RATE'
            write(ios,2008) 'END BATCH   ! Initial conditions'
            call uinitld(revx, hr(np(42)), nf,nrv, ndf, ios)
          endif

        endif

c       Output solution include label

        write(ios,2000)  str(1:13+nstr)
        write(ios,2000) 'STOP'
        close(ios)

        nrv = nrv + revp(2,d)
      end do ! d

c     Final diagnostic information on how elements are assigned

      write(  *,*) ' '
      write(  *,*) ' Total Ghost Stress Elements =',neg
      write(  *,*) ' Total Partn Stress Elements =',nep
      write(  *,*) ' Total all   Stress Elements =',numel
      write(iow,*) ' '
      write(iow,*) ' Total Ghost Stress Elements =',neg
      write(iow,*) ' Total Partn Stress Elements =',nep
      write(iow,*) ' Total all   Stress Elements =',numel

      if(debug .and. ndebug.gt.1) then
        call iprint(domn,1,numel,1,'ELEMENT DOMAIN NUMBER')
      endif

c     Check for an existing solution file:
c     Creates interactive solve file as default

      str           = ' '
      str(1:6)      = 'solve.'
      str(7:nstr+4) = fnamr(2:nstr-1)
      inquire(file = str(1:nstr+4), exist = setval)
      if(.not.setval) then
        open( unit = ios, file = str(1:nstr+4) )
        write(ios,'(a)') ' '
        if(ntasks.le.1) then
          write(ios,'(a)') 'BATCh'
          write(ios,'(a)') '  PETSc ON'
          write(ios,'(a)') 'END'
        endif
        write(ios,'(a)') 'INTERactive'
        write(ios,'(a)') ' '
        close(unit = ios, status = 'keep')
        write(*, 2009) str(1:nstr+4)
      endif

c     Formats

2000  format(/a:,i4:,' Nodes in partition = ',i10/6i10)
2001  format( i9,i2,1p,10e16.8)
2002  format( i9,i4,14i4/(16i4))
2003  format( i9,i2,1p,14e16.8/(1p,16e16.8))
2004  format( 8i9)
2005  format(/a:,i4:,': Total number of domains = ',i10/3i10)
2006  format(16i9)
2007  format(/a,1p,3e14.6)
2008  format( a:,1p,3e14.6)
2009  format(/'  *INFORMATION*  File ',a,' written with INTERactive',
     &        ' solution only.')
2010  format('  ',a,' ',a:,1p,3e16.8)

      end

      subroutine ueqn1(eq,partn,domd,jc,ir, domeq, d,ndf,numnp,neq)

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Set number of entries in matrix for each node in domain

c     Inputs:
c       eq(*)      - Equation number data
c       partn(*)   - Nodal partitionlist
c       domd(2,*)  - Range of equation numbers for each domain
c       jc(*)      - Pointer array for sparse storage of matrix
c       ir(*)      - Row values for sparse matrix
c       d          - Domain number
c       ndf        - Number dof/node
c       numnp      - Number of problem nodes
c       neq        - Number of problem equations

c     Outputs:
c       domeq(2,*) - Number of equations (1) In current domain
c                                        (2) In other domains
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none
      integer    d, ndf,numnp,neq, m,n,nn, i,j,k
      integer    eq(ndf,*),partn(*),domd(2,*),jc(*),ir(*),domeq(2,*)

      do n = 1,neq
        domeq(1,n) = 0
        domeq(2,n) = 0
      end do ! n
      do n = 1,numnp
        if(partn(n).eq.d) then
          do m = 1,ndf
            nn = eq(m,n)
            if(nn.gt.0) then
              domeq(1,nn) = domeq(1,nn) + 1
              do k = 1,ndf
                i = eq(k,n)
                if(i.gt.0 .and. i.ne.nn) then
                  if(i.ge.domd(1,d) .and. i.le.domd(2,d)) then
                    domeq(1,nn) = domeq(1,nn) + 1
                  else
                    domeq(2,nn) = domeq(2,nn) + 1
                  end if
                endif
              end do ! k
              do j = jc(n),jc(n+1)-1
                do k = 1,ndf
                  i = eq(k,ir(j))
                  if(i.gt.0) then
                    if(i.ge.domd(1,d) .and. i.le.domd(2,d)) then
                      domeq(1,nn) = domeq(1,nn) + 1
                    else
                      domeq(2,nn) = domeq(2,nn) + 1
                    end if
                  endif
                end do ! k
              end do ! j
            endif
          end do ! m
        endif
      end do ! n

      end

      subroutine ueqn2(domd,jc,ir, domeq, d,neq)

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Set number of entries in matrix for single domain problem

c     Inputs:
c       domd(2,*)  - Range of equation numbers for each domain
c       jc(*)      - Pointer array for sparse storage of matrix
c       ir(*)      - Row values for sparse matrix
c       d          - Domain number
c       neq        - Number of problem equations

c     Outputs:
c       domeq(2,*) - Number of equations (1) In current domain
c                                        (2) In other domains
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none
      integer    d, neq, n, i,j,jb
      integer    domd(2,*),jc(*),ir(*),domeq(2,*)

      jb = 0
      do n = 1,neq
        domeq(1,n) = 0
        do j = jb+1,jc(n)
          i = ir(j)
            if(i.ge.domd(1,d) .and. i.le.domd(2,d)) then
              domeq(1,n) = domeq(1,n) + 1
            end if
        end do ! j
        domeq(2,n) = jc(n) - jb - domeq(1,n)
        jb         = jc(n)
      end do ! n

      end

      subroutine uoutd1(ix,partn,nodx,revp, domains, nrv)

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Compute number of domain and total nodes for all domains

c     Inputs:
c       ix(*,*)    - Element connection data
c       partn(*)   - Nodal partitionlist
c       domains    - Number of domains

c     Scratch:
c       nodx(*)    - List of local nodes for domain

c     Outputs:
c       revp(4,*)  - Number of (1) domain nodes (2) total nodes
c       nrv        - Total number of nodes in all domains (including
c                    ghost nodes)
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

#     include   "cdata.h"
#     include   "debugs.h"
#     include   "sdata.h"

#     include   "c_0.h"
#     include   "c_contac.h"
#     include   "c_ccp.h"

#     include   "pointer.h"
#     include   "comblk.h"

      logical    eflag, masterfl, slavefl
      integer    domains, d, e, n,nf,nn,nq, nrv, pairn
      integer    ne
      integer    ix(nen1,*), partn(*),nodx(0:*),revp(4,*)

      save

c     Loop over domains

      nrv = 0
      do d = 1,domains

c       Set node numbers for this domain: nodx(n) = 0 indicates a set
c       to exclude nodes.

        do n = 0,numnp
          nodx(n) = 0
        end do ! n

c       Check each element to determine if it has any node in current
c       partition 'd'  - eflag = .true. for any element to be included.

c       ne = 0
        do n = 1,numel
          eflag   = .false.
          do e = 1,nen
            if(ix(e,n).gt.0) then
              if(partn(ix(e,n)).eq.d) then
                eflag         = .true.
              endif
            endif
          end do ! e

c         Count number of elements and set ghost nodes
c         nodx(n) = 1 for any node in current partition (including
c         ghost nodes).

          if(eflag) then
            do e = 1,nen
              if(ix(e,n).gt.0) then
                nodx(ix(e,n)) = 1
              endif
            end do ! e
c           ne = ne + 1
          endif
        end do ! n

        if(ifct) then
          do pairn = 1,cck(3)

c           Check if slave node on pair is in partition
            slavefl = .false.
            do n = mr(np(143)+pairn-1),mr(np(143)+pairn)-1
              if(partn(mr(np(141)+n-1)).eq.d) slavefl = .true.
            end do

c           Add any slave nodes and the whole master list
            if(slavefl) then
              do n = mr(np(143)+pairn-1),mr(np(143)+pairn)-1
                ne = mr(np(141)+n-1)
                if(partn(ne).eq.d) nodx(ne)=1
              end do
              do n = mr(np(144)+pairn-1),mr(np(144)+pairn)-1
                nodx(mr(np(142)+n-1)) = 1
              end do
            endif

c          Check if master node on pair is in partition
            masterfl = .false.
            do n = mr(np(144)+pairn-1),mr(np(144)+pairn)-1
              if(partn(mr(np(142)+n-1)).eq.d) masterfl = .true.
            end do

c           Add any master nodes and the whole slave list
            if(masterfl) then
              do n = mr(np(144)+pairn-1),mr(np(144)+pairn)-1
                ne = mr(np(142)+n-1)
                if(partn(ne).eq.d) nodx(ne)=1
              end do
              do n = mr(np(143)+pairn-1),mr(np(143)+pairn)-1
                nodx(mr(np(141)+n-1)) = 1
              end do
            endif
          end do ! pairn
        endif

c       First count number of nodes belonging to partition using
c       the 'partn' partition array from (par)metis.

        nn = 0
        do n = 1,numnp
          if(partn(n).eq.d) then
            nn      = nn + 1
            nodx(n) = nn
          endif
        end do ! n

        nq = nn       ! Number of nodes in partition

c       Next determine number of ghost nodes to be added using 'nodx'
c       array values that are unity (1).  The 'nodx' array is reset to
c       indicate local node number of ghost nodes.

        do n = 1,numnp
          if(partn(n).ne.d .and. nodx(n).eq.1) then
            nn      = nn + 1
            nodx(n) = nn
          endif
        end do ! n

c       Compute total number partition + ghost nodes in 'd'

        nf = 0
        do n = 1,numnp
          if(nodx(n).gt.0) then
            nf = nf + 1
          endif
        end do ! n
        revp(1,d) = nq   ! Number nodes in partition 'd'
        revp(2,d) = nf   ! Number total (partition + ghost) nodes in 'd'
        nrv       = nrv + nf  ! Accumulate total nodes
        if(debug .and. ndebug.gt.1) then
          call iprint(nodx(1),1,numnp,1,'NODX-domains')
        endif
      end do ! d
      if(debug .and. ndebug.gt.1) then
        call iprint(revp,2,domains,4,'REVP-nodes_in_domain')
      endif

      end

      subroutine uoutd2(ix,partn,nodx,revp,revx, domains)

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Assemble list of global and local nodes for all domains

c     Inputs:
c       ix(*,*)    - Element connection data
c       partn(*)   - Nodal partition list
c       revp(4,*)  - Number of (1) domain nodes (2) total nodes
c       domains    - Number of domains

c     Scratch:
c       nodx(*)    - List of local nodes for domain

c     Outputs:
c       revx(2,*)  - List of (1) global nodes for all domains
c                            (2) local  nodes for all domains
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

#     include   "cdata.h"
#     include   "sdata.h"

#     include   "c_0.h"
#     include   "c_contac.h"
#     include   "c_ccp.h"

#     include   "pointer.h"
#     include   "comblk.h"

      logical    masterfl, slavefl

      logical    eflag
      integer    domains, d, e, n,ne,nn,nq, nrv, pairn
      integer    ix(nen1,*), partn(*),nodx(0:*),revp(4,*),revx(2,*)

      save

c     Loop over domains

      nrv = 0
      do d = 1,domains

c       Set node numbers for this domain

        do n = 0,numnp
          nodx(n) = 0
        end do ! n

        ne = 0
        do n = 1,numel
          eflag   = .false.
          do e = 1,nen
            if(ix(e,n).gt.0) then
              if(partn(ix(e,n)).eq.d) then
                eflag         = .true.
              endif
            endif
          end do ! e

c         Count number of elements and set ghost nodes

          if(eflag) then
            do e = 1,nen
              if(ix(e,n).gt.0) then
                nodx(ix(e,n)) = 1
              endif
            end do ! e
            ne = ne + 1
          endif
        end do ! n

        if(ifct) then
          do pairn = 1,cck(3)

c           Check if slave node on pair is in partition
            slavefl = .false.
            do n = mr(np(143)+pairn-1),mr(np(143)+pairn)-1
              if(partn(mr(np(141)+n-1)).eq.d) slavefl = .true.
            end do

c           Add any slave nodes and the whole master list
            if(slavefl) then
              do n = mr(np(143)+pairn-1),mr(np(143)+pairn)-1
                ne = mr(np(141)+n-1)
                if(partn(ne).eq.d) nodx(ne)=1
              end do
              do n = mr(np(144)+pairn-1),mr(np(144)+pairn)-1
                nodx(mr(np(142)+n-1)) = 1
              end do
            endif

c          Check if master node on pair is in partition
            masterfl = .false.
            do n = mr(np(144)+pairn-1),mr(np(144)+pairn)-1
              if(partn(mr(np(142)+n-1)).eq.d) masterfl = .true.
            end do

c           Add any master nodes and the whole slave list
            if(masterfl) then
              do n = mr(np(144)+pairn-1),mr(np(144)+pairn)-1
                ne = mr(np(142)+n-1)
                if(partn(ne).eq.d) nodx(ne)=1
              end do
              do n = mr(np(143)+pairn-1),mr(np(143)+pairn)-1
                nodx(mr(np(141)+n-1)) = 1
              end do
            endif
          end do ! pairn
        endif

c       Count number of nodes

        nn = 0
        do n = 1,numnp
          if(partn(n).eq.d) then
            nn      = nn + 1
            nodx(n) = nn
          endif
        end do ! n

        nq = nn

        do n = 1,numnp
          if(partn(n).ne.d .and. nodx(n).eq.1) then
            nn      = nn + 1
            nodx(n) = nn
          endif
        end do ! n

c       Get ordered list for nodx

        do n = 1,numnp
          if(nodx(n).gt.0) then
            revx(1,nrv+nodx(n)) = n
            revx(2,nrv+nodx(n)) = nodx(n)
          endif
        end do ! n
        nrv = nrv + revp(2,d)
      end do ! d

      end

      subroutine ugetsend(d,nrv,nodx,revp,revx,igets,dact,domains)

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Compute and output GET and SEND list data

c     Inputs:
c       d          - Domain number
c       nrv        - Pointer to revx list for domain 'd'
c       nodx(*)    - Global nodes in current domain
c       partn(*)   - Nodal partitionlist
c       revp(4,*)  - Number of (1) domain nodes (2) total nodes
c       revx(2,*)  - List of (1) global node number for all domains
c                            (2) local  node number for all domains
c       dact(*,*)  - Active domain get/send array
c       domains    - Number of domains

c     Scratch:
c       igets(*)   - List of values for get/send

c     Outputs:
c       revp(4,*)  - Number of (3) Pointer for GET  data
c                              (4) Pointer for SEND data
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

#     include   "cdata.h"
#     include   "iodata.h"

      integer    domains, d,dd, n,ns,nrv, m,ms,mrv, tot
      integer    nodx(0:*), revp(4,*), revx(2,*), igets(*)
      integer    dact(domains,domains)

c     Get data list

      do dd = 1,domains
        revp(3,dd) = 0
        revp(4,dd) = 0
      end do ! n

      mrv = 0
      tot = 0
      do dd = 1,domains                      ! Compare domain d to dd
        revp(3,dd) = tot
        if(dact(d,dd).gt.0) then
          do m = 1,numnp
            nodx(m) = 0
          end do ! m

          do m = 1,revp(1,dd)
            ms = revx(1,mrv+m)
            nodx(ms) = -1
          end do ! m

          do m = revp(1,d)+1,revp(2,d)
            ns = revx(1,nrv+m)
            if(nodx(ns).eq.-1) then
              nodx(ns) = revx(2,nrv+m)
            endif
          end do

          do n = revp(1,d)+1,revp(2,d)       ! Check ghost nodes for d
            ns = revx(1,nrv+n)               ! Global node of domain d
            ms = nodx(ns)
            if(ms.gt.0) then
              tot        = tot + 1
              igets(tot) = ms
            endif
          end do ! n
        end if
        revp(3,dd) = tot - revp(3,dd)
        mrv        = mrv + revp(2,dd)
      end do ! dd

c     Output get pointer table

      if(tot.gt.0) then
        write(ios,2000) 'GETData','POINter',domains
        write(ios,2001) (revp(3,dd), dd=1,domains)

c       Output list of values

        write(ios,2000) 'GETData','VALUes',tot
        write(ios,2001) (igets(n), n=1,tot)
      endif

c     Send data list

      mrv = 0
      tot = 0
      do dd = 1,domains                      ! Compare domain d to dd
        revp(4,dd) = tot
        if(dact(dd,d).gt.0) then
          do m = 1,numnp
            nodx(m) = 0
          end do ! m

          do m = revp(1,dd)+1,revp(2,dd)
            ms = revx(1,mrv+m)
            nodx(ms) = -1
          end do ! m

          do m = 1,revp(1,d)
            ns = revx(1,nrv+m)
            if(nodx(ns).eq.-1) then
              nodx(ns) = revx(2,nrv+m)
            endif
          end do

          do n = 1,revp(1,d)                 ! Check local nodes for d
            ns = revx(1,nrv+n)               ! Global node of domain d
            ms = nodx(ns)
            if(ms.gt.0) then
              tot        = tot + 1
              igets(tot) = ms
            endif
          end do ! n
        end if
        revp(4,dd) = tot - revp(4,dd)
        mrv        = mrv + revp(2,dd)
      end do ! dd

c     Output send pointer table

      if(tot.gt.0) then
        write(ios,2000) 'SENDdata','POINter',domains
        write(ios,2001) (revp(4,dd), dd=1,domains)

c       Output list of values

        write(ios,2000) 'SENDdata','VALUes',tot
        write(ios,2001) (igets(n), n=1,tot)
      endif

c     formats

2000  format(/a,1x,a,i10)
2001  format(i10)

      end

      subroutine udact(partn,revp,revx,dact, domains)

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Compute domain lists requiring get/send data

c     Inputs:
c       partn(*)   - Nodal partitionlist
c       revp(4,*)  - Number of (1) domain nodes (2) total nodes
c       revx(2,*)  - List of (1) global node number for all domains
c                            (2) local  node number for all domains
c       domains    - Number domains

c     Outputs:
c       dact(*,*)  - Number get/send items in domains
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

#     include   "debugs.h"

      integer    domains, d,dd, n,nrv
      integer    partn(*),revp(4,*),revx(2,*), dact(domains,domains)


c     Get data list

      nrv = 0
      do d = 1,domains
        do dd = 1,domains
          dact(d,dd) = 0
        end do ! dd

c       Check for active get domains

        do n = revp(1,d)+1,revp(2,d)
          dd = partn(revx(1,nrv+n))
          if(dd.ne.d) then
            dact(d,dd) = 1
          endif
        end do ! n
        nrv = nrv + revp(2,d)
      end do ! d

      if(debug) then
        call iprint(dact,domains,domains,domains,'DACT')
      endif

      end

      subroutine node_graph_out(domains,numnp,numedg,point,nodeg)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output node graph to disk for Windows version

c      Inputs:
c        domains   - Number of domains
c        numnp     - Number of nodal points
c        numedg    - Number of graph edges
c        point(*)  - Pointer array for edges
c        nodeg(*)  - Values for nodal edges

c      Outputs:
c        To disk file 'metisgr'
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

#     include   "debugs.h"
#     include   "iodata.h"
#     include   "iofile.h"

      integer    domains,numnp,numedg, n,i
      integer    point(*),nodeg(*)

      open(unit = ios,file = 'metisgr',status = 'unknown')
      rewind ios

      if(debug) then
        write(  *,*) ' No. Nodes =',numnp,' No. edges =',numedg
        write(iow,*) ' No. Nodes =',numnp,' No. edges =',numedg
      endif

      write(ios,2000) numnp,numedg
      do n = 1,numnp
        write(ios,2000) (nodeg(i),i=point(n),point(n+1)-1)
      end do ! n

      close (ios)

c     Format

 2000 format(36i8)

      end

      subroutine u_init_out(u)

      implicit   none

#     include   "cdata.h"
#     include   "comfil.h"
#     include   "iodata.h"
#     include   "sdata.h"

      character  outfile*128
      integer    n
      real*8     u(ndf*numnp)

      outfile      = finp
      outfile(1:1) = 'U'
      open(unit = ios,file = outfile)
      rewind ios
      write(ios,2000) (u(n),n=1,ndf*numnp)
      close(ios,status = 'keep')

c     Format

2000  format(1p,3d25.15)

      end

      subroutine u_init_inp(u)

      implicit   none

#     include   "cdata.h"
#     include   "comfil.h"
#     include   "iodata.h"
#     include   "sdata.h"

      logical    exst
      character  outfile*128
      integer    n
      real*8     u(ndf*numnp)

      outfile      = finp
      outfile(1:1) = 'U'
      inquire(file = outfile, exist=exst)
      if(exst) then
        open(unit = ios,file = outfile)
        rewind ios
        read (ios,2000) (u(n),n=1,ndf*numnp)
        close(ios,status = 'keep')
      else
        write(  *,2001) 'File:',outfile,' does not exist'
      endif

c     Format

2000  format(1p,3d25.15)
2001  format(3a)

      end

      subroutine getgfile(ios,nodpt,numnp)

c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Read array of integers from unit 'ios'

c      Inputs:
c        ios      - Unit number for read
c        numnp    - Length of record

c      Output:
c        nodpt(*) - Record of input data
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    ios,numnp
      integer    nodpt(numnp)

      read(ios,*) nodpt

      end

      subroutine uinitld(revx, u, nf,nrv, ndf, ios)

c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Output initial conditions

c      Inputs:
c        revx(2,*) - Reverse nodal numbering
c        nf        - Number domain nodes
c        nrv       - Offset for current domain
c        ndf       - DOF's/node
c        ios       - Unit number for read

c      Output:
c        u(ndf,*)  - Initial conditions
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      logical    pflag
      integer    n,nf,ng,nrv,ndf,e,ios, revx(2,*)
      real*8     u(ndf,*)

      do n = 1,nf
        ng    = revx(1,nrv+n)
        pflag = .false.
        do e = 1,ndf
          if(u(e,ng).ne.0.0d0) then
            pflag = .true.
            exit
          endif
        end do ! e
        if(pflag) then
          write(ios,2001)  revx(2,nrv+n),0,(u(e,ng),e=1,ndf)
        endif
      end do ! n

c     Formats

2001  format( i9,i2,1p,10e16.8)

      end
