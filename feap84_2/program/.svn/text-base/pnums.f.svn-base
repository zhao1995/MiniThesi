c$Id:$
      subroutine pnums()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add blend option for element types by text       21/01/2007
c       2. Divide computed number of 11 node tets by 8      08/02/2007
c       3. Add 'sweep' option for 2-d to 3-d meshes         17/03/2007
c       4. Set length on assignment to cc and c2            17/04/2007
c       5. Set lengths for 14/15 node tets                  06/09/2007
c       6. Modify user block inputs                         02/09/2008
c       7. Add cubic block type (19); remove format 3004    06/02/2009
c       8. Change 'unumnd' to 'unumel' for user mesh sets   21/11/2010
c       9. Add 'mate' to set of material number             30/11/2010
c      10. Correct 'blen' coding to determine elements      18/07/2010
c      11. Add 'tx' to argument on call to umshlib          26/09/2011
c      12. Change tx to tx(1)                               09/01/2012
c      13. Dimension di(1) and du(1)                        01/05/2012
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Determine maximum node, element and material set
c               number in mesh.

c      Inputs:
c         From read of data file(s)

c      Outputs:
c         numnp     - Maximum node         number input of generated
c         numel     - Maximum element      number input of generated
c         nummat    - Maximum material set number input of generated
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cblend.h'
      include  'cdata.h'
      include  'chdata.h'
      include  'codat.h'
      include  'comfil.h'
      include  'dstars.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'iosave.h'
      include  'nblend.h'
      include  'pload1.h'
      include  'sdata.h'
      include  'umac1.h'

      logical   pcomp, errck, pinput,tinput,vinput, lmesh
      logical   flag,btyp,readfl,savefl,bernfl,eltype,pblktyp,eflg
      character cc*4,c2*8, tx(8)*15
      integer   numn,nume,numm,num,nt,nn,n0,n1,n2,n3,n4,n5,n6,n7,n8,mab
      integer   numn2,nume2,numm2,nseg, iorold
      integer   di(1),dnumb,dnurb,dlayer,nlay,j,loop_l,loop_t,loop_d(20)
      real*8    tb(7),td(16),du(1), loop_v

      save

c     Use input data to determine total number of nodes or elmts

      ldtot     = 0
      loop_t    = 1
      loop_l    = 1
      loop_d(1) = 1
      sptot     = 0
      numn      = 0
      nume      = 0
      numm      = 0

      numsn     = 0
      numsd     = 0
      numfc     = 0
      numbd     = 0
      dnumb     = 0
      mxilr     = 0

      nursn     = 0
      nursd     = 0
      nurfc     = 0
      nurbd     = 0
      dnurb     = 0
      mxnlr     = 0

      lmesh = .false.

  1   errck = tinput(tx,2,td,0)
      if(errck) go to 1
      cc = tx(1)(1:4)
      c2 = tx(2)(1:4)

c     [read] - set read flag for reinputs

      if( pcomp(cc,'read',4) ) then
        lmesh = readfl(tx(2))
        go to 1
      endif

c     [save] - set save flag for saving data

      if(pcomp(cc,'save',4)) then
        lsave = savefl(tx(2))
        go to 1
      endif

c     [coor] - Nodal coordinate inputs

      if(pcomp(cc,'coor',4)) then

  100   errck = pinput(td,1)
        if(errck) go to 100
        num  = nint(td(1))
        numn = max(numn,(num+starnd)*loop_t)
        if(num.gt.0) go to 100

c     [snod] - Super node inputs

      elseif(pcomp(cc,'snod',4)) then

  105   errck = pinput(td,1)
        num   = nint(td(1))
        numsn = max(numsn,num*loop_t)
        if(num.gt.0) go to 105

c     [elem] - Element data inputs

      elseif(pcomp(cc,'elem',4)) then

  110   errck = pinput(td,1)
        if(errck) go to 110
        num  = nint(td(1))
        nume = max(nume,(num+starel)*loop_t)

c       Input additional records

        if(num.gt.0) then
          do n1 = 14,nen,16
            n2 = min(nen+1-n1,16)
  111       errck = pinput(td,n2)
            if(errck) go to 111
          end do ! n1
          go to 110
        endif

c     [side] - Side data inputs

      elseif(pcomp(cc,'side',4)) then

        tx(1) = 'initial'
        do while (.not.pcomp(tx(1),'    ',4))
          errck = tinput(tx,1,td,0)
          if(pcomp(tx(1),'pola',4) .or. pcomp(tx(1),'cart',4) .or.
     &       pcomp(tx(1),'segm',4) .or. pcomp(tx(1),'elip',4) .or.
     &       pcomp(tx(1),'bern',4)) then
            numsd = numsd + loop_t
          endif
        end do ! while

c     [face] - Face data inputs

      elseif(pcomp(cc,'face',4)) then

  114   errck = pinput(td,1)
        if(errck) go to 114
        num = nint(td(1))
        numfc = max(numfc,num)
        if(num.gt.0) go to 114

c     [trib] - Triangular block inputs

      elseif(pcomp(cc,'trib',4)) then

  115   errck = tinput(tx,1,tb,5)
        if(errck) go to 115
        n2  = nint(tb(1))
        if(n2.le.0) then
          write(ilg,3002) n2
          write(iow,3002) n2
          call plstop()
        endif
        n3  = nint(tb(2))
        n4  = nint(tb(3))
        n5  = nint(tb(4))
        n6  = nint(tb(5))
        if(abs(n6).le.1) then
          n7 = n2
        elseif(n6 .eq.2) then
          n7 = mod(n2,2)
          if(n7.ne.0) then
            write(ilg,3007) n2
            write(iow,3007) n2
          else
            n7 = n2/2
          endif
        else
          write(ilg,3008) n6
          write(iow,3008) n6
          call plstop()
        endif

        n1 = ((n2+1)*(n2+2))/2*loop_t       ! Number of block nodes
        if(n3.gt.0) then
          numn = max(numn,n3+n1-1)
        else
          numn = numn+n1
        endif
        if(n4.gt.0) then
          nume = max(nume,n4+n7*n7*loop_t-1)
        elseif(n4.eq.0) then
          nume = nume+n7*n7*loop_t
        endif

c     [bloc] - block inputs

      elseif(pcomp(cc,'bloc',4)) then

  120   errck = tinput(tx,1,tb,7)
        if(errck) go to 120
        n2  = nint(tb(1))
        if(n2.le.0) then
          write(ilg,3002) n2
          write(iow,3002) n2
          call plstop()
        endif
        n3  = nint(tb(2))
        n4  = nint(tb(3))
        n5  = nint(tb(4))
        n6  = nint(tb(5))
        n7  = nint(tb(6))
        n8  = nint(tb(7))

        bernfl = .false.
        nn  = -1
  121   errck = tinput(tx,1,td,7)
        if(errck) go to 121
        if(pcomp(tx(1),'laye',4)) then
          dlayer = nint(tb(1))
          if(dlayer.eq.1) then
            nlay = n2
          elseif(dlayer.eq.2) then
            nlay = n3
          elseif(dlayer.eq.3) then
            nlay = n4
          endif

          j = 1
          do while(j.le.nlay)
  123       errck = tinput(tx,0,td,1)
            if(errck) go to 123
            j = j + 16
          end do ! while
          go to 121
c       elseif(pcomp(tx(1),'bern',4)) then
          bernfl = .true.
          nn     = nn + 1
          go to 121
c       No layers set a coordinate number
        else
          mab    = 0
          eltype = .false.
          do while (.not.eltype)
            eltype = pblktyp(tx(1),td, n8,n3,mab)
            if(.not.eltype) then
              errck = tinput(tx,1,td,2)
            endif
          end do !
          if(eltype) then
            nn = nn + 1
          else
            go to 121
          endif
        end if ! layers

        n7 = max(0,n7)
  124   errck = pinput(td,1)
        nn = nn + 1
        if(nint(td(1)).ne.0) go to 124
        if(n8.eq.0 .and. n4.gt.0 .and. ndm.eq.3 .and. nn.ge.7) then
          n8 = 10
        endif

        if(n8.lt.10 .or. abs(n8).eq.16) then
          nt = 1
          if(bernfl .or. nn.le.3) then       ! Line elements
            n1 = n2+1
            if(n3.le.1) then
              n0 =  n2
            elseif(n3.le.9) then
              n0 =  n2/n3
            else
              n0 =  n2/2
            endif
          else                   ! Block of elements
            if(n3.le.0) then
              write(ilg,3002) n2,n3
              write(iow,3002) n2,n3
              call plstop()
            endif
            n0 =  n2*n3
            n1 = (n2+1)*(n3+1) + n7*n3
            if(abs(n8).eq.16) then
              n0 = n0/9
            elseif(abs(n8).eq.7) then
              n0 = n0/2
              if(n8.eq.-7) n1 = n1 + n0
            elseif(n8.ge.8) then
              n0 = n0/4
              if(n8.eq.8) n1 = n1 - n0
            elseif(n8.gt.0) then
              n0 = n0*2
            elseif(n8.eq.-1) then
              n0 = n0*4
              n1 = n1 + n0
            endif
          endif
          if(n4.gt.0) then
            numn = max(numn,n4+n1*loop_t-1)
          else
            numn = numn+n1*loop_t
          endif
          if(n5.gt.0) then
            nume = max(nume,n5+n0*loop_t-1)
          elseif(n5.eq.0) then
            nume = nume+n0*loop_t
          endif
        elseif(n8.lt.20) then
          nt = n4
          if(n4.le.0) then
            write(ilg,3002) n2,n3,n4
            write(iow,3002) n2,n3,n4
            call plstop()
          endif
          if(n5.gt.0) then
            if(n8.eq.15) then                     ! 11-node tetrahedron
              numn = max(numn,n5 + ((n2+1)*(n3+1)*(n4+1)
     &                           + (6*n2*n3*n4)/8)*loop_t-1)
            elseif(n8.eq.17) then                 ! 14-node tetrahedron
              numn = max(numn,n5 + ((n2+1)*(n3+1)*(n4+1)
     &                           + (n2*n3 + n3*n4 + n4*n1)/2
     &                           + (6*n2*n3*n4)/4)*loop_t-1)
            elseif(n8.eq.18) then                 ! 15-node tetrahedron
              numn = max(numn,n5 + ((n2+1)*(n3+1)*(n4+1)
     &                           + (n2*n3 + n3*n4 + n4*n1)/2
     &                           + (9*n2*n3*n4)/4)*loop_t-1)
            else                                  ! Brick elements
              numn = max(numn,n5 + (n2+1)*(n3+1)*(n4+1)*loop_t-1)
            endif
          else
            if(n8.eq.15) then                     ! 11-node tetrahedron
              numn = numn + ((6*n2*n3*n4)/8)*loop_t
            elseif(n8.eq.17) then                 ! 14-node tetrahedron
              numn = numn + ((6*n2*n3*n4)/4
     &                     + (n2*n3+n3*n4+n4*n2)/2)*loop_t
            elseif(n8.eq.18) then                 ! 15-node tetrahedron
              numn = numn + ((9*n2*n3*n4)/4
     &                     + (n2*n3+n3*n4+n4*n2)/2)*loop_t
            else                                  ! Brick elements
              numn = numn + (n2+1)*(n3+1)*(n4+1)*loop_t
            endif
          endif
          if(n6.gt.0) then
            if(n8.eq.10) then                   !  8-node hexahedron
              nume = max(nume,n6 + n2*n3*n4*loop_t-1)
            elseif(n8.eq.11) then               !  4-node tetrahedron
              nume = max(nume,n6 + 6*n2*n3*n4*loop_t-1)
            elseif(n8.eq.12 .or. n8.eq.14) then ! 20/27-node hexahedron
              nume = max(nume,n6 + n2*n3*n4/8*loop_t-1)
            elseif(n8.eq.13 .or. n8.eq.15) then ! 10/11-node tetrahedron
              nume = max(nume,n6 + 6*(n2*n3*n4)/8*loop_t-1)
            elseif(n8.eq.17 .or. n8.eq.15) then ! 14/15-node tetrahedron
              nume = max(nume,n6 + 6*(n2*n3*n4)/8*loop_t-1)
            elseif(n8.eq.19) then               ! 64-node hexahedron
              nume = max(nume,n6 + n2*n3*n4/27*loop_t-1)
            endif
          elseif(n6.eq.0) then
            if(n8.eq.10) then                   !  8-node hexahedron
              nume = nume + n2*n3*n4*loop_t
            elseif(n8.eq.11) then               !  4-node tetrahedron
              nume = nume + n2*6*n3*n4*loop_t
            elseif(n8.eq.12 .or. n8.eq.14) then ! 20/27-node hexahedron
              nume = nume + n2*n3*n4/8*loop_t
            elseif(n8.eq.13 .or. n8.eq.15) then ! 10/11-node tetrahedron
              nume = nume + 6*(n2*n3*n4)/8*loop_t
            elseif(n8.eq.17 .or. n8.eq.18) then ! 14/15-node tetrahedron
              nume = nume + 6*(n2*n3*n4)/8*loop_t
            elseif(n8.eq.19) then               ! 64-node hexahedron
              nume = nume + n2*n3*n4/27*loop_t
            endif
          endif
        elseif(n8.lt.30) then
          nt = 1
          if(n3.le.0) then
            write(ilg,3002) n2,n3
            write(iow,3002) n2,n3
            call plstop()
          endif
          n0 =  n2*n3
          n1 = (n2+1)*(n3+1) + n7*n3
          if(abs(n8).eq.27) then
            n0 = n0/2
            if(n8.eq.-7) n1 = n1 + n0
          elseif(n8.ge.28) then
            n0 = n0/4
            if(n8.eq.8) n1 = n1 - n0
          elseif(n8.gt.20) then
            n0 = n0*2
          endif
          if(n4.gt.0) then
            numn = max(numn,n4+n1*loop_t-1)
          else
            numn = numn+n1*loop_t
          endif
          if(n5.gt.0) then
            nume = max(nume,n5+n0*loop_t-1)
          elseif(n5.eq.0) then
            nume = nume+n0*loop_t
          endif
        else
          if(ndm.le.2) then
            nt = 1
            n4 = 1
            n0 = nint(tb(3))
            n1 = nint(tb(4))
            n5 = nint(tb(5))
          else
            nt = n4
            n4 = nint(tb(3))
            n0 = nint(tb(4))
            n1 = nint(tb(5))
            n5 = nint(tb(6))
          endif
          n6 = n0
          n7 = n1
          n8 = n8 - 30
          call ublk(n8,nn,n2,n3,n4,du,du,di,di,du(1),du(1),du(1),
     &              n1,n0,ndm,nen1,n5,tx(1),.false., 1)
          if(n6.eq.0) then
            nume = nume + n0*loop_t
          else
            nume = max(nume,n6+(n0-n6)*loop_t)
          endif
          if(n7.eq.0) then
            numn = numn + n1*loop_t
          else
            numn = max(numn,n7+(n1-n7)*loop_t)
          endif
        endif

c     [blen] - blending function inputs

      elseif(pcomp(cc,'blen',4)) then

        dnumb = dnumb + 1
        errck = tinput(tx(2),1,td,7)

        n0 = nint(td(1))
        n1 = max(1,nint(td(2)))

        if(pcomp(tx(2),'surf',4) .or. pcomp(tx(2),'cap',3)) then
          n2   = 1
          n3   = nint(td(3))
          n4   = nint(td(4))
          n5   = nint(td(5))
          n6   = nint(td(6))
          btyp = .true.
        elseif(pcomp(tx(2),'soli',4)) then
          n2   = max(1,nint(td(3)))
          n3   = nint(td(4))
          n4   = nint(td(5))
          n5   = nint(td(6))
          n6   = max(10,nint(td(7)))
          btyp = .false.
        elseif(pcomp(tx(2),'line',4)) then
          n2   = 1
          n3   = nint(td(3))
          n4   = nint(td(4))
          n5   = nint(td(5))
          n6   = nint(td(6)) + 30
          btyp = .true.
        else
          write(ilg,3003) tx(2)
          write(iow,3003) tx(2)
          call plstop()
        endif
        if(n5.lt.0) then
          mxilr = max(mxilr,n0,n1,n2)
        endif

        errck  = tinput (tx(2),1,tb,7)
        di(1)  = n6
        mab    = 0
        eltype = .false.
        do while (.not.eltype)
          eltype = pblktyp(tx(2),tb, di(1),n0,mab)
          if(di(1).ne.n6) then
            n6   = di(1)
          endif
          if(.not.eltype) then
            errck = tinput(tx(2),1,tb,7)
          endif
        end do ! while

        if(btyp) then
          flag = nint(tb(3)).eq.0
        else
          flag = .false.
        endif
        call pnumbl(ndm,n0,n1,n2,n6, n7,n8, flag)

c       Set maximum node and element numbers

        if(n3.gt.0) then
          numn = max(numn,n3+n8*loop_t-1)
        else
          numn = numn+n8*loop_t
        endif

        if(n4.gt.0) then
          nume = max(nume,n5+n7*loop_t-1)
        elseif(n4.eq.0) then
          nume = nume+n7*loop_t
        endif

c     [mate] - material data set inputs

      elseif(pcomp(cc,'mate',4)) then

        errck = vinput(yyy(16:30),15,td,1)
        n0    = nint(td(1))
        if(n0.eq.0 .and. numm.ne.0) then
          write(ilg,3006)
          write(  *,3006)
          call plstop()
        else
          numm = max(numm,n0,1)
        endif

c     [para] - parameter inputs

      elseif(pcomp(cc,'cons',4).or.pcomp(cc,'para',4)) then

        coflg = .true.
        call pconst(.false.)

c     [swee]p - 2-d to 3-d mesh

      elseif(pcomp(cc,'swee',4)) then

        errck = tinput(tx(1),1,td,0)
        inquire(file=tx(1),exist=eflg)

        if(eflg) then

          open(unit=99,file=tx(1))

c         Get size of 2-d mesh describing sweeep

          iorold = ior
          ior    = 99
          tx(1)   = 'start'
          do while( .not.pcomp(tx(1),'end',3))

            errck = tinput(tx(1),1,td,6)

            if(pcomp(tx(1),'feap',4)) then
              errck = tinput(tx(1),0,td,3)
              numn2 = nint(td(1))
              nume2 = nint(td(2))
              numm2 = nint(td(3))
            endif
          end do ! while
          close(99,status = 'keep')
          ior = iorold

c         Read number of segments & compute number of nodes & elements

          errck = tinput(tx(1),1,td,1)
          nseg  = max(1,nint(td(1)))
          numn  = numn + (numn2*nseg + numn2)*loop_t
          nume  = nume + (nume2*nseg)*loop_t
        else
          write(iow,*) 'FILE: ',tx(1),' does not exist.'
          call plstop()
        endif

c     [spin] - Spin inputs data

      elseif(pcomp(cc,'spin',4)) then

        sptot = sptot + 1

c     [nurb] - NURB coordinate inputs

      elseif(pcomp(cc,'nurb',4)) then

  125   errck = pinput(td,1)
        num   = nint(td(1))
        nursn = max(nursn,num*loop_t)
        if(num.gt.0) go to 125

c     [nsid] - NURB side data inputs

      elseif(pcomp(cc,'nsid',4)) then

        tx(1) = 'initial'
        do while (.not.pcomp(tx(1),'    ',4))
          errck = tinput(tx,1,td,0)
          if(pcomp(tx(1),'pola',4) .or. pcomp(tx(1),'cart',4) .or.
     &       pcomp(tx(1),'segm',4) .or. pcomp(tx(1),'elip',4) .or.
     &       pcomp(tx(1),'bern',4) .or. pcomp(tx(1),'side',4)) then
            nursd = nursd + loop_t
          endif
        end do ! while

c     [nble] - NURB blending function inputs

      elseif(pcomp(cc,'nble',4)) then

        dnurb = dnurb + 1
        errck = tinput(tx(2),1,td,7)

        n0 = nint(td(1))
        n1 = max(1,nint(td(2)))

        if(pcomp(tx(2),'surf',4)) then
          n2   = 1
          n3   = nint(td(3))
          n4   = nint(td(4))
          n5   = nint(td(5))
          n6   = nint(td(6))
          btyp = .true.
        elseif(pcomp(tx(2),'soli',4)) then
          n2   = max(1,nint(td(3)))
          n3   = nint(td(4))
          n4   = nint(td(5))
          n5   = nint(td(6))
          n6   = max(10,nint(td(7)))
          btyp = .false.
        elseif(pcomp(tx(2),'line',4)) then
          n2   = 1
          n3   = nint(td(3))
          n4   = nint(td(4))
          n5   = nint(td(5))
          n6   = nint(td(6)) + 30
          btyp = .true.
        else
          write(ilg,3003) tx(2)
          write(iow,3003) tx(2)
          call plstop()
        endif
        if(n5.lt.0) then
          mxnlr = max(mxnlr,n0,n1,n2)
        endif

        errck  = tinput (tx(2),1,tb,7)
        eltype = pblktyp(tx(2),tb, di(1),n0,mab)
        if(eltype) then
          flag = nint(tb(3)).eq.0
        else
          n6   = di(1)
          if(btyp) then
            errck= pinput(td,4)
            flag = nint(td(4)).eq.0
          else
            flag = .false.
          endif
        endif
        call pnumbl(ndm,n0,n1,n2,n6, n7,n8, flag)

c       Set maximum node and element numbers

        if(n3.gt.0) then
          numn = max(numn,n3+n8*loop_t-1)
        else
          numn = numn+n8*loop_t
        endif

        if(n4.gt.0) then
          nume = max(nume,n5+n7*loop_t-1)
        elseif(n4.eq.0) then
          nume = nume+n7*loop_t
        endif

c     [*nod] - star node number set

      elseif(pcomp(cc,'*nod',4)) then

        errck  = vinput(yyy(16:30),15,td,1)
        starnd = nint(td(1))

c     [*ele] - star element number set

      elseif(pcomp(cc,'*ele',4)) then

        errck  = vinput(yyy(16:30),15,td,1)
        starel = nint(td(1))

c     [loop] - start looping

      elseif(pcomp(cc,'loop',4)) then

c       Initialize any data that needs to be counted during loops

        if(loop_l.eq.1) then
          dnumb = 0
        endif

c       Set loop data

        call setval(tx(2),15,loop_v)
        loop_l = loop_l + 1
        loop_d(loop_l) = max(1,nint(loop_v))
        loop_t         = loop_t*loop_d(loop_l)

c     [next] - end looping

      elseif(pcomp(cc,'next',4)) then

c       Multiply anything that is counted through outer loop

        if(loop_l.eq.2) then
          dnumb = dnumb*loop_t
          numbd = numbd + dnumb
          dnumb = 0

          dnurb = dnurb*loop_t
          nurbd = nurbd + dnurb
          dnurb = 0
        endif

c       Close loop set counters for next loop

        loop_t = loop_t/loop_d(loop_l)
        loop_l = loop_l - 1

c     [load] - load table

      elseif(pcomp(cc,'load',4)) then

        if(pcomp(c2,'prop',4) .or. pcomp(c2,'    ',4)) then
          ldtot = ldtot + 1
        endif

c     [end] - end of mesh input data

      elseif(pcomp(cc,'end',3)) then

        if(lsave) then
          write(ilg,3000)
          write(iow,3000)
          call plstop()
        elseif(loop_l.ne.1) then
          write(ilg,3005)
          write(iow,3005)
          call plstop()
        endif

c       Set final number of blends

        numbd = numbd + dnumb
        nurbd = nurbd + dnurb

        go to 200

c     User functions for node and element counting

      else

        do j = 1,10
          if(pcomp(cc,umshc(j),4)) then
            uct    = umshc(j)
            unumnd = 0
            unumel = 0
            call umshlib(j,tx,.false.)
            numn = max(numn,(unumnd+starnd)*loop_t)
            nume = max(nume,(unumel+starel)*loop_t)
          endif
        end do ! j

      endif
      go to 1

 200  numnp  = max(numnp ,numn)
      numel  = max(numel ,nume)
      nummat = max(nummat,numm)

c     Check that nodes, elements and material sets exist

      if(numnp.le.0 .or. numel.le.0 .or. nummat.le.0) then
        write(ilg,3001) numnp,numel,nummat
        write(iow,3001) numnp,numel,nummat
        call plstop()
      endif

c     Find start of problem for input of mesh data

      rewind ior
 300  read(ior,1000) cc
      if(.not.pcomp(cc,'feap',4)) go to 300

c     Advance past control record and return to 'PCONTR'

      read(ior,1000) cc

c     Input/Output Formats

1000  format(a4,11x,i15)

3000  format(' *ERROR* PNUMS: No closing SAVE,END statement for input.')

3001  format(' *ERROR* PNUMS: Problem size does not permit a solution:'/
     &       '         Number of nodes         = ',i9/
     &       '         Number of elements      = ',i9/
     &       '         Number of material sets = ',i9/)

3002  format(' *ERROR* PNUMS: Number of block increments incorrect:'/
     &       '         Number of 1-increments = ',i6/:
     &       '         Number of 2-increments = ',i6/:
     &       '         Number of 3-increments = ',i6/)

3003  format(' *ERROR* PNUMS: Blending option: ',a,' not supported')

3005  format(' *ERROR* PNUMS: Incorrect LOOP-NEXT structure.')

3006  format(' *ERROR* PNUMS: Zero Material Number Specified')

3007  format(' *ERROR* PNUMS: Number of increments =',i5/
     &       '                Must be an even number')

3008  format(' *ERROR* PNUMS: Wrong block type: Input = ',i5)

      end
