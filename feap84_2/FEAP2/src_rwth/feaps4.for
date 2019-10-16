      subroutine plshowm(arcf,rlnew)
c-----------------------------------------------------------------------
c
c      Purpose: show macro informations                   
c
c      Inputs:
c         arcf     - flag arc-length method
c         rlnew    - load factor
c         futher parameters via common
c
c      Outputs:
c
c      Open:                                    
c      arcl: steplength,direction                                 
c
c-----------------------------------------------------------------------
      USE ext2
      USE fdata
      USE iofile
      USE pdata2
      USE plodf
      USE plodfa
      USE plong
      USE prlod
      USE psize
      USE rdata
      USE slid3
      USE soltyp
      USE subdt
      USE tdata
      implicit double precision (a-h,o-z)
      character y(2)*3
      character*55 sol(10)
      logical arcf
      data y /'yes','no '/
      data sol /'standard profil solver',
     +          'SM sparse matrix solver (CSR, without minimum degree)',
     +          'SM sparse matrix solver (CSR, with minimum degree)   ',
     +          'SuperLU  solver with CSR storage',
     +          'Pardiso  solver with CSR storage',
     +          'PCG      solver with CSR storage',
     +          'PGMRES   solver with CSR storage',
     +          'PPGMRES2 solver with CSR storage',
     +          'Pardiso  solver with CSR storage',
     +          'simplex optimization algorithm'/
      xmb   = 4.0/1024
      used  = kmax*xmb
      avail = maxm*xmb
      perc  = 100.0*used/avail
      if(ior.lt.0) then
cww        if(idev.eq.4) call clwopen('Show Values MACRO',20,10,600,460,1,2)
        if(idev.eq.4) call clwopen('Show Values MACRO',1,2)
        if(idev.eq.1.or.idev.eq.2.or.idev.eq.3) write(*,1000)
        write(*,2001) ttim
        write(*,2002) dt
        if(     arcf) write(*,2003) prop*rlnew
        if(.not.arcf) write(*,2003) prop
        write(*,2004) mfmax
        if(     arcf  ) write(*,2005) y(1)
        if(.not.arcf  ) write(*,2005) y(2)
        if(     extflg) write(*,2006) y(1)
        if(.not.extflg) write(*,2006) y(2)
        if(     fl(9) ) write(*,2007) y(1)
        if(.not.fl(9) ) write(*,2007) y(2)
        if(     contfl) write(*,2008) y(1)
        if(.not.contfl) write(*,2008) y(2)
        if(     fl(1) ) write(*,2009) y(1)
        if(.not.fl(1) ) write(*,2009) y(2)
        if(     fl(2) ) write(*,2010) y(1)
        if(.not.fl(2) ) write(*,2010) y(2)
        write(*,2011) tol
        if(     pfr   ) write(*,2012) 'prin'
        if(.not.pfr   ) write(*,2012) 'nopr'
        write(*,2013) npldf
        write(*,2014) ipl(2,1),ipl(2,2)
        write(*,2021) (noden(ii),ii=1,10)
        write(*,2015) ipola
        write(*,2016)
        write(*,2017) a(3,1),a(4,1),a(5,1),a(6,1),iexp(1)
        write(*,2018) a(1,1),a(2,1)
        write(*,2019)
        write(*,2020) avail,used,perc
        write(*,2022) sol(istyp+1)
        if(idev.eq.4) call clwclose(1,2,1)
      end if
      return
c.... formats
1000  format(/3x,'M A C R O   P A R A M E T E R S'/)
2001  format(5x,'Actual time........................... ',f12.5)
2002  format(5x,'Actual time increment................. ',f12.5)
2003  format(5x,'Actual load parameter................. ',e12.5)
2004  format(5x,'No. of calculated Eigenvectors (pairs) ',i5)
2005  format(5x,'Arclength calculation ................ ',a3)
2006  format(5x,'Extended System calc. ................ ',a3)
2007  format(5x,'Dynamic   calculation ................ ',a3)
2008  format(5x,'Contact   calculation ................ ',a3)
2009  format(5x,'Computation of consist. mass matrix... ',a3)
2010  format(5x,'Computation of lumped   mass matrix... ',a3)
2011  format(5x,'Tolerance for calculation............. ',e14.5)
2012  format(5x,'Print option set to................... ',a4)
2013  format(5x,'Tplo (load-disp. curve)(0=off 1=on)... ',i4)
2014  format(5x,'Tplo acts on.... node = ',i4,'.... dof = ',i4)
2021  format(5x,'avail. nodes:',10(1x,i4) )
2015  format(5x,'Displ. for polar coord. (cart = 0)...  ',i4)
2016  format(5x,'Proportional load table 1:')
2017  format(5x,'prop = ',f5.2,'+',f5.2,'*t+',f5.2,'*[sin ',f5.2,'*t]**'
     1          ,i3)
2018  format(5x,'tmin  = ',e10.5,'         tmax  = ',e10.5)
2019  format(5x,'Total Memory in blank common in KB:')
2020  format(5x,'avail: ',f9.1,' used: ',f9.1,' percent used: ',f5.2)
2022  format(5x,'using : ',1x,a55)
      end
c
      subroutine qloada(u,f,dr,id,numel,nneq)
c-----------------------------------------------------------------
c
c     not active (only for constant loads)  
c     active is SR PLOADQ via PLOADS
c
c      Purpose: calculate element load vectors based on values defined 
c               in SR qloadi (MESH>QLOA) and add to load vector f at the
c               beginning of process. Thus only linear load terms are 
c               taken into account. 
c               Within SR PLOADQ linear and nonlinear terms for 
c               load vectors are calculated 
c               -e.g. Temperature-or follower loads
c
c      Inputs:
c         u          - Nodal solution vector
c         id(*)      - Equation numbers for each active dof
c         numel      - Number of elements in mesh
c         nneq       - No of unknowns in problem 
c
c      Scratch:
c         dr(*)       
c
c      Outputs:
c
c         f(*)       - Nodal load vector for mesh
c
c
c.... ww bs uka 02/06
c-----------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical fa,tr
      dimension u(*),f(*),dr(*)
      fa=.false.
      tr=.true.
c.... calculate load terms in dr       
      call formfe(u,dr,dr,dr,fa,tr,fa,fa,22,1,numel,1)
c.... add in vector f taking into account b.c. 
      call pmoveca(id,dr,f,nneq)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine qloadi(q,ix,nen1,ie,nie,numel,prt)
c-----------------------------------------------------------------------
c
c      Purpose: set load values  for element load vectors which are 
c               calculated in SR qloada at the beginning of execution 
c               under isw=22
c
c      Inputs:
c         q(numel,10) - array of loads on each element
c         ix(nen1,*)  - Element nodal connections of mesh
c         ie(nie,*)   - Assembly information for material set
c         nen1        - Dimension for ix array
c         nie         - Dimension for ie array
c         numel       - Number of elements in mesh
c         prt         - Print Flag
c
c      Outputs:
c
c      Comments:
c
c      input is used within element formulation under isw=22           |
c      input values depend on element: q(1)-q(10)                      |
c      see manual on elements                                          |
c      elmt01: mat,-                                                   |
c      elmt02: mat,q1,q2,n1,n2                                         |
c      elmt03: mat,q1,q2,n1,n2                                         |
c      elmt04: mat,qxL,qyL,qzL,qxG,qyG,qzG                             |
c      elmt05: mat,b1,b2                                               |
c      elmt06: mat,-                                                   |
c      elmt07: mat,q                                                   |
c      elmt08: mat,to do                                               |
c      elmt09: mat,b1,b2,b3,bx,by,bz                                   |
c      elmt10: mat,-                                                   |
c      elmt11: mat,qx,qy,qz                                            |
c      elmt12: mat,-                                                   |
c      elmt13: mat,q                                                   |
c      elmt14: mat,-                                                   |
c      elmt15: mat,q1,q2,ifol                                          |
c      elmt16: mat,to do                                               |
c      elmt17: mat,q                                                   |
c      elmt18: mat,q                                                   |
c      elmt19: mat,q                                                   |
c      elmt20: mat,q1,q2                                               |
c      elmt21: mat,qx,qy,qz,dT                                         |
c      elmt22: mat,qx,qy                                               |
c      elmt30: mat,qx,qy,qz,ltyp                                               |
c
c      Call to QLOA can be done several times, values are added!
c                                                                      |
c-----------------------------------------------------------------------
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      dimension td(11),q(numel,10),q1(10)
      integer ix(nen1,*),ie(nie,*)
c
      if(ior.lt.0) write(  *,2000) 
                   write(iow,2000) 
c.... read input of load data, end with q=0  
100   call pzero(td,11) 
      call pzero(q1,10) 
      if(ior.lt.0) write(*,3001)
      call dinput(td,11)
      if(errck) go to 100
      qii = dsqrt(dot(td,td,11))
      if(qii.le.0) return ! end for zero loading
c.... Mat.Number      
      ma = td(1)
c.... Associated Element-typ
      iel = ie(nie-1,ma)
      if(iel.eq.0) 
     +   stop '***  error: element no. 0 calculated using macro QLOA.'
c.... Associated Element-loads
      do i = 1,10
        q1(i) = td(i+1)
      end do
      if(ior.lt.0) then 
        write(  *,2001)
        write(  *,2002) ma,iel,(q1(i),i=1,10)      
      end if 
      write(iow,2001)
      write(iow,2002)   ma,iel,(q1(i),i=1,10)      
c
c.... read associated elements end with e=0
101   if(ior.lt.0) write(*,3002)
102   call dinput(td,3)
      if(errck) go to 101
      kk = dsqrt(dot(td,td,3))
      if(kk.eq.0) goto 100 ! all elements read
      ia=td(1)
      ib=td(2)
      incr=td(3)
      if(ib.eq.0) ib=ia 
      if(incr.eq.0) incr=1
      if(ia.gt.numel.or.ia.lt.1) then
        write(*,*) 'Wrong element number in Macro QLOA  ',ia  
        stop
      end if     
      if(ib.gt.numel.or.ib.lt.1) then
        write(*,*) 'Wrong element number in Macro QLOA  ',ib  
        stop
      end if     
c      
c.... add element loading         
      do i = ia,ib,incr
        if(ma.ne.ix(nen1,i)) then
          write(  *,2004) i,ma
          write(iow,2004) i,ma
          stop
        end if
        do k = 1,10
          q(i,k) = q(i,k) + q1(k)
        end do 
      end do 
      if(ior.lt.0) write(  *,2003) (ii,ii=ia,ib,incr)
                   write(iow,2003) (ii,ii=ia,ib,incr)

c.... new record 
      goto 102   
c.... formats 
2000  format(/,'  QLOA element loads',$)
2001  format(/,2x,'Mat.No',1x,'EL.typ',
     +      '    q01     ','    q02     ','    q03     ','    q04     ',
     +      '    q05     ','    q06     ','    q07     ','    q08     ',
     +      '    q09     ','    q10     ')
     +          
2002  format(1x,i5,1x,i5,1x,10f12.5)
2003  format('  Elements',12(1x,i5)/(10x,12(1x,i5)))
2004  format('  Element ',i5,' belongs not to Material  ',i5)   

3001  format('Input: q(1)-q(10) due to element >',$)
3002  format('Input: (elem1,elem2,...., >',$)
      end 
c
      subroutine showcoor(x,ndm,numnp,n1,n2,n3)
c----------------------------------------------------------------------
c
c      Purpose: print nodal coordinates
c
c      Inputs:
c         x(ndm,*)   - Nodal coordinates of mesh
c         ndm        - Spatial dimension of mesh
c         numnp      - Number of nodes in mesh
c         n1         - first node
c         n2         -  last node
c         n3         - increment
c
c      Outputs:
c
c      Comments:
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension x(ndm,*)
c
      if(n2.eq.0) n2 = n1
      if(n3.eq.0) n3 = 1
        n1 = max(1,min(n1,numnp))
        n2 = max(1,min(numnp,n2))
        if(n2-n1.ne.0) n3 = isign(n3,n2-n1)
      write(*,2000) (i,'coord',i=1,ndm)
      do i =n1,n2,n3
        write(*,2001) i,(x(j,i),j=1,ndm)
      enddo
2000  format('N o d a l   C o o r d i n a t e s',/,
     1'  node',6(i6,a6):/(6x,6(i6,a6)))
2001  format(i6,1p6e12.5:/(6x,1p6e12.5))
      end
c
      subroutine showboun(id,ndf,numnp,n1,n2,n3)
c----------------------------------------------------------------------
c
c      Purpose: print nodal boundaries
c
c      Inputs:
c         id(ndf,*)  - Equation numbers for each active dof
c         ndf        - Number dof/node
c         numnp      - Number of nodes in mesh
c         n1         - first node
c         n2         -  last node
c         n3         - increment
c
c      Outputs:
c
c-----------------------------------------------------------------------
      integer id(ndf,*),idl(20) 
c
      if(ndf.gt.20) call drawmess('idl too small in showboun',1,0)
      if(n2.eq.0) n2 = n1
      if(n3.eq.0) n3 = 1
        n1 = max(1,min(n1,numnp))
        n2 = max(1,min(numnp,n2))
        if(n2-n1.ne.0) n3 = isign(n3,n2-n1)
      write(*,2000) (i,i=1,ndf)
      do i =n1,n2,n3
        call pzeroi(idl,20)
        do j=1,ndf
          if(id(j,i).lt.0) idl(j)=1 
        enddo
      write(*,2001) i,(idl(j),j=1,ndf)
      enddo
2000  format('N o d a l   b o u n d a r y  c o n d i t i o n s',/,
     1       6x,'node',9(i2,'-b.c.')/(10x,9(i2,'-b.c.')))
2001  format(i10,9i7/(10x,9i7))
      end
c
      subroutine showforc(f,ndf,numnp,n1,n2,n3)
c----------------------------------------------------------------------
c
c      Purpose: print nodal forces/prescribed displacements
c
c      Inputs:
c         f(ndf,*)   - Equation numbers for each active dof
c         ndf        - Number dof/node
c         numnp      - Number of nodes in mesh
c         n1         - first node
c         n2         -  last node
c         n3         - increment
c
c      Outputs:
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension f(ndf,*)
c
      if(n2.eq.0) n2 = n1
      if(n3.eq.0) n3 = 1
        n1 = max(1,min(n1,numnp))
        n2 = max(1,min(numnp,n2))
        if(n2-n1.ne.0) n3 = isign(n3,n2-n1)
      write(*,2000) (i,'force',i=1,ndf)
      do i =n1,n2,n3
        write(*,2001) i,(f(j,i),j=1,ndf)
      enddo
2000  format('N o d a l   F o r c e s',/,
     1'  node',6(i6,a6):/(6x,6(i6,a6)))
2001  format(i6,1p6e12.5:/(6x,1p6e12.5))
      end
c
      subroutine showelem(ix,nen,nen1,numel,n1,n2,n3)
c----------------------------------------------------------------------
c
c      Purpose: print nodes at elements
c
c      Inputs:
c         ix(nen1,*) - Element nodal connections of mesh
c         nen1       - Dimension for ix array
c         nen        - Max. Number of nodes per element
c         numel      - Number of elements in mesh
c         n1         - first element
c         n2         -  last element
c         n3         - increment
c
c      Outputs:
c
c-----------------------------------------------------------------------
      integer ix(nen1,*)
c
      if(n2.eq.0) n2 = n1
      if(n3.eq.0) n3 = 1
        n1 = max(1,min(n1,numel))
        n2 = max(1,min(numel,n2))
        if(n2-n1.ne.0) n3 = isign(n3,n2-n1)
        write(*,2000) (k,k=1,nen)
      do i =n1,n2,n3
          write(*,2001) i, ix(nen1,i), (ix(k,i),k=1,nen)
      enddo
2000  format('E l e m e n t s'//4x,'elmt    matl',
     1   8(i3,' node')/(16x,8(i3,' node')))
2001  format(10i8/(16x,8i8))
      end
c
      subroutine shownode(x,id,f,ndm,ndf,numnp,n1)
c----------------------------------------------------------------------
c
c      Purpose: print nodal values: coordinates, forces, Boundary cond.
c
c      Inputs:
c         x(ndm,*)   - Nodal coordinates of mesh
c         ix(nen1,*) - Element nodal connections of mesh
c         f(ndf,*)   - Equation numbers for each active dof
c         ndm        - Spatial dimension of mesh
c         ndf        - Number dof/node
c         numnp      - Number of nodes in mesh
c         n1         - Node number to show 
c
c      Outputs:
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension x(ndm,*),f(ndf,*)
      integer id(ndf,*),idl(20)
c
      if(ndf.gt.20) call drawmess('idl too small in showboun',1,0)
        n1 = max(1,min(n1,numnp))
c...  coordinates
      write(*,2000) n1,(i,'coord',i=1,ndm)
      write(*,2001) (x(j,n1),j=1,ndm)
c.... forces
      write(*,2002) (i,'force',i=1,ndf)
      write(*,2001) (f(j,n1),j=1,ndf)
c...  boundary conditions 
      write(*,2004) (i,i=1,ndf)
      call pzeroi(idl,20)
      do j=1,ndf
        if(id(j,n1).lt.0) idl(j)=1 
      enddo
      write(*,2005) (idl(j),j=1,ndf)
c
2000  format('N o d a l  V a l u e s  for node ',i6,/,
     +       'Coordinates',/,6(i6,a6):/(6(i6,a6)))
2001  format(1p6e12.5:/(1p6e12.5))
2002  format('Forces',/,6(i6,a6):/(6(i6,a6)))
2004  format('Boundary  conditions',/,
     1      9(i2,'-b.c.')/(10x,9(i2,'-b.c.')))
2005  format(9i7/(9i7))
      end
c
      subroutine pmacio (jct,lct,ct,wd,nwd,nlp,nif,ll,ncmds,asc1,asc2)
c----------------------------------------------------------------------
c
c      Purpose: Macro instruction input subprogram
c
c      Inputs:
c         wd(*)     - List of command language names
c         nwd       - Number of commands in wd
c         nlp       - Location of 'loop' command
c         nif       - Location of 'if'   command
c         ncmds     - number of command cards
c
c      Scratch:
c         ct(3,*)   - Values of command parameters
c
c      Outputs:
c         jct(ll)   - Command numbers to execute
c         lct(ll)   - Command options
c         asc1      - only for macro exit,,n1 
c                     n1=asc1>0: restart file is in ASCII-mode 
c         asc2      - only for updh and shared memory 
c
c      Comments 
c
c       Macro  |  V1           |  V2                  |  V3  
c-----------------------------------------------------------------------
c       LOOP   |  It.max       |  Pos. 1. after LOOP  |
c       NEXT   |  actual It.No |  Pos. 1. in    LOOP  |
c              |  or Itmax     |                      |
c
c-----------------------------------------------------------------------
c       IF     |  -            |  Pos. next ELSE      |  Pos. ENDI
c       ELSE   |  -            |  Pos. next ELSE      |  Pos. ENDI
c       ENDI   |  -            |  -                   |  -
c
c-----------------------------------------------------------------------
c     special cases for exit
c     ll=-2 QUIT,q,RINP(irfeap=2), RINP,new(irfeap=1)
c     ll=-2 REME,new(irfeap=3), REME,old(irfeap=3) in PMACR
c     ll=-1 exit,e
c
c----------------------------------------------------------------------
c
      USE bdata
      USE fdata
      USE idata
      USE iofile
      USE rfeap
      USE yydata 
      implicit real*8 (a-h,o-z)
      logical pcomp,flg
      character*4 lct(*),wd(nwd),clab1,clab2,yy*1
      integer jct(*)
      real*8  ct(3,*)
cww   real    tary

      dimension yy(80)

      save nnn
      data nnn/0/
      asc1 = 0.0d0
      asc2 = 0.0d0
      

c.... initiate the read of macro statements
cww   if(ior.gt.0) write(iow,2001) o,head
      if(ior.gt.0) write(iow,2001)

c.... read macro cards
      ll = 1
      nn = nnn
      nnx = nlp + 1 ! pos next = loop+1 
      nel = nif + 1 ! pos else = if+1
      nei = nif + 2 ! pos endi = if+2
      jct(1) = nlp
      ct(1,1) = 1.0

      ct(2,1) = 1.0 ! from IF ??

      flg = .true.
cww100   if(ior.lt.0.and.flg.and.pfr) write(*,2002)

100   continue
      if(nn.gt.ncmds) nn = 0
      if(ll.gt.ncmds-5.and.ior.lt.0) write(*,3001) ncmds,ncmds
cww   if(ior.lt.0) write(*,2003) nn+1,ll
      if(ior.lt.0) write(*,2003) 
      ll = ll + 1
      if(ll.gt.ncmds) then
        write(yyy,3002) ncmds
        call drawmess(yyy,1,0)
cww       stop
        return
      end if

c.... input command
      call parse(clab1,clab2,ct(1,ll))

c...  copy yyy->yy (USE yydata) 
      do i = 1,80
        yy(i) = yyy(i:i)
      end do


      if(ior.lt.0) then
      
         if(pcomp(clab1,'help',4)) then
c           show macro commands
            jflag = 0
            do 111 jj = 1,nwd
              if(pcomp(clab2,wd(jj),4)) then
              call pman(2,wd(jj))
              jflag = 1
              end if
111         continue
            if(pcomp(clab2,'end ',4).or.pcomp(clab2,'exit',4).or.
     +         pcomp(clab2,'quit',4).or.pcomp(clab2,'hist',4).or.
     +         pcomp(clab2,'help',4).or.pcomp(clab2,'proc',4)) then
               call pman(2,clab2)
               jflag = 1
            end if
            if(jflag.eq.0) call phelpm(wd,nwd,clab2)
            ll = ll - 1
            go to 100

         else if(pcomp(clab1,'!   ',1)) then
c           repeat previous or match command
            call histex(wd,clab1,jct,lct,ct,nwd,nlp,nnx,ll,is)
            if(is.eq.0) go to 140
            ll = ll - 1
            go to 100
         end if
      end if

      if(pcomp(clab1,'hist',4)) then
c       Hist: Perform set of commands from history
        call phist(wd,clab2,jct,lct,ct,ll,is)
        nnn = nn
        if(is.eq.0) then
          go to 150
        else
          ll = ll - 1
          go to 100
        end if

      else if(pcomp(clab1,'proc',4)) then
c       Proc: Enter or execute procedure
        call proced(yy(16),yy(31),wd,nwd,ll,jct,lct,ct,flg,1)
        ll = ll - 1
        go to 100
      end if

c     End: end    =with rest in batch  rinp(new)
c          exit(e)=with rest           quit(q)=without rest     
      if(ior.gt.0.and.pcomp(clab1,'end ',4)) go to 150 
      if(ior.lt.0) then
        if(pcomp(clab1,'quit',4).or.pcomp(clab1,'q   ',4)) then
          go to 249
        else if(pcomp(clab1,'exit',4).or.pcomp(clab1,'e   ',4)) then
          go to 250
        else if(pcomp(clab1,'rinp',4)) then    ! restart program
          irfeap = 2 
          call deallocall
          if(pcomp(clab2,'new',3))irfeap = 1 ! restart program with new file
          go to 249
        end if
      else if(pcomp(clab1,'end ',4)) then
        go to 150
      end if

c.... set execution flag
      lo      = ll
      lct(ll) = clab2

c.... look in list
      do 110 j = 1,nwd
110   if(pcomp(clab1,wd(j),4)) go to 130

c.... look at existing procedure
      call proced(yy,yy(31),wd,nwd,ll,jct,lct,ct,flg,2)
      if(flg) go to 135
      call errclr('PMACIO')
      ll = ll - 1
      go to 100

130   jct(ll) = j

c.... save information for history
135   if (hadd) then
        do lc = lo,ll
          nn = nn + 1
          if(nn.le.size(js)) then
            js(nn) = jct(lc)
            ljs(nn) = lct(lc)
            vjs(1,nn) = ct(1,lc)
            vjs(2,nn) = ct(2,lc)
            vjs(3,nn) = ct(3,lc)
          else
            write(*,*) 'max number of macros reached'       
          end if
        end do ! lc
      end if
      nnn = nn
      if(ior.gt.0) write(iow,2000) clab1,clab2,ct(1,ll),ct(2,ll)
     1     ,ct(3,ll)
      if (ior.lt.0) go to 140
      go to 100

140   ll = ll + 1
150   jct(ll)= nnx

c...  set ascii/bin for rest via end
      asc1=ct(1,ll)
      asc2=ct(2,ll)

c     Check loop/next pairs
      j = 0
      do l = 1,ll
        if(jct(l).eq.nlp) j = j + 1
        if(jct(l).eq.nnx) j = j - 1
        if(j.lt.0) then
          if(ior.gt.0) then
            go to 402
          else
            ll  = ll - 2
            flg = .false.
            write(*,4002)
            go to 100
          end if
        end if
      end do ! l

      if(j.ne.0) then
        if(ior.gt.0) then
          go to 400
        else
          ll = ll - 1
          flg = .false.
          go to 100
        end if
      end if

      flg = .true.
c     Set loop/next markers
      do l = 1,ll-1
        if(jct(l).eq.nlp) then 
          j = 1
          k = l + 1
          sv = ct(2,l) ! from F83

          do i = k,ll
            if(jct(i).eq.nlp) j = j + 1
            if(j.gt.9) go to 401
            if(jct(i).eq.nnx) j = j - 1
            if(j.eq.0) go to 200
          end do ! i 
          go to 400

200       ct(1,i) = 0.0d0 ! from F83
          ct(2,i) = l
          ct(3,i) = sv    ! from F83
          
          ct(2,l) = i

        end if 
      end do ! l

c     Check if/endif pairs

      j = 0
      do l = 1,ll
        if(jct(l).eq.nif) j = j + 1
        if(jct(l).eq.nei) j = j - 1
        if(j.lt.0) then
          if(ior.gt.0) then
            go to 403
          else
            ll  = ll - 2
            flg = .false.
            write(*,4003)
            go to 100
          end if ! ior
        end if ! j
      end do ! l
      
      if(j.ne.0) then
        if(ior.gt.0) then
          go to 404
        else
          ll = ll - 1
          flg = .false.
          go to 100
        end if ! ior
      end if ! j

c     Set if/else/endif markers

      flg = .true.
      do l = 1,ll-1
        if(jct(l).eq.nif .or.jct(l).eq.nel) then

c         Locate the matching endif command

          j  = 1
          k  = 0
          do i = l+1,ll
            if(jct(i).eq.nel .and. k.eq.0) then
              k = i
            end if
            if(jct(i).eq.nif) j = j + 1
            if(j.gt.9) go to 401
            if(jct(i).eq.nei) j = j - 1
            if(j.eq.0) go to 205
          end do ! i
          go to 404

205       if(k.eq.0) then
            ct(2,l) = i - 1
          else
            ct(2,l) = k - 1
          end if
          ct(3,l) = i - 1

        end if
      end do ! l

      return

249   ll = -2  ! QUIT
      return
250   asc1 = ct(1,ll)
      asc2 = ct(2,ll)
      ll = -1  ! EXIT
      return

c.... error messages
400   write(iow,4000)
      if(ior.lt.0) write(*,4000)
cww   stop
      return
401   write(iow,4001)
      if(ior.lt.0) write(*,4001)
cww   stop
      return
402   write(iow,4002)
      if(ior.lt.0) write(*,4002)
cww    stop

403   write(iow,4003)
      if(ior.lt.0) write(*,4003)
cww   stop
      return

404   write(iow,4004)
      if(ior.lt.0) write(*,4004)
cww   stop
      return

2000  format(7x,a4,1x,a4,1x,3g12.5)

cww2001  format(a1,19a4,a3//'  m a c r o   i n s t r u c t i o n s'//
2001  format(/'  m a c r o   i n s t r u c t i o n s'/
     1  2x,'macro statement',2x,'variable 1',2x,'variable 2',
     2      2x,'variable 3')

cww2002  format(' Input a macro instruction: Enter "help" for list of ',
cww  1 'commands.'/' Enter "exit" to end with restart save, "quit" to ',
cww  2 'quit without restart save.')
cww2003  format('   List',i3,'  Macro',i3,'> ',$)

2003  format('Macro>',$)

cww2004  format(' *End of macro execution*',40x,'t=',f9.2)

3001  format(' *WARNING* Maximum number of command statements =',i4
     &     ,/'    Use history edit to reduce or program will stop'
     &     ,' when',i4,' is reached')

3002  format(' *ERROR* PMACIO: Maximum number of command instructions',
     &       '                  is limited to',i4)

4000  format(' error in pmacio ** unbalanced loop/next macros')

4001  format(' error in pmacio ** loops nested deeper than 8')

4002  format(' *ERROR* PMACIO: "loop" must precede "next" instruction')

4003  format(' *ERROR* PMACIO: "if" must precede "endi"f instruction')

4004  format(' *ERROR* PMACIO: Unbalanced if/endif commands')

      end
c
      subroutine histex(wd,clab1,jct,lct,ct,nwd,nlp,nnx,ll,is)
c----------------------------------------------------------------------
c
c      Purpose: Control command language execution by history inputs
c
c      Inputs:
c         wd(*)    - List of solution commands
c         clab1    - Command to execute
c         nwd      - Number of commands in wd
c         nlp      - Loop command number
c         nnx      - Next command number
c         ll       - Command number to execute

c      Outputs:
c         jct(ll)  - Command number
c         lct(ll)  - Command option
c         ct(3,ll) - Values n1,n2,n3
c         is       - Error if 0, otherwise 1
c    
c----------------------------------------------------------------------
      USE idata
      implicit real*8 (a-h,o-z)
      logical pcomp
      character*1 clab1(4),wd(4,nwd),y
      character*4 lct(*)
      integer jct(*)
      real*8  ct(3,*)
      is = 0
      if(nn.le.0) then
        write(*,2000)
        is = 1
      else
        if(clab1(2).ne.'!') then
          nc = 3
          if(clab1(4).eq.' ') nc = 2
          if(clab1(3).eq.' ') nc = 1
          if(clab1(2).eq.' ') nc = 0
        else
          nc = 0
        end if
        if(nc.gt.0) then
          do 101 n = nn,1,-1
            if(pcomp(clab1(2),wd(1,js(n)),nc)) go to 102
101       continue
          write(*,2001)
          is = 1
          return
        else
          n = nn
        end if
102     continue
        if(js(n).ne.nlp .and. js(n).ne.nnx) then
          jct(ll) = js(n)
          lct(ll) = ljs(n)
          ct(1,ll) = vjs(1,n)
          ct(2,ll) = vjs(2,n)
          ct(3,ll) = vjs(3,n)
103       write(*,2002) (wd(i,jct(ll)),i=1,4),lct(ll),(ct(i,ll),i=1,3)
          read(*,1000,err=103,end=900) y
cww104  if(y.ne.' ' .and. y.ne.'y' .and. y.ne.'Y') is = 1
          if(y.ne.' ' .and. y.ne.'y' .and. y.ne.'Y') is = 1
        else
          write(*,2003)
          is = 1
        end if
      end if
      return
c.... eof encountered
900   call  endclr ('PMACIO',y)
cww   goto 104
cww>
        if(y.ne.' ' .and. y.ne.'y' .and. y.ne.'Y') is = 1
      return
cww<
c
1000  format(a1)
2000  format(/' *** No previous instruction of this type executed')
2001  format(/' *** Error *** no match of this macro name')
2002  format(/' Macro to be executed.'/10x,'--> ',4a1,1x,a4,3e12.5
     1       /' Enter y or <CR> to accept.-->',$)
2003  format(/' *** Error *** loop/next execution not permitted')
      end
c
      integer function ipos(file,nn)
c----------------------------------------------------------------------
c
c      Purpose: determine length of a character string
c
c      Inputs:
c         file(nn) - name
c
c      Outputs:
c         ipos     - length of character string
c    
c----------------------------------------------------------------------
      character*1 file(nn)
      do 100 n = nn,1,-1
        if(file(n).ne.' ') go to 200
100   continue
      ipos = 0
      return
200   ipos = n
      return
      end
c
      integer function ipos1(file,nn)
c----------------------------------------------------------------------
c
c      Purpose: determine length of a character string until '\'
c
c      Inputs:
c         file(nn) - name
c
c      Outputs:
c         ipos     - length of character string
c    
c----------------------------------------------------------------------
      USE pdata2
      character*1 file(nn),st
      if(idev.le.2) st = '/'
      if(idev.eq.3) st = '\'
      if(idev.eq.4) st = '\'
      do 100 n = nn,1,-1
        if(file(n).eq.st) go to 200
100   continue
      ipos1 = 0
      return
200   ipos1 = n
      return
      end
c
      subroutine plotdf(plo,node,ndir,nstedf,mkflg,mmc,mmst,incmk,isw,
     1                  u,v,a,st,id,nneq,numnp)
c----------------------------------------------------------------------
c
c      Purpose: plot load deflection curve, time history etc.
c               main program
c
c      Input:
c        plo(10,*)  - array of plot values
c        node       - node to plot
c        ndir       - dof  to plot 
c        nstedf     - number of points of curve
c        mkflg      - mark curve if 1
c        mmc        - color of mark
c        mmst       - type pf marker (0-6)
c        incmk      - increment of mark
c        isw        - 1-11 type of curve, see PMACR
c        u(*)       - displacement vector
c        v(nneq)    - velocity     vector
c        a(nneq)    - acceleration vector
c        st(numnp,*)- stress vector
c        id(ndf,*)  - Equation numbers for each active dof
c        nneq       - No of unknowns in problem 
c        numnp      - No of nodes in problem 
c
c      Output:
c
c----------------------------------------------------------------------
      USE plodfs
      USE tplomax
      implicit double precision(a-h,o-z)
      dimension plo(10,*),u(*),v(nneq,*),a(nneq,*),st(numnp,*),id(*)
      character gchar*4,xname*12,yname*12
      data scalma,scalmi /1.00d0,1.00d0/
      if(isw.eq.1) then
        xname = 'Displacement'
        yname = 'Loadfactor'
        n1    = 2
        n2    = 1
      else if(isw.eq.2) then
        xname = 'Time'
        yname = 'Displacement'
        n1    = 3
        n2    = 2
      else if(isw.eq.3) then
        xname = 'Displacement'
        yname = 'Velocity'
        n1    = 2
        n2    = 4
      else if(isw.eq.4) then
        xname = 'Displacement'
        yname = 'Reaction'
        n1    = 2
        n2    = 5
      else if(isw.eq.5) then
        xname = 'Displacement'
        yname = 'Determinant'
        n1    = 2
        n2    = 6
      else if(isw.eq.6) then
        xname = 'Time'
        yname = 'Velocity'
        n1    = 3
        n2    = 4
      else if(isw.eq.7) then
        xname = 'Time'
        yname = 'Acceleration'
        n1    = 3
        n2    = 7
      else if(isw.eq.8) then
        xname = 'Time'
        yname = 'Loadfactor'
        n1    = 3
        n2    = 1
      else if(isw.eq.9) then
        xname = 'Time'
        yname = 'Reaction'
        n1    = 3
        n2    = 5
      else if(isw.eq.10) then
        xname = 'Time'
        yname = 'Determinant'
        n1    = 3
        n2    = 6
      else if(isw.eq.11) then
        xname = 'Displacement'
        yname = 'Time'
        n1    = 2
        n2    = 3
      else if(isw.eq.12) then
        xname = 'Displacement'
        yname = 'Stress'
        write(yname,'(a4,i1,a3,i4)') 'Stre',nstri,' No',nstrno 
        n1    = 2
        n2    = 8
      else if(isw.eq.13) then
        xname = 'Time'
        write(yname,'(a4,i1,a3,i4)') 'Stre',nstri,' No',nstrno 
        n1    = 3
        n2    = 8
      else if(isw.eq.14) then
        xname = 'Displacement'
        yname = 'Us1d'
        n1    = 2
        n2    = 9
      else if(isw.eq.15) then
        xname = 'Displacement'
        yname = 'Us2d'
        n1    = 2
        n2    = 10
      else if(isw.eq.16) then
        xname = 'Time'
        yname = 'Us1t'
        n1    = 3
        n2    = 9
      else if(isw.eq.17) then
        xname = 'Time'
        yname = 'Us2t'
        n1    = 3
        n2    = 10
      end if
c.... add values of actual iteration (only for plot)
      call plotdf1(u,v,a,st,id,plo,nneq,numnp,1)
c.... plots load deflection curve with different graphic routines
      if(nstedf.eq.1) return
c.... compute xmax,xmin,ymax,ymin
      if(imaxx.eq.0) then
        xmin  = 0.d0
        xmax  = plo(n1,1)
        do i = 2,nstedf
          if(plo(n1,i).gt.xmax)  xmax  = plo(n1,i)
          if(plo(n1,i).lt.xmin)  xmin  = plo(n1,i)
        enddo
      else
        xmin  = xmint
        xmax  = xmaxt     
      end if

      if(imaxy.eq.0) then
        ymin  = 0.d0
        ymax  = plo(n2,1)
        do i = 2,nstedf
          if(plo(n2,i).gt.ymax)  ymax  = plo(n2,i)
          if(plo(n2,i).lt.ymin)  ymin  = plo(n2,i)
        enddo
      else
        ymin  = ymint 
        ymax  = ymaxt     
      end if

      xmax  = xmax*scalma
      xmin  = xmin*scalmi
      ymax  = ymax*scalma
      ymin  = ymin*scalmi
      call gopen(xmin,xmax,-1.D0,ymin,ymax,-1.D0,xname,yname,1,1)
      call gmove1(1,plo(n1,1),plo(n2,1),1,1)
      call gwrite(1,0.15d0,0.10d0,'Node : ',1)
      write(gchar,'(i4)') node
      call gwrite(1,0.30d0,0.10d0,gchar,1)
      call gwrite(1,0.43d0,0.10d0,', Dir.: ',1)
      write(gchar,'(i4)') ndir
      call gwrite(1,0.58d0,0.10d0,gchar,1)
      inc = -1
      imc = mmc

      do 300 i=1,nstedf
        inc = inc+1
        call pppcol(7)
        call gdraw(1,plo(n1,i),plo(n2,i))
        if (mkflg.eq.1.and.inc.eq.incmk) then
         if(i.eq.nstedf) imc = 2                 ! last point in red
           call gmark(plo(n1,i),plo(n2,i),imc,mmst)
           inc = 0
        end if
300   continue
      call ginit(-1)
c.... set to original state
      call plotdf1(u,v,a,st,id,plo,nneq,numnp,2)
      return
      end
c
      subroutine plotdf1(u,v,a,st,id,plo,nneq,numnp,isw)
c----------------------------------------------------------------------
c
c      Purpose: store actual values in plo-array for plotting including
c               actual solution (only local)  
c
c      Input:
c        u(*)       - displacement vector
c        v(nneq)    - velocity     vector
c        a(nneq)    - acceleration vector
c        st(numnp,*)- stress vector
c        id(ndf,*)  - Equation numbers for each active dof
c        plo(10,*)  - array of plot values
c        nneq       - No of unknowns in problem 
c        numnp      - No of nodes in problem 
c        isw        - =1 add actual values
c                   - =2 set back to original values
c
c      Output:
c        plo(10,*)  - array of plot values
c
c----------------------------------------------------------------------
      USE arcl
      USE endata
      USE ext2
      USE fdata
      USE plodf
      USE plodfa
      USE plodfs
      USE plodfu
      USE rdata
      USE tdata
      implicit double precision (a-h,o-z)
      dimension u(*),v(nneq,*),a(nneq,*),plo(10,*),st(numnp,*),id(*)
c.... only in case of convergence
cww      if(abs(aengy).gt.tol*rnmax*1.d9.and.nploc.eq.1) return
c
      if(isw.eq.1) then
        nstedf = nstedf + 1
        if(nstedf.gt.1) then
          if(arcf.or.extflg) then
            plo(1,nstedf) = propld(ttim,0)*rlnew*pf
          else
            plo(1,nstedf) = propld(ttim,0)*pf
          end if
          plo(6,nstedf) = detc
          j = ipl(1,1)
          plo(2,nstedf) =  u(j)
          if(ipl(2,2).lt.0) plo(2,nstedf) = -u(j)
          plo(3,nstedf) = ttim
          if(fl(9)) then
            n = id(j)
            if(n.gt.0) then
              plo(4,nstedf) = v(n,1)
              plo(7,nstedf) = a(n,1)
            end if
          end if
          if(flreac) then
            plo(5,nstedf) = react*pf
          end if
          if(nstri.ne.0) then 
            if(.not.fl(11)) then
              call drawmess('Calculate first nodal stresses',1,0)
              return
            end if  
            plo(8,nstedf)=st(nstrno,nstri)
          end if
          plo( 9,nstedf) = valuse1
          plo(10,nstedf) = valuse2
        end if
      else
        do i =1,10
          plo(i,nstedf) = 0.0d0
        enddo
      nstedf=nstedf-1
      end if
      return
      end      
c
      subroutine storpl(v,a,plo,nstedf,nneq,id,j)
c----------------------------------------------------------------------
c
c      Purpose: store velocity/acceleration component in plo-array
c
c      Input:
c        v(nneq)    - velocity     vector
c        a(nneq)    - acceleration vector
c        plo(10,*)  - array of plot values
c        nstedf     - number of points of curve
c        nneq       - No of unknowns in problem 
c        id(ndf,*)  - Equation numbers for each active dof
c        j          - Dof to plot 
c
c      Output:
c        plo(10,*)  - array of plot values
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension v(nneq,*),a(nneq,*),plo(10,nstedf),id(*)
      n = id(j)
      n = iabs(n)     !  then for all  nodes, else only for free nodes
      if(n.gt.0) then
         plo(4,nstedf) = v(n,1)
         plo(7,nstedf) = a(n,1)
      end if
      return
      end
c
      subroutine storps(st,plo,numnp,nstedf,nstri,nstrno)
c----------------------------------------------------------------------
c
c      Purpose: store stress nstri at node nstrno in plo-array
c
c      Input:
c        st(numnp,*) - stress vector
c        plo(10,*)   - array of plot values
c        nstedf      - number of points of curve
c        numnp       - No of nodes in problem 
c
c      Output:
c        plo(10,*)  - array of plot values
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension st(numnp,*),plo(10,nstedf)
      plo(8,nstedf) = st(nstrno,nstri)

      write( *,*) nstri,nstrno,st(nstrno,nstri)
      write(16,*) nstri,nstrno,st(nstrno,nstri)
   
      return
      end
c
      subroutine storpl1(u,v,a,plo,nneq,nstedf,id,ndf,pf)
c----------------------------------------------------------------------
c
c      Purpose: store selected node values to file for use with [tplo]
c               without n2 (sign of disp), n3 (factor for load/reac)   
c
c      Input:
c        u(*)        - displacement vector
c        v(nneq)     - velocity     vector
c        a(nneq)     - acceleration vector
c        plo(10,*)   - array of plot values
c        nneq        - No of unknowns in problem 
c        nstedf      - number of points of curve
c        id(ndf,*)   - Equation numbers for each active dof
c        ndf         - Number dof/node
c        pf          - multiplier for symmetry
c
c      Output:
c
c----------------------------------------------------------------------
      USE fdata
      USE plodfa
      implicit double precision (a-h,o-z)
      dimension u(*),v(*),a(*),plo(10,*),id(*)
      dimension uu(6),vv(6),aa(6),rr(6)
C
      write(25,'(a,i6,6(1x,g10.4))') 't',0,plo(3,nstedf),
     + plo(1,nstedf)/pf,plo(6,nstedf),plo(8,nstedf),plo(9,nstedf),
     + plo(10,nstedf) ! t dum,time,lambda,det,stress,valuse1,valuse2 

      do 100 i =1,10
        inode = noden(i)
        if(inode .eq.0) goto 100
        in1 = ndf*(inode-1)        ! field-index in u of node inode
        call pzero (uu,6)  
        call pzero (aa,6)  
        call pzero (vv,6)  
        call pzero (rr,6)  
        ndfmax = min(ndf,6)        ! only up to 6 dofs possible
        do ii = 1,ndfmax
          uu(ii) = u(in1+ii)       ! displacements of node inode
            if(flreac)then 
            rr(ii) = reacc(i,ii)   ! reactions     of node inode
          end if
            if(fl(9)) then
            vv(ii) = v(in1+ii)     ! velocities    of node inode
            aa(ii) = a(in1+ii)     ! acceleration  of node inode
          end if
        end do
          write(25,1000) 'd',inode,(uu(iii),iii=1,ndfmax) ! displ. to disk
          if(flreac) 
     +    write(25,1000) 'r',inode,(rr(iii),iii=1,ndfmax) ! reactions  
          if(fl(9)) then
          write(25,1000) 'v',inode,(vv(iii),iii=1,ndfmax) ! velocities
          write(25,1000) 'a',inode,(aa(iii),iii=1,ndfmax) ! acceleration
        end if                                      

100   continue
1000  format(a1,i6,6(1x,e10.4))
      return
      end
c
      subroutine vazero(urate,xm,neq,nneq,nrt)
c----------------------------------------------------------------------
c
c      Purpose: look for zero mass d.o.f. - set the rate terms zero
c
c      Input:
c        urate(*)   - velocity vector
c         neq       - Number of active equations
c         nneq      - Total number of terms in solution arrays numnp * ndf
c         nrt       -
c
c      Output:
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xm(neq),urate(nneq,nrt)
      do 100 n = 1,neq
        if(xm(n).eq.0.0d0) then
          do 50 i = 1,nrt
           urate(n,i) = 0.0d0
50        continue
        end if
100   continue
      return
      end
c
      subroutine pmesh(idl,ie,d,id,x,ix,f,t,f0,ndd,nie,ndf,ndm,nen1,iii,
     +                 prt)
c----------------------------------------------------------------------
c
c      Purpose: Data input routine for mesh description
c
c      Inputs:
c         idl(nst)    - scratch array
c         ie(nie,*)   - Assembly information for material set
c         d(ndd,*)    - Material set parameters
c         id(ndf,*)   - Equation numbers for each active dof
c         x(ndm,*)    - Nodal coordinates of mesh
c         ix(nen1,*)  - Element nodal connections of mesh
c         f(ndf,*)    - Nodal force and displacement values
c         t(*)        - Nodal temperature values
c         f0(ndf,*)   - Nodal force values F0(=without t)
c         ndd         - Dimension for d array
c         nie         - Dimension for ie array
c         ndf         - Number dof/node
c         ndm         - Spatial dimension of mesh
c         nen1        - Dimension for ix array
c         iii         - Initialization indicator
c         prt         - Flag, print input data if true
c
c      Outputs:
c         iii         - Initialization indicator
c         Depends on commands specified
c
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE codat
      USE comfil
      USE debugs
      USE dirdat
      USE eldata
      USE errblk
      USE errchk
      USE hdata
      USE inptc
      USE iodata
      USE iofile
      USE iosave
      USE isogeo
      USE jinteg
      USE mdat2
      USE pcrit
      USE pdata2
      USE plong
      USE prisdat
      USE qload
      USE rsum
      USE sldata
      USE slid3
      USE smpak
      USE soltyp
      USE ximp
      USE yydata
      USE doalloc
      implicit real*8 (a-h,o-z)
      logical prt,error,pcomp,lexist,lopen,fldir,flqloa,flrsum
      logical flnmp, flkv1, flkv2
      character*4 wd(92),cc,c2,yyy1*100
cww   character*4 an(3)
c     character*4 va(2)
      character fnamr*229,fnams*229,fext*4
      integer ie(nie,*),id(ndf,*),ix(nen1,*),idl(*)
      integer ice(nen),lce(nen)
      real*8  d(ndd,*),x(ndm,*),f(ndf,*),f0(ndf,*),t(*),td(16)
      real*8,  allocatable, dimension(:) :: rtemp
      integer, allocatable, dimension(:) :: itemp
      save   isfile,llo
      data wd/'coor','elem','mate','boun','forc','temp','end ','prin',
     1      'nopr','page','bloc','pola','ebou','angl','sloa','para',
     2      'sphe','btem','glob','icon','pars','nopa','trib','jint',
     3        'vbou','gbou','rbou','base','el3b','link','opti','solv',
     4        'curv','edge','debu','load','eloa','tran','rot ','impf',
     5        'cmod','nege','neco','back','geom','segm','regi','pres',
     6        'fixe','poin','aloa','eang','blox','rsum','disp','gmes',
     7        'gele','gcor','ndvi','rndm','qloa','isec','icor','elfr',
     8        'vang','inte','stop','blco','bsys','ynod','yedg','ybou',
     9        'yloa','mes1','mes2','mes3','mes4','mes5','edis','elec',
     +        'epsq','cons','knv1','knv2','nmpq','loa0','feap','back',
     1        'cylt','macr','tie ','ycon'/ 
      data list/92/
cww   data an/'  an','gles','    '/
c     data va/' val','ue  '
cww   data ldir/0/,mdir/1/,xdir/4*0.d0/,knode/0/,fldir/.false./
      ldir =0
      mdir =1
      xdir=0.d0
      knode=0
      fldir=.false.
cwd   initialization of NURBS values
      NURn=0;NURm=0;NURp=0;NURq=0      


c.... macros which are in WD as dummys
c     30 link
c     32 opti
c     66 inte
c     67 stop
c     87 feap
c     88 back
c     89 cylt
c     90 macr
c     91 tie
c     92 ycon 


c.... macros which are not documented
c     BSYS ok read COOR,ELEM,MATE from file Bname-binary  
c     PAGE ok set printer page eject code
c     ELEC ok modify node numbering on element

c.... macros which are obsolete
c     CONS obsolete but still ok 
c     FORC obsolete but still ok


c.... initialize arrays and set error detection values
      error = .false.
      ll = 1
      jflgu = .false.
      if(iii.ge.0) then ! 1. call, not for call from pmacr! 
        istyp = 0
        mqloa = 1 
        flqloa= .true.
        coflg = .false.
        lread = .false.
        lsave = .false.
        lfile = ios
        clfl  = .true.
        contfl = .false.
        prt    = .true.
        nums   = 0
        iblk   = 0
        debug  = 0
        dbgtxt = ' '
c...    set array for principal stresses to default values
        nprip(1) = 1 
        nprip(2) = 2 
        nprip(3) = 3 
        nprip(4) = 4 
        nprip(5) = 5 
        nprip(6) = 6 
        nprip(7) = 7 
        nprip(8) = 8 
        nptyp = 1 
c
        bang(1:size(bang)) = 0.d0
c....   set boundary code/forced values to zero
        do 101 n = 1,numnp
          do 100 i = 1,ndf
            id(i,n) = 0
            f(i,n) = 0.0
100       continue
          if(iii.eq.0) then
c....   set nodal coordinate/temperatures values to specified values
            x(1,n) = -999.
            t(n) = 0.0
        end if
101     continue
        if(iii.eq.0) then
c.... set material number on each element to zero
          do 102 n = 1,numel
            ix(nen1,n) = 0
102       continue
        end if
      end if
103   if(ior.lt.0) write(*,2008) ll
      call pintio(yyy,10)
      read(yyy,1000,err=110,end=900) cc,c2
      if(inptctrl.eq.1) call perform1(2,cc)
      if(pcomp(cc,'read',4)) then
        if(pcomp(c2,'end',3)) then
          inquire(file=fnamr,opened=lopen,exist=lexist)
          if(lexist.and.lopen) close(lfile)
          if(lread) then
            ior   = isfile
            lread = .false.
            ll    = llo
          end if
        else
          fnamr =  finp
          fext  =  c2
          call addext(fnamr,fext)
        inquire(file=fnamr,exist=lread)
        if(lread) then
            open(ios,file=fnamr,status='old')
        else
            write(yyy1,2014) fnamr(1:55)
          call drawmess(yyy1,1,0)
        end if
 
          if(.not.lread) go to 103
          llo    = ll
          isfile = ior
          ior    = lfile
        end if
        go to 103
      end if
      if(pcomp(cc,'save',4)) then
        if(pcomp(c2,'end',3)) then
          inquire(file=fnams,opened=lopen,exist=lexist)
          if(lexist.and.lopen.and.lsave) then
            write(lfile,2009)
            close(lfile)
          end if
          lsave = .false.
        else
          fnams =  finp
          fext  =  c2
          call addext(fnams,fext)
        inquire(file=fnams,exist=lsave)
        if(lsave) then
            open(ios,file=fnams,status='old')
        else
          open(ios,file=fnams,status='new')
        end if
          lsave =  .true.
        end if
        go to 103
      end if
      if(ior.lt.0.and.pcomp(cc,'help',4)) then
        jflag = 0
        do 104 jj = 1,list
          if(pcomp(c2,wd(jj),4)) then
            call pman(1,wd(jj))
            jflag = 1
          end if
104     continue
        if(jflag.eq.0) call phelpmp(wd,list,'MESH',c2)
        go to 103
      end if
      go to 120
110   write(*,2018) yyy
      call  errclr ('PMESH ')
      go to 103
120   do 130 i = 1,list
130   if(pcomp(cc,wd(i),4)) go to 140
      if(ior.lt.0)  call errclr('PMESH ')
      go to 103
140   ll = ll + 1
      if(inptctrl.eq.1) call perform1(1,cc) ! show actual position of input
c
c----------------------------------------------------------------------
c             c  e  m  b  f  t  e  p  n  p  b  p  e  a  s  p  s       |
c             o  l  a  o  o  e  n  r  o  a  l  o  b  n  l  a  p       |
c             o  e  t  u  r  m  d  i  p  g  o  l  o  g  o  r  h       |
c             r  m  e  n  c  p     n  r  e  c  a  u  l  a  a  e       |
      go to ( 1, 2, 3, 4, 5, 6, 7, 8, 8,10,11,12,13,14,15,16,17,
c----------------------------------------------------------------------
c             t  g  i  p  n  t  j  v  g  r  d  e  l  o  s  c  e       |
c             b  l  c  a  o  r  i  b  b  b  r  l  i  p  o  u  d       |
c             l  o  o  r  p  i  n  o  o  o  e  3  n  t  l  r  g       |
c             o  b  n  s  a  b  t  u  u  u  c  b  k  i  v  v  e       |
     1       18,19,20,21,21,23,24,25,26,27,28, 2,30,31,32,33,34,
c----------------------------------------------------------------------
c             d  l  e  t  r  i  x  n  n  b  g  s  r  p  f  p  a       |
c             e  o  l  r  o  m  m  e  e  a  e  e  e  r  i  o  l       |
c             b  a  o  a  t  p  o  g  c  c  o  g  g  e  x  i  o       |
c             u  d  a  n     f  d  e  o  k  m  m  i  s  e  n  a       |
     2       35,36,37,38,39,40,41,42,42,42,42,42,42,42,42,50,51,  
c----------------------------------------------------------------------
c             e  b  r  d  g  g  g  n  r  q  i  i  e  v  i  s  b       |
c             a  l  s  i  m  e  c  d  n  l  s  c  l  a  n  t  l       |
c             n  o  u  s  e  l  o  v  d  o  e  o  f  n  t  o  c       |
c             g  x  m  p  s  e  r  i  m  a  c  r  r  g  e  p  o       |
     3       52,53,54,55,59,59,59,59,60,61,62,63,64,65,66,67,68,     
c----------------------------------------------------------------------
c             b  y  y  y  y  m  m  m  m  m  e  e  e  c  k  k  n       |
c             s  n  e  b  l  e  e  e  e  e  d  e  p  o  n  n  m       |
c             y  o  d  o  o  s  s  s  s  s  i  l  s  n  v  v  p       |
c             s  d  g  u  a  1  2  3  4  5  s  c  q  s  1  1  q       |
     4       69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,             
c----------------------------------------------------------------------
c             l  f  b  c  m  t  y       |
c             o  e  a  y  a  i  c       |
c             a  a  c  l  c  e  o       |
c             0  p  k  t  r     n       |
     5       86,87,88,89,90,91,92),i             
c----------------------------------------------------------------------


c     87 feap
c     88 back
c     89 cylt
c     90 macr
c     91 tie
c     92 ycon 


c
c.... [COOR] nodal coordinate data input
1     if(numnp.gt.numnpic) inptctrl=1
      call genvec(ndm,x,' coordinates',prt,error,.true.)
      go to 103
c.... [ELEM] element data input and EL3B
2     continue
      if(numel.gt.numelic) inptctrl=1
      ityp = i
      l    = 0
      ilx  = 0
      if(ior.lt.0) write(*,2010)
      if(inptctrl.eq.1.and.numel.ne.0) call perform2('elem',1,0,0)
      ipos = 0
      do 220 i = 1, numel, 50 

        ipos = ipos+1
        if(inptctrl.eq.1.and.numel.ne.0.and.ipos.eq.4) then
          ipos = 0
          call perform2('elem',2,i,numel) ! plot actual state every 200 elements
        end if

        if (prt) then
cww                    write(iow,2001) o, head, (k,k=1,nen)
cww       if(ior.lt.0) write(*  ,2001) o, head, (k,k=1,nen)
                       write(iow,2001)          (k,k=1,nen)
          if(ior.lt.0) write(*  ,2001)          (k,k=1,nen)
        end if
        j = min(numel,i+49)
        do 215 n = i, j           ! loop over 50 el to print
          if (l .lt. n) then
            if (.false.) then
              call  errclr ('PMESH ')
            end if
            llx = ilx
c....       input the element records - N.B. limit is 16 nos. / record
            il = min(nen+3,16)
200         call dinput(td,il)
            if(errck) go to 200
            l  = td(1)
            lk = td(2)
            do 201 k = 1,min(nen,14)
              idl(k) = td(k+2)
201         continue
            if(nen.gt.13) then
              inc = 1
              do 204 ii = 1,(nen+3)/16
                if(ii.eq.(nen+3)/16) inc = 0
                is = il+1
                il = min(is+15,nen+3)
202             call dinput(td,il-is+1)
                if(errck) go to 202
                do 203 k = 1,il-is + inc
                      idl(k+is-3) = td(k)
203             continue
204           continue
              il = il -is + 1
            end if
            lx  = td(il)
            ilx = lx
            if (l .eq. 0) then
              l = numel + 1
            end if
            if (lx .eq. 0) then
              lx = 1
            end if
          end if
          if (l .lt. n) then  !  element l appears after element n
cww                      write(iow,3001) l,n
cww                      write(*  ,3001) l,n
cww         if(ior.lt.0) write(*  ,3001) l,n
            write(yyy,3001) l,n
            call drawmess(yyy,1,0)
            error = .true.
          else if ((l .gt. n) .and. (llx .ne. 0)) then !  generate elements
            do 210 k = 1, nen
c....         do not increment for node 3 = 3d beam = EL3B
              ix(k,n) = ix(k,n-1) + nx
              if(ityp.eq.29.and.k.eq.3) ix(k,n) = ix(k,n-1)
              if (ix(k,n-1) .eq. 0) then
                ix(k,n) = 0
              end if
              if ((ix(k,n) .gt. numnp) .or. (ix(k,n) .lt. 0)) then ! wrong node
cww                          write(iow,3002) n
cww                          write(*  ,3002) n
cww             if(ior.lt.0) write(*  ,3002) n
                write(yyy,3002) n
                call drawmess(yyy,1,0)
                error = .true.
              end if
210         continue
            ix(nen1,n) = ix(nen1,n-1)
          else
            nx = lx
            do 205 k = 1, nen
              if ((idl(k) .gt. numnp) .or. (idl(k) .lt. 0)) then ! wrong node
                          write(iow,3002) n
                          write(*  ,3002) n
cww          if(ior.lt.0) write(*  ,3002) n
cww                       write(yyy,3002) n
cww                  call drawmess(yyy,1,0)
                error = .true.
              end if
cww                          ix(k,l) = idl(k)    ! set nodes for element l
              if(l.le.numel) ix(k,l) = idl(k)    ! set nodes for element l
205         continue
cww                        ix(nen1,l) = lk
            if(l.le.numel) ix(nen1,l) = lk
          end if
          if ((prt) .and. (.not. error)) then
            write(iow,2002) n, ix(nen1,n), (ix(k,n),k=1,nen)
            if(ior.lt.0) then
              write(*,2002) n, ix(nen1,n), (ix(k,n),k=1,nen)
            end if
          end if
215     continue
220   continue
        if(inptctrl.eq.1.and.numel.ne.0) call perform2('elem',3,0,0)
      go to 103
c.... [MATE] material data input
cww3  if(prt) write(iow,2004) o,head
3     if(prt) write(iow,2004) 
      do 307 n = 1,nummat
300   if(ior.lt.0) write(*,2011)
      ma =0
      iel=0
c.... input the mate records - limit is 16 nos. / record
      il = min(ndf+2,16)
      call dinput(td,il)
      if(errck) go to 300
      ma  = td(1)
      if(ior.lt.0 .and. (ma.lt.0 .or.ma.gt.nummat)) go to 300
      if(ma.le.0 .or. ma.gt.nummat) go to 103
      iel = td(2)
      if(ie(nie-1,ma).ne.0 .and. iel.le.0) iel = ie(nie-1,ma)
      do 302 k = 1,min(ndf,14)
        idl(k) = td(k+2)
302   continue
      if(ndf.gt.14) then
        do 305 ii = 1,(ndf+2)/16
        is = il+1
        il = min(is+15,ndf+2)
303     call dinput(td,il-is+1)
        if(errck) go to 303
        do 304 k = 1,il-is+1
          idl(k+is-3) = td(k)
304     continue
305     continue
      end if
c
      do 313 i = 1,ndf
        ie(i,ma) = idl(i)
313   continue
      do 314 i = 1,ndf
        if(idl(i).ne.0) go to 316
314   continue
c.... reset all zero inputs
      do 315 i = 1,ndf
        ie(i,ma) = i
315   continue
316   ie(nie-1,ma) = iel
c.... set flags for number of history terms
      mct = 0
      nh1 = 0
      nh3 = 0
c.... obtain inputs from the element routine
cww   if(prt) write(iow,2003) ma,iel,(i,i=1,ndf),(ie(i,ma),i=1,ndf)
      if(prt) write(iow,2003) ma,iel,(i,i=1,ndf)                     
      if(prt) write(iow,2015)        (ie(i,ma),i=1,ndf)
      call elmlib(d(1,ma),idl,x,ix,idl,idl,idl,idl,idl,idl,
     1            ndf,ndm,ndf,iel,1)

c.... set number(here=nh1,nh3) of history terms for material
      if(nh1.ne.0) then
        ie(nie,ma) = nh1
      else
        ie(nie,ma) = mct
      end if
      if(nh3.ne.0) then
        ie(nie-2,ma) = nh3
      else
        ie(nie-2,ma) = mct
      end if
307   continue
      go to 103
c.... [BOUN] Read in the restraint conditions for each node
cww4     if(prt) write(iow,2000) o,head,(i,i=1,ndf)
4     if(prt) write(iow,2000) (i,i=1,ndf)
c.... reset +0/1 in case of input from macro, then profil(1)=setq necessary!
      if(ior.lt.0) then
        do in=1,numnp
          do idf = 1,ndf
            if(id(idf,in).gt.0) then
              id(idf,in)=0
            else
              id(idf,in)=1
            end if
          end do
        end do                                
        write(*,*)'Profile destroyed by new B.C.:Use Macr Nsys to reset'
      end if
c----------
      iii= 1
      n  = 0
      ng = 0
400   l  = n
      lg = ng
401   if(ior.lt.0) write(*,2012)
c.... input the restraint records - limit is 16 nos. / record
      il = min(ndf+2,16)
      call dinput(td,il)
      if(errck) go to 401
      n  = td(1)
      ng = td(2)
      do 402 k = 1,min(ndf,14)
        idl(k) = td(k+2)
402   continue
      if(ndf.gt.14) then
          do 405 ii = 1,(ndf+2)/16
            is = il+1
            il = min(is+15,ndf+2)
403         call dinput(td,il-is+1)
            if(errck) go to 403
            do 404 k = 1,il-is+1
              idl(k+is-3) = td(k)
404         continue
405       continue
      end if
cww   if(n.le.0.or.n.gt.numnp) go to 410
      if(n.eq.0) goto 410
      if(n.lt.0.or.n.gt.numnp) then
cww     write(iow,3005) n
cww     write(*  ,3005) n
        write(yyy,3005) n
        call drawmess(yyy,1,0)
        error = .true.
cww      goto 7
      end if
      do 407 i = 1,ndf
        idli = idl(i)
        id(i,n) = idli
        if (l.eq.0) then
          idil = 0
        else
          idil = id(i,l)
        end if
cww 407 if(l.ne.0.and.idl(i).eq.0.and.id(i,l).lt.0) id(i,n) = -1
407   if(l.ne.0.and.idli.eq.0.and.idil.lt.0) id(i,n) = -1
      lg = isign(lg,n-l)
408   l = l + lg
      if((n-l)*lg.le.0) go to 400
      do 409 i = 1,ndf
409   if(id(i,l-lg).lt.0) id(i,l) = -1
      go to 408
410   do 413 n = 1,numnp
      do 411 i = 1,ndf
411   if(id(i,n).ne.0) go to 412
      go to 413
412   if(prt) write(iow,2007) n,(id(i,n),i=1,ndf)
413   continue
      go to 103
c.... [FORC/LOAD] force/load data input same macro for forc,load,disp
cww5  call genvec(ndf,f,' force/displ',prt,error,.false.)
5     call genvec(ndf,f,' ext. loads ',prt,error,.false.)
      go to 103
c.... [TEMP] temperature data input
6     call genvec(1,t,' temperature',prt,error,.false.)
      go to 103
c.... [END] end of mesh data inputs
cww7  if(error) stop
7     if(error)
     +  call drawmess(' Error in Mesh data, modify input file !!',1,0)
      if(lread) then
        close(lfile)
        ior   = isfile
        lread = .false.
        ll    = llo
      end if
      if(inptctrl.eq.1) call perform1(2,cc)
      return
c.... [PRIN] [NOPR] print/noprint of data
8     prt = i.eq.8
      go to 103
c.... [PAGE] set printer page eject code
10    if(ior.gt.0) read(ior,1000,err=111) o
      if(ior.lt.0) read(*,1000,err=111,end=111) o
      go to 103
111   call  errclr ('PMESH ')
      go to 10
c.... [BLOC] generate block of nodes and elements
11    if(iii.lt.0) write(iow,3003)
      call blkgen(ndm,ndf,nen,nen1,x,ix,id,prt,iblk)
      if(iblk.eq.1) then
        error=.true.
cww        goto 7
      end if
c>>wd   quick'n'dirty for perfect circle
!      n_changed=0
!      n=100
!      do i=n+1,(n+1)*(n+1),(n+1)
!        xroot=sqrt(x(1,i)**2+x(2,i)**2)
!       if (xroot.ne.0.0d0) then
!            if (abs(x(1,i)).gt.abs(x(2,i))) x(1,i)=-sqrt(9-x(2,i)**2)
!            if (abs(x(1,i)).lt.abs(x(2,i))) x(2,i)=sqrt(9-x(1,i)**2)
!!            x(1,i-1)=0.5*(x(1,i-2)+x(1,i))
!!            x(2,i-1)=0.5*(x(2,i-2)+x(2,i))
!            n_changed=n_changed+1
!        end if
!      end do
!      do i=(n+1)*(n+1)+(n+1),2*(n+1)*(n+1),(n+1)
!        xroot=sqrt(x(1,i)**2+x(2,i)**2)
!        if (xroot.ne.0.0d0) then
!            if (abs(x(1,i)).gt.abs(x(2,i))) x(1,i)=-sqrt(9-x(2,i)**2)
!            if (abs(x(1,i)).lt.abs(x(2,i))) x(2,i)=sqrt(9-x(1,i)**2)
!!            x(1,i-1)=0.5*(x(1,i-2)+x(1,i))
!!            x(2,i-1)=0.5*(x(2,i-2)+x(2,i))
!            n_changed=n_changed+1
!        end if
!      end do
c<<wd      
      go to 103
c.... [POLA] convert polar to cartesian coordinates
12    call polar(x,ndm,prt)
      go to 103
c.... [EBOU] set edge boundary constraints
c.... [EBOU,GAP]  set edge boundary constraints with GAP=val
13    call pedges(x,id,idl,ndm,ndf,numnp,c2,prt)
      go to 103
c.... [ANGL] set local coordinates system at single node, only 1 choice possible for ang+eang+vang!
14    call dinput(td,3)
      iang1 = td(1)
      iang2 = td(2)
      if(iang1.eq.0) iang1=1 
      if(iang2.eq.0) iang2=2
      itrot = td(3)
      ia(1) = iang1
      ia(2) = iang2
      if(ndf.eq.2) then
        if(iang1.ne.1) stop 'error input angl1 ANGL'
        if(iang2.ne.2) stop 'error input angl2 ANGL'
      else if(ndf.gt.2) then
        if(iang1.lt.1.or.iang1.gt.3) stop 'error input angl1 ANGL'
        if(iang2.lt.1.or.iang2.gt.3) stop 'error input angl2 ANGL'
      end if
      call genvec(1,bang,'  angles    ',prt,error,.false.)
      if(prt) write(iow,2013) iang1,iang2
      go to 103
c.... [SLOA] surface loadings
15    j = 0
      call ploadi
      go to 103
c.... [PARA] set parameter variables
16    coflg = .true.
      call pconst(prt)
      go to 103
c.... [SPHE] convert polar to cartesian coordinates
17    call sphere(x,ndm,prt)
      go to 103
c.... [BTEM] input block of interpolated temperatures
18    call blktem(ndm,t,prt)
      go to 103
c.... [GLOB] specify global nodes (if any) mot used
c     must be deleted
19    continue
      go to 103
c.... [ICON] input and check contact data
20    call incon(x,ndm,prt)
      contfl = .true.
      go to 103
c.... [PARS/NOPA] parsing/noparsing of statements
21    coflg = i.eq.21
      go to 103
c 22  NOPA under 21 
c
c.... [TRIB] triangular block generator
23    if(iii.lt.0) write(iow,3003)
      call blktri(ndm,nen,nen1,x,ix,prt,iblk)
      if(iblk.eq.1) then
        error=.true.
cww        goto 7
      end if
      go to 103
c.... [JINT] input elementnumbers for J-integral
24    call dinput(td,1)
      njint = td(1)
      njint = max(njint,1)
      call ialloc(ajint,njint,'JINT',jflgu)
      call injint(njint,ajint,prt)
      go to 103
c.... [VBOU] set boundary constraints for a region
25    call pblock(x,id,idl,ndm,ndf,numnp,prt)
      go to 103
c.... [GBOU] set edge boundary constraints, generate for idf > 3
26    call pedgeg(x,id,idl,ndm,ndf,numnp,prt)
      go to 103
c.... [RBOU] set boundary constraints for a region in r,phi,z
27    call pblockr(x,id,idl,ndm,ndf,numnp,prt)
      go to 103
c.... [BASE] define nodal director field
28    ldir = 1
      call dinput(xdir,4)
      idtyp = xdir(1)
cwd   second parameter
      ibasepara = xdir(2)
      if(idtyp.eq.1.or.idtyp.eq.2) then
        knode=nen*numel
      else
        knode=numnp
      end if
      call ralloc(basea,10*2*knode,'BASE',fldir)
      mdir  = 10
cwd   for isogeometric analysis: director for every control point
c       is needed to determine number of rotational DOFS in function ro56      
      if (idtyp.eq.21) call pdireciga5(ibasepara,x,basea,numnp,ndm,
     1        1,numel,AIninc,AInien,AInipa,ngetnen(),prt,ipr)
      if (idtyp.eq.20) call pdireciga6(idtyp,x,basea,numnp,ndm,1,
     +             numel,prt)
cwd
c>>wd   for isogeometric analysis: director for every control point
c     is calculated globally, isw=18 is dummy routine, in isw=3
c     director vector is loaded with pdirec1
c     n1=100: director vector is averaged from linear control point mesh
c             tangential vectors from connection
c     n1=101: director vector is averaged from linear control point mesh
c             tangential vectors from neighbouring "element"  
c     n1=102: Normal vector from control point to shell lamina is used
!      if (idtyp.eq.100.or.idtyp.eq.101) then
!        call pdireciga1(idtyp,x,basea,numnp,ndm,1)
!        call pdireciga1(idtyp,x,basea,numnp,ndm,2)
!        test123=basea
!        test234=m(mdir+9)
!      end if
c<<wd
      go to 103
c29.. [EL3B] elem for 3d - 3 node beam  -> 2
c30/31[LINK, OPTI] dummy only for wd field
30    goto 103
31    goto 103
c
c     [SOLV], istyp
c
c      istyp = 0  => standard   solver FEAP
c            = 1  => SM         solver (CSR sparse matrix without minimum degree)
c            = 2  => SM         solver (CSR sparse matrix with    minimum degree)
c            = 3  => SuperLU    solver CSR 
c            = 4  => Pardiso    solver CSR  direct
c            = 5  => PBCG       solver CSR 
c            = 6  => PGMRES     solver CSR 
c            = 7  => PGMRES2    solver CSR 
c            = 8  => Pardiso/ML solver CSR  iterative
c            = 9  => simplex optimization process
c
c     Solver  typ    na(d)             nau          nal Pointer-arrays    
c-----------------------------------------------------------------------
c     0    Standard   K_D [neq]         K_U[jd(neq)] K_L jd
c     1-2  SM         K   [jd(neq+1)-1] .......          jd=ia kptr=ja   
c     3    SuperLU    K   [jd(neq+1)-1] .......       -  jd=ia ja ka(Diag)               
c     4    Pardiso    K   [jd(neq+1)-1] .......       -  jd=ia ja ka(Diag)               
c     5-7  PCG/PGMRES K   [jd(neq+1)-1] .......       -  jd=ia ja ka(Diag)               
c     8    Pardiso/ML K   [jd(neq+1)-1] .......       -  jd=ia ja ka(Diag)               
c
c     Comments:
c     length solver 1-8: sym: jd(neq+1)-1 = jd(neq) ->ok
c                           u-sym: jd(neq+1)-1 = jd(neq)+...(but <neq)->ok  
c     
c
32    if(ior.lt.0) write(*,2017)
      call dinput(td,7)
      istyp=td(1)
      do i = 1,6
        ctis(i) = td(i+1) 
      end do

      if(istyp.gt.9) then
         stop 'Solver only 0-9 possible'
      end if  

      if(idev.eq.3) then 
       if(istyp.eq.3) then
        if (prt) then
         write(  *,*)  
     +   'Solver 03 not possible for INTEL, solver 02 will be used!' 
         write(iow,*)
     +   'Solver 03 not possible for INTEL, solver 02 will be used!' 
        end if  
        istyp=2
       end if 
      end if

      if(idev.eq.4) then 
       if(istyp.eq.4) then
        if (prt) then
         write(  *,*)  
     +   'Solver 04 not possible for SALFORD, solver 02 will be used!' 
         write(iow,*)
     +   'Solver 04 not possible for SALFORD, solver 02 will be used!' 
        end if  
        istyp=2
       end if 
      end if

      goto 103
c
c.... [CURVE] input parameters for curves 
33    call curvin(m)
      goto 103
c
c.... [EDGE] set boundary constraints along line
34    call pedgeb(x,id,idl,ndm,ndf,numnp,prt)
      go to 103
c
c     [DEBU]  0: debug = false, 1: debug = true
35    call dinput(td,1)
      idbg = td(1)
      if(idbg.eq.1) debug =  1
      goto 103
c     [LOAD]  same as force
36    goto 5
c
c.... [ELOA]D set loads along lines
37    continue
c     temporary
      allocate(rtemp(numnp))
      allocate(itemp(numnp))
      call eload(x,f,id,rtemp,itemp,ndm,ndf,numnp,prt)
      deallocate(rtemp)
      deallocate(itemp)
      go to 103
c
c.... [TRAN]s define new coordinates xnew= xold+dx
38    call transx(x,ndm,numnp,prt)
      go to 103
c
c.... [ROT ] 
39    call rot(x,ndm,numnp,prt)
      go to 103
c
c.... [IMPF] nodal imperfection data input: X=X+a*PHI
c.... use only after EBOU, ELOAD etc.
40    continue
      call ralloc(uimp,numnp*ndm,'IMPF',flimp)  
      call genvec1(ndm,x,uimp,'x=x+u       ','u=a*u0      ',
     +prt,error,.true.)
      mimp=2 
      go to 103
c
c.... [cmod] modify  cartesian coordinates  due to type
c     [cmod, 1,na,ne,r   ] : Type1: sphere       z = sqrt(r*r-x*x-y*y)
c     [cmod, 2,na,ne,r,f ] : Type2: parabola     z = [-f/(r*r)](x*x+y*y)+f
c     [cmod, 3,na,ne,a   ] : Type3: hypar        z = a*x*y
c     [cmod, 4,na,ne,a,c ] : Type4: hyperboloid  r = a/c*sqrt(c**2+z**2)
c     [cmod, 5,na,ne,r,a ] : Type5: hyperboloid  z = r*(sqrt(1+x**2/a**2)-1)
c     [cmod, 6,na,ne,r,a ] : Type6: sphere       z = z + sqrt(r*r-y*y)-a
c     [cmod, 7,na,ne,r,a ] : Type7: twisted beam p = pi/2/a, y=y*cos p, z=y*sin p
c.... use only after EBOU, ELOAD etc.
41    call xmodif(x,ndm,numnp,prt)
      go to 103
c
c42-49[NEGE, NECO, BACK, GEOM, SEGM, REGI, PRES, FIXE] dummy only for wd field
42    goto 103
c
c.... [POIN]T data input(load,b.c.) for a point(x,y), node must be there!
50    call ppload(x,f,id,ndm,ndf,numnp,prt)
      go to 103
c.... [ALOA]D constant loads per area                                             
51    call aload(x,f,ix,id,ndm,ndf,nen1,prt)
      go to 103
c
c.... [EANG] set local coordinates system along line, only 1 choice possible for ang+eang+vang!
52    call dinput(td,3)
      iang1 = td(1)
      iang2 = td(2)
      if(iang1.eq.0) iang1=1 
      if(iang2.eq.0) iang2=2
      itrot = td(3)
      if(iang1.eq.0) iang1=1 
      if(iang2.eq.0) iang2=2
      ia(1) = iang1
      ia(2) = iang2
      if(ndf.eq.2) then
        if(iang1.ne.1) stop 'error input angl1 EANG'
        if(iang2.ne.2) stop 'error input angl2 EANG'
      else if(ndf.gt.2) then
        if(iang1.lt.1.or.iang1.gt.3) stop 'error input angl1 EANG'
        if(iang2.lt.1.or.iang2.gt.3) stop 'error input angl2 EANG'
      end if
      call panglee(x,bang,ndm,numnp,prt,iang1,iang2)
      go to 103
c
c.... [BLOX] generate block of nodes and elements with delta x
53    if(iii.lt.0) write(iow,3003)
      call blkgendx(ndm,nen,nen1,x,ix,prt,iblk)
      if(iblk.eq.1) then
        error=.true.
      end if
      go to 103
c
c.... [RSUM] generate nodes for reaction force in TPLO       
54    call ialloc(irpt,numnp,'RSUM',flrsum)
      call genrsum(irpt,x,prt)
      go to 103
c
c.... [DISP] displacement data input same macro for forc,load,disp
55    call genvec(ndf,f,' presc.displ',prt,error,.false.)
      go to 103
c56-59[GMESH, GELE, GCOR, NDVI] dummy only for wd field
59    goto 103
c
c.... [RNDM] nodal imperfection due to random parameter(0<rndm<1): data input: X=X+rndm*fact
c.... for nodes ia,ib,inc
c.... use only after EBOU, ELOAD etc.
60     call genvec2(ndm,x,'x=x+rnd*fact',prt,error,.true.)
      go to 103
c
c.... [QLOA]D 10 load values for element loads to be calculated on element level under isw=22                                            
cww61    mqloa = 1
61    if(flqloa) then
        call ralloc(aqloa,10*numel,'QLOA',flqloa)
        mqloa = 10
      end if 
      call qloadi(aqloa,ix,nen1,ie,nie,numel,prt)
      go to 103
c.... [ISEC] set global dof 6=1 and release at intersections
62    call isect(x,id,ndm,ndf,numnp,prt)
      go to 103
c
c.... [ICOR] imperfect nodal coordinates input: X_perf output: X_imperf
c.... use only after final definition of coordinates, bcs and loads
63    call genvec3(ndm,x,prt)
      go to 103
c      
c.... [ELFR]EE element input free of number
c     read numel elements without generation and without increasing sequence 
64    continue
      if(numel.gt.numelic) inptctrl=1
      if(inptctrl.eq.1.and.numel.ne.0) call perform2('elfr',1,0,0)
      ipos = 0
      iread= 0   
c.... input the element records - N.B. limit is 16 nos. / record
      il = min(nen+3,16)
640   call dinput(td,il)
      if(errck) go to 640
      l  = td(1)
      lk = td(2)
      iread = iread+1
      ipos  = ipos+1
      if(inptctrl.eq.1.and.numel.ne.0.and.ipos.eq.200) then
        ipos = 0
        call perform2('elfr',2,iread,numel) ! plot actual state every 200 elements
      end if
      if(l.eq.0) goto 642 ! finish
      do k = 1,min(nen,14)
        idl(k) = td(k+2)
      end do   
      if(nen.gt.13) then
        inc = 1
        do ii = 1,(nen+3)/16
          if(ii.eq.(nen+3)/16) inc = 0
          is = il+1
          il = min(is+15,nen+3)
641       call dinput(td,il-is+1)
          if(errck) go to 641
          do k = 1,il-is + inc
            idl(k+is-3) = td(k)
          end do  
        end do   
        il = il -is + 1
      end if
      do k = 1, nen
        if ((idl(k) .gt. numnp) .or. (idl(k) .lt. 0)) then ! wrong node
          write(yyy,3002) n
          call drawmess(yyy,1,0)
        end if
        ix(k,l) = idl(k)    ! set nodes for element l
      end do   
      ix(nen1,l) = lk       ! set increment (not necessary) 
      goto 640
c...  print
642     do i = 1,numel,50
        if (prt) then
                       write(iow,2001)          (k,k=1,nen)
          if(ior.lt.0) write(*  ,2001)          (k,k=1,nen)
        end if
        j = min(numel,i+49)
        do n = i, j           ! loop over 50 el to print
          if (prt) then
            if(ix(1,n).eq.0) then ! not used            
                 write(iow,2016) n
              if(ior.lt.0) then
                 write(  *,2016) n
              end if
            else 
                 write(iow,2002) n, ix(nen1,n), (ix(k,n),k=1,nen)
              if(ior.lt.0) then
                 write(  *,2002) n, ix(nen1,n), (ix(k,n),k=1,nen)
              end if
            end if
          end if
        end do
      end do  
      if(inptctrl.eq.1.and.numel.ne.0) call perform2('elfr',3,0,0)
      goto 103
c
c.... [VANG] set local coordinates system at volume, only 1 choice possible for ang+eang+vang!
65    call dinput(td,4)
      iang1 = td(1)
      iang2 = td(2)
      itrot = td(3)
      inpc  = td(4)
      if(iang1.eq.0) iang1=1 
      if(iang2.eq.0) iang2=2
      if(inpc.eq.0)  inpc=1
      ia(1) = iang1
      ia(2) = iang2
      if(ndf.eq.2) then
        if(iang1.ne.1) stop 'error input angl1 VANG'
        if(iang2.ne.2) stop 'error input angl2 VANG'
      else if(ndf.gt.2) then
        if(iang1.lt.1.or.iang1.gt.3) stop 'error input angl1 VANG'
        if(iang2.lt.1.or.iang2.gt.3) stop 'error input angl2 VANG'
      end if
      call panglev(x,bang,ndm,numnp,prt,iang1,iang2,inpc)
      go to 103
c
c66/67[INTE, STOP] dummy only for wd field
66    go to 103
67    go to 103
c
c.... [BLCO] generate block of nodes and elements
68    if(iii.lt.0) write(iow,3003)
      call blkgenco(ndm,nen,nen1,x,ix,prt,iblk)
      if(iblk.eq.1) then
        error=.true.
      end if
      go to 103
c
c.... [BSYS] read COORdinates and ELEMents from file Bname-binary 
c     NUMNP and NUMEL have to be specified for 1. Macro FEAP!!
69    call readerb(x,ix,ndm,numnp,nen,nen1,numel)
      go to 103
c
c---- [YNOD] dummy only for wd field (ylt-nodes of cylt)
70    goto 103
c
c---- [YEDG] dummy only for wd field (ylt-edges of cylt)
71    continue
      goto 103
c
c---- [YBOU] read boundary conditions for yield-line element
72    continue
      call yltbou(d(1,ma),x,ix,id,bang,ndf,ndm,prt,iang1,iang2)      
      goto 103
c
c---- [YLOA] read area/line/point load for yield-line procedure
73    continue
      call yltloa(d(1,ma),x,ix,id,f,ndf,ndm,prt)            
      goto 103      
c
c.... 'MES1' user input command
74    call umesh1(idl,ie,d,id,x,ix,f,t,ndd,nie,ndf,ndm,nen1,iii,prt)
      goto 103
c
c.... 'MES2' user input command
75    call umesh2(idl,ie,d,id,x,ix,f,t,ndd,nie,ndf,ndm,nen1,iii,prt)
      goto 103
c
c.... 'MES3' user input command
76    call umesh3(idl,ie,d,id,x,ix,f,t,ndd,nie,ndf,ndm,nen1,iii,prt)
      goto 103
c
c.... 'MES4' user input command
77    call umesh4(idl,ie,d,id,x,ix,f,t,ndd,nie,ndf,ndm,nen1,iii,prt)
      goto 103
c
c.... 'MES5' user input command
78    call umesh5(idl,ie,d,id,x,ix,f,t,ndd,nie,ndf,ndm,nen1,iii,prt)
      goto 103
c
c.... [EDIS] set prescribed values for edge boundary constraints
79    call pedgesd(x,f,ndm,ndf,numnp,prt)
      go to 103
c
c.... [ELEC] modify node numbering on element! ELEM has to be used before!
c     Pos1,Pos2,Pos3,Pos4,....,Posnn
c
c     not documented, necessary for input from other programs, 
c     typical element numbering clockwise, but in FEAP anti clockwise
c     solution: read elements with ELEM, modify position via ELEC
c
c     for mixed meshes 3/4node do not change node4
c     Node 4 is=0 for triangles in ELEM and
c     Node 4 is set to node 3 with ELEC, 
c     this is ok for e.g.B/D element, but not for DKQ  
c
80    call dinput(td,nen)
      do n=1,nen
        ice(n)=td(n) 
      end do 
      write(iow,2019) (ice(n),n=1,4)
      do n = 1,numel
        do nn = 1,nen 
          lce(nn) = ix(nn,n) ! copy 
        end do  
        do nn = 1,nen
          kk = ice(nn) 
          ix(nn,n) = lce(kk) ! new storage 
c         this is valid only for input mesh with 3/4 nodes
c         in case of node 4=0: node 4=node3 and mat no =2
          if(ix(nn,n).eq.0) then 
            ix(nn,n)   = ix(nn-1,n) 
            ix(nen1,n) = 2    
          end if
        end do  
      write(iow,*) n,(ix(nn,n),nn=1,nen)
      end do
      goto 103  
c.... [EPSQ] read strains for prescribed displ. for MACRO EPSQ
81    call epsq_me(prt)
      go to 103
c
c.... [CONS] set parameter variables = PARA 16
82    go to 16
c
c.... [KNV1] sets value of knot vector in xsi1-direction for isogeoFEAP
83    lenKnv1=0
c     max accuracy is 10+1 chars
      do 8301 i=1,NURnpatch
        lenKnv1 = lenKnv1 + (nNURnmpq(i,1)+nNURnmpq(i,3)+1)
8301  continue
      call setKnotVectLength(lenKnv1,1)
      call ralloc(AInkv1,lenKnv1,'KV1',flkv1)
      do 831 i=1,lenKnv1
        call dinput(td,1)
        call setKv(i,td(1),AInkv1,lenKnv1)
831   continue
      go to 103
c
c.... [KNV2] sets value of knot vector in xsi2-direction for isogeoFEAP
84    lenKnv2 = 0
c     max accuracy is 10+1 chars
      do 8401 i=1,NURnpatch
        lenKnv2 = lenKnv2 + (nNURnmpq(i,2)+nNURnmpq(i,4)+1)
8401  continue
      call setKnotVectLength(lenKnv2,2)
      call ralloc(AInkv2,lenKnv2,'KV2',flkv2)
      do 841 i=1,lenKnv2
        call dinput(td,1)
        call setKv(i,td(1),AInkv2,lenKnv2)
841   continue
      go to 103
c
c.... [NMPQ] sets order and number of control points for isogeoFEAP
c     has to be read before KNV1 and KNV1 and after coord!     
85    call dinput(td,1)
      NURnpatch = td(1)
      call ialloc(AInmpq,NURnpatch*6,'NMPQ',flnmp)
      do 8501 i=1,NURnpatch
        call dinput(td,6) ! for NURBS-Surface set to 6, for Solids to 8 
        call initIGA2(td,AInmpq,NURnpatch,i)
8501  continue
      call initIGA(ix,nen1,numnp,numel,nen,ipr)
c.... computing element data from NURBS input
c     data of last element into scratch array (seems rather nonsense)
c      do i = 1,nen
c        idl(i)=ix(i,numel)
c      end do
      go to 103
c
c.... [LOA0] load data for vector F0 (without t, e.g. dead load)
86    call genvec(ndf,f0,' loads F0   ',prt,error,.false.)
      go to 103
c
c     [FEAP, BACK, CYLT, MACR, TIE, YCON] dummy only for wd field
87    goto 103
88    goto 103
89    goto 103
90    goto 103
91    goto 103
92    goto 103

c
900   call  endclr ('PMESH ',cc)
cww   stop
      return
1000  format(a4,6x,a4)
cww2000  format(a1,19a4,a3//'     nodal b.c.'//
2000  format(/'     nodal b.c.'/
     1       6x,'node',9(i2,'-b.c.')/(10x,9(i2,'-b.c.')))
cww2001  format(a1,19a4,a3//5x,'elements'//4x,'elmt    matl',
2001  format(/5x,'elements'//6x,'elmt    matl',
     1   8(i3,' node')/(18x,8(i3,' node')))
2002  format(2x,10i8/(18x,8i8))
cww2003  format(/5x,'material set',i3,' for element type',i2,5x,//
cww     1        10x,'degree of freedom assignments    local    global' /
cww     2        42x, 'number',4x,'number'/(36x,2i10))
2003  format(5x,'material set',i3,' for element type',i3,/
     1       5x,'dofs relation: local  number:',20i3/)
2015  format(5x,'               global number:',20i3)
cww2004  format(a1,19a4,a3//'     material properties')
2004  format(/2x,'m a t e r i a l  p r o p e r t i e s')
2007  format(i10,9i7/(10x,9i7))
2008  format(' Enter "help" for list of commands, "end" to exit'/
     1       '     Mesh ',i3,'> ',$)
2009  format('read,end')
2010  format(' Input: elmt#, matl#, (ix(i),i=1,nen), inc'/3x,'>')
2011  format(' Input: matl#, elmt type'/3x,'>',$)
2012  format(' Input: node#, inc., (b. codes, i=1,ndf)'/3x,'>',$)
2013  format(//'    Angles set for directions ',i1,i1)
2014  format(' No file named ',a55,' exists.                      ')
2016  format(i10,' has not been input or generated')
2017  format(' Input: 1-7 Parameter for solver'/3x,'>',$)
2018  format(a80)
2019  format(' Element numbering from ELEM has been modified with ELEC',
     +/, 'Pos1,Pos2,Pos3,Pos4 = ',4i3,/)

cww3001  format(' **ERROR** element',i5,' appears after element',i5)
cww3002  format(' **ERROR** element',i5,' has illegal nodes')
3003  format(' **WARNING** element connections necessary to use '
     1      ,'block in macro program')
cww3004  format(' **ERROR**  no global nodes specified on the control')
3001  format(' Element',i7,' appears after element',i5)
3002  format(' Element',i7,' has illegal nodes')
3005  format(' Node number ',i7,' for BOUN not possible ')
      end
c
      subroutine addext(fnam,fext)
c----------------------------------------------------------------------
c
c      Purpose: adds character string to file fnam->fnam.fext
c
c      Input:
c         fnam(229)  -  character string without extension fnam
c         fext(4)    -  extension to add
c
c      Output:
c         fnam(229)  -  character string with extension    fnam.fext
c
c----------------------------------------------------------------------
      character fnam*229,fext*4    
      character*1 fnam1(229),fext1(4)
      do i = 1,229
        fnam1(i) = fnam(i:i)
      enddo
      do i = 1,4
        fext1(i) = fext(i:i)
      enddo
      iposl = ipos(fnam1,229)
      iposx = ipos(fext1,4)
      ii = iposl + 1
      do i = ii,229  
          fnam1(i) = ' '
      enddo    
      fnam1(ii) = '.'
      if((ii+iposx).gt.229) then
        call drawmess('Filename + Extension is to long(<229!)',1,0)
        return
      end if
      do i = 1,iposx
          fnam1(ii+i) = fext1(i)
      enddo    
      do i = 1,229
        fnam(i:i) = fnam1(i)
      enddo
      end
c
c
      subroutine pmove(a,b,nn)
c----------------------------------------------------------------------
c
c      Purpose: Move real array a into b
c
c      Inputs:
c         a(*)      - Array to move
c         nn        - Length of array to move
c
c      Outputs:
c         b(*)      - Moved array
c
c----------------------------------------------------------------------
      double precision a(nn),b(nn)
      do 100 n = 1,nn
100   b(n) = a(n)
      return
      end
c
      subroutine pmovec(id,a,b,nn)
c----------------------------------------------------------------------
c
c      Purpose:  Move compressed array a into uncompressed array b
c
c      Inputs:
c         id(*)     - Equation numbers for each active dof
c         a(*)      - Compressed array to move
c         nn        - Length of uncompressed array
c
c      Outputs:
c         b(*)    - Uncompressed move of a (zero undefined values)
c
c----------------------------------------------------------------------
      double precision a(nn),b(nn)
      integer id(*)
      call pzero(b,nn)
      do 100 n = 1,nn
        j = id(n)
        if (j.gt.0) b(n) = a(j)
100   continue
      return
      end
c
      subroutine pmoveca(id,a,b,nn)
c----------------------------------------------------------------------
c
c      Purpose:  ADD compressed array a into uncompressed array b
c
c      Inputs:
c         id(*)     - Equation numbers for each active dof
c         a(*)      - Compressed array to add
c         nn        - Length of uncompressed array
c
c      Outputs:
c         b(*)    - Uncompressed add of a (zero undefined values)
c
c----------------------------------------------------------------------
      double precision a(nn),b(nn)
      integer id(*)
      do 100 n = 1,nn
        j = id(n)
        if (j.gt.0) b(n) = b(n) + a(j)
100   continue
      return
      end
c
      subroutine pmovei(ia,ib,nn)
c----------------------------------------------------------------------
c
c      Purpose: Copy integer array ia into ib
c
c      Inputs:
c         ia(nn)    - Array to copy
c         nn        - Length of array to copy
c
c      Outputs:
c         ib(nn)    - Moved array
c
c----------------------------------------------------------------------
      integer ia(nn),ib(nn)
      do 100 n = 1,nn
        ib(n) = ia(n)
100   continue
      end
c
      subroutine pmoves(a,b,nn,fac)
c----------------------------------------------------------------------
c
c      Purpose: Copy real array a x fac into b 
c
c      Inputs:
c         a(nn)     - Array to copy
c         nn        - Length of array to copy
c         fac       - Factor to multiply
c
c      Outputs:
c         b(nn)     - Moved array
c
c----------------------------------------------------------------------
      double precision a(nn),b(nn),fac
      do 100 n = 1,nn
100   b(n) = a(n)*fac
      return
      end
c
      subroutine pmovedt(ip,b,numnp,ndf)
c----------------------------------------------------------------------
c
c      Purpose: add displacement values for tied nodes
c
c      Inputs:
c         ip(*)     - Node number list tied nodes
c         b(ndf,*)  - Nodal displacement values without tied nodes
c         numnp     - Number of nodes in mesh
c         ndf       - Number dof/node
c
c      Outputs:
c         b(ndf,*)  - Nodal displacement values for all nodes
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ip(*),b(ndf,*)
      do ni = 1,numnp
        n = iprttie(ip,ni)
        if(ni.ne.n) then
          do i = 1,ndf
            b(i,ni) = b(i,n)
          enddo
        end if
      enddo
      return
      end
c
c
      subroutine pmovest(ip,s,numnp,nstv)
c----------------------------------------------------------------------
c
c      Purpose: add stress values for tied nodes
c
c      Inputs:
c         ip(*)        - Node number list tied nodes
c         s(numnp,nstv)- Nodal stress values without tied nodes
c         numnp        - Number of nodes in mesh
c         nstv         - Number stress values/node
c
c      Outputs:
c         s(numnp,nstv)- Nodal stress values for all nodes
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ip(numnp),s(numnp,nstv)
      do ni = 1,numnp
        n = iprttie(ip,ni)
        if(ni.ne.n) then
          do i = 1,nstv
            s(ni,i) = s(n,i)
          enddo
        end if
      end do
      return
      end
c
c
      subroutine pnums(numnp,numel,ndm,nen,nummat,prt)
c----------------------------------------------------------------------+
c                                                                      |
c      Purpose: Determine maximum node and element number in mesh.     |
c               Determine number of materials (only if nummat=0)       |                                             |
c      Inputs:                                                         |
c         From read of data file(s)                                    |
c         prt    - Print input values if true
c
c      Outputs:                                                        |
c         numnp     - Maximum node    number                           |
c         numel     - Maximum element number                           |
c         nummat    - Maximum number of materials                      |
c                                                                      |
c      Comment:       read only necessary values                       |
c                     calculates always all 3 values                   | 
c                     for MATE only correct if each material           |
c                     has its own MATE-card                            | 
c                                                                      |
c       OPEN:  blox,btem,trib, total number to use e.g. for pola       |
c                                                                      |
c----------------------------------------------------------------------|  
      USE codat 
      USE iofile
      implicit  none
      logical   pcomp,prt
      character cc*4, yyy*80
      integer   numnp,numel,nummat           
      integer   numn,nume,num,n2,n3,n4,n5,n6,n7,n8,n9,
     +          ndm,nen,numna,numea, numm, nummi
      real*8    td(15)
      coflg = .false.
 
c     Use input data to determine total number of nodes or elmts
 
      numn = 0
      nume = 0
      numm = nummat

1     call pintio(yyy,8)
      read(yyy,1000,err=1) cc
 
      if(pcomp(cc,'cons',4).or.pcomp(cc,'para',4)) then      

        coflg = .true.
        call pconst(prt)

      else if(pcomp(cc,'coor',4)) then      
 
  100   call dinput(td,1)
        num = td(1)
        numn = max(numn,num)
        if(num.gt.0) go to 100

      else if(pcomp(cc,'elem',4)) then  

  110   call dinput(td,1)
        num = td(1)
        nume = max(nume,num)
        if((3+nen).gt.16) call dinput(td,1) ! dummy, max=29 nodes 
        if(num.gt.0) go to 110

      else if(pcomp(cc,'el3b',4)) then  

  115   call dinput(td,1)
        num = td(1)
        nume = max(nume,num)
        if((3+nen).gt.16) call dinput(td,1) ! dummy, max=29 nodes  
        if(num.gt.0) go to 115

      else if(pcomp(cc,'elfr',4)) then  

  118   call dinput(td,3+nen)
        num = td(1)
        nume = max(nume,num)
        if((3+nen).gt.16) call dinput(td,1) ! dummy, max=29 nodes  
        if(num.gt.0) go to 118

      else if(pcomp(cc,'bloc',4).or.pcomp(cc,'blox',4)) then  

        call dinput(td,8)
        n2  = td(2) ! nr
        n3  = td(3) ! ns 
        n4  = td(4) ! ni   nt
        n5  = td(5) ! ne   ni
        n6  = td(6) ! mat  ne
        n7  = td(7) ! inc  mat
        n8  = td(8) ! ntyp
c...... n8: 0-9,16 2D n8: 10-15,19 3D
        if(n8.lt.10.or.n8.eq.16) then   ! 2D
          if(n2.le.0.or.n3.le.0) then
            write(iow,3002) n2,n3
            return
          end if
          numna = (n2+1)*(n3+1)
          numea = n2*n3
          if(n7.ne.0) numna = numna+n3*n7      ! increments
          if(n8.ge.1.and.n8.le.6) then         ! lin.   triangles
            numea = numea*2
          else if(n8.eq.8.or.n8.eq.9) then     ! quadr. quadrilaterals
            numea = numea/4
          else if(n8.eq.16) then               ! cubic quadrilaterals
            numea = numea/9
          end if

          numn = max(numn,numna+n4-1)
          nume = max(nume,numea+n5-1)

        else if(n8.le.15.or.n8.eq.19) then  ! 3D
          if(n2.le.0.or.n3.le.0.or.n4.le.0) then
            write(iow,3002) n2,n3,n4
            return
          end if
          numna = (n2+1)*(n3+1)*(n4+1)
          numea = n2*n3*n4
          if(n8.eq.12.or.n8.eq.13.or.n8.eq.14) then  ! quadratic elements
            numea = numea / 8
          else if(n8.eq.19) then  ! cubic elements
            numea = numea / 27
          else if(n8.eq.15) then         ! biquadratic/linear
            numea = numea / 4
c            if(n4.ne.1) then
c              write(  *,*) 'Error Bloc tinc must be 1!'
c              write(iow,*) 'Error Bloc tinc must be 1!'
c              stop  
c            end if 
          else if(n8.eq.11) then         ! tetraheder
            numea = numea *6
          end if
          numn = max(numn,numna+n5-1)
          nume = max(nume,numea+n6-1)
        else
                       write(*  ,*)'ERROR - USER BLOCK NOT SUPPORTED'
          if(ior.lt.0) write(iow,*)'ERROR - USER BLOCK NOT SUPPORTED'
          return
        end if

      else if(pcomp(cc,'blco',4)) then  

        call dinput(td,9)
        n2 = td(2)  ! nr 
        n3 = td(3)  ! ns
        n4 = td(4)  ! ni
        n5 = td(5)  ! ne
        n6 = td(6)
        n7 = td(7)
        n8 = td(8)
        n9 = td(9)  ! 8/18/16 node elmts

        numna = (n2+1)*(n3+1)*2        ! 2 times nodes
        numea = n2*n3*3                ! 3 times elements 
        if(n9.eq.2) numea = n2*n3*2.25 ! 2+1/4times elements 
        if(n9.eq.3) numea = n2*n3*2.25 ! 2+1/4times elements 
        numn = max(numn,numna+n4-1)
        nume = max(nume,numea+n5-1)

      else if(pcomp(cc,'blox',4)) then
                     write(*  ,*) 'Pnums: blox not implemented'
        if(ior.lt.0) write(iow,*) 'Pnums: blox not implemented' 

      else if(pcomp(cc,'btem',4)) then
                     write(*  ,*) 'Pnums: btem not implemented'
        if(ior.lt.0) write(iow,*) 'Pnums: btem not implemented' 

      else if(pcomp(cc,'trib',4)) then
                     write(*  ,*) 'Pnums: trib not implemented'
        if(ior.lt.0) write(iow,*) 'Pnums: trib not implemented'

      else if( pcomp(cc,'mate',4).and.nummat.eq.0) then   ! only for nummat=0
        call dinput(td,1)
        nummi = td(1)
        numm  = max(nummi,numm)

      else if(pcomp(cc,'end',3)) then
        go to 200
      end if
      go to 1

 200  numnp  = numn
      numel  = nume
      nummat = numm 
      if(numnp.le.0 .or. numel.le.0 .or. nummat.le.0) then
                     write(iow,3001) numnp,numel,nummat
        if(ior.lt.0) write(*  ,3001) numnp,numel,nummat
        return
      end if
c     Find start of problem for input of mesh data

      rewind ior
      read(ior,1000) cc  ! first  line again
      read(ior,1000) cc  ! second line again


c     Input/Output Formats

1000  format(a4)

3001  format(' *ERROR* Problem does not permit a solution:'/
     1       '         Number of nodes         = ',i9/
     2       '         Number of elements      = ',i9/ 
     3       '         Number of materials     = ',i9/)

3002  format(' *ERROR* Number of block increments incorrect:'/
     1       '         Number of 1-increments = ',i6/:
     2       '         Number of 2-increments = ',i6/:
     3       '         Number of 3-increments = ',i6/)

      end
c
      subroutine polar(x,ndm,prt)
c----------------------------------------------------------------------
c
c      Purpose: convert polar to cartesian coordinates
c
c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh in polar coordinates
c         ndm       - Spatial dimension of mesh
c         prt       - Print flag
c
c      Outputs:
c         x(ndm,*)  - Nodal coordinates of mesh in cartesian coordinates
c
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      dimension x(ndm,*),td(6)
      character xnam*2,ynam*2,yyy*80
      if(ndm.eq.1) return
      mct = 0
      th = datan(1.0d0)/45.0
100   if(ior.lt.0) write(*,3001)
      call dinput(td,6)
      if(errck) go to 100
      ni  = td(1)
      ne  = td(2)
      inc = td(3)
      x0  = td(4)
      y0  = td(5)
      itype = td(6)
      if(ni.le.0) return
      if(ni.gt.numnp.or.ne.gt.numnp) go to 300
      inc = isign(max(iabs(inc),1),ne-ni)
      if(ne.eq.0) ne = ni
      n = ni
200   if(itype.eq.0 .or. itype. eq. 12) then
        r = x(1,n)
        x(1,n) = x0 + r*cos(x(2,n)*th)
        x(2,n) = y0 + r*sin(x(2,n)*th)
        xnam = 'x0'
        ynam = 'y0'
      else if(itype.eq.23 .and. ndm.eq.3) then
        r = x(2,n)
        x(2,n) = x0 + r*cos(x(3,n)*th)
        x(3,n) = y0 + r*sin(x(3,n)*th)
        xnam = 'y0'
        ynam = 'z0'
      else if(itype.eq.13 .and. ndm.eq.3) then
        r = x(1,n)
        x(1,n) = x0 + r*cos(x(3,n)*th)
        x(3,n) = y0 + r*sin(x(3,n)*th)
        xnam = 'x0'
        ynam = 'z0'
      end if
      if(mct.gt.0) go to 250
cww   if(prt) write(iow,2000) o,head,xnam,x0,ynam,y0,(i,i=1,ndm)
      if(prt) write(iow,2000)        xnam,x0,ynam,y0,(i,i=1,ndm)
      if(prt.and.ior.lt.0)
cww  1        write(*,2000)   o,head,xnam,x0,ynam,y0,(i,i=1,ndm)
     1        write(*,2000)          xnam,x0,ynam,y0,(i,i=1,ndm)
      mct = 50
250   if(prt) write(iow,2001) n,(x(i,n),i=1,ndm)
      if(prt.and.ior.lt.0) write(*,2001) n,(x(i,n),i=1,ndm)
      mct = mct - 1
      n = n + inc
      if((ne-n)*inc.ge.0) go to 200
      if(mod(ne-ni,inc).eq.0) go to 100
      ni = ne
      n = ne
      go to 200
c.... error
cww300write(iow,3000) ni,ne
cww   if(ior.lt.0) write(*,3000) ni,ne
cww   stop
300   write(yyy,3000) ni,ne
      call drawmess(yyy,1,0)
      return
c.... formats
cww2000  format(a1,19a4,a3//
2000  format(/
     1 '  cartesian coordinates computed from polar input with ',
     2 a2,' = ',g12.5,', ',a2,' = ',g12.5/2x,'node',6(i6,'-coord'))
2001  format(i8,6f12.4)
cww3000  format('  **fatal error 16** attempt to convert ndes ni = ',i6,
cww     1 ' to ne = ',i6)
3000  format(
     +' Try to convert nodes in POLA from ',i4,' to ',i4,'(> max node)')
3001  format(' Input: node-1,node-2,inc, x0, y0'/'   >',$)
      end
c
      subroutine prefl(id,x,ix, ndm,ndf,nen,nen1,
     1                 numnp,numel, nnp1,nel1, neq,nt,xt)
c----------------------------------------------------------------------
c.... reflect mesh about x(nt,n) = xt
c     not used in actual version

c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical error
      real*8  x(ndm,*)
      integer id(ndf,*),ix(nen1,*)
c.... set conversion parameters
      error = .false.
      xt2 = xt + xt
      nnn = numnp + nnp1 + 1
      nne = numel + nel1 + 1
c.... compute maximum equation number in current mesh
      neq = 0
      nqm = 0
      do 50 n1 = 1,nnp1
      do 50 i  = 1,ndf
        neq = max(id(i,n1),neq)
        nqm = min(id(i,n1),nqm)
50    continue
      nq2 = neq + neq  + 1
c.... loop over nodes to set coordinates and boundary codes
      do 130 n1 = 1,nnp1
        n2 = nnn - n1
c.... set up coordinates on new nodes
        do 100 i = 1,ndm
          if(i.ne.nt) then
            x(i,n2) = x(i,n1)
          else
            x(i,n2) = xt2 - x(i,n1)
          end if
100     continue
c.... assign boundary condition/equation number
        do 110 i = 1,ndf
          if    (id(i,n1).gt.0) then
            id(i,n2) = nq2 - id(i,n1)
          else if(id(i,n1).lt.0) then
            id(i,n2) = nqm + id(i,n1)
          else
            id(i,n2) = 0
          end if
110     continue
130   continue
c.... set up new element lists
      do 230 n1 = 1,nel1
        n2 = nne - n1
        ix(nen1,n2) = ix(nen1,n1)
c.... determine maximum number of nodes on this element
        do 200 i = nen,1,-1
          i2 = i + 1
          if(ix(i,n1).ne.0) go to 210
200     continue
c.... assign numbers to new element
210     do 220 i = 1,i2-1
          if(ix(i,n1).gt.numnp) then
            write(*,2000) i,n1,ix(i,n1)
            error = .true.
          else if(ix(i,n1).ne.0) then
            ix(i2-i,n2) = nnn - abs(ix(i,n1))
          else
            ix(i2-i,n2) = 0
          end if
220     continue
230   continue
      if(error) then
cww       stop
        return
      end if
      return
2000  format(' **ERROR** dof',i2,' for element',i5,' is',i6)
      end
c
      subroutine proced (name,v,wd,nwd,ll,jct,lct,ct,flg,iopl1)
c----------------------------------------------------------------------
c
c      Purpose: Define or execute a procedure for command language
c               solution
c
c      Inputs:
c         name      - Procedure name (1-8 characters)
c         vv        - Character form of procedure parameters
c         wd(*)     - List of possible commands
c         nwd       - Number of possible commands in wd
c         ll        - Number of procedure command
c         ct(3,*)   - Parameters for procedure execution
c         flg       - Flag, true if procedure already exists
c         iopl1     - Switch: Define procedure if = 1; else execute.
c
c      Outputs:
c         jct(*)    - List of commands to execute
c         lct(*)    - List of command options
c
c.... read,write and define a procedure
c----------------------------------------------------------------------
      USE feapprog
      USE iodata
      USE iofile
      USE pdata2
      USE proc
      implicit double precision (a-h,o-z)
      character name(15)*1,v(3)*15,wd(nwd)*4,lct(*)*4,yyy*80
      character fnam*229,x1*75,macc*229
      character*1 xx(75),yy(75),name1*4
      logical pcomp,flg
      dimension jct(*),ct(3,*),va(3)

c.... case: proc,--  1,2  return   3,4 go into editor
      if(idev.eq.1 .and. name(1).eq.' ') return  
      if(idev.eq.2 .and. name(1).eq.' ') return  

c.... set path
      do i = 1,4
        name1(i:i) = name(i)
      end do
      if(name1(1:4).eq.'path') then
        call readstr('Set new path for PCD-Files:',procpath)
        iposl = ipos(procpath,229)
        if(iposl.eq.0)return
        procpath(iposl+1:iposl+1) = '\'
        return
      end if
c.... set file name to store procedure
      do j = 1,229
        macc(j:j) = ' ' 
        fnam(j:j) = ' ' 
      end do
      iposl = ipos(procpath,229)
      do j = 1,iposl
        macc(j:j) = procpath(j:j)
        fnam(j:j) = procpath(j:j)
      end do 
      do i = 1,8
          if(name(i).eq.' ') go to 120
      end do
      i = 9
120   name(i  ) = '.'
      name(i+1) = 'p'
      name(i+2) = 'c'
      name(i+3) = 'd'
c.... move name to fnam/macc to open file
      do j = 1,i-1
        i1 = iposl+j
        macc(i1:i1) = name(j)
      end do   
      do j = 1,i+3
        i1 = iposl+j
        fnam(i1:i1) = name(j)
      end do   
c.... check if file exists and open file
      if((idev.eq.3.or.idev.eq.4).and.iopl1.eq.1) goto 131
      call opnfil(macc,fnam,iopl1,flg,ir)
      if(ir.eq.1) return  
c.... write a new procedure
131   if(iopl1.eq.1) then
        if(idev.eq.3.or.idev.eq.4) then
          call efeap
        else
c....     write variable names onto file
cww       write(ios,2000) v
cww>>     write with commas
          yyy = ' '
          yyy( 1:15) = v(1)
          yyy(16:30) = v(2)
          yyy(31:45) = v(3)
          call trans2(yyy)
          write(ios,2001) yyy
cww<<
c....     input procedure commands and add to file
          if(ior.lt.0) write(*,3000)
          n = 0
140       n = n + 1
          if(ior.lt.0) write(*,3001) n
          call pintio(yyy,15)
          if(pcomp(yyy,'end ',4)) go to 170
          do 150 i = 1,nwd
            if(pcomp(yyy,wd(i),4)) go to 160
150       continue
          if(ior.lt.0) then
            write(*,3002) yyy
          else
            write(iow,3002) yyy
cww         stop
            return
          end if
          n = n - 1
          go to 140
c....     modify yyy
160       call trans2(yyy)
          write(ios,2001) yyy
          go to 140
170       close(ios)
        end if
      else if(flg) then
c....   read an existing procedure
        va(1) = ct(1,ll)
        va(2) = ct(2,ll)
        va(3) = ct(3,ll)
        ll = ll - 1
c       procedures  without formats
        read(ios,'(75a)') xx        ! first line in procedure file
        call acheck(xx,yy,15,75,75)
        do i=1,45
          yyy(i:i)=yy(i)
        enddo
        read(yyy,'(3a15)') v
200     read(ios,'(a)',end=300) x1  ! read  line in procedure file
c....   remove leading blanks
        do  i=1,75
          if( x1(i:i).ne.' ') then
            ii=1
            do j=i,75
              xx(ii) = x1(j:j)
              ii=ii+1
            end do 
            goto 210
          end if
        end do
210     call acheck(xx,yy,15,75,75)
        yyy=' '
        do i=1,75
          yyy(i:i)=yy(i)
        end do
        if(.not.pcomp(v(1),'    ',4)) call setpcd(yyy,v(1),va(1))
        if(.not.pcomp(v(2),'    ',4)) call setpcd(yyy,v(2),va(2))
        if(.not.pcomp(v(3),'    ',4)) call setpcd(yyy,v(3),va(3))
        call setmac(yyy,wd,nwd,ll,jct,lct,ct)
        go to 200
300     close(ios)
      end if
      return
2001  format(a80)
3000  format('  Procedure Definition - Terminate with "end". ')
3001  format('  Procedure Statement',i3,' >',$)
3002  format('  **WARNING** Illegal command ',a4,' in procedure.')
      end
c
      subroutine trans2 (yyy)
c----------------------------------------------------------------------
c
c       Purpose: remove blanks and add commas  of string yyy
c
c       Input:
c          yyy  - original string with blanks
c
c       Output:
c          yyy  - modified string with commas
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character*1 yyy(80),x(80),xx*15
      character*15 yy(5)
      do i = 1,80
         x(i) = ' '
      enddo
c.... copy values
      do  i = 1,5
        do k = 1,15
          kk = (i-1)*15 + k
          xx(k:k) = yyy(kk)
        enddo
        yy(i) = xx
      enddo
c.... number of terms 
      do i = 5,1,-1
          if(yy(i).ne.'               ') goto 10
      enddo
c.... remove blanks and add commas
10    kk = 1
      do 20 j = 1,i
        if(yy(j).eq.'              0') yy(j) = ' ' ! empty = 0
        xx = yy(j)
        do 30 k = 1,15
          if(xx(k:k).eq.' ') goto 30
          x(kk) = xx(k:k)
          kk = kk + 1
30      continue
        if(j.lt.i) then    ! add comma, not for last term
          x(kk) = ','
          kk = kk + 1
        end if
20    continue
c.... copy back
      do i = 1,80
        yyy(i) = x(i)
      enddo
      return
      end
c
      subroutine opnfil(macc,name,iopld,exst,ir)
c----------------------------------------------------------------------
c
c      Purpose: Open file ios for FEAP I/O operations
c
c      Inputs:
c         macc  - Name for report
c         name  - Name of file to open
c         iopld - Indicator on type of file to open or report
c         ir    - 
c
c      Outputs:
c         exst - Flag, true if file exists
c
c----------------------------------------------------------------------
      USE iodata
      USE pdata2
      implicit double precision (a-h,o-z)
c.... iopl in commom = idummy
      character*229 macc,name
      character*280 yyy   !>229+45 from text
      character y*1
      logical exst
      ir = 0
      nc = ipos(macc,229)
c.... test if name exist
      inquire(file=name,exist=exst)
      if(exst) then
          if(iopld.eq.1) then
c....     file exist, decide what to do
            write(*,2000) macc(1:nc)
10          read (*,1000,err=11,end=12) y
            goto        13
11          call        errclr ('OPNFIL')
            goto        10
12          call        endclr ('OPNFIL',y)
13          if(y.ne.'y'.and.y.ne.'Y') then
            ir = 1
            return
          end if
          end if
          open(ios,file=name,status='old')
      else
          if(iopld.eq.2) then
            write(yyy,2001) macc(1:nc)
          call drawmess(yyy,1,0)
            return
          else if(iopld.eq.-2) then
            write(yyy,2002) macc(1:nc)
          call drawmess(yyy,1,0)
            return
          end if
        open(ios,file=name,status='new')
      end if
      return
1000  format(a1)
2000  format(' A procedure named ',a,' exist. Continue?  (y/n) >',$)
2001  format(' No procedure or macro command named ',a,' exists.')
2002  format(' No file named ',a,' exists.                      ')
      end
c
      subroutine setpcd(yy,v,va)
c----------------------------------------------------------------------
c
c      Purpose: Put string into string in widths of 15
c
c      Inputs:
c         yy(*)     - String of input data
c         v         - Character string for compare
c         va        - String to insert for character string
c
c      Outputs:
c         yyy(*)    - String after substitution
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character*(*) v,yy
      logical pcomp
      if(pcomp(yy(31:34),v,4)) call readv1(yy(31:45),va)
      if(pcomp(yy(46:49),v,4)) call readv1(yy(46:60),va)
      if(pcomp(yy(61:64),v,4)) call readv1(yy(61:75),va)
      return
      end
c
      subroutine readv1(yyy,va)
c----------------------------------------------------------------------
c
c      Purpose: Read string into string 
c
c      Inputs:
c         yyy(*)    - String of input data
c         va        - String to insert 
c
c      Outputs:
c         yyy(*)    - String after substitution
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character*15 yyy
      write(yyy,2000) va
      return
2000  format(e15.8)
      end
c
      subroutine setmac(yyy,wd,nwd,ll,jct,lct,ct)
c----------------------------------------------------------------------
c
c      Purpose: Set macro commands for procedures use
c
c      Inputs:
c         yyy      - String with command name
c         wd(nwd)  - Active command names
c         nwd      - Number of active commands

c      Outputs:
c         ll       - Execution command numbers
c         jct(*)   - List of command numbers
c         lct(*)   - String for command option
c         ct(3)    - command parameters
c----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      logical pcomp
      character*4 clab1,clab2,wd(nwd),lct(*)
      character*(*) yyy
      dimension jct(*),ct(3,*)
      ll = ll + 1
      read(yyy,1000) clab1,clab2,(ct(i,ll),i=1,3)
      do 100 i = 1,nwd
        if(pcomp(clab1,wd(i),4)) go to 110
100   continue
        if(ior.lt.0) write(*,3000) clab1
      ll = ll - 1
      return
110   jct(ll) = i
      lct(ll) = clab2
      return
1000  format(2(a4,11x),3f15.0)
3000  format(' **WARNING** Illegal command ',a4,' in procedure.')
      end
c
      subroutine profil (jd,idl,id,ix,ien,inn,ieb,iop,prt)
c----------------------------------------------------------------------
c
c      Purpose: Compute profile of global arrays
c
c      Inputs:
c        ix(*)  - Element nodal connection list
c        ien(*) - edge
c        inn(*) - edge
c        ieb(*) - edge
c        iop    - Switch to control operation
c                  = 1 to set up equation numbers of dof's
c                  = 2 to compute the column/row lengths and profile.
c        prt    - Flag, print solution properties if true
c
c      Scratch:
c        idl(*) - Array to hold temporary information
c
c      Outputs:
c        id(*)  - Equation numbers for degree of freedoms     (iop = 1)
c        jd(*)  - Pointer array to row/column ends of profile (iop = 2)
c
c----------------------------------------------------------------------
c
c.... type declaration for variables
      USE cdata
      USE edgdat
      USE iofile
      USE iscsr
      USE sdata
      USE slid1
      USE slid3
      USE smpak
      USE soltyp
      use plong
      logical prt,edfl
      character*55 sol(10),un*3
      integer iop
      integer mm,jd_len
      real*8  dops,stime
c.... type declaration for arrays
      integer jd(*),idl(*),id(*),ix(*),ien(*),inn(*),ieb(*)
      common /jdlen/ jd_len
      data sol /'standard profil solver',
     +         'SM sparse  matrix solver (CSR, without minimum degree)',
     +         'SM sparse  matrix solver (CSR, with minimum degree)   ',
     +         'SuperLU    solver with CSR storage',
     +         'Pardiso    solver with CSR storage (direct)',
     +         'PCG        solver with CSR storage',
     +         'PGMRES     solver with CSR storage',
     +         'PPGMRES2   solver with CSR storage',
     +         'Pardiso/ML solver with CSR storage (iterative)',
     +         'simplex optimization algorithm'/
      data un /'   '/
 
      edfl  = nde*ned.gt.0
      nde1  = nde
      if(nde1.eq.0) nde1=1 
      if(iop.eq.1) then
c....   set up the equation numbers
        call seteq(id,ien,ieb,jd,ndf,nde1,numnp, neq, edfl)

      else if(iop.eq.2) then
c....   compute column heights
        if ( istyp.eq.0) then  
c....     standard solver
          call pzeroi(jd,neq)
          call rstprf(jd,idl,id,ix,ien,inn,ieb,
     1                nde1,ndf,nen1,nen,neq,numnp,numel,edfl)
c....     modify for contact slidelines
          if(contfl) call adjchf(cl00,cl01,cl06,cl07,id,jd,ndf,nsl)
c....     compute diagonal pointers for profile
          call nwprof(jd,neq)

        else if ( istyp.eq.1 .or. istyp.eq.2) then
c....     SM solver
          call mkptr2(isymcsr)
          call mkptr_csr(ix,id,jd,nen1,ndf,isymcsr)          

        else if ( istyp.eq.3) then  
c....     SuperLU  
          call mkptr3(isymcsr)
          call mkptr_csr(ix,id,jd,nen1,ndf,isymcsr)
          if(isymcsr.eq.2) un = 'un-'  

        else if ( istyp.eq.4) then  
c....     Pardiso direct 
          call mkptr4(isymcsr)    
          call mkptr_csr(ix,id,jd,nen1,ndf,isymcsr)
          if(isymcsr.eq.2) un = 'un-'  

        else if ( istyp.eq.5 ) then  
c....     PBCG  
          call mkptr5(isymcsr)
          call mkptr_csr(ix,id,jd,nen1,ndf,isymcsr)          
          if(isymcsr.eq.2) un = 'un-'

        else if ( istyp.eq.6 ) then  
c....     PGMRES  
          call mkptr6(isymcsr)
          call mkptr_csr(ix,id,jd,nen1,ndf,isymcsr)
          if(isymcsr.eq.2) un = 'un-'

        else if ( istyp.eq.7) then  
c....     PGMRES2  
          call mkptr7(jd,isymcsr)          
          call mkptr_csr(ix,id,jd,nen1,ndf,isymcsr)
          if(isymcsr.eq.2) un = 'un-'

        else if ( istyp.eq.8) then  
c....     Pardiso/ML iterative 
          call mkptr8(isymcsr)    
          call mkptr_csr(ix,id,jd,nen1,ndf,isymcsr)
          if(isymcsr.eq.2) un = 'un-'  

        else if ( istyp.eq.15) then
c....     simplex optimization
          continue

        end if
c....   estimate solve time stime (dops = est. no. ops/sec)
c....   estimate col.height mm  
        if(neq.gt.0) then
          jd_len = 0
          if(istyp.eq.0) then
             jd_len = jd(neq)+neq ! Standard
          else if(istyp.ge.1.and.istyp.le.9) then 
            jd_len = jd(neq+1)-1!CSR=SM/PCG/PGMRES/SuperLU/Pardiso
          end if  
          mm    = (jd_len+neq)/neq
cww       dops  =  160000  ! org TAYLOR estimation 
          dops  = 5000000  ! new estimation 
          stime = neq*mm
          stime = stime*mm/dops
c...      print results
          write(iow,2001) ndm,ndf,neq,nde,mm,numnp,jd_len,numel,stime,
     1                    nummat
          write(iow,2002) un,'symmetric',sol(istyp+1)

          if(prt) then
            write(*,2001) ndm,ndf,neq,nde,mm,numnp,jd_len,numel,
     1                    stime,nummat
            write(*,2002) un,'symmetric',sol(istyp+1)
          end if
        end if
      end if
c
2001  format(/4x,'E q u a t i o n / P r o b l e m   S u m m a r y:'/
     1 5x,'Space dimension (ndm) =',i9,   4x,'Number dof (ndf) =',i7/
     2 5x,'Number of equations   =',i9,   4x,'Number dof (nde) =',i7/
     3 5x,'Average col. height   =',i9,   4x,'Number nodes     =',i7/
     4 5x,'No. terms in profile  =',i9,   4x,'Number elements  =',i7/
     5 5x,'Est. factor SOLV-sec  =',e11.4,2x,'Number materials =',i7)
2002  format(5x,'using : ',a4,a9,1x,a55)
      end
c
c-----------------------------------------------------------------------
c
      subroutine profupd(prt)
c----------------------------------------------------------------------
c
c      Purpose: update of profile after system modifications
c               or change of solver 
c
c      Inputs:
c       n2    ld(nst)         - eq. numbers of dofs of actual element
c       n7    id(ndf,numnp)   - equation numbers for each active dof
c       n8    x(ndm,numnp)    - nodal coordinates of mesh
c       n9    ix(nen1,numel)  - element nodal connections of mesh
c       n12   jd(*)           - pointer array for row/columns of tangent
c       ne1   ien(*)          - edge
c       ne2   inn(*)          - edge
c       ne3   ieb(*)          - edge
c
c      Outputs:
c       n7    id(ndf,numnp)   - equation numbers for each active dof
c       n12   jd(*)           - pointer array for row/columns of tangent
c                             - and csrja, csrka 
c       nli1  iplk1(numnp)    - List of linked nodes to plot
c       nli2  iplk2(ndf,numnp)- List of bc for linked nodes to plot
c
c      Comments:
c       up to now only for solver 0+12   
c       LINK should be ok - check only for simple test   

c      Open:
c      input ctis for solver is missing.
c      tie =??
c      bei großen Problemen sehr langsam!! also möglichst vermeiden.
c      Idee Einbau: Teile in PCONTR als Unterprogramm, Daten dann heraus-
c      schreiben und hier wieder Einlesen-analog LINK      
c      wenn numnp, numel sich ändert, werden alle pseta-Aufrufe falsch!
c      Lösung A: dynamische Speicherverwaltung
c      Lösung B: alle Felder größer anlegen; numnp_max,numel_max 
c      Weder A noch B ist realisiert. Daher geht nur: Reduktion der RB
c 
c----------------------------------------------------------------------
c
      USE cdata
      USE edgdat
      USE mdata
      USE nolink
      USE psize
      USE sdata
      implicit double precision (a-h,o-z)
      logical prt

c.... set up equation numbers of dof's new
      call profil(jdt12,eeqn,psid,econ,edge1,edge2,edge3,1,prt)

c.... modify with respect to link ( no new definitions allowed)
      if(flnk) then
        call pzeroi (link1,numnp)   
        call pconsi2(link2,numnp*ndf,1)   
        call plink(eeqn,psid,coor,ndm,ndf,numnp,neq,
     +                     link1,link2,nlk,prt)
      end if 

c.... set profile new
      call profil(jdt12,eeqn,psid,econ,edge1,edge2,edge3,2,prt)
 
      return
      end
c
      subroutine seteq(id,ien,ieb,jd,ndf,nde,numnp, neq, edfl)
c----------------------------------------------------------------------
c.... set the equation numbers for all types of points
c----------------------------------------------------------------------
c.... type declaration for variables
      USE mxasz
      logical edfl
      integer ndf,nde,numnp,nad,neq,n,nn
      integer i,ii,i1,i2, j
c..... type declaration for arrays
      integer id(ndf,numnp),ien(numnp),ieb(nde,*),jd(ndf,numnp)
      neq = 0
      nad = 0
      do 40 n = 1,numnp
        if(nren.eq.0) then
          nn = n
        else
          nn = optin(n)
        end if
        if(edfl) then
          i1 = ien(nn) + 1
          if(nn.lt.numnp) i2 = ien(nn+1)
          if(i1.le.i2) then
            do 20 ii = i1,i2
            do 10 i = 1,nde
             j = ieb(i,ii)         
             if(j.eq.0) then
               neq = neq + 1
               ieb(i,ii) = neq 
             else
               nad  = nad - 1
               ieb(i,ii) = nad
             end if
10            continue
20          continue
          end if
        end if
        if(nn.eq.0) go to 40 
        do 30 i = 1,ndf
          j = id(i,nn)
          if(j.eq.0) then
            neq     = neq + 1
            id(i,nn) = neq
          else
            nad     = nad - 1
            id(i,nn) = nad
          end if
30      continue
40    continue
c
      end
