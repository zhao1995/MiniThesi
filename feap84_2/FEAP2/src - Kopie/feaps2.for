      subroutine dosets(inp,outp,res,sav,plt)
c----------------------------------------------------------------------
c
c      Purpose: Add initial character to default files using input
c               filename
c
c      Inputs:
c         inp  - Input   filename without inital character
c         otp  - Output  filename without inital character
c         res  - Restart filename without inital character
c         sav  - Save    filename without inital character
c         plt  - Plot    filename without inital character
c
c      Outputs:
c         inp  - Input   filename with inital character
c         otp  - Output  filename with inital character
c         res  - Restart filename with inital character
c         sav  - Save    filename with inital character
c         plt  - Plot    filename with inital character
c
c----------------------------------------------------------------------
      character*229 inp,outp,res,sav,plt
      integer ipos
      call dochar2(inp,ipos)
      if(outp.eq.' ') then
        outp = inp
        call dochar1(outp,'O',ipos)
      end if
      if(res.eq.' ') then
        res = inp
        call dochar1(res,'R',ipos)
      end if
      if(sav.eq.' ') then
        sav = inp
        call dochar1(sav,'R',ipos)
      end if
      if(plt.eq.' ') then
        plt = inp
        call dochar1(plt,'P',ipos)
      end if
      return
      end
c
      subroutine dochar1(nam,char,ipos)
c----------------------------------------------------------------------
c
c      Purpose: set first character of filename without path
c
c      Inputs:
c         nam   - filename without path and no special inital character
c         char  - inital character
c         ipos  - position of first character of filename without path
c
c      Outputs:
c         nam   - filename without path with special inital character
c
c----------------------------------------------------------------------
      character*1 char,nam*229
      integer*4 ipos
c.... set first character
      nam(ipos:ipos) = char
      return
      end
c
      subroutine dochar2(newf,ipos)
c----------------------------------------------------------------------
c
c      Purpose: find first character of filename without path, 
c               scan from back for '\'
c
c      Inputs:
c         newf  - filename 
c
c      Outputs:
c         ipos  - position of first character of filename without path
c
c----------------------------------------------------------------------
      character   newf*229
      character   newf1(229)*1
      integer*4   ipos
c.... copy name
      do i = 1,229
          newf1(i) = newf(i:i)
      end do
      newf= ' '
c.... find filename
      do i=229,1,-1
#ifdef _LINUX_
      if (newf1(i).eq.'/') go to 20
#else
      if (newf1(i).eq.'\') go to 20
#endif
      end do
      ipos = 0
      goto 30
20    ipos = i
c.... path
      do i = 1,ipos
          newf(i:i) = newf1(i)
      end do
c.... name
30    do i = ipos+1,229
          newf(i:i) = newf1(i)
      end do

      ipos = ipos+1

      return
      end
c
c
      integer function cvtc2i(cdev)
c----------------------------------------------------------------------
c      Purpose: Convert a character to an integer
c               ASCII character codes assumed

c      Inputs:
c         cdev    - Character array to convert

c      Outputs:
c         cvtc2i  - Integer value of character
c
c----------------------------------------------------------------------
      character*1 cdev(12)
      nz= ichar('0')
      n = 0
      do i = 1,12
        if(cdev(i).eq.char(0) .or. cdev(i).eq.' ') go to 200
        n = 10*n + (ichar(cdev(i)) - nz)
      end do ! i
200   cvtc2i = n

      end
c
      function cvti2c(iint)
c----------------------------------------------------------------------
c      Purpose: Convert an integer to a character
c               ASCII character codes assumed

c      Inputs:
c         cdev    - Integer to convert

c      Outputs:
c         cvti2c  -  value of character
c
c      J.Hebel TUD 12/10 
c----------------------------------------------------------------------
      implicit none
      character(len=11) :: cvti2c  ! String, Ausgabe
      integer, intent(in) :: iint  ! Integer, Eingabe
      character(len=11) :: buff    ! interne Datei
      character(len=5)  :: fmt     ! Ausgabeformat

      if (iint.GT.0) then
        fmt="(I1)"
        if (iint/10.GE.1) fmt="(I2)"
        if (iint/100.GE.1) fmt="(I3)"
        if (iint/1000.GE.1) fmt="(I4)"
        if (iint/10000.GE.1) fmt="(I5)"
        if (iint/100000.GE.1) fmt="(I6)"
        if (iint/1000000.GE.1) fmt="(I7)"
        if (iint/10000000.GE.1) fmt="(I8)"
        if (iint/100000000.GE.1) fmt="(I9)"
        if (iint/1000000000.GE.1) fmt="(I10)"
      else if (iint.LT.0) then
        fmt="(I2)"
        if (-iint/10.GE.1) fmt="(I3)"
        if (-iint/100.GE.1) fmt="(I4)"
        if (-iint/1000.GE.1) fmt="(I5)"
        if (-iint/10000.GE.1) fmt="(I6)"
        if (-iint/100000.GE.1) fmt="(I7)"
        if (-iint/1000000.GE.1) fmt="(I8)"
        if (-iint/10000000.GE.1) fmt="(I9)"
        if (-iint/100000000.GE.1) fmt="(I10)"
        if (-iint/1000000000.GE.1) fmt="(I11)"
      end if

      write(buff,fmt) iint
      cvti2c = buff

      end
c
      double precision function dot(a,b,n)
c----------------------------------------------------------------------
c      Purpose: dot (scalar) product of two vectors
c
c      Inputs:
c         a(*)  - Vector 1
c         b(*)  - Vector 2
c         n     - length of vectors
c
c      Outputs:
c         dot   - Scalar product
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension a(*),b(*)
      dot = 0.0
      do 100 i = 1,n
        dot = dot + a(i)*b(i)
100   continue
      return
      end
c
c      double precision function ddot(n,a,ia,b,ib)
c----------------------------------------------------------------------
c      INTEL:   use MKL LIB
c      Salford: see plot_sal_cw.for
c
      double precision function dotid(a,b,id,n)
c----------------------------------------------------------------------
c
c      Purpose: dot (scalar) product of two vectors using id array.
c
c      Inputs:
c         a(*)  - Vector 1
c         b(*)  - Vector 2 to be accessed using id array
c         id(*) - Equation pointer array
c         nn    - length of vectors
c
c      Outputs:
c         dotid - Scalar product
c
c      Comment: Typical use:  vector(neq) dot vector(nneq)
c
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      integer id(*)
      real*8  a(*),b(*)
      temp = 0.d0
      do 100 i = 1,n
        j = id(i)
        if(j.gt.0) temp = temp + a(i)*b(j)
  100 continue
      dotid = temp
      return
      end
c
      subroutine dpacku ( u, id, nneq )
c----------------------------------------------------------------------
c
c      Purpose: Expand solution vector to fill entries
c
c      Inputs:
c         u(*)  - Vector
c         id(*) - Equation pointer array
c         nneq  - length of vectors
c
c      Outputs:
c         u(*)  - Vector
c
c----------------------------------------------------------------------
      integer i, j, nneq, id(nneq)
      real*8  u(nneq)
c
      do 275 i = nneq,1,-1
        j = id(i)
        if(j.gt.0) then
          u(i) = u(j)
          u(j) = 0.0d0
        end if
275   continue
c
      end
c
      subroutine dredu(al,au,ad,jh,flg,dj)
c----------------------------------------------------------------------
c
c      Purpose: Reduce diagonal in symmetric/unsymmetric triangular 
c               decomposition
c
c      Inputs:
c         al(*)  - Lower terms in row
c         au(*)  - Upper terms in column
c         ad(*)  - Reduced diagonals of previous equations
c         jh     - Length of row/column
c         flg    - .false.=sym  .true.=unsym
c
c      Outputs:
c         dj     - reduced diagonal for current equation

c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical flg
      dimension al(jh),au(jh),ad(jh)
      if(.not.flg) then
c.... computation for symmetric matrices
        do 100 j = 1,jh
          ud    = au(j)*ad(j)
          dj    = dj - au(j)*ud
          au(j) = ud
100     continue
c.... computation for unsymmetric matrices
      else
        do 200 j = 1,jh
          ud    = au(j)*ad(j)
          dj    = dj - al(j)*ud
          au(j) = ud
          al(j) = al(j)*ad(j)
200     continue
      end if
      return
      end
c
      subroutine edgupd(du,ien,ieb,f,u,nde,neq,numnp)
c----------------------------------------------------------------------
c
c      Purpose: update the edge degree-of-freedoms
c
c----------------------------------------------------------------------
c..... type declaration for variables
      USE ddata
      USE prlod
      USE tdata
      integer nde,neq,numnp
      integer i,ii,i1,i2, j, n,nl2,nl3
      real*8  cc1, ub
c..... type declaration for arrays
      integer ien(*), ieb(nde,*)
      real*8  du(*), f(nde,*), u(nde,*)
c
      nl2 = ien(numnp)
      nl3 = nl2+nl2
      cc1 = 1.0d0
cww      if (nop.eq.3)  HHT
cww      if (nop.eq.4) cc1 = c1 !Generalized-alpha ??
      do 200 n = 1,numnp-1
        i1 = ien(n) + 1
        i2 = ien(n+1)
        do 150 ii = i1,i2
          do 140 i = 1,nde
            j = ieb(i,ii)
            if (j.gt.0) then
c.... for the active degrees-of-freedom compute values from solution
c.... where 'du(j)' is and increment of 'u' for active dof 'j'.
              u(i,ii)     = u(i,ii)     + cc1*du(j)
              u(i,ii+nl2) = u(i,ii+nl2) + cc1*du(j)
              u(i,ii+nl3) = cc1*du(j)
            else
c.... for fixed degrees-of-freedom compute values from forced inputs
              ub          = f(i,ii)*prop
              du(neq-j)   = ub - u(i,ii)
              u(i,ii+nl3) = ub - u(i,ii)
              u(i,ii+nl2) = u(i,ii+nl2) + ub - u(i,ii)
              u(i,ii)     = ub
            end if
140       continue
150     continue
200   continue
      end
c
c------------------------------------------------------------------------
c
      subroutine eigi (mnal,mnau,mna,jd,neq,phi,x1,nstep,tol,evi)  
c------------------------------------------------------------------------
c
c      Purpose: compute the lowest eigenvalue evi and eigenvector phi
c               using inverse iteration 
C               coded according to habil. Wagner pp. (knebel jan. 94)
c
c      Inputs:
c        mnal(*)   - factored lower triangular terms
c        mnau(*)   - factored upper triangular terms
c        mna (*)   - factored diagonal terms
c         jd(*)    - Pointer to row/column ends in 'al' and 'au'.
c        neq       - Number of unknowns
C        nstep     -  max. number of iter. steps ( default=50 )
C        tol       -  used iter. tolerance       ( default=1.d-05 )
c
c      Scratch:
C        x1(neq)   -  temporary used vector = dr
c 
c      Outputs:
C        phi(neq)  -  eigenvector
C        evi       -  eigenvalue
c
c      Comment:
c        use [tang,,1] before calling this subroutine
c
c------------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension phi(*),x1(*),mnal(*),mnau(*),mna(*),jd(*)
C
      istep = 0
      evold = 1.d+10
      do  i=1,neq                     ! set starting vector to 1.0
c        phi(i) = 1.0d0 ! problems in finding non symmetric evs
c        x1(i)  = 1.0d0
        call random_number(rnum)
        phi(i) = rnum        
        call random_number(rnum)
        x1(i)  = rnum

      end do
100   istep = istep + 1
      call dasol (mnal,mnau,mna,x1,jd,neq,dummy)   ! solve for x1
c
      q1 = ddot(neq,x1,1,phi,1)        ! compute rayleigh-quotient rq
      q2 = ddot(neq,x1,1,x1,1)
      rq = q1/q2
c
      q2 = sqrt(q2)            
      do i=1,neq                  ! scale x1
        x1(i) = x1(i)/q2
      end do
c
      evi = rq                    ! set eigenvalue  'evi'  
      call pmove(x1,phi,neq)      ! set eigenvector 'phi'
c
      tole = evold-evi            ! tolerance 
      if( tole .lt. tol) then     ! check convergence
        if(ior.lt.0) write(*  ,1000) istep,evi,tole
                     write(iow,1000) istep,evi,tole
        call scalev(phi,neq)  ! scale to maximum element of 1.0
        return
      else
        evold = evi
        if(istep.le.nstep) goto 100
        if(ior.lt.0) write(*  ,2000) nstep,evi,tole
                     write(iow,2000) nstep,evi,tole
        call scalev(phi,neq)  ! scale to maximum element of 1.0
        return
      end if
c
1000  format(/'Inverse Iteration: ',i2,' iterations',/,
     +     3x,'eigenvalue         ',g12.5,' tol = ',g12.5,/)
2000  format(/'no convergence for eigenvalue after',i2,' iterations',/,
     +     3x,'last eigenvalue was ',g12.5,' tol = ',g12.5,/)

      return
      end
c
      subroutine elmterr(ix,xl,ndm,numel,e_omex,e_ome)
c---------------------------------------------------------------------+
c
c      Purpose: calculate and plot/print error distribution on element 
c               level
c
c      Inputs:
c         ix(nen1,*)    - Element nodal connections of mesh
c         xl(ndm,*)     - Nodal coordinates of element
c         ndm           - Spatial dimension of mesh
c         numel         - Number of elements in mesh
c
c      Outputs:
c        e_omex(numerr)      - element errors for print/plot
c        e_ome(numel,numerr) - element errors for remesh 
c
c      Comments: 
c        e_name(numerr)      - name of error
c        e_omex(numerr)      - element errors for print/plot
c        e_ome(numel,numerr) - element errors for remesh
c        numerr=2
c
c        iet(1)= 1    plot  errors     
c        iet(1)= 2    print errors     
c        iet(1)= 3    remesh     
c
c        iet(2)= 1-3  type of error norm
c
c---------------------------------------------------------------------+
      USE eldata
      USE errin1
      USE errin2
      USE errnam
      USE fdata
      USE iofile
      implicit double precision (a-h,o-z)
      dimension ix(*),xl(ndm,*), e_omex(*),e_ome(numel,*) 
c.... domain error = sum of element errors
      do i=1,numerr
        e_om(i) = e_om(i) + e_omex(i) 
      end do

      if (iet(1).eq.3) then   ! reme(sh)
c....   create error field for refinement action
c....   already set into relation of the global energy 'e_bar'
        do i=1,numerr
          if (abs(e_bar(i)).gt. 1e-10)
     +      e_ome(n,i) = sqrt(e_omex(i))/e_bar(i) 
        end do
      else ! print/plot
c....   output of errors
        do i=1,numerr
          e_omex(i) = sqrt(e_omex(i))/e_bar(i) 
        end do
      end if
      if (iet(1).eq.2) then     
c....   output on  'ofile'/output on  screen:  macro>erro
        if(pfr) then
          if (mct .eq. 0) then
            ival = eval
                            write (iow,2000) ival,(e_name(i),i=1,numerr)
            if (ior .lt. 0) write (  *,2000) ival,(e_name(i),i=1,numerr)
            mct = 50
          end if
          mct = mct - 1
                          write(iow,2001) n,(e_omex(i),i=1,numerr)
          if (ior .lt. 0) write(*  ,2001) n,(e_omex(i),i=1,numerr)
        end if
c....   output on file:  macro>erro,save
        if(ioerr.eq.1) write(22,2001) n,(e_omex(i),i=1,numerr)
c....   output into array e_ome, Konflikt zu oben??!! ww
        do i=1,numerr
          e_ome(n,i) = e_omex(i)  
        end do

      else if (iet(1) .eq. 1) then                                          
c.... output on  screen:  plot>erro
        icol = 1 + e_omex(iet(2))
        icol = min(14,max(1,icol))
        icol = 31 - icol
        call pppcol (icol)
        call plot9 (iel,ix,xl,ndm,nel,1)
      end if
2000  format(/'    M e s h   R e f i n e m e n t s   F o r',i4,'%',
     +        '    E r r o r'//'    elmt',5(1x,a15) /)
2001  format(i8,1p5e12.4)
      end
c
      subroutine eload(x,f,id,dist,idis,ndm,ndf,numnp,prt)
c-----------------------------------------------------------------------
c
c      Purpose: calculate loads along a line
c
c      Inputs:
c         x(ndm,*)    - Nodal coordinates of mesh
c         f(ndf,*)    - Nodal force of mesh
c         ix(nen1,*)  - Element nodal connections of mesh
c         ndm         - Spatial dimension of mesh
c         ndf         - Number dof/node
c         numnp       - Number of nodes in mesh
c         prt         - Print flag

c
c      Outputs:
c         f(ndf,*)    - Nodal force of mesh
c
c      Scratch
c         dist(numnp) - distance of node to ref. point 
c         idis(numnp) - list of nodes
c
c      Comments: 
c       1=line(P_3=0),2=parabola(P_3>0,rad=0),3=circle(P_3>0,rad>0)
c       input:                                                     
c       iload,P_1(x,y,[z]),P_2(x,y,[z]),[P_3(x,y,[z])],tol,ibou,rad
c       iload: 1 (linear) or 2 (quadratic) elements                
c       ibou : 0 add if no boun  1 set for boun                    
c                                                                  
c       Load case 1: const. loads in global 1-3 direction          
c       Load case 2: linear loads in global 1-3 direction          
c                                                                  
c       loads act in global direction on the element projection     
c       defined by 0: ds, 1=dx, 2= dy, 3=dz                         
c                                                                  
c                                                                  
c       not documented:                                            
c       load case 3: constant normal loads                          
c       load case 4: constant normal loads  axisymmetric           
c       Load case 5: boundary terms for warping function           
c       Load case 6: boundary terms for stress  function
c                                                       
c-----------------------------------------------------------------------
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      character*11 cload
      character*8  ltyp(3)
      character*20 ltyp1(3)
      dimension x(ndm,numnp),f(ndf,numnp),td(16),dist(*),idis(*),
     +          dl(3),q0(6),qa(3),qb(3),dq(3),q1(3),q2(3),ils(6),
     +          id(ndf,*),td1(14),ddl(3),dl1(3),dl2(3),dlsf(6),ddlsf(3),
     +          dsf1(3),dsf2(3)
      dimension dp(3),p(3),p1(3),p2(3),p3(3),a(3),dx(3),dn(3),dh1(3),
     +          dh2(3),dh3(3)
      dimension shp(2,5),pg(5),wg(5),xl(2,3)
      dimension fm(ndf),fma(ndf)
      data blank/-999.d0/,tol/1000.d0/ q0(3) /0.d0/
      data ltyp/'line    ','parabola','circle  '/
      data ltyp1/' intermediate point ',' intermediate point ',
     +           ' origin of circle   '/
c.... read input of line data  - limit is 16 nos. / record
100   if(ior.lt.0) write(*,3001)
c.... input coordinates
      call dinput(td,3+ndm*3)
      if(errck) go to 100
c.... all coordinates 0
      xii = dsqrt(ddot(ndm*3,td,1,td,1))
      if(xii.le.0) return
c...  clear all arrays
      call pzero(dist,numnp)
      call pzeroi(idis,numnp)
      call pzero(dp,3)
      call pzero(p ,3)
      call pzero(p1,3)
      call pzero(p2,3)
      call pzero(p3,3)
      call pzero(a ,3)
      call pzero(dx,3)
      call pzero(dn,3)
      call pzero(dh1,3)
      call pzero(dh2,3)
      call pzero(dh3,3)
c.... middle node only within 0.25 < P_3 < 0.75
c.... to do only for parabola!
c
c.... points of line
      toluser= td(1+ndm*3)
      if(dabs(toluser).gt.1.e-10) tol = toluser
      ibou   = td(2+ndm*3)
      radius = td(3+ndm*3)
      do i = 1,ndm
        p1(i) = td(i      )
        p2(i) = td(i+  ndm)
        p3(i) = td(i+2*ndm)
      end do
c.... calculate distance tol, vector a = p2-p1, aa = a*a
      do i = 1,ndm
        dx(i) = pdiff(x(i,1),ndm,numnp)
        a(i)  = p2(i) - p1(i)
      end do
c.... circle tolerance
      tolc = dsqrt(dot(dx,dx,3))/tol
      aa   = dot(a,a,3)
c.....find type of line
      dtyp = dot(p3,p3,3)
      if(radius.ne.0.d0) then
c....   circle
        ityp = 3
cww        if(ndm.eq.3) stop 'eload for circle not in 3d in SR elaod'
c....   calculate boundary point 1 of circle line
        dxp1 = p1(1)-p3(1)
        dyp1 = p1(2)-p3(2)
        fr1 = dabs(dxp1*dxp1 + dyp1*dyp1 - radius**2)
        fr1 = dsqrt(fr1)
        if(fr1.ge.tolc) then
          write(iow,2007)
          return
        end if
        phi1 = datan2(dyp1,dxp1)
        if(phi1.lt.0.d0) phi1 = phi1+pi2
c....   calculate boundary point 2 of circle line
        dxp2 = p2(1)-p3(1)
        dyp2 = p2(2)-p3(2)
        fr2 = dabs(dxp2*dxp2 + dyp2*dyp2 - radius**2)
        fr2 = dsqrt(fr2)
        if(fr2.ge.tolc) then
          write(iow,2007)
          return
        end if
        phi2 = datan2(dyp2,dxp2)
        if(phi2.le.0.d0) phi2 = phi2+pi2
      else if(dtyp.gt.0.0d0.and.radius.eq.0.d0) then
c....   parabola
        ityp = 2
c....   define parabola equation
c....   vector p3-p1
        do i = 1,3  
          dp(i) = p3(i) - p1(i)
        end do
        dlamb = dot(dp,a,3)/aa
        if(dlamb.gt.1.d0) then
          call drawmess('error in defining parabola',1,0)
          return
        end if
c....   find vector d=dp normal to a = p2-p1
        do i = 1,3  
          dp(i) = p3(i) - (p1(i) + dlamb*a(i))
        end do
        ddp = dsqrt(dot(dp,dp,3))
c....   set abbreviation  dn
        do i = 1,3  
          dn(i) = dp(i)/(dlamb*dlamb-dlamb)
        end do
c....   -> parabola x=p1+dlamb*a+dlamb*(dlamb-1)*dn
c....   test
        if (ddp.lt.1.e-5) then
c....     point is on line
          ityp = 1
        end if
      else
        ityp = 1
      end if
c.... print base data
      if(prt)              write(iow,2002) ltyp (ityp),p1,p2,
     +                                     ltyp1(ityp),p3
      if(prt.and.ior.lt.0) write(*  ,2002) ltyp (ityp),p1,p2,
     +                                     ltyp1(ityp),p3
      if(prt)              write(iow,2006) tol,ibou,radius
      if(prt.and.ior.lt.0) write(*  ,2006) tol,ibou,radius         
c.... find nodes 
      ii = 0
c.... loop over all nodes
      do 200 n = 1,numnp
        if(x(1,n).eq.blank) goto 200
        if(ityp.eq.1) then
c....     line
c....     test if node is near line p=x  dp=p-p1
          do i = 1,ndm
            p(i) = x(i,n)
            dp(i) = p(i) - p1(i)
          end do
c....     point on a = p2-p1:  x = p1 + dlamb*a
          dlamb = dot(dp,a,3)/aa
          if(dlamb.lt.-1.d0/tol.or.dlamb.gt.(1.d0+1.d0/tol)) goto 200
c.....    vector d = dp = p - x
          do i=1,3  
            dp(i) = dp(i) - dlamb*a(i)
          end do
c.....    length of d
          dd = dsqrt(dot(dp,dp,3))
        else if(ityp.eq.2) then
c....     parabola
c....     condition (p-x)*t = 0  t = x'
c....     ->cubic eq. asc*dlamb**3+bsc*dlamb**2+csc*dlamb+dsc = 0
c....     abbreviations
          do i = 1,ndm
            p(i)   = x(i,n)
            dh1(i) = x(i,n) - p1(i)
            dh2(i) = a(i)   - dn(i)
          end do
c.....    test, if node is in parabola plane
          call vecp(a(1),dn(1),dh3(1))
          hsc =  dabs(dot(dh1,dh3,3))
          if(hsc.gt.tolc) goto 200
c.....    test, if node is on parabola
          esc =       dot(dn ,dn ,3)
          fsc =       dot(a  ,a  ,3)
c.....    cubic equation
          asc = -2.d0*esc
          bsc =  3.d0*esc
          csc =  2.d0*dot(dn,dh1,3) - esc - fsc
          dsc =       dot(dh1,dh2,3)
          qsc=(bsc/(3.d0*asc))**3-bsc*csc/(6.d0*asc**2)+dsc/(2.d0*asc)
          psc=(3.0d0*asc*csc-bsc*bsc)/(9.d0*asc*asc)
          diskr = qsc*qsc + psc*psc*psc
c.....    test diskriminante
          if(diskr.le.0.d0) goto 200
          wdiskr = dsqrt(diskr)
c.....    solution (value xlamb for parabola x)
          y11 = (-qsc+wdiskr)
          if(y11.lt.0.d0) y11 = -dabs(y11)**(1.d0/3.d0)
          if(y11.gt.0.d0) y11 =      (y11)**(1.d0/3.d0)
          y12 = (-qsc-wdiskr)
          if(y12.lt.0.d0) y12 = -dabs(y12)**(1.d0/3.d0)
          if(y12.gt.0.d0) y12 =      (y12)**(1.d0/3.d0)
          y1 = y11+y12
          xlamb = y1 - bsc/(3.d0*asc)
          if(xlamb.lt.-1.d0/tol.or.xlamb.gt.1.d0+1.d0/tol) goto 200
c....     node on parabola
          do i = 1,3
            dp(i) = p1(i) + xlamb*a(i) + xlamb*(xlamb-1.d0)*dn(i)
          end do
c.....    difference vector dp = p-x(xlamb)
          dp(1) = x(1,n) - dp(1)
          dp(2) = x(2,n) - dp(2)
          if(ndm.eq.3) dp(3) = x(3,n) - dp(3)
          if(ndm.eq.2) dp(3) =        - dp(3)
c.....    length of d
          dd = dsqrt(dot(dp,dp,3))
        else if(ityp.eq.3) then
c....     node on circle
          dxpn = x(1,n)-p3(1)
          dypn = x(2,n)-p3(2)
          frn = dabs(dxpn*dxpn + dypn*dypn - radius**2)
          dd  = dsqrt(frn)
          if(dd.ge.tolc) goto 200
          phin = datan2(dypn,dxpn)
          if(phin.lt.0.d0) phin = phin+pi2
c.....    special case:  phin in 4.quadrant and eq. 360 deg. or .lt.eps
cww       if(phin.eq.0.d0.and.phi1.gt.0.d0) phin = phin+pi2
          if(phin.lt.eps .and.phi1.gt.0.d0) phin = phin+pi2
c
          if(phin.lt.phi1-eps) goto 200
          if(phin.gt.phi2+eps) goto 200
        end if
        if(dd.lt.tolc) then
c....     node is on curve
          ii = ii+1
c....     calculate distance to point P_1
          if(ityp.le.2) then
            d1   = p(1)-p1(1)
            d2   = p(2)-p1(2)
            d3   = p(3)-p1(3)
            dist(ii) = sqrt(d1*d1+d2*d2+d3*d3)
          else
            dist(ii) = (phin - phi1)*radius
          end if
          idis(ii) = n
        end if
200   continue
c.... sort nodes due to distance to point P_1
        call eload1(dist,idis,ii,tolc)
c.... set load vector
101   if(ior.lt.0) write(*,3002)
3002  format('Input:iload,icase,loads due to icase >',$)
      call dinput(td1,14)
      if(errck) go to 101
      iload = td1(1)
      if(iload.eq.1) cload='  linear   '
      if(iload.eq.2) cload=' quadratic '
      if(prt)              write(iow,2009) cload
      if(prt.and.ior.lt.0) write(*  ,2009) cload
c.....copy td1->td to have old numbers
      do i = 1,13
        td(i) = td1(i+1)
      end do
      icase = td(1)
      if(icase.eq.1) then
c....   load case 1: const. loads in global 1-3 direction 
c       on the element projection il defined by 0= ds, 1=dx, 2= dy, 3=dz 
c....   1, q0_x, il_x, q0_y, il_y, q0_z, il_z
c....      m0_x, im_x, m0_y, im_y, m0_z, im_z
         q0(1) = td(2)
        ils(1) = td(3) 
         q0(2) = td(4)
        ils(2) = td(5) 
         q0(3) = td(6)
        ils(3) = td(7) 
         q0(4) = td(8)
        ils(4) = td(9) 
         q0(5) = td(10)
        ils(5) = td(11) 
         q0(6) = td(12)
        ils(6) = td(13) 

        call eload2(ils)
c        
        if(prt)              write(iow,2003) q0,ils      
        if(prt.and.ior.lt.0) write(*  ,2003) q0,ils        
c....   loop over all elements
        if(iload.eq.1) then       ! linear
           inc = 1
           c1  = 0.5d0 
        else if(iload.eq.2) then   ! quadratic
           inc = 2
           c1  = 1.d0/6.d0 
        end if
        do i = 1,ii-inc,inc
          n1 = idis(i)
          n2 = idis(i+inc)
          if(iload.eq.2) n3 = idis(i+1)
          dl(1) = dabs( x(1,n2)-x(1,n1) )
          dl(2) = dabs( x(2,n2)-x(2,n1) )
          dl(3) = 0.d0
          if(ndm.eq.3) dl(3) = dabs( x(3,n2)-x(3,n1) )
          ndf1 = min(ndf,6)
c....     load projection length
          call eload3(ils,dlsf,dl)
c############### ilm
c.....    nodal loads
cww       do j = 1,ndm 
          do j = 1,ndf1 
            if(ibou.eq.0) then
c....         load projection 
              dls = dlsf(j)
c
c....         from input file(ior.ge.0)  id=0/1            
c....         interactive    (ior.lt.0)  id=-i/k via macro>mesh         
              if(ior.ge.0) then          
                if(id(j,n1).eq.0) f(j,n1) = f(j,n1) + c1*q0(j)*dls
                if(id(j,n2).eq.0) f(j,n2) = f(j,n2) + c1*q0(j)*dls
              else if(ior.lt.0) then          
                if(id(j,n1).gt.0) f(j,n1) = f(j,n1) + c1*q0(j)*dls
                if(id(j,n2).gt.0) f(j,n2) = f(j,n2) + c1*q0(j)*dls
              end if
              if(iload.eq.2) then
                if(ior.ge.0) then          
                  if(id(j,n3).eq.0) f(j,n3)=f(j,n3) + 4.d0*c1*q0(j)*dls
                else if(ior.lt.0) then         
                  if(id(j,n3).gt.0) f(j,n3)=f(j,n3) + 4.d0*c1*q0(j)*dls
                end if 
              end if
            else
              f(j,n1) = q0(j)
              f(j,n2) = q0(j)
              if(iload.eq.2) f(j,n3) = q0(j)
            end if
          end do
        end do
      else if(icase.eq.2) then
c....   load case 2: linear  loads in global 1-3 direction  
c       on the element projection il defined by 0: ds, 1=dx, 2= dy, 3=dz 
c       q_x = q_ax + (q_bx-q_ax)*xsi
c       q_y = q_ay + (q_by-q_ay)*xsi
c       q_z = q_az + (q_bz-q_az)*xsi
c....   2, q_xa, q_xb, il_x, q_ya, q_yb, il_y, q_za, q_zb, il_z 
         qa(1) = td(2)
         qb(1) = td(3)
        ils(1) = td(4)
         qa(2) = td(5)
         qb(2) = td(6)
        ils(2) = td(7)
         qa(3) = td(8)
         qb(3) = td(9)
        ils(3) = td(10)
        dq(1) = (qb(1)-qa(1)) 
        dq(2) = (qb(2)-qa(2)) 
        dq(3) = (qb(3)-qa(3)) 
        ils(4) = 0
        ils(5) = 0
        ils(6) = 0
        if(prt) write(iow,2004) qa(1),qb(1),ils(1),
     +                          qa(2),qb(2),ils(2),qa(3),qb(3),ils(3)
        if(prt.and.ior.lt.0) 
     +          write(*  ,2004) qa(1),qb(1),ils(1),
     +                          qa(2),qb(2),ils(2),qa(3),qb(3),ils(3)
        call eload2(ils)

        dl(1) = p2(1) - p1(1)
        dl(2) = p2(2) - p1(2)
        dl(3) = 0.d0
        if(ndm.eq.3) dl(3) = p2(3) - p1(3)
        call eload3(ils,dlsf,dl)
c....   loop over all elements
        if(iload.eq.1) then       ! linear
           inc = 1
           c11 = 2.d0/6.d0           
           c12 = 1.d0/6.d0
           c21 = 1.d0/6.d0           
           c22 = 2.d0/6.d0
        else if(iload.eq.2) then   ! quadratic
           inc = 2
           c11 = 1.d0/6.d0           
           c12 = 0.d0     
           c21 = 0.d0                
           c22 = 1.d0/6.d0
           c31 = 2.d0/6.d0
           c32 = 2.d0/6.d0
        end if
        do i = 1,ii-inc,inc
          n1 = idis(i)
          n2 = idis(i+inc)
          if(iload.eq.2) n3 = idis(i+1)
          ddl(1) = dabs( x(1,n2)-x(1,n1) )
          ddl(2) = dabs( x(2,n2)-x(2,n1) )
          ddl(3) = 0.d0
          if(ndm.eq.3) ddl(3) = dabs( x(3,n2)-x(3,n1) )
          call eload3(ils,ddlsf,ddl)
          dl1(1) = x(1,n1) - p1(1)
          dl1(2) = x(2,n1) - p1(2)
          dl1(3) = 0.d0
          if(ndm.eq.3) dl1(3) = x(3,n1) - p1(3)
          call eload3(ils,dsf1,dl1)
          dl2(1) = x(1,n2) - p1(1)
          dl2(2) = x(2,n2) - p1(2)
          dl2(3) = 0.d0
          if(ndm.eq.3) dl2(3) = x(3,n2) - p1(3)
          call eload3(ils,dsf2,dl2)
c.....    nodal loads
cww       do j = 1,2   
          ndf1 = min(ndf,3)  ! in case of linear loads only for 1-3 dofs
          do j = 1,ndf1
c....       load projection dls 
            dls = dlsf(j)
c....       load projection ddls
            ddls = ddlsf(j)
c...        q at point 1
c....       load projection ds1
            ds1 = dsf1(j)
            q1(j) = qa(j)+dq(j)/dls*ds1
c...        q at point 2
c....       load projection ds2
            ds2 = dsf2(j)
            q2(j) = qa(j)+dq(j)/dls*ds2                
c
            if(ibou.eq.0) then
             if(ior.ge.0) then          
              if(id(j,n1).eq.0) 
     +        f(j,n1) = f(j,n1) + ( c11*q1(j)+c12*q2(j) )*ddls
              if(id(j,n2).eq.0) 
     +        f(j,n2) = f(j,n2) + ( c21*q1(j)+c22*q2(j) )*ddls
             else if(ior.lt.0) then          
              if(id(j,n1).gt.0) 
     +        f(j,n1) = f(j,n1) + ( c11*q1(j)+c12*q2(j) )*ddls
              if(id(j,n2).gt.0) 
     +        f(j,n2) = f(j,n2) + ( c21*q1(j)+c22*q2(j) )*ddls
             end if
             if(iload.eq.2) then
              if(ior.ge.0) then          
               if(id(j,n3).eq.0) 
     +         f(j,n3) = f(j,n3) + ( c31*q1(j)+c32*q2(j) )*ddls
              else if(ior.lt.0) then          
               if(id(j,n3).gt.0) 
     +         f(j,n3) = f(j,n3) + ( c31*q1(j)+c32*q2(j) )*ddls
              end if 
             end if
            else
             f(j,n1) = q0(j)
             f(j,n2) = q0(j)
             if(iload.eq.2)  f(j,n3) = q0(j)
            end if
          end do
        end do
      else if(abs(icase).eq.3.or.abs(icase).eq.4) then
c....   load case 3: constant normal loads
c....   load case 4: constant normal loads  axisymmetric
        q0(1) = td(2)
        if(icase.eq.3) then
          if(prt)               write(iow,2010) q0(1)
          if(prt.and.ior.lt.0)  write(*  ,2010) q0(1)
        else if(icase.eq.4) then
          if(prt)               write(iow,2011) q0(1)
          if(prt.and.ior.lt.0)  write(*  ,2011) q0(1)
        end if   
c....   loop over all elements
        if(iload.eq.1) then       ! linear
          inc  = 1
          nel  = 2
          lint = 2
        else if(iload.eq.2) then   ! quadratic
          inc  = 2
          nel  = 3
          lint = 3
        end if
        do i = 1,ii-inc,inc
          n1 = idis(i)
          n2 = idis(i+inc)
          if(iload.eq.2) n3 = idis(i+1)
          xl(1,1) = x(1,n1)
          xl(2,1) = x(2,n1)
          xl(1,2) = x(1,n2)
          xl(2,2) = x(2,n2)
          if(iload.eq.2) then
            xl(1,2) = x(1,n3)
            xl(2,2) = x(2,n3)
            xl(1,3) = x(1,n2)
            xl(2,3) = x(2,n2)
          end if
          call gaus1D(lint,pg,wg)
          do l = 1,lint
            xsi = pg(l)
            call shape1D(shp,xsi,xl,ds,ndm,nel)
            dxl = wg(l)*ds
            xx  = 0.d0
            yy  = 0.d0
            xt  = 0.d0
            yt  = 0.d0
            do k = 1,nel
              xx = xx + shp(2,k)*xl(1,k)
              yy = yy + shp(2,k)*xl(2,k)
              xt = xt + shp(1,k)*xl(1,k)
              yt = yt + shp(1,k)*xl(2,k)
            end do      
            xn =  yt
            yn = -xt
c           load in x and y direction
            qxn = xn*q0(1)
            qyn = yn*q0(1)
            if(icase.eq.4) then
              qxn = qxn*pi2*xx
              qyn = qyn*pi2*xx
            end if        
c
            if(iload.eq.1) then
              if(ior.ge.0) then          
               if(id(1,n1).eq.0) f(1,n1) = f(1,n1) + qxn*shp(2,1)*dxl
               if(id(1,n2).eq.0) f(1,n2) = f(1,n2) + qxn*shp(2,2)*dxl
               if(id(2,n1).eq.0) f(2,n1) = f(2,n1) + qyn*shp(2,1)*dxl
               if(id(2,n2).eq.0) f(2,n2) = f(2,n2) + qyn*shp(2,2)*dxl
              else if(ior.lt.0) then          
               if(id(1,n1).gt.0) f(1,n1) = f(1,n1) + qxn*shp(2,1)*dxl
               if(id(1,n2).gt.0) f(1,n2) = f(1,n2) + qxn*shp(2,2)*dxl
               if(id(2,n1).gt.0) f(2,n1) = f(2,n1) + qyn*shp(2,1)*dxl
               if(id(2,n2).gt.0) f(2,n2) = f(2,n2) + qyn*shp(2,2)*dxl
              end if
            else if(iload.eq.2) then
              if(ior.ge.0) then          
               if(id(1,n1).eq.0) f(1,n1) = f(1,n1) + qxn*shp(2,1)*dxl
               if(id(1,n2).eq.0) f(1,n2) = f(1,n2) + qxn*shp(2,3)*dxl
               if(id(1,n3).eq.0) f(1,n3) = f(1,n3) + qxn*shp(2,2)*dxl
               if(id(2,n1).eq.0) f(2,n1) = f(2,n1) + qyn*shp(2,1)*dxl
               if(id(2,n2).eq.0) f(2,n2) = f(2,n2) + qyn*shp(2,3)*dxl
               if(id(2,n3).eq.0) f(2,n3) = f(2,n3) + qyn*shp(2,2)*dxl
              else if(ior.lt.0) then          
               if(id(1,n1).gt.0) f(1,n1) = f(1,n1) + qxn*shp(2,1)*dxl
               if(id(1,n2).gt.0) f(1,n2) = f(1,n2) + qxn*shp(2,3)*dxl
               if(id(1,n3).gt.0) f(1,n3) = f(1,n3) + qxn*shp(2,2)*dxl
               if(id(2,n1).gt.0) f(2,n1) = f(2,n1) + qyn*shp(2,1)*dxl
               if(id(2,n2).gt.0) f(2,n2) = f(2,n2) + qyn*shp(2,3)*dxl
               if(id(2,n3).gt.0) f(2,n3) = f(2,n3) + qyn*shp(2,2)*dxl
              end if
            end if
          end do   
        end do
      else if(abs(icase).eq.5) then
c....   load case 5: boundary terms for warping function
c       Int [(ny*z - nz*y) N_k] ds
        if(prt)               write(iow,2005)
        if(prt.and.ior.lt.0)  write(*  ,2005)
c....   loop over all elements
        if(iload.eq.1) then       ! linear
          inc  = 1
          nel  = 2
          lint = 2
        else if(iload.eq.2) then   ! quadratic
          inc  = 2
          nel  = 3
          lint = 3
        end if
        do i = 1,ii-inc,inc
          n1 = idis(i)
          n2 = idis(i+inc)
          if(iload.eq.2) n3 = idis(i+1)
          xl(1,1) = x(1,n1)
          xl(2,1) = x(2,n1)
          xl(1,2) = x(1,n2)
          xl(2,2) = x(2,n2)
          if(iload.eq.2) then
            xl(1,2) = x(1,n3)
            xl(2,2) = x(2,n3)
            xl(1,3) = x(1,n2)
            xl(2,3) = x(2,n2)
          end if
          call gaus1D(lint,pg,wg)
          do l = 1,lint
            xsi = pg(l)
            call shape1D(shp,xsi,xl,ds,ndm,nel)
            dxl = wg(l)*ds
            xx  = 0.d0
            yy  = 0.d0 
            xt  = 0.d0
            yt  = 0.d0
            do k = 1,nel
              xx = xx + shp(2,k)*xl(1,k)
              yy = yy + shp(2,k)*xl(2,k)
              xt = xt + shp(1,k)*xl(1,k)
              yt = yt + shp(1,k)*xl(2,k)
            end do      
            xn =  yt
            yn = -xt
c
            q = xn*yy - yn*xx
c
            if(icase.lt.0) q = -q
            if(iload.eq.1) then
              if(ior.ge.0) then          
               if(id(1,n1).eq.0) f(1,n1) = f(1,n1) + q*shp(2,1)*dxl
               if(id(1,n2).eq.0) f(1,n2) = f(1,n2) + q*shp(2,2)*dxl
              else if(ior.lt.0) then          
               if(id(1,n1).gt.0) f(1,n1) = f(1,n1) + q*shp(2,1)*dxl
               if(id(1,n2).gt.0) f(1,n2) = f(1,n2) + q*shp(2,2)*dxl
              end if
            else if(iload.eq.2) then
              if(ior.ge.0) then          
               if(id(1,n1).eq.0) f(1,n1) = f(1,n1) + q*shp(2,1)*dxl
               if(id(1,n2).eq.0) f(1,n2) = f(1,n2) + q*shp(2,3)*dxl
               if(id(1,n3).eq.0) f(1,n3) = f(1,n3) + q*shp(2,2)*dxl
              else if(ior.lt.0) then          
               if(id(1,n1).gt.0) f(1,n1) = f(1,n1) + q*shp(2,1)*dxl
               if(id(1,n2).gt.0) f(1,n2) = f(1,n2) + q*shp(2,3)*dxl
               if(id(1,n3).gt.0) f(1,n3) = f(1,n3) + q*shp(2,2)*dxl
              end if
            end if
          end do   
        end do
      else if(abs(icase).eq.6) then
c       Int [(-n_y*cy*zz**2 - n_z*cz*yy**2) N] ds 
        azz = td(2)
        ayy = td(3)
        ayz = td(4)
        qy  = td(5)
        qz  = td(6)
        xnue= td(7)
        y0  = td(8)
        z0  = td(9)
        if(prt)             write(iow,2008) azz,ayy,ayz,y0,qy,qz,xnue,z0
        if(prt.and.ior.lt.0)write(*  ,2008) azz,ayy,ayz,y0,qy,qz,xnue,z0
c....   loop over all elements
        if(iload.eq.1) then       ! linear
          inc  = 1
          nel  = 2
          lint = 2
        else if(iload.eq.2) then   ! quadratic
          inc  = 2
          nel  = 3
          lint = 3
        end if
        do i = 1,ii-inc,inc
          n1 = idis(i)
          n2 = idis(i+inc)
          if(iload.eq.2) n3 = idis(i+1)
          xl(1,1) = x(1,n1)
          xl(2,1) = x(2,n1)
          xl(1,2) = x(1,n2)
          xl(2,2) = x(2,n2)
          if(iload.eq.2) then
            xl(1,2) = x(1,n3)
            xl(2,2) = x(2,n3)
            xl(1,3) = x(1,n2)
            xl(2,3) = x(2,n2)
          end if
          call gaus1D(lint,pg,wg)
          do l = 1,lint
            xsi = pg(l)
            call shape1D(shp,xsi,xl,ds,ndm,nel)
            dxl = wg(l)*ds
            yy  = -y0
            zz  = -z0
            xt  = 0.d0
            yt  = 0.d0
            do k = 1,nel
              yy = yy + shp(2,k)*xl(1,k)
              zz = zz + shp(2,k)*xl(2,k)
              xt = xt + shp(1,k)*xl(1,k)
              yt = yt + shp(1,k)*xl(2,k)
            end do      
            xn =  yt
            yn = -xt
c
            det = ayy*azz-ayz*ayz             
            a1  = (qy*azz-qz*ayz)/det
            a2  = (qz*ayy-qy*ayz)/det
            dnue=0.5d0*xnue/(1.d0+xnue)
c
            q = dnue*(-xn*a1*zz*zz - yn*a2*yy*yy)
c
            if(icase.lt.0) q = -q
            if(iload.eq.1) then
              if(ior.ge.0) then          
               if(id(1,n1).eq.0) f(1,n1) = f(1,n1) + q*shp(2,1)*dxl
               if(id(1,n2).eq.0) f(1,n2) = f(1,n2) + q*shp(2,2)*dxl
              else if(ior.lt.0) then          
               if(id(1,n1).gt.0) f(1,n1) = f(1,n1) + q*shp(2,1)*dxl
               if(id(1,n2).gt.0) f(1,n2) = f(1,n2) + q*shp(2,2)*dxl
              end if  
            else if(iload.eq.2) then
              if(ior.ge.0) then          
               if(id(1,n1).eq.0) f(1,n1) = f(1,n1) + q*shp(2,1)*dxl
               if(id(1,n2).eq.0) f(1,n2) = f(1,n2) + q*shp(2,3)*dxl
               if(id(1,n3).eq.0) f(1,n3) = f(1,n3) + q*shp(2,2)*dxl
              else if(ior.lt.0) then          
               if(id(1,n1).gt.0) f(1,n1) = f(1,n1) + q*shp(2,1)*dxl
               if(id(1,n2).gt.0) f(1,n2) = f(1,n2) + q*shp(2,3)*dxl
               if(id(1,n3).gt.0) f(1,n3) = f(1,n3) + q*shp(2,2)*dxl
              end if
            end if
          end do   
        end do
      end if
c
      if(prt)              write(iow,2000) (k,k=1,ndf)
      if(prt.and.ior.lt.0) write(*  ,2000) (k,k=1,ndf)
      call pzero(fm,ndf)  
      call pzero(fma,ndf)  
      do 410 n = 1,ii
        in = idis(n)
        do nn=1,ndf
          fm(nn)=fm(nn)+f(nn,in)
          fma(nn)=fma(nn)+abs(f(nn,in))
        end do 
        if(prt)                write(iow,2001) in,(f(k,in),k=1,ndf)
          if(prt.and.ior.lt.0) write(*  ,2001) in,(f(k,in),k=1,ndf)
410   continue
       
        if(prt) then                
          write(iow,2012) ( fm(k),k=1,ndf)
          write(iow,2013) (fma(k),k=1,ndf)
          if(ior.lt.0) then 
            write(  *,2012) ( fm(k),k=1,ndf)
            write(  *,2013) (fma(k),k=1,ndf)
          end if
        end if
      go to 100
2000  format(/' Actual sum of n o d a l  forces along line'/,
     +        ' node',6(i3,'- force  ') )
2001  format(i5,6(e12.5))
2002  format(/,' Load generation with ELOAD on ',a8,/,
     +         ' From         Point ',3(1x,e12.5),/,
     +         ' To           Point ',3(1x,e12.5),/,
     +           a20                 ,3(1x,e12.5),/)
2003  format(' Load case 1: const. loads in global 1-3 direction:',/,
     +         6(1x,e12.5),/,' load projection to ',
     +         i2,'/',i2,'/',i2,'/',i2,'/',i2,'/',i2,' direction.')
2004  format(' Load case 2: linear loads in global 1-3 direction:',/,
     +       ' q1x=:',e12.5,' q2x=:',e12.5,' ilx=:',i2,/ 
     +       ' q1y=:',e12.5,' q2y=:',e12.5,' ily=:',i2,/  
     +       ' q1z=:',e12.5,' q2z=:',e12.5,' ilz=:',i2)
2005  format(' Load case 5: boundary terms for warping function:',/,
     +       ' Int [(n_y*z - n_z*y) N] ds') 
2006  format(' tol = ',e12.5,' for boundary(=1) =',i5,/,
     +         ' Radius(only ityp=3)',e12.5)
2007  format(' ELOAD: Point is not on defined circle')
2008  format(' Load case 6: boundary terms for stress function:',/,
     +       ' Int [(n_y*cy*zz**2 - n_z*cz*yy**2) N] ds',/,
     +       ' azz = ',e12.5,' ayy = ',e12.5,' ayz = ',e12.5,
     +       ' y0  = ',e12.5,/,
     +       ' qy  = ',e12.5,' qz  = ',e12.5,' xnue= ',e12.5,
     +       ' z0  = ',e12.5)
2009  format(' for ',a11,' Elements',/)
2010  format(' Load case 3: const. normal loads in 1-2 plane:',
     +         1x,e12.5 )
2011  format(' Load case 4: const. normal loads axisymmetric:',
     +         1x,e12.5 )
2012  format(' SUM ',6(e12.5))
2013  format(' SUMA',6(e12.5))
3001  format('Input:P_1(x,y,[z]),P_2(x,y,[z]),[P_3(x,y,[z])],',/,
     +        'tol,ibou,rad >',$)
      end
c
      subroutine eload1(dist,idis,numel,tol)
c-----------------------------------------------------------------------
c
c      Purpose: sort the node numbers due to shortest distance from 
c               reference point
c
c      Inputs:
c         dist(*)   - Distance of Node to reference point
c         idis(*)   - Node number
c         numel     - Number of nodes in list
c         tol       - Tolerance for 'same' node 
c
c      Outputs:
c         dist(*)   - Distance of Node to reference point
c         idis(*)   - Node number
c
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension dist(*),idis(*)
c.... sort on first value
      do n = 1,numel-1
          do i = n+1,numel
            if(dist(i).lt.dist(n)) then
              sn = dist(n)
              dist(n) = dist(i)
              dist(i) = sn
              in      = idis(n)
              idis(n) = idis(i)
              idis(i) = in
            end if
          end do
      end do
c...  test on same distance
      do n = 1, numel-1
        d1 = dist(n)
        i1 = idis(n)
        do m = n+1,numel
          d2 = dist(m)
          i2 = idis(m)
c...      delete values on nodes with higher number
          if(dabs(d2-d1).lt.tol) then
            if(i1.lt.i2 .and. i1.ne.0) then
              dist(m) = 0.d0
              idis(m) = 0
            else if(i2.lt.i1 .and. i2.ne.0) then
              dist(n) = 0.d0
              idis(n) = 0
            end if
          end if
        end do
      end do      
c...  delete zero values
      numel1 = numel
      n = 0
10    continue
      n = n+1
      i1 = idis(n)
        if(i1.eq.0) then
          do m = n+1,numel
             dist(m-1) = dist(m)
             idis(m-1) = idis(m)
          end do
          numel1 = numel1 - 1
          n = n - 1
        end if
      if(n.lt.numel1) goto 10
c.... actual number of nodes        
      numel = numel1
      return
      end
c
      subroutine eload2(ils)
c-----------------------------------------------------------------------
c
c      Purpose: check for errors in choosing load projection, 1-3,4-6                        |

c      Inputs:
c         ils(*)   - input planes for loads
c
c      Outputs:    - errors
c
c-----------------------------------------------------------------------
      dimension ils(6)

      do i=1,4,3 
        if(ils(i).eq.1)  stop 'error ELOA: input dir for load q_x'
        if(ils(i).eq.12) stop 'error ELOA: input dir for load q_x'
        if(ils(i).eq.13) stop 'error ELOA: input dir for load q_x'

        if(ils(i+1).eq.2)  stop 'error ELOA: input dir for load q_y'
        if(ils(i+1).eq.12) stop 'error ELOA: input dir for load q_y'
        if(ils(i+1).eq.23) stop 'error ELOA: input dir for load q_y'

      if(ils(i+2).eq.3)  stop 'error ELOA: input dir for load q_z'
      if(ils(i+2).eq.13) stop 'error ELOA: input dir for load q_z'
      if(ils(i+2).eq.23) stop 'error ELOA: input dir for load q_z'
      end do       

      return
      end
c
      subroutine eload3(ils,dlsf,dl)
c-----------------------------------------------------------------------
c
c      Purpose: calculate element length for load projection

c      Inputs:
c         ils(*)   - input planes for loads
c         dl(3)    - length dx dy dz
c
c      Outputs:   
c        dlsf      - projected length ds
c
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension ils(6),dlsf(6),dl(3)
      do j = 1,6
        dlsf(j) = dsqrt( dot(dl,dl,3) )     
        if(ils(j).eq. 1) dlsf(j) = dabs(dl(1))
        if(ils(j).eq. 2) dlsf(j) = dabs(dl(2))
        if(ils(j).eq. 3) dlsf(j) = dabs(dl(3))
        if(ils(j).eq.12) dlsf(j) = dsqrt(dl(1)*dl(1)+dl(2)*dl(2))
        if(ils(j).eq.13) dlsf(j) = dsqrt(dl(1)*dl(1)+dl(3)*dl(3))
        if(ils(j).eq.23) dlsf(j) = dsqrt(dl(2)*dl(2)+dl(3)*dl(3))
      end do
      return
      end
c      
      subroutine endclr (subnam,chr)
c----------------------------------------------------------------------
c
c      Purpose: end-of-file clearing routine
c
c      Inputs:
c         subnam   - name of Subroutine with error
c         chr      - error data
c
c      Outputs:    - error message 
c
c----------------------------------------------------------------------
      USE iofile
      character     subnam*(*),chr*(*)
c
      if (ior.gt.0)  then
         write(iow,2000) subnam
cww      stop
       return
      else
         backspace  5
         chr = ' '
      end if
      return
c.... error message, only for batch mode
 2000 format (' ** ERROR in ',a6,' ** end of file encountered')
      end
c
      subroutine    errclr (subnam)
c----------------------------------------------------------------------
c
c      Purpose: input error clearing routine
c
c      Inputs:
c         subnam   - name of Subroutine with error
c
c      Outputs:    - error message 
c
c----------------------------------------------------------------------
      USE iofile
      character     subnam*(*), yyy*80
c
      write(iow,2002) subnam
      write(  *,2002) subnam
 
cww   write(yyy,2002) subnam
cww   call drawmess(yyy,1,-2)
      return
2002  format (' ** ERROR in ',a6,' ** reinput last record')
      end
c
      subroutine esideb(x,ien,inn,ieb,ef,ndm,nde,numnp,prt)
c----------------------------------------------------------------------
c      Purpose: input and set the b.c. for edge nodes
c
c     Inputs:
c         x(i,n)      - nodal coordinates: i=1,ndm ; n=1,numnp
c         ien(n)      - pointer array to edge sides: n=1,numnp
c         inn(i)      - list of edge sides: i = ien(n)+1,ien(n+1)
c         ndm         - space dimension of finite element mesh
c         nde         - number of edge d.o.f.
c         numnp       - number of nodal points
c
c     Outputs:
c
c         ieb(i,n)      - b.c.  array to edge nodes: i=1,nde
c                         n is defined by ien and inn.
c         ef (i,n)      - force array to edge nodes: i=1,nde
c                         n is defined by ien and inn.
c
c----------------------------------------------------------------------
c....  Declare variable types
      USE iofile
      logical pcomp,prt
      character wrd*4,labl*8,lab13*13
      integer ndm, nde, numnp, i, i1, i2, idr
      integer j,   n,   n1,  n2
      real*8  x1, x2, tol
c..... Declare array types
      integer ien(*),       inn(*),      ieb(nde,*), ilb(9)
      real*8  x(ndm,numnp), ef(nde,*),   elf(9),     td(10)
      data labl/'--------'/, lab13/'-------------'/, tol/1.e-8/
1     if(ior.lt.0) then
        write(*,2000)
        read(*,1000) wrd
      else
        read(ior,1000) wrd
      end if
      if(pcomp(wrd,'ende',4)) return
      if(pcomp(wrd,'edbc',4)) then
c.... input nodal pairs and boundary conditions
        if(prt) then
          if(ior.lt.0) then
            write(*,2010) (i,i=1,nde)
            write(*,2020) (labl,i=1,nde)
          end if
          write(iow,2010) (i,i=1,nde)
          write(iow,2020) (labl,i=1,nde)
        end if
100     if(ior.lt.0) write(*,2001)
        call dinput(td,nde+2)
        i1 = min(td(1),td(2))
        i2 = max(td(1),td(2))
        if(i1.le.0.or.i1.ge.numnp) go to 400
        do 105 n = 1,nde
          ilb(n) = td(n+2)
105     continue
        n1 = ien(i1) + 1
        n2 = ien(i1+1)
        do 110 n = n1,n2
          if(inn(n).eq.i2) go to 120
110     continue
        write(*,3000) i1,i2
        go to 100
120     do 130 i = 1,nde
          ieb(i,n) = ilb(i)
130     continue
        if(prt) then
          if(ior.lt.0) write(*,2030) i1,i2,(ilb(i),i=1,nde)
          write(iow,2030) i1,i2,(ilb(i),i=1,nde)
        end if
        go to 100
c.... set the b.c. for this node
      else if(pcomp(wrd,'edgb',4)) then
200     if(ior.lt.0) write(*,2002)
        call dinput(td,nde+2)
        idr = td(1)
        idr = min(ndm,idr)
        if(idr.le.0) go to 400
        do 205 n = 1,nde
          ilb(n) = td(n+2)
205     continue
        if(prt) then
          if(ior.lt.0) then
            write(*,2015) idr,td(2),(i,i=1,nde)
            write(*,2020) (labl,i=1,nde)
          end if
          write(iow,2015) idr,td(2),(i,i=1,nde)
          write(iow,2020) (labl,i=1,nde)
        end if
        do 230 n1 = 1,numnp-1
          x1 = x(idr,n1)
          i1 = ien(n1) + 1
          i2 = ien(n1+1)
          do 220 i = i1,i2
            n2 = inn(i)
            x2 = x(idr,n2)
            if(abs(td(2)-x1).lt.tol .and. abs(td(2)-x2).lt.tol) then
              do 210 j = 1,nde
                ieb(j,i) = ilb(j)
210           continue
              if(prt) then
                if(ior.lt.0) write(*,2030) n1,n2,(ilb(j),j=1,nde)
                write(iow,2030) n1,n2,(ilb(j),j=1,nde)
              end if
            end if
220       continue
230     continue
        go to 200
      else if(pcomp(wrd,'edgf',4)) then
c.... input nodal pairs and boundary force/displacements
        if(prt) then
          if(ior.lt.0) then
            write(*,2040) (i,i=1,nde)
            write(*,2050) (lab13,i=1,nde)
          end if
          write(iow,2040) (i,i=1,nde)
          write(iow,2050) (lab13,i=1,nde)
        end if
300     if(ior.lt.0) write(*,2001)
        call dinput(td,nde+2)
        i1 = min(td(1),td(2))
        i2 = max(td(1),td(2))
        if(i1.le.0.or.i1.ge.numnp) go to 400
        do 305 n = 1,nde
          elf(n) = td(n+2)
305     continue
        n1 = ien(i1) + 1
        n2 = ien(i1+1)
        do 310 n = n1,n2
          if(inn(n).eq.i2) go to 320
310     continue
        write(*,3000) i1,i2
        go to 300
320     do 330 i = 1,nde
          ef(i,n) = elf(i)
330     continue
        if(prt) then
          if(ior.lt.0) write(*,2060) i1,i2,(elf(i),i=1,nde)
          write(iow,2060) i1,i2,(elf(i),i=1,nde)
        end if
        go to 300
      end if
400   go to 1
c
1000  format(a4)
2000  format('  Input: edbc or edgb'/'  ->',$)
2001  format('  Input:  n1, n2, b.c.-codes'/'  ->',$)
2002  format('  Input: dir, xx, b.c.-codes'/'  ->',$)
2010  format(/10x,'  Edge Node Boundary Conditions'
     1      //10x,'  Nodal Pair |  ',8(i3,'-B.C.'))
2015  format(/10x,'  For x(',i1,') = ',e12.5
     1      //10x,'  Edge Node Boundary Conditions'
     2      //10x,'  Nodal Pair |  ',8(i3,'-B.C.'))
2020  format(10x,'  --n1---n2--+--',8a8)
2030  format(11x,2i5,'  |',8i8)
2040  format(/10x,'  Edge Node Forced Conditions'
     1      //10x,'  Nodal Pair |  ',8(i3,'-Forc/Disp'))
2050  format(10x,'  --n1---n2--+--',8a13)
2060  format(11x,2i5,'  |',1p8e13.4)
3000  format(10x,' *ERROR* No pair:',2i5)
      end
c
      subroutine esiden(ix,ie,ien,inn,ine,nie,nen1,numel,numnp)
c----------------------------------------------------------------------
c
c     Purpose: compute the connection list for edge nodes
c
c     Inputs:
c
c         ix(i,numel) - nodal connection lists: i = 1,nen1
c                     - ma = ix(nen1,numel) is material set
c         ie(nie-1,i) - element type for ma
c         iedg(j,i)   - edge sequence for element type i: j = 1,nedg
c         nen1        - dimension for ix-array
c         nedg(i)     - number of edges on element type i
c         numel       - number of elements
c         numnp       - number of nodes
c
c     Outputs:
c
c         ien(n)      - pointer array to edge sides: n=1,numnp-1
c         inn(i)      - list of edge sides: i = ien(n-1)+1,ien(n)
c
c     Temporary:
c
c         ine(i)      - temporary array: i = 1 to <= nedg*numnp
c
c----------------------------------------------------------------------
c..... Declare variable types
      USE edgdat
      integer nie, nen1, numel, numnp, i, i1, i2, j, j1, j2, k, k1, k2
      integer iel, n, ma, ne
c..... Declare array types
      integer ix(nen1,numel), ie(nie,*), ien(numnp), inn(*), ine(*)
c
c  1. count number of times a node is on any element
c
      call pzeroi(ien, numnp)
      do 110 n = 1,numel
        ma  = ix(nen1,n)
        iel = ie(nie-1,ma)
        do 100 i = 1,nedg(iel)
          i1 = iedg(i,iel)
          i1 = abs(ix(i1,n))
          if(i1.gt.0) then
            ien(i1) = ien(i1) + 1
          end if
100     continue
110   continue
      do 120 n = 2,numnp
        ien(n) = ien(n) + ien(n-1)
120   continue
c
c  2. find elements for each node
c
      call pzeroi(ine, ien(numnp))
      do 230 n = 1,numel
        ma  = ix(nen1,n)
        iel = ie(nie-1,ma)
        do 220 i = 1,nedg(iel)
          i1 = iedg(i,iel)
          i1 = abs(ix(i1,n))
          if(i1.gt.0) then
            j1 = 1
            if(i1.gt.1) j1 = ien(i1-1) + 1
            do 200 j = j1,ien(i1)
              if(ine(j).eq.0) go to 210
200         continue
210         ine(j) = n
          end if
220     continue
230   continue
c
c  3. compute the node list profile
c
      k1 = 1
      k2 = 0
      j1 = 1
      do 340 n = 1,numnp
c.... find elements associated with node
        j2 = ien(n)
        if(j1.le.j2) then
          do 330 j = j1,j2
            ne  = ine(j)
            ma  = ix(nen1,ne)
            iel = ie(nie-1,ma)
            ned = nedg(iel)
            do 320 i = 1,ned
              i1 = iedg(i,iel)
              i1 = abs(ix(i1,ne))
              if(i1.eq.n) then
                if(i.eq.1) then
                  i2 = iedg(ned,iel)
                else
                  i2 = iedg(i-1,iel)
                end if
                i2 = abs(ix(i2,ne))
                if(i2.gt.n) then
c.... check to see if the edge pair already exists
                  if(k1.le.k2) then
                    do 300 k = k1,k2
                      if(i2.eq.inn(k)) go to 305
300                 continue
                  end if
c.... add new edge pair
                  k2 = k2 + 1
                  inn(k2) = i2
                end if
305             if(i.eq.ned) then
                  i2 = iedg(  1,iel)
                else
                  i2 = iedg(i+1,iel)
                end if
                i2 = abs(ix(i2,ne))
                if(i2.gt.n) then
c.... check to see if the edge pair already exists
                  if(k1.le.k2) then
                    do 310 k = k1,k2
                      if(i2.eq.inn(k)) go to 320
310                 continue
                  end if
c.... add new edge pair
                  k2 = k2 + 1
                  inn(k2) = i2
                end if
              end if
320         continue
330       continue
        end if
        j1 = j2 + 1
        ien(n) = k2
        k1 = k2 + 1
340   continue
c.... shift the list up for later uses
      do 350 n = numnp,2,-1
        ien(n) = ien(n-1)
350   continue
      ien(1) = 0
c
      end
c
      subroutine eprint(ien,inn,numnp)
c----------------------------------------------------------------------
c
c     Purpose: print the nodal lists for edge nodes
c
c     Inputs:
c
c         ien(n)      - pointer array to edge sides: n=1,numnp-1
c         inn(i)      - list of edge sides: i = ien(n-1)+1,ien(n)
c         numnp       - number of nodes
c
c     Outputs:
c
c----------------------------------------------------------------------
c..... Declare variable types
      USE iofile
      integer numnp, n, j, j1, j2
c..... Declare array types
      integer ien(numnp),inn(*)
      if(ior.lt.0) write(*,2000)
      write(iow,2000)
      do 100 n = 1,numnp-1
        j1 = ien(n) + 1
        j2 = ien(n+1)
        if(ior.lt.0) write(*,2001) n,(inn(j),j=j1,j2)
        write(iow,2001) n,(inn(j),j=j1,j2)
100   continue
c
2000  format(/10x,'  Edge Connection List'
     1      //10x,'  Node | Side Connections'
     1       /10x,' ------+-----------------')
2001  format(i15,'  |',9i5)
      end
c
      subroutine evalex(xs,v,val,ns,error)
c----------------------------------------------------------------------
c     Purpose: Identify expression in character string and evaluate

c     Inputs:
c         xs(*) - Character string of input data
c         v(*)  - Array of real values
c         ns    - Length of character string

c     Outputs:
c        val   - Value of expression
c        error - Error indicator if true
c
c.... evaluate an expression: (+) add,   (-) subtract, (*) multiply,
c                             (/) divide,(^) power    
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical error
      character*1 xs(ns),x(80),y,op(25)
      dimension v(*)
c.... first pack the expression by removing any ' ' characters
      call pevpak(xs,ns, x,num)
      k  = 0
      do 100 i = 1,num
        if(k.eq.0 .and. x(i).eq.' ') go to 100
        if((x(i).eq.'+') .or. (x(i).eq.'-')) then
          if(i.gt.2) then
            y = x(i-1)
            if(y.eq.'e'.or.y.eq.'d'.or.y.eq.'E'.or.y.eq.'D') then
              if((x(i-2).ge.'0'.and.x(i-2).le.'9').or.
     1            x(i-2).eq.'.') go to 100
            end if
          end if
          k     = k + 1
          op(k) = x(i)
          x(i)  = ','
        end if
        if((x(i).eq.'*') .or. (x(i).eq.'/') .or. (x(i).eq.'^')) then
          k     = k + 1
          op(k) = x(i)
          x(i)  = ','
        end if
100   continue
      call dcheck(x,v,num,error)
      if(error) return
c.... compute the value of the expression
      val = v(1)
      if(k.ge.1) then
c..... 1. evaluate all exponentiations
        i = 1
120     continue
       if(op(i).eq.'^') then
          v(i) = v(i) ** v(i+1)
          k = k - 1
          do 110 j = i,k
            v(j+1) = v(j+2)
            op(j) = op(j+1)
110         continue
       end if
          i = i + 1
        if(i.le.k) go to 120
c..... 2. evaluate all multiplications and divisions
        i = 1
140     continue
          if(op(i).eq.'*') v(i) = v(i) * v(i+1)
          if(op(i).eq.'/') v(i) = v(i) / v(i+1)
          if(op(i).eq.'*' .or. op(i).eq.'/') then
            k = k - 1
            do 130 j = i,k
              v(j+1) = v(j+2)
              op(j) = op(j+1)
130         continue
          else
            i = i + 1
          end if
        if(i.le.k) go to 140
c..... 3. evaluate all additions and subtractions
        val = v(1)
        if(k.gt.0) then
          do 160 i = 1,k
            if(op(i).eq.'+') val = val + v(i+1)
            if(op(i).eq.'-') val = val - v(i+1)
160       continue
        end if
      end if
      return
      end
c
      subroutine pevpak(xs,ns, x,n)
c----------------------------------------------------------------------
c
c      Purpose: Remove unwanted blanks from character string
c      Inputs:
c         xs(*)   - Unpacked character string
c         ns      - length of unpacked string

c      Outputs:
c         x(*)    - Packed character string
c         n       - length of packed string
c
c----------------------------------------------------------------------
      character*1 xs(ns), x(ns)
      n = 0
      do 100 i = 1,ns
        if(xs(i).ne.' ') then
          n = n + 1
          x(n) = xs(i)
        end if
100   continue
      return
      end
c
      subroutine f3d(shp,ul,f,detf,ndf,nel,isw)
c----------------------------------------------------------------------
c
c     Purpose: calculate F for 3-D problems
c
c     Inputs:
c         shp(3,*)  - Shape function array
c         ul(ndf,*) - element displacements
c         ndf       - Number dof/node
c         nel       - Number nodes/element
c         isw       = 1 Ref.konfiguration, bekannt X_A: F    = 1+Grad u
c         isw       = 2 Mom.konfiguration, bekannt x_i: F^-1 = 1-grad u
c
c     Outputs:
c         f(3,*)    - F
c         detf      - det F
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      real*8 shp(4,*),ul(ndf,*),f(3,*),finv(3,3),af(3,3)
c
c.... isw = 1 ==> Referenzkonfiguration, bekannt X_A: F = 1 + Grad u
c.... isw = 2 ==> Momentankonfiguration, bekannt x_i: F^-1 = 1 - grad u
c
      goto(1,2),isw
c
c.... Berechnung von F
c
1     continue
      do 12 i = 1,3
        do 11 j = 1,3
          f(i,j) = 0.0d0
          do 10 k = 1,nel
            f(i,j) = f(i,j) + ul(i,k)*shp(j,k)
  10      continue
  11    continue
        f(i,i) = f(i,i) + 1.0d0
  12  continue
c.... compute adjoint to jacobian
      af(1,1) = f(2,2)*f(3,3) - f(2,3)*f(3,2)
      af(2,1) = f(2,3)*f(3,1) - f(2,1)*f(3,3)
      af(3,1) = f(2,1)*f(3,2) - f(2,2)*f(3,1)
c.... compute determinant of jacobian
      detf    = f(1,1)*af(1,1)+f(1,2)*af(2,1)+f(1,3)*af(3,1)
      return
2     continue
c
c.... Berechnung von F invers
c
      do 22 i = 1,3
        do 21 j = 1,3
          finv(i,j) = 0.0d0
          do 20 k = 1,nel
            finv(i,j) = finv(i,j) - ul(i,k)*shp(j,k)
  20      continue
  21    continue
        finv(i,i) = finv(i,i) + 1.0d0
  22  continue
c.... compute adjoint to jacobian
      af(1,1) = finv(2,2)*finv(3,3) - finv(2,3)*finv(3,2)
      af(1,2) = finv(1,3)*finv(3,2) - finv(1,2)*finv(3,3)
      af(1,3) = finv(1,2)*finv(2,3) - finv(1,3)*finv(2,2)
      af(2,1) = finv(2,3)*finv(3,1) - finv(2,1)*finv(3,3)
      af(2,2) = finv(1,1)*finv(3,3) - finv(1,3)*finv(3,1)
      af(2,3) = finv(1,3)*finv(2,1) - finv(1,1)*finv(2,3)
      af(3,1) = finv(2,1)*finv(3,2) - finv(2,2)*finv(3,1)
      af(3,2) = finv(1,2)*finv(3,1) - finv(1,1)*finv(3,2)
      af(3,3) = finv(1,1)*finv(2,2) - finv(1,2)*finv(2,1)
c.... Berechne Determinante von F invers
      detinv  = finv(1,1)*af(1,1)+finv(1,2)*af(2,1)+finv(1,3)*af(3,1)
c
c.... Berechnung von F
c
      do 121 i = 1,3
        do 120 j = 1,3
          f(i,j) = af(i,j)/detinv
120     continue
121   continue
      detf    = 1.d0/detinv
      return
      end
c
      subroutine b3d(shp,ul,ndf,nel,b,detf)
c----------------------------------------------------------------------
c
c      Purpose: calculate left Cauchy-Green tensor b
c
c      Inputs:
c         shp(3,*)  - Shape function array
c         ul(ndf,*) - element displacements
c         ndf       - Number dof/node
c         nel       - Number nodes/element
c         detf      - det F
c
c      Outputs:
c         b(6)      - b
c
c      Comments: 

c           b_11    b(1)
c           b_22    b(2)
c      b =  b_33  = b(3)
c           b_12    b(4)
c           b_23    b(5)
c           b_13    b(6)
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c
      real*8 shp(4,*),ul(ndf,*),f(3,3),b(*),bi(3,3)
c
c.... Berechne Deformationsgradienten F
c
      call f3d(shp,ul,f,detf,ndf,nel,2)
c
c.... Berechne linken Cauchy-Green Tensor b
c
      do 339 i = 1,3
        do 340 j = i,3
          bi(i,j) = 0.d0
          do 341 k = 1,3
            bi(i,j) = bi(i,j) + f(i,k)*f(j,k)
 341      continue
 340    continue
 339  continue
      do 350 i = 1,3
        b(i) = bi(i,i)
 350  continue
      b(4) = bi(1,2)
      b(5) = bi(2,3)
      b(6) = bi(1,3)
c
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine formfe(u,b,a,c,afl,bfl,cfl,dfl,isw,nl1,nl2,nl3)
c----------------------------------------------------------------------
c
c      Purpose: Forms finite element arrays as required

c      Inputs:
c         u      - nodal solution vectors
c         b      - FEM Vector
c         a      - FEM diagonal and upper part of array
c         c      - FEM lower part of array
c         afl    - If true assemble 'a' array (which includes 'au')
c         bfl    - If true assemble 'b' array(vector)
c         cfl    - If true assemble 'c' array (array a is unsymmetric)
c         dfl    - If true assemble 'b' uncompressed (e.g.reactions)
c         isw    - Solution switch controlling action to be taken
c         nl1    - First element to be processed
c         nl2    - Last  element to be processed
c         nl3    - Increment to 'nl1'

c      Outputs:
c         b      - Values for FEM Vector
c         a      - Values for FEM diagonal and upper array part
c         c      - Values for FEM lower array part
c
c----------------------------------------------------------------------
c
c.... Declare variable types
      USE cdat1
      USE cdata
      USE edgdat
      USE hdata
      USE hdatam
      USE iofile
      USE mdata
      USE mdat2
      USE ndata
      USE sdata
      USE slid1
      USE slid2
      USE slid3
      USE slid4
      implicit real*8 (a-h,o-z)
      logical afl,bfl,cfl,dfl,rfl,edfl
      character yyy*80
      integer nh1pf,nh2pf,nh3pf
c.... Declare array types
      integer m
      real*8  u(*), b(*), a(*), c(*) 
c
c.... check for contact: reset profile and compute penetrations
c
      rfl = (bfl .and. (.not. dfl)) .and. (isw .ne. 5)
      if (contfl .and. (isw.eq.3 .or. rfl)) then
          edfl = nde*ned.gt.0
          nde1 = nde
          if(nde1.eq.0) nde1=1
          call rstprf(jdt12,eeqn,psid,econ,edge1,edge2,edge3,
     +            nde1,ndf,nen1,nen,neq,numnp,numel,edfl)
          if (cl29(1) .le. 1) then
            call contas(jdt12,psid,u,cl00,cl01,cl02,cl03,
     1            cl33,cl04,cl05,cl06,cl07,cl08,cl34,
     2            cl09,cl35,cl10,cl11,cl12,cl13,cl14,
     3            cl15,cl16,cl17,cl18,cl19,cl20,cl21,
     4            cl22,cl23,cl24,cl25,cl26,cl27,cl28,
     5            cl29,cl30,ndf)
          else
cww                    write(iow,2002) m(l29)
cww         if(ior.lt.0) write(*  ,2002) m(l29)
cww         stop
          write(yyy,2002) cl29(1)
          call drawmess(yyy,1,0)
          return
          end if
      end if
      
!csck removed because local h-arrays are in modules now
!c.... set local h1,h2,h3-arrays
!c.... default  
!      nh1pf = 1
!      nh2pf = 1
!      nh3pf = 1
!c.... in case of existance
!      if(nhmax.gt.0) then
!        nh1pf = (nh1-1)*ipr+1  ! nh1m 
!        nh2pf = (nh2-1)*ipr+1  ! nh2m 
!      end if
!      if(nh3max.gt.0) then
!        nh3pf = (nh3-1)*ipr+1  ! nh3m 
!      end if
      
c.... form appropriate finite element arrays
c
      call pform(edis,ecor,etem,eeqn,epve,ekma,eh1,eh2,
     1 eh3,nmat,edma,psid,coor,econ,gloa,gtem,jdt12,
     2 glo0,u,trans,b,a,c,bang,aang,ndd,nie,ndf,ndm,nen1,nst,
     3 afl,bfl,cfl,dfl,isw,nl1,nl2,nl3)
c
c.... calculate and add the contact tangent stiffness
c
      if (contfl .and. (isw .eq. 3 .or. rfl)) then
          call contat(cl24,cl25,cl08,cl34,cl09,cl35,
     1                cl26,cl27,jdt12,a,b,psid,
     2                cl00,cl01,cl02,cl28,cl30,
     3                cl31,nsl,ndf,neq,afl,rfl)
      end if
c
c.... reset update flag for history variables
c
      hflgu  = .false.
      h3flgu = .false.
c
c.... formats
c
 2002 format('*** INCON *** Unknown option -> ',i3)
c
      end
c
      double precision function gamma1(f0,f,id,u,dr,du,t,s,nqe,nneq)
      USE cdata
      USE fdata
      USE ndata
      USE prlod
      USE tdata
      implicit double precision (a-h,o-z)
      logical fa,tr  
      dimension f0(*),f(*),id(*),u(*),dr(*),du(*),t(*)
      data fa,tr/.false.,.true./
c.... get a search displacement
      nneq2 = nneq + nneq
      do 100 n = 1,nneq
      j = id(n)
      if(j.gt.0) then
        t(n)      = u(n) + s*du(j)
        t(n+nneq) = u(n+nneq) + s*du(j)
        t(n+nneq2)= s*du(j)
      else
        db        = f0(n) + f(n)*prop
        t(n+nneq2)= db - u(n)
        t(n+nneq) = u(n+nneq) -u(n) + db
        t(n)      = db
      end if
100   continue
c.... compute a residual
      call ploads(t,dr,prop,.false.,.false.,.false.)
      call pload(id,f,f0,dr,nneq,prop)
c.... update the residual for lumped mass inertial effects
cww   hflgu/h3flgu for history terms->.false.?
      call formfe(t,dr,dr,dr,fa,tr,fa,fa,6,1,numel,1)
      if(fl(9)) then
        call pmove (trans(nw),t,nneq)
        call pmove (massm,t(nneq+1),nqe)
        ctem = s*c2
        do 110 n = 1,nneq
          j = id(n)
          if(j.gt.0) dr(j) = dr(j) - t(nneq+j)*(t(n)+ctem*du(j))
110     continue
      end if
c.... compute the value of gamma
      gamma1 = ddot(nqe,du,1,dr,1)
      return
      end
c
      subroutine genrsum(nodesrf,x,prt)
c--------------------------------------------------------------------
c
c      Purpose: generate list of nodes for macro rsum
c               node first, node last, incr

c      Inputs:
c         nodesrf     - array for node numbers with reactions 
c         x(ndm,numnp)- nodal coodinates
c         prt         - Print generated data if true
c
c      Outputs:        - nodesrf: array of node numbers for sum of reac
c                      - nfs1:    max. number of nodes  for sum of reac     
c
c      Comments generation of tied nodes is allowed (no influnece)
c--------------------------------------------------------------------
      USE cdata
      USE iofile
      USE rsum
      USE sdata
      implicit double precision (a-h,o-z)
      logical prt
      dimension td(3),nodesrf(*),x(ndm,numnp)
      data blank/-999.d0/
      if(ior.lt.0) write(*,2000)
c ... read dof,type   
      call dinput(td,2)
      ndfrs  = td(1)
      irstyp = td(2)
      if(irstyp.eq.0) irstyp=1
      if (irstyp.gt.2) then
        call drawmess('RSTyp not defined !',1,0)
        STOP
      end if

      if(irstyp.le.1) then 
c ...   read values
101     if(ior.lt.0) write(*,2001)
        call dinput(td,3)
        n1 = td(1)
        n2 = td(2)
        ni = td(3)
        if(n1.eq.0) goto 110
c...    set nodes
        if(ni.eq.0) then  ! only node n1
          nfs1=nfs1+1  
          nodesrf(nfs1) = n1
        else              ! generate n1,n2,ni
          ninc = (n2-n1)/ni + 1   
          do i = 1,ninc 
            n = n1 + (i-1)*ni
            if(n.le.n2) then
              nfs1=nfs1+1 
              nodesrf(nfs1) = n
            end if 
          end do   
        end if
        goto 101
      else
c ...   read values
        if(ior.lt.0) write(*,2002)
        call dinput(td,2)
        idir = td(1)
        x0   = td(2)
c...    find nodes
c....   search circle        
        dx = pdiff(x(idir,1),ndm,numnp)/1000.
        nfs1 = 0
        do i = 1,numnp,1
          if(x(1,i).ne.blank.and.abs(x(idir,i)-x0).le.dx) then
            nfs1 = nfs1 + 1
            nodesrf(nfs1) = i
          end if
        enddo
        goto 110
      end if
        
110   continue
c.... print results
      if(.not.prt) return
            if(ior.le.0 )write(*  ,2010) ndfrs
                   write(iow,2010) ndfrs
        if(ior.le.0) then
          write(*  ,2020) (nodesrf(k),k=1,nfs1)
        else
          write(iow,2020) (nodesrf(k),k=1,nfs1)
        end if
      return
2000  format(' Input: node#, inc.'/3x,'>',$)
2001  format(' Input: dof, type'/3x,'>',$)
2002  format(' Input: direction, xi-value.'/3x,'>',$)
2010  format(/5x,'Nodes for sum of reaction forces for dof  ',i3)       
2020  format(20(5x,'nodes ',10(i8,1x)/))
      end
c
      subroutine genvec(ndm,x,cdd,prt,err,prtz)
c--------------------------------------------------------------------
c
c      Purpose: Generate real data arrays by linear interpolation
c
c      Inputs:
c         ndm       - Number of data items/node
c         cdd(3)    - Header identifier
c         prt       - Print generated data if true
c         prtz      - Print zero data if true
c
c      Outputs:
c         x(ndm,*)  - Generated data
c         err       - Error flag, true if error occurs
c--------------------------------------------------------------------
      USE bdata
      USE cdata
      USE errchk
      USE inptc
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt,err,prtz
      character*4 cdd*12,cd(3),yyy*160
      dimension x(ndm,*),xl(100),td(16)
      data bl/-999.d0/
      save ipos

      if(inptctrl.eq.1.and.numnp.ne.0.and.cdd.eq.' coordinates') 
     +   call perform2('coor',1,0,0)
      
c.... xl(ndm) is necessary
      if(ndm.gt.100) then
        call drawmess('xl in genvec is set only to 100 values',1,0)
        return
      end if
      mct = 0
      n   = 0
      ng  = 0
      do 50 i = 1,3
  50  cd(i) = cdd(4*(i-1)+1:4*(i-1)+4)
100   l   = n
      lg  = ng

      if(inptctrl.eq.1.and.numnp.ne.0.and.cdd.eq.' coordinates') then
        if(l.gt.ipos*200) then ! plot actual state every 200 coordinates 
          ipos = ipos+1 
          call perform2('coor',2,l,numnp)
        end if
      end if
      
c.... call input routine - values returned in td and then moved
101   if(ior.lt.0) write(*,2010)
      il = min(ndm+2,16)
      call dinput(td,il)
      if(errck) go to 101
      n   = td(1)
      ng  = td(2)
      do 102 i = 1,min(ndm,14)
        xl(i) = td(i+2)
102   continue
      if(ndm.gt.14) then
        do 205 ii = 1,(ndm+2)/16
          is = il+1
          il = min(is+15,ndm+2)
203       call dinput(td,il-is+1)
          if(errck) go to 203
          do 204 k = 1,il-is+1
            xl(k+is-3) = td(k)
204       continue
205     continue
      end if
cww   if(n.gt.numnp.and.ior.gt.0) write(iow,3001) n,(cd(i),i=1,3)
cww   if(n.gt.numnp.and.ior.lt.0) write(  *,3001) n,(cd(i),i=1,3)
cww   if(n.gt.numnp) then
cww      write(yyy,3001) n,(cd(i),i=1,3)
cww      call drawmess(yyy,1,0)
cww      err = .true.
cww   end if
      if(n.gt.numnp) then
        write(yyy,3001) n,(cd(i),i=1,3),numnp
        call drawmess(yyy,1,0)
        n = numnp
      end if
      if(n.le.0.or.n.gt.numnp) go to 109
      do 103 i = 1,ndm
103   x(i,n) = xl(i)
      if(lg) 104,100,104
104   lg = isign(lg,n-l)
      xli =(iabs(n-l+lg)-1)/iabs(lg)
      do 105 i = 1,ndm
105   xl(i) = (x(i,n)-x(i,l))/xli
106   l = l + lg
      if((n-l)*lg.le.0) go to 100
      if(l.le.0.or.l.gt.numnp) go to 108
      do 107 i = 1,ndm
107   x(i,l) = x(i,l-lg) + xl(i)
      go to 106
108   continue
cww   if(ior.gt.0) write(iow,3000) l,(cd(i),i=1,3)
cww   if(ior.lt.0) write(  *,3000) l,(cd(i),i=1,3)
      write(yyy,3000) l,(cd(i),i=1,3)
      call drawmess(yyy,1,0)
      err = .true.
      go to 100

109   if(inptctrl.eq.1.and.numnp.ne.0.and.cdd.eq.' coordinates') 
     +   call perform2('coor',3,0,0)
   
      if(.not.prt) return

      do 113 j = 1,numnp
        if(prtz) go to 111
        do 110 l = 1,ndm
          if(x(l,j).ne.0.0d0) go to 111
110     continue
        go to 113
111     mct = mct - 1
        if(mct.gt.0) go to 112
        mct = 50
        if(ior.gt.0) then
cww       write(iow,2000) o,head,(cd(l),l=1,3),(l,cd(1),cd(2),l=1,ndm)
          write(iow,2000)        (cd(l),l=1,3),(l,cd(1),cd(2),l=1,ndm)
        else
cww       write(  *,2000) o,head,(cd(l),l=1,3),(l,cd(1),cd(2),l=1,ndm)
          write(  *,2000)        (cd(l),l=1,3),(l,cd(1),cd(2),l=1,ndm)
        end if
112     if(ior.gt.0) then
          if(x(1,j).eq.bl) write(iow,2008) j
          if(x(1,j).ne.bl) write(iow,2009) j,(x(l,j),l=1,ndm)
        else
          if(x(1,j).eq.bl) write(  *,2008) j
          if(x(1,j).ne.bl) write(  *,2009) j,(x(l,j),l=1,ndm)
        end if
113   continue
      return
cww2000  format(a1,19a4,a3//5x,'nodal',3a4//6x,'node',6(i7,a4,a2)/
2000  format(/5x,'nodal',3a4//6x,'node',6(i7,a4,a2)/
     1       (10x,6(i7,a4,a2)))
2008  format(i10,' has not been input or generated')
2009  format(i10,6f13.5:/(10x,6f13.5))
2010  format(' Input: node#, inc., values'/3x,'>',$)
cww3000  format(' **ERROR** attempt to generate node',i5,' in '
cww     1       ,3a4)
3000  format(' Try to generate node',i5,' for macro ',3a4)
cww3001  format(' **ERROR** attempt to input node',i5,', terminate',
cww     1 ' input of nodes in ',3a4)
3001  format(' Try to make input values for node',i5,' for macro ',3a4,
     +       ' !   Input continued with final node',i5)
      end
c
      subroutine genvec1(ndm,x,u,cdd,cdi,prt,err,prtz)
c--------------------------------------------------------------------
c
c      Purpose: Generate real data arrays by linear interpolation
c               same as genvec! for imperfections X = X +(!) alpha*phi
c               furtermore: store u=alpha*phi 
c
c      Inputs:
c         ndm       - Number of data items/node
c         cdd(3)    - Header identifier 1 
c         cdi(3)    - Header identifier 2
c         prt       - Print generated data if true
c         prtz      - Print zero data if true
c
c      Outputs:
c         x(ndm,*)  - Generated data  x=x+alpha*phi
c         u(ndm,*)  - Generated data  u=  alpha*phi
c         err       - Error flag, true if error occurs
c
c--------------------------------------------------------------------
      USE bdata
      USE cdata
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt,err,prtz
      character*4 cdd*12,cdi*12,cd(3),cd1(3),yyy*80
      dimension x(ndm,*),u(ndm,*),xl(ndm),ul(ndm),td(16)
      data bl/-999.d0/
      mct = 0
      n   = 0
      ng  = 0
      do 50 i = 1,3
      cd1(i) = cdi(4*(i-1)+1:4*(i-1)+4)
  50   cd(i) = cdd(4*(i-1)+1:4*(i-1)+4)
c.....BEGIN NEW
c.... read factor
      call dinput(td,1)
      alpha = td(1)
      if(alpha.eq.0.d0) return
c.... END NEW
100   l   = n
      lg  = ng
c.... call input routine - values returned in td and then moved
101   if(ior.lt.0) write(*,2010)
      il = min(ndm+2,16)
      call dinput(td,il)
      if(errck) go to 101
      n   = td(1)
      ng  = td(2)
      do 102 i = 1,min(ndm,14)
        xl(i) = td(i+2)
102   continue
      if(ndm.gt.14) then
        do 205 ii = 1,(ndm+2)/16
        is = il+1
        il = min(is+15,ndm+2)
203     call dinput(td,il-is+1)
        if(errck) go to 203
        do 204 k = 1,il-is+1
           xl(k+is-3) = td(k)
204     continue
205     continue
      end if
cww   if(n.gt.numnp.and.ior.gt.0) write(iow,3001) n,(cd(i),i=1,3)
cww   if(n.gt.numnp.and.ior.lt.0) write(  *,3001) n,(cd(i),i=1,3)
      if(n.gt.numnp) then
        write(yyy,3001) n,(cd(i),i=1,3)
        call drawmess(yyy,1,0)
        err = .true.
      end if
      if(n.le.0.or.n.gt.numnp) go to 109
      do 103 i = 1,ndm
cww103   x(i,n) = xl(i) change to genvec
      u(i,n) =        + alpha*xl(i) ! new
103   x(i,n) = x(i,n) + alpha*xl(i) ! new
      if(lg) 104,100,104
104   lg = isign(lg,n-l)
      xli =(iabs(n-l+lg)-1)/iabs(lg)
      do 105 i = 1,ndm
      ul(i) = (u(i,n)-u(i,l))/xli
105   xl(i) = (x(i,n)-x(i,l))/xli
106   l = l + lg
      if((n-l)*lg.le.0) go to 100
      if(l.le.0.or.l.gt.numnp) go to 108
      do 107 i = 1,ndm
      u(i,l) = u(i,l-lg) + ul(i)
107   x(i,l) = x(i,l-lg) + xl(i)
      go to 106
108   continue
cww   if(ior.gt.0) write(iow,3000) l,(cd(i),i=1,3)
cww   if(ior.lt.0) write(  *,3000) l,(cd(i),i=1,3)
      write(yyy,3000) l,(cd(i),i=1,3)
      call drawmess(yyy,1,0)
      err = .true.
      go to 100
109   if(.not.prt) return
      do 113 j = 1,numnp
      if(prtz) go to 111
      do 110 l = 1,ndm
      if(x(l,j).ne.0.0d+00) go to 111
110   continue
      go to 113
111   mct = mct - 1
      if(mct.gt.0) go to 112
      mct = 50
      if(ior.gt.0) then
cww    write(iow,2000) o,head,(cd(l),l=1,3),(l,cd(1),cd(2),l=1,ndm)
       write(iow,2000) (cd(l),l=1,3),(cd1(l),l=1,3),
     +               (l, cd(1), cd(2),l=1,ndm),(l,cd1(1),cd1(2),l=1,ndm)
      else
cww    write(  *,2000) o,head,(cd(l),l=1,3),(l,cd(1),cd(2),l=1,ndm)
       write(  *,2000)(cd(l),l=1,3),(cd1(l),l=1,3),
     +               (l, cd(1), cd(2),l=1,ndm),(l,cd1(1),cd1(2),l=1,ndm)
      end if
cww   print not phi but X+alpha*phi !!
112   if(ior.gt.0) then
          if(x(1,j).eq.bl) write(iow,2008) j
          if(x(1,j).ne.bl) write(iow,2009) j,(x(l,j),l=1,ndm),
     +                                       (u(l,j),l=1,ndm)
      else
          if(x(1,j).eq.bl) write(  *,2008) j
          if(x(1,j).ne.bl) write(  *,2009) j,(x(l,j),l=1,ndm),
     +                                       (u(l,j),l=1,ndm)
      end if
113   continue
      return
cww2000  format(a1,19a4,a3//5x,'nodal',3a4//6x,'node',6(i7,a4,a2)/
2000  format(/5x,'nodal',3a4,3a4//6x,'node',12(i7,a4,a2)/
     1       (10x,12(i7,a4,a2)))
2008  format(i10,' has not been input or generated')
2009  format(i10,12f13.4:/(10x,6f13.4))
2010  format(' Input: node#, inc., values'/3x,'>',$)
cww3000  format(' **ERROR** attempt to generate node',i5,' in '
cww     1       ,3a4)
3000  format(' Try to generate node',i5,' for macro ',3a4)
cww3001  format(' **ERROR** attempt to input node',i5,', terminate',
cww     1 ' input of nodes in ',3a4)
3001  format(' Try  to input values for node',i5,' for macro ',3a4)
      end
c
c
      subroutine genvec2(ndm,x,cdd,prt,err,prtz)
c--------------------------------------------------------------------
c
c      Purpose: Generate real data arrays by linear interpolation
c            same as genvec1! for random imperfections X = X+fact*random
c
c      Inputs:
c         ndm       - Number of data items/node
c         cdd(3)    - Header identifier
c         prt       - Print generated data if true
c         prtz      - Print zero data if true
c
c      Outputs:
c         x(ndm,*)  - Generated data
c         err       - Error flag, true if error occurs
c
c--------------------------------------------------------------------
      USE bdata
      USE cdata
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt,err,prtz
      character*4 cdd*12,cd(3),yyy*80
      dimension x(ndm,*),td(16)
      data bl/-999.d0/
      mct = 0
      do i = 1,3
        cd(i) = cdd(4*(i-1)+1:4*(i-1)+4)
      end do
c.... read factor
      call dinput(td,2)
      fact  = td(1)
      ntype = td(2)
      if(ior.gt.0) then
        write(iow,2001) fact,ntype
      else
        write(  *,2001) fact,ntype
      end if
      if(fact.eq.0.d0) return
c
c.... call input routine - values returned in td and then moved
100   if(ior.lt.0) write(*,2010)
      call dinput(td,3)
      if(errck) go to 100
      ia  = td(1)
      ie  = td(2)
      ng  = td(3)
c      
      if(ia.lt.0.or.ia.gt.numnp.or.ie.lt.0.or.ie.gt.numnp) then
        write(yyy,3001) n,(cd(i),i=1,3)
        call drawmess(yyy,1,0)
        return
      end if
      if(ia.le.0.or.ia.gt.numnp) go to 109
c.... modify coordinates
      if(ntype.ne.2) call random_seed() 
      do  n = ia,ie,ng
        do  i = 1,ndm
          call random_number(rnum)
           x(i,n)=x(i,n)+fact*(rnum-0.5d0)*2.d0  
        end do
      end do
      go to 100
109   if(.not.prt) return
      do 113 j = 1,numnp
        if(prtz) go to 111
        do 110 l = 1,ndm
          if(x(l,j).ne.0.0d+00) go to 111
110     continue
      go to 113
111   mct = mct - 1
      if(mct.gt.0) go to 112
      mct = 50
      if(ior.gt.0) then
       write(iow,2000)        (cd(l),l=1,3),(l,cd(1),cd(2),l=1,ndm)
      else
       write(  *,2000)        (cd(l),l=1,3),(l,cd(1),cd(2),l=1,ndm)
      end if
c     print  X+fact*random
112   if(ior.gt.0) then
          if(x(1,j).eq.bl) write(iow,2008) j
          if(x(1,j).ne.bl) write(iow,2009) j,(x(l,j),l=1,ndm)
      else
          if(x(1,j).eq.bl) write(  *,2008) j
          if(x(1,j).ne.bl) write(  *,2009) j,(x(l,j),l=1,ndm)
      end if
113   continue
      return
2000  format(/5x,'nodal',3a4//6x,'node',6(i7,a4,a2)/
     1       (10x,6(i7,a4,a2)))
2001  format(/5x,'Modifaction of coordinates with random numbers',/
     1        5x,' factor ',g12.5,/,5x,'  ntype ',i2)  
2008  format(i10,' has not been input or generated')
2009  format(i10,6f13.4:/(10x,6f13.4))
2010  format(' Input: node1, node2, inc.'/3x,'>',$)
3001  format(' Try  to input values for node',i5,' for macro ',3a4)
      end
c
      subroutine genvec3(ndm,x,prt)
c--------------------------------------------------------------------
c
c      Purpose: Generate imperfect coordinates arrays
c                tolc can be set here
c                tie: tolc=1000
c
c      Inputs:
c         ndm       - Number of data items/node
c         prt       - Print generated data if true
c
c      Outputs:
c         x(ndm,*)  - Generated data
c
c--------------------------------------------------------------------
      USE cdata
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      dimension x(ndm,*),td(16),xp(3),xip(3),dx(3),tol(3)
      data tolc /250/ 
c
      np  = 0 
c
c.... set tolerance (from tie)
      toln = numnp
      toln = sqrt(toln)*tolc
c.... calculate distance tol of mesh
      call pzero(tol,3)
      do i = 1,ndm
        tol(i) = pdiff(x(i,1),ndm,numnp)/toln
      end do
c
      if(prt)              write(iow,2000)  
      if(prt.and.ior.lt.0) write(  *,2000)
c        
c...  read data
100   if(ior.lt.0) write(*,2010)
      call dinput(td,2*ndm)
      xp(1)   = td(1)
      xp(2)   = td(2)
      xp(3)   = 0.d0
      if(ndm.eq.3) xp(3) = td(3)
      xip(1)  = td(ndm+1)
      xip(2)  = td(ndm+2)
      xip(3)  = 0.d0
      if(ndm.eq.3) xip(3) = td(ndm+3)
      xpq  = dot(xp,xp,3) 
      xipq = dot(xip,xip,3) 
      if(xpq.eq.0.d0.and.xipq.eq.0.d0) then
                     write(iow,2002) np
        if(ior.lt.0) write(  *,2002) np 
        return
      end if
c      
c.... find nodal point     
      ip = 0
      do n=1,numnp
c....   check point          
        dx(1) = dabs(x(1,n)-xp(1))
        dx(2) = dabs(x(2,n)-xp(2))
        dx(3) = 0.d0
        if(ndm.eq.3) dx(3) = dabs(x(3,n)-xp(3))
        if(dx(1).le.tol(1).and.dx(2).le.tol(2).and.dx(3).le.tol(3)) then 
c....     modify and print coordinates
          np = np+1 
          if(prt) then
            if(ior.lt.0) then
              write(  *,2009) n,(x(l,n),l=1,ndm),(xip(l),l=1,ndm)
            else
              write(iow,2009) n,(x(l,n),l=1,ndm),(xip(l),l=1,ndm)
            end if
          end if
          x(1,n)   = xip(1)
          x(2,n)   = xip(2)
          if(ndm.eq.3) x(3,n) = xip(3)
          ip= 1
        end if
      end do
      if(ip.ne.1) then
                     write(iow,2001) (xp(l),l=1,ndm)
        if(ior.lt.0) write(  *,2001) (xp(l),l=1,ndm)
      end if
      goto 100
c

2000  format(/5x,'Imperfect nodal coordinates',/,
     +        7x,'Node','  old coordinates       new coordinates')
2001  format(5x,' Point not found, coordinates=:  ',3e14.6)   
2002  format(5x,' No of points found:  ',i10)   
2009  format(i10,6e14.6)
2010  format(' Input: x,y,z perfect, x,y,z imperf.'/3x,'>',$)
      end
c
      subroutine ini
c--------------------------------------------------------------------
c
c      Purpose: set initial values (e.g. start adresses for common m)
c
c      Inputs:
c
c      Outputs:
c      flg( 1)  - F 
c      flg( 2)  - F
c      flg( 3)  - T
c      flg( 4)  - T
c      flg( 5)  - T
c      flg( 6)  - T
c      flg( 7)  - T
c      flg( 8)  - F
c      flg( 9)  - F
c      flg(10)  - F
c      flg(11)  - F
c      flg(12)  - F
c
c
c
c--------------------------------------------------------------------
      USE arcl
      USE comfil
      USE ddata
      USE dspos
      USE dtauto
      USE edgdat
      USE epsdh
      USE errnam
      USE ext2
      USE fdata
      USE fe2dat
      USE feapmpi
      USE fileno
      USE fornam
      USE hdatam
      USE iosave
      USE isprec
      USE mate
      USE mxasz
      USE ndata
      USE nolink
      USE pdata1
      USE pdata12
      USE pindex
      USE plodf
      USE plodfa
      USE plodfs
      USE plodfu
      USE plotter
      USE pltran
      USE pnodn
      USE ppers
      USE prlod
      USE proc
      USE psethis
      USE ptext
      USE qload
      USE rsum
      USE smpak
      USE stepc
      USE strnam
      USE subdt
      USE ximp
      implicit real*8(a-h,o-z)
c
      common /colini/ icolold
      common /curvedat/cpar(20,8),nrt1(20,3),ic,nbe,nn3  ! for curves
cww      common /easdata/ ieas   
cww      common /el86a/  nk1  !only c.kne for sdis
c
c.... set initial values defined by user
      call iniuser
c
      do 50 i = 1,11
          if(i.ge.3 .and. i.le.7) then
            fl(i)  = .true.
          else
            fl(i)  = .false.
          end if
50    continue
c
      fl(1)   = .false.
      fl(9)   = .false.
      lsave   = .false.
      lodrv   = .false.    ! sm-solver
c.... history
      hflgu   = .false.
      h3flgu  = .false.
c
cww      nk1  = 1    ! added knebel for sdis
      nrk  = 0    ! dynamic
      nrc  = 0
      nrm  = 0
      nrt  = 0
c.... set parameter for arclength step control
      arcfs   = .false.
      itd  = 6       
      cmax = 1.d0   
      sp   = 0.5d0    
      rm   = 0.5d0    
c.... preconditioner
      lpreco = .true. 
      ippc   = 3
      tolpc  = 1.e-3
      lfil   = 5
c.... extended system
      extflg = .false.
      nmode  = 0
      eps    = 1.d-7
c.... micro problem
      flgfe2 = .true.
c
      nren = 0
      mf   = 0
      detc = 0.d0
      prop = 1.d0
      nv   = 1
      na   = 1
      nal  = 1
      nau  = 1
      nl   = 1   
      nm   = 1   
c.... set parameter for epsd.h
      factx = 0.d0 
      facty = 0.d0 
c.... for graphic
      iso = .false.
      kpers = 0
      ipla  = 0
      call pzero(rotang,3)
      call pzero(tra,9)
      call pzero(vr,3)
c.... for plotting file
      call pzero(xxp,200)
      call pzero(yyp,200)
      icolold = 0
c.... forc and stress description
      do i = 1,11
        forsus(i)=' '
      end do
      do i = 1,25
        strsus(i)=' '
      end do
c.... for curves
      ic  = 0
      nbe = 1
      call pzero ( cpar,20*8 ) 
      call pzeroi( nrt1,20*3 )
c.... for tplo
      npldf1=0
      call pzeroi( noden,10 )
      nploc =1      ! factor for tplo-convergence  0= plot 1=no plot
      valuse1 = 0.d0
      valuse2 = 0.d0
c.... for perform only from pmacr(dos)
      ijump = 0
c.... procedure path (default = path of input file
      procpath = ' '
      ipath = ipos1(finp,229)
      procpath(1:ipath) = finp(1:ipath)
c.... initial number for save file.isno
      isno = 0
cwwc.... initial value for EAS usage 0 = not used
cww      ieas = 0
c.... set text for drawmess
      text1=' '
      text2=' '
c.... set parameter for AUTO
      htol =0.d0
      dtmax=0.d0
      dtdo  = 0.85d0
      dtup1 = 0.80d0
      dtup2 = 1.25d0
      iaback= 0
c...  for MPI
      parform = .true.
c...  set all materials to be plotted
c      do i = 1,40
c        ipma(i) = i
c      end do
c
c.... set parameter for LINK
      nli1=0
c.... for rsum
      nfs1=0
c.... newmesh for POST
      lindex = .true. 
c     /psethis/
      kpset=0

c.... set parameter for ixtie: double ix array for parview
      ixtie=0

c.... set parameter for velo+acce array for parview
      flparv = .true.

c.... set parameter for impf
      flimp = .true.
      mimp  = 1

c.... set parameter for plodfs Tplo, stress nstri at node nstrno 
      nstri  = 0
      nstrno = 0      
      
c.... set parameter for dspos        
      dscor = 0
      
c.... set parameter for mate
      flmat = .true.

c.... description of erro indicators
      e_name(1)  = 'Energy norm    '
      e_name(2)  = 'L_2-norm S     '
      e_name(3)  = 'Y0             '

      return
      end
c
      subroutine injint(njint,nj,prt)
c--------------------------------------------------------------------
c
c      Purpose: read element numbers for J-integral
c
c      Inputs:
c
c      Outputs:
c
c--------------------------------------------------------------------
      USE iofile
      integer nj(njint)
      real*8 td(8)
      logical prt
c.... split input into rows of 8
      nimax = njint/8 + 1
      if(prt) write(iow,2000)
      do 100 i = 1, nimax
        call dinput(td,8)
        do 110 j = 1,8
          nj(8*(i-1) + j) = td(j)
110     continue
        if(prt) write(iow,2001) (nj(8*(i-1) + j),j=1,8)
100   continue
2000  format(//,5x,'Element numbers for J-integral',/)
2001  format(5x,8i6)
      return
      end
c
      subroutine cjint(u,dr,nj,njint)
c--------------------------------------------------------------------
c
c      Purpose: compute J-integral
c
c      Inputs:
c
c      Outputs:
c
c--------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      integer nj(njint)
      dimension u(*),dr(*)
      logical fa
      data fa /.false./
      do 100 i=1,njint
        nelem = nj(i)
cww   hflgu/h3flgu for history terms?
        call formfe(u,dr,dr,dr,fa,fa,fa,fa,16,nelem,nelem,1)
100   continue
      return
      end
c
      subroutine int2d(l,lint,sg)
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form Gauss points and weights for two dimensions

c      Inputs:
c         l       - Number of points/direction

c      Outputs:
c         lint    - Total number of points
c         sg(3,*) - Array of points and weights
c-----[--.----+----.----+----.-----------------------------------------]
      USE eldata
      implicit  none

      integer   i,j,k,l,lint
      integer   lr(9),lz(9),lw(9)

      real*8    g,h, third
      real*8    sg(3,*),ss(5),ww(5)

      save

      data      lr/-1,1,1,-1,0,1,0,-1,0/,lz/-1,-1,1,1,-1,0,1,0,0/
      data      lw/4*25,4*40,64/
      data      third / 0.3333333333333333d0 /

c     Set number of total points

      lint = l*l

c     5 pt. integration

      if(l.eq.0) then

        lint = 5
        g    = sqrt(0.6d0)
        do i = 1,4
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = 5.d0/9.d0
        end do

        sg(1,5) = 0.0d0
        sg(2,5) = 0.0d0
        sg(3,5) = 16.d0/9.d0

c     1x1 integration

      else if(l.eq.1) then
        sg(1,1) = 0.d0
        sg(2,1) = 0.d0
        if(nel.eq.3) sg(2,1) = -third
        sg(3,1) = 4.d0

c     2x2 integration

      else if(l.eq.2) then
        g = sqrt(third)
        do i = 1,4
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = 1.d0
        end do

c     3x3 integration

      else if(l.eq.3) then
        g = sqrt(0.6d0)
        h = 1.d0/81.d0
        do i = 1,9
          sg(1,i) = g*lr(i)
          sg(2,i) = g*lz(i)
          sg(3,i) = h*lw(i)
        end do

c     4x4 integration

      else if(l.eq.4) then
        g     = sqrt(4.8d0)
        h     = third/g
        ss(1) = sqrt((3.d0+g)/7.d0)
        ss(4) = - ss(1)
        ss(2) = sqrt((3.d0-g)/7.d0)
        ss(3) = -ss(2)
        ww(1) = 0.5d0 - h
        ww(2) = 0.5d0 + h
        ww(3) = 0.5d0 + h
        ww(4) = 0.5d0 - h
        i = 0
        do j = 1,4
          do k = 1,4
            i = i + 1
            sg(1,i) = ss(k)
            sg(2,i) = ss(j)
            sg(3,i) = ww(j)*ww(k)
          end do
        end do

c     5x5 integration

      else if(l.eq.5) then

        g     =  sqrt(1120.d0)
        ss(1) =  sqrt((70.d0 + g)/126.d0)
        ss(2) =  sqrt((70.d0 - g)/126.d0)
        ss(3) =  0.0d0
        ss(4) = -ss(2)
        ss(5) = -ss(1)

        ww(1) =  (21.d0*g + 117.6d0)/(g*(70.d0 + g))
        ww(2) =  (21.d0*g - 117.6d0)/(g*(70.d0 - g))
        ww(3) =  2.d0*(1.d0 - ww(1) - ww(2))
        ww(4) =  ww(2)
        ww(5) =  ww(1)

        i = 0
        do j = 1,5
          do k = 1,5
            i = i + 1
            sg(1,i) = ss(k)
            sg(2,i) = ss(j)
            sg(3,i) = ww(j)*ww(k)
          end do
        end do

      endif

      end
c
      subroutine int3d(ll,lint,s)
c--------------------------------------------------------------------
c      Purpose: Gauss quadrature for 3-d element

c      Inputs:
c         ll     - Number of points/direction

c      Outputs:
c         lint   - Total number of quadrature points
c         s(4,*) - Gauss points (1-3) and weights (4)
c
c      Comments: 
c       ll=1:  1         pt. quadrature
c       ll=2:  2 x 2 x 2 pt. quadrature
c       ll=3:  3 x 3 x 3 pt. quadrature
c       ll=4:  4 x 4 x 4 pt. quadrature
c       ll=5:  5 x 5 x 5 pt. quadrature
c       ll=9:  special 9 pt. quadrature Simo, Armero, Taylor 1993 
c       ll=10: special 4 pt. quadrature Ref=? 
c
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension s(4,*),ig(4),jg(4),sw(2,5)
cww      dimension xi3(3),wg3(3)
      data ig/-1,1,1,-1/,jg/-1,-1,1,1/
cww      xi3(1) = -dsqrt(0.6d0)
cww      xi3(2) =  0.d0
cww      xi3(3) =  dsqrt(0.6d0) 
c

      sqt13  = 1.d0/sqrt(3.0d0)
      sqtp6  = sqrt(0.6d0)
      five9  = 5.d0/9.d0
      thty29 = 32.d0/9.d0


c     if(ll.eq.1) then
cc....   1 pt. quadrature
c       lint = 1
c       do i = 1,3
c         s(i,1) = 0.0d0
c       end do ! i
c       s(4,1) = 8.0d0
c
c     else if(ll.eq.2) then
cc....   2 x 2 x 2 pt. quadrature
cc....   first bottom(as nodes), then top (as nodes)
c       lint = 8
c       g    = sqt13
c       do i = 1,4
c         s(1,i)   = ig(i)*g ! bot
c         s(1,i+4) = s(1,i)  ! top
c         s(2,i)   = jg(i)*g ! bot
c         s(2,i+4) = s(2,i)  ! top
c         s(3,i)   =  g      ! bot
c         s(3,i+4) = -g      ! top
c         s(4,i)   = 1.d0
c         s(4,i+4) = 1.d0
c       end do ! i
c
c      else if(ll.eq.3) then
cc....   3 x 3 x 3 pt. quadrature
cc....   first xsi, then eta, then zeta
c        lint = 27
c        wg3(1) = 5.d0/9.d0
c        wg3(2) = 8.d0/9.d0
c        wg3(3) = 5.d0/9.d0
c        kg = 0
c        do i=1,3
c          do j=1,3
c            do k=1,3
c              kg=kg+1
c              s(1,kg)=xi3(k)
c              s(2,kg)=xi3(j)
c              s(3,kg)=xi3(i)
c              s(4,kg)=wg3(i)*wg3(j)*wg3(k)
c            end do  
c          end do
c        end do  
c

      if(ll.le.5) then
c....   ll x ll x ll pt. quadrature
        call int1d(ll,sw)
        lint = 0
        do k = 1,ll
          do j = 1,ll
            do i = 1,ll
              lint = lint + 1
              s(1,lint) = sw(1,i)
              s(2,lint) = sw(1,j)
              s(3,lint) = sw(1,k)
              s(4,lint) = sw(2,i)*sw(2,j)*sw(2,k)
            end do ! i
          end do ! j
        end do ! k

      else if(ll.eq.9) then
c....   9 pt. quadrature Simo, Armero, Taylor 1993 
c....   first bottom(as nodes), then top (as nodes)
c....   outer points like 3x3x3 + central point  
        lint = 9
        g    = sqtp6
        do i = 1,4
          s(1,i)   = ig(i)*g ! bot
          s(1,i+4) = s(1,i)  ! top
          s(2,i)   = jg(i)*g ! bot
          s(2,i+4) = s(2,i)  ! top
          s(3,i)   =  g      ! bot
          s(3,i+4) = -g      ! top
          s(4,i)   = five9
          s(4,i+4) = five9
        end do ! i
        s(1,9)     =  0.d0
        s(2,9)     =  0.d0
        s(3,9)     =  0.d0
        s(4,9)     =  thty29

      else if(ll.eq.10) then
c....   special 4 pt. quadrature
c....   bottom: on diagonal 1-3, top on diagonal 4-2,  
        lint = 4
        g    = sqt13
        do i = 1,4
          s(1,i) = ig(i)*g
          s(2,i) = s(1,i)
          s(3,i) = jg(i)*g
          s(4,i) = 2.0d0
        end do ! i
        s(2,3) = -g
        s(2,4) =  g

      else
        stop 'Error in SR INT3D: Int.Order unknown'
      end if
      return
      end
c
      subroutine int1d(l,sw)
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Gauss quadrature for 1-d element

c      Inputs:
c         l     - Number of points

c      Outputs:
c         sw(1,*) - Gauss points
c         sw(2,*) - Gauss weights
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

cww      include  'pconstant.h'

      integer   l
      real*8    sw(2,*), t
      real*8    sqt13,sqtp6,sqt48,five9,eight9,one3

      save

      sqt13  = 1.d0/sqrt(3.d0)  
      sqtp6  = sqrt(0.6d0)  
      sqt48  = sqrt(4.8d0)  
      five9  = 5.d0/9.d0 
      eight9 = 8.d0/9.d0 
      one3   = 1.d0/3.d0 

      if(l.eq.1) then

        sw(1,1) = 0.0d0
        sw(2,1) = 2.0d0

      else if(l.eq.2) then

        sw(1,1) = -sqt13                     ! sqrt(1/3)
        sw(1,2) = -sw(1,1)
        sw(2,1) = 1.0d0
        sw(2,2) = 1.0d0

      else if(l.eq.3) then

        sw(1,1) = -sqtp6                    ! sqrt(0.6)
        sw(1,2) = 0.0d0
        sw(1,3) = -sw(1,1)
        sw(2,1) = five9
        sw(2,2) = eight9
        sw(2,3) = sw(2,1)

      else if(l.eq.4) then

        t       =  sqt48                    ! sqrt(4.8)
        sw(1,1) = -sqrt((3.d0+t)/7.d0)
        sw(1,2) = -sqrt((3.d0-t)/7.d0)
        sw(1,3) = -sw(1,2)
        sw(1,4) = -sw(1,1)
        t       =  one3/t                   ! 1/3/t
        sw(2,1) =  0.5d0 - t
        sw(2,2) =  0.5d0 + t
        sw(2,3) =  sw(2,2)
        sw(2,4) =  sw(2,1)

      else if(l.eq.5) then

        t       =  1120.0d0
        t       =  sqrt(t)

        sw(1,1) = (70.d0+t)/126.d0
        sw(1,2) = (70.d0-t)/126.d0

        t       =  1.d0/(15.d0 * (sw(1,2) - sw(1,1)))

        sw(2,1) = (5.0d0*sw(1,2) - 3.0d0)*t/sw(1,1)
        sw(2,2) = (3.0d0 - 5.0d0*sw(1,1))*t/sw(1,2)
        sw(2,3) =  2.0d0*(1.d0 - sw(2,1) - sw(2,2))
        sw(2,4) =  sw(2,2)
        sw(2,5) =  sw(2,1)

        sw(1,1) = -sqrt(sw(1,1))
        sw(1,2) = -sqrt(sw(1,2))
        sw(1,3) =  0.0d0
        sw(1,4) = -sw(1,2)
        sw(1,5) = -sw(1,1)

c     Compute points and weights

c      else

c        call gausspw(l,sw)

      endif

      end
c
      subroutine int3dt(ll,lint,s)
c----------------------------------------------------------------------
c      Purpose: Gauss quadrature for 3-d tetrahedon

c      Inputs:
c         ll     - Number of points/direction

c      Outputs:
c         lint   - Total number of quadrature points
c         s(4,*) - Gauss points (1-3) and weights (4)
c
c      Comments: 
c       ll=1:  1 pt. quadrature
c       ll=4:  4 pt. quadrature
c       ll=5:  5 pt. quadrature   neg. weight !!
c       ll=6:  6 pt. quadrature   added ww from Smith/Griffiths
c
c     P. Wriggers, Feb. 1992, THD
c     Checked WW BS KIT 01/2010 
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension s(4,*),a1(4),a2(4),a3(4),
     1                 b1(4),b2(4),b3(4)
      data a1/1.,1.,1.,0./,a2/1.,1.,0.,1./,a3/1.,0.,1.,1./
      data b1/0.,0.,0.,1./,b2/0.,0.,1.,0./,b3/0.,1.,0.,0./
      if(ll.eq.1) then
c....   1 pt. quadrature
        lint = 1
        do i = 1,3
          s(i,1) = 0.25d0
        end do
        s(4,1) = 1.d0/6.d0
      else if(ll.eq.4) then
c....   4 pt. quadrature
        lint = 4
        a = (5.d0 - sqrt(5.d0))/20.d0
cww     b = (5.d0 - 3.d0*sqrt(5.d0))/20.d0
        b = (5.d0 + 3.d0*sqrt(5.d0))/20.d0
        c = 1.d0/24.d0
        s(1,1) = a
        s(2,1) = a
        s(3,1) = a
        s(4,1) = c
        do i = 2,4
          s(1,i) = a*a1(i) + b*b1(i)
          s(2,i) = a*a2(i) + b*b2(i)
          s(3,i) = a*a3(i) + b*b3(i)
          s(4,i) = c
        end do
      else if(ll.eq.5) then
c....   5 pt. quadrature
        lint = 5
        d = 0.25d0
        a = 1.d0/6.d0
        b = 0.5d0
        c = 3.d0/40.d0
        s(1,1) = d
        s(2,1) = d
        s(3,1) = d
        s(4,1) = -2.d0/15.d0  ! negative weight for Point 1!!
        do i = 1,4
          s(1,i+1) = a*a1(i) + b*b1(i)
          s(2,i+1) = a*a2(i) + b*b2(i)
          s(3,i+1) = a*a3(i) + b*b3(i)
          s(4,i+1) = c
        end do
      else if(ll.eq.6) then
c....   6 pt. quadrature
        lint = 6
        call pzero(s,24)
        s(1,1) = -1.d0
        s(1,2) =  1.d0
        s(2,3) = -1.d0
        s(2,4) =  1.d0
        s(3,5) = -1.d0
        s(3,6) =  1.d0
        c = 1.d0/36.d0
        do i = 1,6
          s(4,i) = c
        end do  
      else
        stop 'Error in SR INT3DT: Int.Order unknown'
      end if
      return
      end
c
      subroutine intio(y,wd1,wd2,v,error)
c----------------------------------------------------------------------
c
c      Purpose: read a FEAP input macro for PMACIO>PARSE

c      Inputs:
c         y(*)    - input array to check

c      Outputs:
c         wd1(4)  - macro 1    
c         wd1(4)  - macro 2
c         v(3)    - macro parameters
c         error   - error on input
c
c----------------------------------------------------------------------
      logical error,pcomp
      double precision v(3)
      character*4 wd1,wd2
      character*75 y
      character*1 yval(45)
      call chkblk(y,15,75)

      if(pcomp(y,'proc',4)) then
        read(y,'(2(a4,11x))',err=400) wd1,wd2
      else
c       org         
cww     read(y,'(2(a4,11x),3f15.0)',err=400) wd1,wd2,v
c       new version getting value from defined cons/param
        read(y,'(2(a4,11x),45(a1))',err=400) wd1,wd2,yval
        call setval(yval(1:15),15, v(1))
        call setval(yval(16:30),15, v(2))
        call setval(yval(31:45),15, v(3))
      end if

      return
400   error = .true.
      return
      end
c
      subroutine invec(ndm,id,x,cd,prt)
c----------------------------------------------------------------------
c
c      Purpose: generate initial values for an array e.g. velocity
c
c      Inputs:
c         id(ndm,*)   - Equation numbers for each active dof
c         ndm         - Number dof/node
c         cd(2)       - Name of array
c         prt         - print flag
c
c      Outputs:
c         x(ndm,*)    - array
c
c      Comments: 
c       new version ww 3.10.03
c       set values only if no b.c. occur!
c       set all values for specified nodes
c       values for other nodes are not changed (e.g. 0 or set before!)
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      character*4 cd(2),yyy*160
      dimension id(ndm,*),x(*),xl(6),td(8)
      n = 0
      ng = 0
100   l  = n
      lg = ng
101   if(ior.lt.0) write(*,3002)
      call dinput(td,8)
      if(errck) go to 101
      n  = td(1)
      ng = td(2)
      if(n.gt.numnp) then
cww                  write(iow,3001) n,(cd(i),i=1,2)
cww     if(ior.lt.0) write(  *,3001) n,(cd(i),i=1,2)
cww     call errclr('INVEC ')
        write(yyy,3001) n,(cd(i),i=1,2),numnp
        call drawmess(yyy,1,0)
        n = numnp+1
      end if
      if(n.le.0.or.n.gt.numnp) go to 109
      do 103 i = 1,ndm
        ii = id(i,n)
        xl(i) = td(i+2)
        if(ii.gt.0) x(ii) = xl(i)
103   continue
      if(lg) 104,100,104
104   lg  = isign(lg,n-l)
      xli =(iabs(n-l+lg)-1)/iabs(lg)
      do 105 i = 1,ndm
        ii = id(i,n)
        ll = id(i,l)
        x1 = 0.d0
        x2 = 0.d0 
        if(ii.gt.0) x1 = x(ii)
        if(ll.gt.0) x2 = x(ll)
105     xl(i) = (x1-x2)/xli
cww105     xl(i) = (x(ii)-x(ll))/xli
106   l = l + lg
      if((n-l)*lg.le.0) go to 100
      if(l.le.0.or.l.gt.numnp) go to 108
      do 107 i = 1,ndm
        ii = id(i,l)
        ll = id(i,l-lg)
        x2 = 0.d0 
        if(ll.gt.0) x2 = x(ll)
        if(ii.gt.0) x(ii) = x2 + xl(i)
cww     if(ii.gt.0) x(ii) = x(ll) + xl(i)
107   continue
      go to 106
108   continue
cww                write(iow,3000) l,(cd(i),i=1,2)
cww   if(ior.lt.0) write(  *,3000) l,(cd(i),i=1,2)
      write(yyy,3000) l,(cd(i),i=1,2)
      call drawmess(yyy,1,0)
      goto 100
109   if(.not.prt) return
                   write(iow,2000) (cd(l),l=1,2),(l,cd(1),cd(2),l=1,ndm)
      if(ior.lt.0) write(  *,2000) (cd(l),l=1,2),(l,cd(1),cd(2),l=1,ndm)
c...  print values
      do 113 j = 1,numnp
        call pzero(xl,6)
        do 110 i = 1,ndm
          ii = id(i,j)
          if(ii.gt.0) xl(i) = x(ii)
110     continue       
                   write(iow,2009) j,(xl(i),i=1,ndm)
      if(ior.lt.0) write(  *,2009) j,(xl(i),i=1,ndm)
113   continue
      return
c
cww2000  format(a1,19a4,a3//5x,'nodal',2a4//6x,'node',6(i7,a4,a2)/
2000  format(/5x,'initial nodal',2a4//7x,'node',6(i7,a4,a2)/
     1       (10x,6(i7,a4,a2)))
2009  format(i10,6f13.4:/(10x,6f13.4))
3000  format('  **ERROR** attempt to generate node',i5,' in',2a4)
cww3001  format(5x,'**ERROR** attempt to input node',i5,', terminate',
cww     1' input of nodes in',2a4)
3001  format(' Try to make input values for node',i5,' for macro init',
     1        2a4,' !   Input continued with final node',i5)
3002  format(' Input: node, inc, v(i),i=1,ndf'/'   >',$)
      end
c
      subroutine invert(a,nmax,ndm)
c----------------------------------------------------------------------
c
c      Purpose: Invert small square matrix
c
c      Inputs:
c         a(ndm,*) - Matrix to be inverted
c         nmax     - Size of upper submatrix to invert
c         ndm      - Dimension of array
c
c      Outputs:
c         a(ndm,*) - Submatrix replaces original terms, others not
c                    changed
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension a(ndm,ndm)
      do 200 n = 1,nmax
      d = a(n,n)
      do 100 j = 1,nmax
100   a(n,j) = -a(n,j)/d
      do 150 i = 1,nmax
      if(n.eq.i) go to 150
      do 140 j = 1,nmax
      if(n.ne.j) a(i,j) = a(i,j) +a(i,n)*a(n,j)
140   continue
150   a(i,n) = a(i,n)/d
      a(n,n) = 1.0/d
200   continue
      return
      end
c
      subroutine iprint(ia,ii,jj,mm,aname)
c----------------------------------------------------------------------
c
c      Purpose: Output array of integer values
c
c      Inputs:
c         ia(mm,*) - Array to output
c         ii       - Number of rows to output
c         jj       - Number of columns to output
c         mm       - Dimension of array
c         name     - Name of array to appear with outputs
c
c      Outputs:
c         none
c
c----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      character*10 aname
      dimension ia(mm,*)
      nn=(jj+5)/6
      jb=0
      do 110 n=1,nn
      ja=jb + 1
      jb=ja+5
      if(jj.lt.(ja+5)) jb=jj
      if(ior.gt.0) then
        write(iow,2000)aname,(j,j=ja,jb)
        do 100 i=1,ii
        write(iow,2001)i,(ia(i,j),j=ja,jb)
  100   continue
      else
        write(*,2000)aname,(j,j=ja,jb)
        do 105 i=1,ii
        write(*,2001)i,(ia(i,j),j=ja,jb)
  105   continue
      end if
  110 continue
      return
 2000 format(5x,'matrix ',a10/1x,'row/col',i4,8i8)
 2001 format(i4,9i8)
      end
c
      subroutine isect(x,id,ndm,ndf,numnp,prt)
c----------------------------------------------------------------------
c
c      Purpose: set global dof 6=1 and release at straight intersections
c               special macro for 5 dof shells                                   |
c               numnp must be known!
c
c      Inputs:
c         x(ndm,*)    - Nodal coordinates of mesh
c         id(ndf,*)   - Equation numbers for each active dof
c         ndm         - Spatial dimension of mesh
c         ndf         - Number dof/node
c         numnp       - Number of nodes in mesh
c         prt         - print flag
c
c      Outputs:
c         id(ndf,*)   - Equation numbers for each active dof
c
c----------------------------------------------------------------------
      USE bdata
      USE errchk
      USE iofile
      USE isecn
      implicit double precision (a-h,o-z)
      logical prt   
      dimension x(ndm,numnp),id(ndf,numnp),td(7)
      dimension dp(3),p(3),p1(3),p2(3),a(3),dx(3)
      data blank/-999.d0/,tol/1000.d0/
      data eps /0.001d0/
      call pzero(dp,3)
      call pzero(p ,3)
      call pzero(p1,3)
      call pzero(p2,3)
      call pzero(a ,3)
      call pzero(dx,3)
      call pzero(td,7)
      if(ndm.ne.3) stop 'error: isect only for ndm=3'
      is=0 
  
      if(prt)              write(iow,2000)
      if(prt.and.ior.lt.0) write(*  ,2000)
c
c.... set dof6=1 for all nodes       
      do n = 1,numnp
        id(6,n)=1 
      end do
c      
c.... release dof 6 at straight intersections      
100   if(ior.lt.0) write(*,3001)
      call dinput(td,7)
      if(errck) go to 100
      dii = dsqrt(dot(td,td,6))
c
      if(dii.le.0) go to 4
c      
      toluser = td(7)
      if(dabs(toluser).gt.1.e-10) tol = toluser
c
      do i = 1,ndm
        p1(i) = td(i)
        p2(i) = td(i+3)
      end do
c
c.... calculate distance tol, vector a = p2-p1, aa = a*a
      do i = 1,ndm
        dx(i) = pdiff(x(i,1),ndm,numnp)
        a(i)  = p2(i) - p1(i)
      end do
c.... circle tolerance
      tolc = dsqrt(dot(dx,dx,3))/tol
      aa   = dot(a,a,3)
c
c.... print base data
      if(prt)              write(iow,2001) p1,p2,toluser
      if(prt.and.ior.lt.0) write(*  ,2001) p1,p2,toluser
c
c.... find nodes 
      do 200 n = 1,numnp
        if(x(1,n).eq.blank) goto 200
c....   test if node is near line p=x  dp=p-p1
        do i = 1,ndm
          p (i) = x(i,n)
          dp(i) = p(i) - p1(i)
        end do
c....   point on a = p2-p1:  x = p1 + dlamb*a
        dlamb = dot(dp,a,3)/aa
        if(dlamb.lt.-1.d0/tol.or.dlamb.gt.(1.d0+1.d0/tol)) goto 200
c.....  vector d = dp = p - x
        do i=1,3  
          dp(i) = dp(i) - dlamb*a(i)
        end do
c.....  length of d
        dd = dsqrt(dot(dp,dp,3))
c...... set dof6=0 for point
        if(dd.lt.tolc) then
          id(6,n) = 0
          is = is+1
          isecno(is) = n ! for plot
          if(is.gt.1000) stop 'field isecno in routine isect too small'
        end if
200   continue
      go to 100

4         if(prt) then
            write(iow,2002) 
            write(iow,2003) (isecno(n),n=1,is)
          end if
          if(prt.and.ior.lt.0) then
            write(iow,2002) 
            write(*  ,2003) (isecno(n),n=1,is)
          end if
2000  format(/,' Intersections generated with ISEC on ')
2001  format(' From         Point ',3(1x,e12.5),/,
     +       ' To           Point ',3(1x,e12.5),/,
     +       ' User defined Tol   ',e12.5)
2002  format(/,' The following nodes lie on intersections:')
2003  format(2x,14i5,/,2x,14i5)
3001  format(' Input: P_1(x,y,z) P_2(x,y,z) tol >',$)
      return
      end
c
      subroutine just(y,k,n0)
c----------------------------------------------------------------------
c
c      Purpose: Justify alphanumeric data in a string:
c               - Numbers are right justified
c               - Alphanumerics remain left justified
c
c      Inputs:
c         y*(*) - Unjustified string of data
c         k     - Length of string
c         n0    - Field width for justification
c
c      Outputs:
c         y*(*) - Justified string of data
c----------------------------------------------------------------------
      character*1 y(*)
      n1 = n0 - 1
      do i = 1,k,n0
c.... find last character in string
        do j = i,i+n1
          if(y(j).ne.' ') go to 100
        end do   
        y(i+n1) = '0'
100     if(y(i+n1).eq.' ') then
c....     identify a number in the field and right justify
          if((y(i).ge.'0'.and.y(i).le.'9') 
     +                   .or. (y(i).eq.'-')
     +                   .or.(y(i).eq.'+') 
     +                   .or. (y(i).eq.'.')) then
            do j = i+n1-1,i,-1
              if(y(j).ne.' ') go to 110
            end do  
110         nl = n1 + i - j
            do l = j,i,-1
              y(l+nl) = y(l)
              y(l) = ' '
            end do
          end if
        end if
      end do
      return

      end
c
      subroutine length
c----------------------------------------------------------------------
c
c      Purpose: calculate length of array m in use and maximal used time
c
c      Inputs:
c
c      Outputs:
c
c----------------------------------------------------------------------
      USE iofile
      USE plong
      USE psize
      implicit double precision(a-h,o-z)
      real*4 tary
      integer nhunds,nsecs,nmins,nhrs
c.... used memory
      xmb   = 4.0/1024
      used  = kmax*xmb
      avail = maxm*xmb
      perc  = 100.0*used/avail
      if(ior.lt.0) write(*  ,1000) maxm,avail,kmax,used,perc
                   write(iow,1000) maxm,avail,kmax,used,perc

c.... detailed usage of m
      call pseta_out(2) 
c.... allocated arraxy
      call pmemo_out(2) 

c.... used time
      call etimef(tary)
      nhunds= mod(tary*100  ,100.)
      nsecs = mod(tary      ,60.)
      nmins = mod(tary/  60.,60.)
      nhrs  =     tary/3600.
                   write(iow,1001) nhrs,nmins,nsecs,nhunds
      if(ior.lt.0) write(*  ,1001) nhrs,nmins,nsecs,nhunds

      
      return
1000  format(1x,'*FEAP: used Memory in blank common in Elements/KB:',/
     1       1x,' avail....................',i15,1x,f15.1,/,
     1       1x,' used.....................',i15,1x,f15.1,/,
     1       1x,' percent used..................',f5.2)
1001  format(/,1x,'*FEAP: used time for calculation  ',
     1       i5,' h ',i2,' m ',i2,' s ',i2,' 1/00 s ')

      end
c
cww      subroutine   memlim (machne,ipr)
c----------------------------------------------------------------------
c
c      Purpose: core-memory limit setting routine
c
c      Inputs:
c
c      Outputs:
c
c----------------------------------------------------------------------
cww      USE psize
c..                maxm = last  useful blank common word
c..                noff = first useful blank common word
cww      integer*4 length
c
cww      goto (100,200), machne
c.... set virtual minimum/maximum for unix or VMS based computers
cww  100 noff   = 1
cww      length = 2*ipr*(maxm-noff)
c                     [bytes]
cww      goto  900
c.... set absolute minimum/maximum for CMS based computers
c         lower limit = offset from first address of common array
c         upper limit = maximum useful core-memory
cww  200 length = 75000
cww      call memory (m,ibias,length,0)
cww      ibias=1 !ww alles dummy
cww      noff   =    1 +  ibias/2*ipr
cww      goto  900
c
cww  900 maxm = noff + length/4
c
cww      return
cww      end
c
c
cww      subroutine   memory (m,ibias,length,n)
cww             ...dummy...
cww      return
cww      end
c
c
      subroutine mateerr(matenew,ix,e_ome,numel,nen1,v1,n2,n3)
c---------------------------------------------------------------------+
c
c      Purpose: calculate new material array with respect to error
c
c      Inputs:
c        ix(nen1,numel)      - element conectivity array
c        e_ome(numel,numerr) - element errors 
c        v1                  - value for modification
c        n2                  - no. of errornorm (1<n1<numerr)
c        n3                  - no. of material (def.=2)
c
c      Outputs:
c        matenew(numel)      - element material number ma
c
c      Comments:
c        change actual material to material n3(def=2) 
c        v1>0: change for value larger  than v1 
c        v1<0: change for value smaller than v1 
c
c      WW IBS KIT 04/15
c---------------------------------------------------------------------+
      USE fdata
      USE iofile
      implicit double precision (a-h,o-z)
      dimension matenew(numel),e_ome(numel,*),ix(nen1,*)
      mct=0
      
      do 10 n = 1,numel 
        mat1 = ix(nen1,n)   
        matenew(n) = mat1 ! no change 
         
        if(v1.gt.0.and.e_ome(n,n2).gt.    v1)  goto 20
        if(v1.lt.0.and.e_ome(n,n2).lt.abs(v1)) goto 20
        goto 10
20      mat2       = n3 ! change
        matenew(n) = mat2 
        if(pfr) then ! print changes
          if (mct .eq. 0) then
                             write (iow,2000)
            if (ior .lt. 0) write (  *,2000)
            mct = 50
          end if
          mct = mct - 1
          if(mat2.ne.mat1) then
                            write(iow,2001) n,e_ome(n,n2),mat1,mat2 
            if (ior .lt. 0) write(*  ,2001) n,e_ome(n,n2),mat1,mat2 
          end if
        end if
10    continue

2000  format(/'    M A T E R I A L  M o d i f i c a t i o n s'//
     +        '    Elmt',' Error       ',' Old MA',' NEW MA'/)
2001  format(i8,1pe12.4,2i7)

      end
c
      subroutine meshck(x,ix,nen,nen1,ndm,numnp,numel,nummat,errs)
c----------------------------------------------------------------------
c
c      Purpose: Perform check on mesh data to ensure nodes/elements input
c
c      Inputs:
c         x(ndm,*)    - Nodal coordinates of mesh
c         ix(nen1,*) - Element nodal connection lists
c         nen        - Maximum number of nodes/element
c         nen1       - Dimension for ix array
c         ndm         - Spatial dimension of mesh
c         numnp      - Number of nodes in mesh
c         numel      - Number of elemenst in mesh
c         nummat     - Number of material sets
c
c      Outputs:
c         errs       - Flag, true if errors detected
c
c----------------------------------------------------------------------
      USE errblk
      USE iofile
      implicit double precision (a-h,o-z)
      logical errs
      dimension ix(nen1,numel),x(ndm,numnp)
      data blk/-999.0d0/
      errs = .false.
      if(iblk.eq.1) return  ! no check if error in blkgen
c.... test on elements
      do 120 n = 1,numel
c....   test if element has nodes and material
        nensum = 0
          do 105 i = 1,nen
           ii = abs(ix(i,n))
           nensum = nensum + ii
105     continue
        imat = ix(nen1,n)
        if((imat.le.0 .or. imat.gt.nummat).and.(nensum.eq.0)) then
c....     element and material not defined 
                       write(iow,2003) n
                       write(*  ,2003) n
cww       if(ior.lt.0) write(*  ,2003) n
          errs = .true.
        else if((imat.gt.0 .and. imat.le.nummat).and.(nensum.eq.0)) then
c....     nodes not defined for element, material defined 
                       write(iow,2005) n
                       write(*  ,2005) n
cww       if(ior.lt.0) write(*  ,2005) n
          errs = .true.
        else if((imat.le.0 .or. imat.gt.nummat).and.(nensum.ne.0)) then
c....     material not defined,  nodes defined for element 
                       write(iow,2000) n
                       write(*  ,2000) n
cww       if(ior.lt.0) write(*  ,2000) n
          errs = .true.
        end if  
120   continue
cww   if(errs) return
c.... test on possible element nodes
      do 121 n = 1,numel
        do 110 i = 1,nen
          ii = abs(ix(i,n))
          if(ii.gt.numnp) then
c....       test node numbers at elements (ii>numnp)
                         write(iow,2001) ii,n
                         write(*  ,2001) ii,n
cww         if(ior.lt.0) write(*  ,2001) ii,n
            errs = .true.
          else if(ii.ne.0) then
            if(x(1,ii).eq.blk) then
c....         test coordinate input of node at elements
                           write(iow,2002) ii,n
                           write(*  ,2002) ii,n
cww           if(ior.lt.0) write(*  ,2002) ii,n
              errs = .true.
            end if
          end if
110     continue
121   continue
cww   if(errs) return
c
c.... test on nodal coordinates
      do 123 ii = 1,numnp
        if(x(1,ii).eq.blk) then
                     write(iow,2004) ii
                     write(*  ,2004) ii
cww       (ior.lt.0) write(*  ,2004) ii
          errs = .true.
        end if
123   continue
      return
c
cww2000  format(5x,' **ERROR** Material type',i6,' for element ',i6)
cww2001  format(5x,' **ERROR** Node number  ',i6,' for element  ',i6, 
cww     1                ' is greater than maximum')
cww2002  format(5x,' **ERROR** Data for node ',i6,'( at element',i6,
cww     1                ' ) not input')
2001  format('  Node number ',i4,' for element ',i4,' > max. node')
2002  format('  Data for node ',i4,'( at element',i4,' ) not input')
2000  format('  Mat./El.  ',i6,' not defined')
2003  format('  Element   ',i6,' not defined')
2004  format('  Node      ',i6,' not defined')
2005  format('  Nodes/El. ',i6,' not defined')
      end
c
      subroutine modify(b,ld,s,dul,nst)
c----------------------------------------------------------------------
c
c      Purpose: Modify element residual for effects of specified
c               boundary values.
c
c               b(i) = b(i) - s(i,j)*dul(j)
c
c      Inputs:
c         b(*)   - Residual 
c         ld(*)  - Array with negative entries where boundary
c                  solution to be imposed
c         s(*,*) - Element tangent array
c         dul(*) - Value of specified solution increments
c         nst    - Dimension of element arrays

c      Outputs:
c         b(*)   - Residual modified for effect of increments
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension b(*),s(nst,nst),dul(nst),ld(nst)
c.... loop over the columns and search for boundary terms
      do 110 j = 1,nst
        if(ld(j).lt.0) then
c.... loop over the rows to modify active equations
          do 100 i = 1,nst
            ii = ld(i)
            if(ii.gt.0) then
              !$OMP ATOMIC
              b(ii) = b(ii) - s(i,j)*dul(j) ! critical OMP
            end if
100       continue
        end if
110   continue
      return
      end
c
      subroutine mprint(a,ii,jj,mm,aname)
c----------------------------------------------------------------------
c
c      Purpose: Output array of integer values
c
c      Inputs:
c         a(mm,*)  - Array to output
c         ii       - Number of rows to output
c         jj       - Number of columns to output
c         mm       - Dimension of array for first subscript
c         aname    - Name of array to appear with outputs
c
c      Outputs:
c         none
c
c----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      character*(*) aname
      dimension a(mm,*)
      nn = (jj+5)/6
      jb = 0
 
      do 110 n=1,nn
      ja = jb + 1
      jb = min(jj,ja+5)
      if(ior.gt.0) then
        write(iow,2000) aname,(j,j=ja,jb)
        do 100 i=1,ii
        write(iow,2001) i,(a(i,j),j=ja,jb)
  100   continue
      else
        write(*,2000) aname,(j,j=ja,jb)
        do 105 i=1,ii
        write(*,2001) i,(a(i,j),j=ja,jb)
  105   continue
        write(iow,2000) aname,(j,j=ja,jb)
        do 106 i=1,ii
        write(iow,2001) i,(a(i,j),j=ja,jb)
  106   continue
      end if
  110 continue

      return
 2000 format(5x,'matrix ',a/1x,'row/col',i6,5i12)
 2001 format(i4,3x,1p6e12.4)
      end
c
      subroutine mprinte(a,ii,jj,mm,eps,aname)
c----------------------------------------------------------------------
c
c      Purpose: Output array of integer values like MPRIN but with EPS 
c
c      Inputs:
c         a(mm,*)  - Array to output
c         ii       - Number of rows to output
c         jj       - Number of columns to output
c         mm       - Dimension of array for first subscript
c         aname    - Name of array to appear with outputs
c         eps      - values less eps are set to zero
c
c      Outputs:
c         none
c
c----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      character*(*) aname
      dimension a(mm,*)
      nn = (jj+5)/6
      jb = 0
 

      do i = 1,ii
        do j = 1,jj
          if(abs(a(i,j)).lt.eps) a(i,j)=0.d0  
        end do
      end do
 
      do 110 n=1,nn
      ja = jb + 1
      jb = min(jj,ja+5)
      if(ior.gt.0) then
        write(iow,2000) aname,(j,j=ja,jb)
        do 100 i=1,ii
        write(iow,2001) i,(a(i,j),j=ja,jb)
  100   continue
      else
        write(*,2000) aname,(j,j=ja,jb)
        do 105 i=1,ii
        write(*,2001) i,(a(i,j),j=ja,jb)
  105   continue
        write(iow,2000) aname,(j,j=ja,jb)
        do 106 i=1,ii
        write(iow,2001) i,(a(i,j),j=ja,jb)
  106   continue
      end if
  110 continue

      return
cww 2000 format(5x,'matrix ',a/1x,'row/col',i6,5i12)
cww 2001 format(i4,3x,1p6e12.4)
 2000 format(5x,'matrix ',a/1x,'row/col',i6,5i14)
 2001 format(i4,3x,1p6e14.6)
      end
c
      subroutine mprintp(a,ii,jj,mm,aname)
c----------------------------------------------------------------------
c
c      Purpose: Output array of integer values 
c               in file with respect to processor/thread number
c
c      Inputs:
c         a(mm,*)  - Array to output
c         ii       - Number of rows to output
c         jj       - Number of columns to output
c         mm       - Dimension of array for first subscript
c         aname    - Name of array to appear with outputs
c
c      Outputs:
c         none
c
c      Comments:  
c           1) Files foutpar_0i,i=1,8 have to be opened with DEBU,par
c           2) direct print write(51+ip,*)  xxxx 
c              with  ip = OMP_GET_THREAD_NUM()
c           3) print matrices with MPRINTP   
c
c           4) THREAD numbers are 0,1,2,3,....
c
c----------------------------------------------------------------------
      USE foutp
      implicit double precision (a-h,o-z)
      integer OMP_GET_THREAD_NUM
      character*(*) aname
      dimension a(mm,*)
  
      ip = OMP_GET_THREAD_NUM()

      nn = (jj+5)/6
      jb = 0
  
      do 110 n=1,nn
      ja = jb + 1
      jb = min(jj,ja+5)
      if(ior.gt.0) then
        write(51+ip,2000) aname,ip,(j,j=ja,jb)
        do 100 i=1,ii
        write(51+ip,2001) i,(a(i,j),j=ja,jb)
  100   continue
      else
        write(*,2000) aname,ip,(j,j=ja,jb)
        do 105 i=1,ii
        write(*,2001) i,(a(i,j),j=ja,jb)
  105   continue
        write(51+ip,2000) aname,ip,(j,j=ja,jb)
        do 106 i=1,ii
        write(51+ip,2001) i,(a(i,j),j=ja,jb)
  106   continue
      end if
  110 continue

      return
 2000 format(5x,'matrix ',a,5x,'THREAD',i4,/1x,'row/col',i6,5i12)
 2001 format(i4,3x,1p6e12.4)
      end
c
      subroutine norm(x,y,n)
c----------------------------------------------------------------------
c
c      Purpose: normalize vector y to unit vector x
c
c      Inputs:
c         y(n)  - vector general length
c         n     - length of vector
c
c      Outputs:
c         x(n)  - vector unit length
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension x(*),y(*)
      scale = dsqrt(dot(y,y,n))
      do 100 i = 1,n
100   x(i) = y(i)/scale
      return
      end
c
      subroutine numass(b,neq,mq,jd,imas)
c----------------------------------------------------------------------
c
c      Purpose: Determine number of diagonal entries in B array of
c               generalized eigen problem:  A*x = B*x*lambda
c
c      Inputs:
c         b(*)  - Diagonals of B array
c         neq   - Number of equations
c         jd(*) - Pointer array for row/columns of tangent
c         imas  - 1=CMAS  2=LMAS
c
c      Outputs:
c         mq    - Number of non-zero entries in B-diagonal
c
c      Comments:
c       subspace size reduction by number of nonzero lumped mass terms
c 
c----------------------------------------------------------------------
      USE iofile
      USE iscsr
      USE soltyp
      implicit double precision (a-h,o-z)
      dimension b(*),jd(*)
      nn = 0
      
      if ((istyp.eq.1 .or. istyp.eq.2).and.imas.eq.1) then ! SM
        do n = 1,neq
          if( b(jd(n)).ne.0.0d0) nn = nn + 1
        end do   

      else if((istyp.ge.3.and.istyp.le.8).and.imas.eq.1) then !  all other CSR
        call dnzero_csr(neq,b,csrka,nn,1,1)

      
      else  ! standard-solver for all and others for diagonals
        do n = 1,neq
          if(b(n).ne.0.0d0) nn = nn + 1
        end do  
      end if
            
                   if(nn.lt.mq) write(iow,2000) nn
      if(ior.lt.0.and.nn.lt.mq) write(*  ,2000) nn
      mq = min0(mq,nn)
      return
2000  format(' subspace size reduced to',i4,' by number of nonzero lumpe
     1d mass terms')
      end
c
      subroutine opnum(ix,nd,ln,ne,ndw,msum,nfrnt,
     .                 numnp,numel,nen,nen1,nbn,prt)
c-----------------------------------------------------------------------
c      Purpose: Calculation of best order to number equations.
c               Numbers equations for minimum front width/profile.

c               Ref: M. Hoit and E.L. Wilson, 'An Equation
c                    Numbering Algorithm Based on a Minimum
c                    Front Criteria,' Computers & Structures,
c                    v 16, No. 1-4, pp225-239, 1983.
c
c      Inputs:
c         ix(nen1,*)     - Element nodal connection list
c         numnp          - Number of nodes in mesh
c         numel          - Number of elements in mesh
c         nen            - Maximum number of nodes/element
c         nen1           - Dimension of ix  array
c         nbn            - max. number of neigbouring nodes at a node
c         prt            - Print results if true
c
c      Scratch:
c         nd         - Nodes currently in front
c         ln         - Pointer array for element number
c         ne         - Element numbers connected to nodes
c         ndw        - Node weights
c         msum       - Element weights
c
c      Outputs:
c         nfrnt(i)   - Original node for new number i
c         msum(i)    - New element number for element i
c
c----------------------------------------------------
c
c..... type declaration for variables
      USE iofile
      logical prt
      integer ie, l, m,minw,ml,mm,mh, n,nf,nn,nsum,
     1        numnp,numel,nen,nen1,node,nume,nstart,numb
c..... type declaration for arrays
      integer ix(nen1,numel),nd(numnp),ln(numnp),ne(*),ndw(numnp),
     1        msum(numel),nfrnt(numnp)
c..... intrinsic declaration
      intrinsic abs
c
c---- initialization --------------------------------
      write(iow,2010)
      write(*  ,2010)
2010  format(/,9x,'Hoit/Wilson Front Optimization (OPTI)')
      if(prt) write (iow,2000)
      node = 0
      nume = 0
c
c---- set the start lists
      call nodel(ix,nd,ln,ne,numnp,numel,nsum,nen,nen1,nbn)
      do 100 m=1,numel
          msum(m) = 1
  100 continue
      call perform(1,3,2,29)
c---- locate starting node --------------------------
  200 nstart  = 0
      call nodew(ndw,msum,ix,nen,nen1,numnp,numel,nstart)
      if(nstart.ne.0) then
          if(prt) write (iow,2004) nstart
          nfrnt(1) = nstart
          numb = 1
c----   find next element to be added to front -------
  350     ie = 0
          minw = 32000
c----   (a) loop over existing nodes on front --------
          do 400 nn=1,numb
            nf = nfrnt(nn)
            mh = ln(nf)
            ml = 1
            if(nf.ne.1) ml = ln(nf-1) + 1
c----     (b) for each node on front check elements ---
            do 380 mm=ml,mh
c
              m = ne(mm)
              if(msum(m).gt.0) then
c----         (c) evaluate increase or decrease in front --
                nsum = 0
                do 370 l= 1,nen
                      n = abs(ix(l,m))
                      if(n.gt.0) then
                        if(ndw(n).eq.msum(m)) nsum = nsum - 1
                        if(nd(n).ge.0)        nsum = nsum + 1
                      end if
  370           continue
c----         (d) compare with previous minimum -----------
                if(nsum.lt.minw) then
                      minw = nsum
                      ie   = m
                end if
              end if
  380       continue
c
  400     continue
c----   substract element sums from nodal values -----
          m = ie
          if(prt) write (iow,2001) m
c
          do 450 l=1,nen
c----     (a) reduce node sums by element weights ------
            n = abs(ix(l,m))
            if(n.gt.0) then
              ndw(n) = ndw(n) - msum(m)
              call front(nfrnt,nd,n,numb,1)
c----       (b) check if equation is to be numbered ------
              if(ndw(n).eq.0) then
                node = node + 1
                nd(n) = node
                call front(nfrnt,nd,n,numb,2)
                if(prt) write (iow,2002) n
              end if
c
            end if
  450     continue
c----   (b) remove element from system --------------
          nume = nume + 1
          msum(m) = - nume
c----   check if front has been reduced to zero -----
          if(numb.eq.0) go to 200
          go to 350
      end if
      call perform(2,3,2,29)
c---- put in final order - add inactive nodes -----
      do 500 n=1,numnp
c
        if(nd(n).eq.0) then
            node = node + 1
            nfrnt(node) = n
            if(prt) write (iow,2003) n
        else
            nn = nd(n)
            nfrnt(nn) = n
        end if
c
  500 continue
c
      do 600 n = 1,numel
          ne(n) = abs(msum(n))
  600 continue
      do 700 n = 1,numel
          m     = ne(n)
          msum(m) = n
  700 continue
      call perform(3,3,2,29)
c.... formats
 2000 format(/1x,'Calculation of Minimum Front'/)
 2001 format(5x,'Next element on front=',i6)
 2002 format(5x,'Next node numbered =',i6)
 2003 format(5x,'Node with no element attached =',i7)
 2004 format(5x,'Starting node number =',i7/)
      end
c
      subroutine front(nfrnt,nd,n,numb,ntag)
c----------------------------------------------------------------------
c
c      Purpose: Compute front width for profile optimization.
c
c      Inputs:
c         n         - Number of node to add or remore
c         numb      - Number of nodes on front
c         ntag      - Indicator to add or remove node from front
c
c      Outputs:
c         nfrnt(*)  - Current nodes on front
c         nd(*)     - Indicator on active nodes
c
c----------------------------------------------------------------------
c..... type declaration for variables
      integer l, m, n,numb,ntag
c..... type declaration for arrays
      integer nfrnt(*),nd(*)
      if(ntag.eq.1) then
c---- add node to front if new node ------------------
        do 100 m=1,numb
          if(n.eq.nfrnt(m)) return
  100   continue
        numb = numb + 1
        nfrnt(numb) = n
        nd(n) = -1
      else
c---- remove node from front --------------------------
        do 200 m=1,numb
          if(n.eq.nfrnt(m)) go to 220
  200   continue
c
  220   numb = numb - 1
        if(m.le.numb) then
c
          do 600 l=m,numb
            nfrnt(l) = nfrnt(l+1)
  600     continue
        end if
      end if
c
      end
c
      subroutine optn(ix,newno,jt,newjt,jmem,memjt,numnp,numel,
     1                nen,nen1,nren,nbn,prt)
c---------------------------------------------------------------------
c
c      Purpose: Calculation of best order to number equations.
c               Bandwidth-optimization cuthill-mckee algorithm
c               (c) w.wagner,hannover feb. 1994
c
c               nodes without elements possible?
c
c      Input:
c         ix(nen1,*)     - Element nodal connection list
c         newno(numnp)   - new node numbers 
c         numnp          - Number of nodes in mesh
c         numel          - Number of elements in mesh
c         nen            - Maximum number of nodes/element
c         nen1           - Dimension of ix  array
c         nren           - ?? not used??
c         nbn            - max. number of neigbouring nodes at a node
c         prt            - Print results if true
c
c      Scratch:
c         jt(numnp)        -  numbering array
c         newjt(numnp)     -  numbering array
c         jmem(numnp)      -  total no of possible neigbour nodes
c         memjt(numnp,nbn) -  node no. of all neigbouring nodes of a node
c
c      Output:
c         newno(numnp)     -  first:  new node number for node i
c                          -  second: original node for new number i
c
c      Comments:
c      nbn = max. number of neigbouring nodes at a node
c      possible no. of nodes
c      nbn =  8      3+4 node-element plate
c      nbn = 11      4   node-element shell with T
c      nbn = 14      4   node-element shell with X
c      nbn = 17      4   node-element I-beam with stiffener
c      nbn = 20      8   node-element plate
c      nbn = 24      9   node-element plate
c      nbn = 26      8   node-volume-element 3d
c
c---------------------------------------------------------------------
      USE iofile
      logical prt
      integer jt(numnp),newjt(numnp),jmem(numnp),newno(numnp),
     1        memjt(numnp,nbn),ix(nen1,numel)
      character yyy*80
      idiff=0
      nbmax = 0
c.... set array memjt
      do 60 j=1,numel
        do 50 i=1,nen
            jnti=ix(i,j)
            do 40 ii=1,nen
              if(ii.eq.i) goto 40
              jjt=ix(ii,j)
              mem1=jmem(jnti)
              if(mem1.eq.0) goto 30
              do 20 iii=1,mem1
                index1=memjt(jnti,iii)
                if(index1.eq.jjt) goto 40
20            continue
30            jmem(jnti)=jmem(jnti)+1
            nbmax=max(nbmax,jmem(jnti))
              if(nbmax.gt.nbn) then
                write(yyy,300)jmem(jnti),nbn
              call drawmess(yyy,1,0)
              stop
              end if
              memjt(jnti,jmem(jnti))=jjt
              if(iabs(jnti-jjt).gt.idiff) idiff=iabs(jnti-jjt)
40          continue
50        continue
          call perform(j,2*numel,2,29)
60      continue
c...  print
      write(iow,100) idiff
        write(*  ,100) idiff
      if(prt) then
        write(iow,55)
        do 57 i = 1,numnp
57      write(iow,56) i,jmem(i),(memjt(i,j),j=1,nbmax)
      end if
      minmax=idiff
c
c.... calculate optimal node numbering
      do 65 i = 1,numnp
        newno(i)= i
65    continue
      do 110 ik=1,numnp
        do 120 j=1,numnp
               jt(j) =0
            newjt(j) =0
120     continue
        maxi=0
        i=1
        newjt(1)=ik
        jt(ik)=1
        k=1
130     index=newjt(i)
        k4=jmem(index) ! sck array bounds should be checked, index can be 0 while arrays start at 1
        do 140 jj=1,k4
            k5=memjt(index,jj)
            if(jt(k5).gt.0) goto 140
            k=k+1
            newjt(k)=k5
            jt(k5)=k
            ndiff=iabs(i-k)
            if(ndiff.ge.minmax) goto 111
            if(ndiff.gt.maxi) maxi=ndiff
140     continue
        if(k.eq.numnp) goto 150
        i=i+1
        goto 130
150     minmax=maxi
c...    newno = newno(old number)
        do 155 j=1,numnp
           newno(j)=jt(j)
155     continue
111     call perform(numnp+ik,2*numnp,2,29)
110   continue
c
c...  print results
      write(iow,105) minmax
      write(*  ,105)  minmax
      if(prt) then
        write(iow,230)
        do j=1,numnp
            write(iow,250) j,newno(j)
        end do
      end if
c...  reformulation for feap: oldno = newno(new number)
c...  modify node numbering only if bandwith is reduced
      if(minmax.lt.idiff) then
          do j = 1,numnp
            jt(j) = newno(j)
          end do
          do j = 1,numnp
            ii        = jt(j)
            newno(ii) = j
          end do
      else
          do j = 1,numnp
            newno(j) = j
          end do
      end if
c
c.... Formats
55    format(' Node    No.assoc. nodes    node numbers')
56    format(13(1x,i5))
100   format(/,9x,'Cuthill-McKee Optimization of nodes (OPTN)',/,
     1         9x,'Node difference before optimization =',i5)
105   format(  9x,'Node difference after  optimization =',i5)
230   format(/,' old node no. new node no.')
250   format(3x,i6,3x,i6)
300   format(
     +'Too few storage in OPTN nbn: required:',i5,' available:',i5)
      return
      end
c
      subroutine shape1D(shp,xsi,xl,ds,ndm,nel)
c-----------------------------------------------------------------------
c
c      Purpose: One dimensional shape functions and natural derivatives
c
c      Inputs:
c         xi        - Isoparametric coordinate: ( -1 < xi < 1 )
c         nel       - Number of nodes / element   : ( 2 or 3 )
c
c      Outputs:
c         shp(2,3)  - Shape functions and spatial derivatives
c                     (natural derivatives only)
c         shp(1,i)  - Shape function spatial derivative: N_i,xi
c                     (natural derivatives only)
c         shp(2,i)  - Shape function                   : N_i
c
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension shp(2,*), xl(ndm,*), xs(3)

      go to (1,2,3,4) nel-1

c.... linear element
1     shp(2,1) = (1 - xsi)*0.5d0
      shp(2,2) = (1 + xsi)*0.5d0
c.... derivatives
      shp(1,1) = -1.d0/2.d0
      shp(1,2) =  1.d0/2.d0
      goto 100

c.... quadratic element
2     one = 1.d0
      hal = 1.d0/2.d0
      shp(2,1) = -hal * xsi * (one-xsi)
      shp(2,2) =  one - xsi*xsi
      shp(2,3) =  hal * xsi * (one+xsi)
c.... derivatives
      shp(1,1) = -(0.5d0 - xsi)
      shp(1,2) = -(2.0d0 * xsi)
      shp(1,3) =  (0.5d0 + xsi)
      goto 100

c.... cubic element
3     one = 1.d0
      dri = 1.d0/3.d0
      sec = 1.d0/16.d0
      xs2 = xsi*xsi
      shp(2,1) = -9.d0*sec * (dri+xsi) * (dri-xsi) * (one-xsi)
      shp(2,2) = 27.d0*sec * (one+xsi) * (dri-xsi) * (one-xsi)
      shp(2,3) = 27.d0*sec * (one+xsi) * (dri+xsi) * (one-xsi)
      shp(2,4) = -9.d0*sec * (one+xsi) * (dri+xsi) * (dri-xsi)
c.... derivatives
      shp(1,1) = sec *(  1.d0 + 18.d0*xsi - 27.d0*xs2)
      shp(1,2) = sec *(-27.d0 - 18.d0*xsi + 81.d0*xs2)
      shp(1,3) = sec *( 27.d0 - 18.d0*xsi - 81.d0*xs2)
      shp(1,4) = sec *( -1.d0 + 18.d0*xsi + 27.d0*xs2)
      goto 100

c.... quart. element
4     one = 1.d0
      hal = 1.d0/2.d0
      dri = 1.d0/3.d0
      xs2 = xsi*xsi
      xs3 = xs2*xsi
      shp(2,1) =  2.d0*dri *(hal+xsi) *xsi *(hal-xsi) *(one-xsi)
      shp(2,2) = -8.d0*dri *(one+xsi) *xsi *(hal-xsi) *(one-xsi)
      shp(2,3) =  4.d0 *(one+xsi) *(hal+xsi) *(hal-xsi) *(one-xsi)
      shp(2,4) =  8.d0*dri *(one+xsi) *(hal+xsi) *xsi *(one-xsi)
      shp(2,5) = -2.d0*dri *(one+xsi) *(hal+xsi) *xsi *(hal-xsi)
c.... derivatives
      shp(1,1) = dri *(  hal -   one*xsi -  6.d0*xs2 +  8.d0*xs3) 
      shp(1,2) = dri *(-4.d0 + 16.d0*xsi + 12.d0*xs2 - 32.d0*xs3) 
      shp(1,3) =      (       -10.d0*xsi +             16.d0*xs3) 
      shp(1,4) = dri *( 4.d0 + 16.d0*xsi - 12.d0*xs2 - 32.d0*xs3) 
      shp(1,5) = dri *( -hal -   one*xsi +  6.d0*xs2 +  8.d0*xs3) 
100   continue

      do 200 i = 1,ndm
         xs(i) = 0.d0
        do 200 j = 1,nel
           xs(i) = xs(i) + xl(i,j)*shp(1,j)
200   continue
      if(ndm.eq.2) xs(3) = 0.d0

      call vnorm(xs,ds)

      do 300 i = 1,nel
        shp(1,i) = shp(1,i)/ds
300   continue

      return
      end
c
      subroutine gaus1D(ngaus,pg,wg)
c-----------------------------------------------------------------------
c
c      Purpose: Form Gauss points and weights for 1D 1-5 points
c
c      Inputs:
c        ngaus     - Number of points
c
c      Outputs:
c         pg(*)    - position Gauss point
c         wg(*)    - weight   Gauss point
c
c-----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension pg(*),wg(*)

      go to (1,2,3,4,5) ngaus
c.... 1 point integration
1     pg(1) = 0.0d0
      wg(1) = 2.0d0
      return

c.... 2 point integration
2     pg(1) = -dsqrt(3.0d0)/3.0d0
      pg(2) = -pg(1)
      wg(1) =  1.0d0
      wg(2) =  1.0d0
      return

c.... 3 point integration
3     pg(1) = -dsqrt(3.d0/5.d0)
      pg(2) = 0.0d0
      pg(3) = -pg(1)
      wg(1) = 5.d0/9.d0
      wg(2) = 8.d0/9.d0
      wg(3) = wg(1)
      return

c.... 4 point integration
4     b = (2.d0/7.d0)*dsqrt(6.d0/5.d0)
      a = 3.d0/7.d0
      pg(1) = -dsqrt(a + b)
      pg(2) = -dsqrt(a - b)
      pg(3) = -pg(2)
      pg(4) = -pg(1)
      wg(1) = 0.5d0 - (dsqrt(5.d0/6.d0))/6.d0
      wg(2) = 1.d0 - wg(1)
      wg(3) = wg(2)
      wg(4) = wg(1)
      return

c.... 5 point integration
5     a = 5.d0/9.d0
      b = (2.d0/9.d0)*dsqrt(10.d0/7.d0)
      pg(1) = -dsqrt(a + b)
      pg(2) = -dsqrt(a - b)
      pg(3) =  0.d0
      pg(4) = -pg(2)
      pg(5) = -pg(1)
      d = (13.d0/90.d0)*dsqrt(7.d0/10.d0)
      c = 161.d0/450.d0
      wg(1) = c - d
      wg(2) = c + d
      wg(3) = 128.d0/225.d0
      wg(4) = wg(2)
      wg(5) = wg(1)
      return
      end
c

