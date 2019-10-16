      subroutine promul(al,au,ad,b,c,jp,ne,ityp,iadd)
c----------------------------------------------------------------------
c
c      Purpose: Routine to form c = c +/- a*b
c
c      Inputs:
c         al(*)  - Lower part of matrix
c         al(*)  - Upper part of matrix
c         ad(*)  - Diagonal of matrix
c         b(*)   - Vector to multiply
c         jp(*)  - Pointer for row/column ends of profile
c         ne     - Number equations
c         ityp   - =1: a is a square matrix  stored in form with
c                        respect to solver
c                - =2: a is a diagonal matrix  stored as vector
c
c         iadd   - =1: add to c,  =2: subtract from c

c      Outputs:
c         c(*)   - Vector with added matrix vector sum
c
c      Comments:
c
c      istyp  = 0      => standard solver         a in profile
c               1 or 2 => SM sparse matrix solver a in profile of SM
c               3      => SuperLU
c               4      => Pardiso
c               5      => PBCG
c               6      => PGMRES
c               7      => PGMRES2
c               8      => Pardiso iterative
c
c----------------------------------------------------------------------
      USE cdata
      USE iscsr
      USE mdata
cww      USE psize
      USE soltyp
      implicit double precision (a-h,o-z)
      dimension al(*),au(*),ad(*),b(*),c(*),jp(*)
      logical add
      if(ityp.eq.1) then
        if(iadd.eq.1) then
          add=.true.
        else
          add=.false.
        end if
c
        if (istyp.eq.0) then
c....     standard solver
          call promul1(al,au,ad,b,c,jp,ne,add)

        else if (istyp.eq.1.or.istyp.eq.2) then
c....     SM sparse matrix solver
          call promul_csr(ad,b,c,jp,csrja,ne,add,isymcsr)

        else if (istyp.ge.3.and.istyp.le.8) then
c....     all CSR solver
          call promul_csr(ad,b,c,jp,csrja,ne,add,isymcsr)

        end if

      else if(ityp.eq.2) then

        if(iadd.eq.1) then
          do n = 1,neq
            c(n) = c(n) + ad(n) * b(n)
          end do

        else if(iadd.eq.2) then
          do n = 1,neq
            c(n) = c(n) - ad(n) * b(n)
          end do

        end if

      end if
      return
      end
c
      subroutine promul1(al,au,ad,b,c,jp,ne,add)
c----------------------------------------------------------------------
c.... routine to form c = c +/- a*b where a is a  square matrix
c.... stored in profile form, b,c are vectors, and jp locates the
c.... bottom of columns or rows in a.
c.... if add .true. : add to c, if add .false. : subtract from c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical add
      dimension al(*),au(*),ad(*),b(*),c(*),jp(*)
      neq = iabs(ne)
      jd  = 0
c.... multiply lower triangular part
      if(ne.gt.0) then
        do 100 j = 1,neq
          js = jd
          jd = jp(j)
          if(jd.gt.js) then
            jh = jd - js
            ab = ddot(jh,al(js+1),1,b(j-jh),1)
            if(add) then
              c(j) = c(j) + ab
            else
              c(j) = c(j) - ab
            end if
          end if
100     continue
c.... do diagonal part
        do 150 j = 1,neq
          if(add) then
            c(j) = c(j) + ad(j)*b(j)
          else
            c(j) = c(j) - ad(j)*b(j)
          end if
150     continue
      end if
c.... multiply the upper triangular part
      do 200 j = 1,neq
        js = jd
        jd = jp(j)
        bj = b(j)
        if(add) bj = -bj
        if(jd.gt.js) then
           jh = jd - js
           call daxpty(jh,-bj,au(js+1),c(j-jh))
        end if
200   continue
      return
      end
c
c
      double precision function propld(t,j)
c----------------------------------------------------------------------
c      Purpose: Proportional load table (j load cards,maximum 10)
c
c      Inputs:
c         t         - Time for lookup
c         j         - Number of proportional load to define.
c                     N.B. Only used for inputs
c
c
c      Input data
c                  ityp = td(1)  type of function  def = 1
c                  iexp = td(2)  exponent          def = 1
c                  a(1) = td(3)  tmin              def = 0
c                  a(2) = td(4)  tmax              def = 10^8
c                  a(3) = td(5)  A1                def = 0
c                  a(4) = td(6)  A2                def = 1
c                  a(5) = td(7)  A3                def = 0
c                  a(6) = td(8)  A4                def = 0
c
c      Load types:
c               1.) general function
c                   prop = A1 + A2*(t-tmin) + A3*sin[A4*(t-tmin)]**iex
c                         for: tmin <= t <= tmax
c
c               2.) Sawtooth - cyclic loading
c                   prop =
c                         for: tmin <= t <= tmax
c
c               3.) Polynomial
c                   prop = A1 + A2*(t-tmin) + A3*(t-tmin)^2 + A4*(t-tmin)^3
c                         for: tmin <= t <= tmax
c
c      Outputs:
c         propld    - Value of total proportional load
c
c      Comments: up to 10 cards are allowed
c
c
c----------------------------------------------------------------------
      USE errchk
      USE iofile
      USE prlod
      implicit double precision (a-h,o-z)
      dimension td(8)
      save nprop
      if(j.gt.0) then
c....   input table of proportional loads
                     write(iow,2000)
        if(ior.lt.0) write(*  ,2000)
        if(ior.lt.0) write(*  ,2003)
        do 100 i=1,j
          iexp(i)  = 1
          do jj=1,6           ! initialize parameters
            a(jj,i) = 0.d0
          end do
101       call dinput(td,8)
          if(errck) go to 101
          ik(i)   = td(1)
          if(ik(i).gt.0) then !
            iexp(i)  = td(2)
            a(1,i)  = td(3)
            a(2,i)  = td(4)
            a(3,i)  = td(5)
            a(4,i)  = td(6)
            a(5,i)  = td(7)
            a(6,i)  = td(8)
            if(ik(i).eq.2 .and. a(4,i).eq.0.0d0) then
                           write(iow,2004)
              if(ior.lt.0) write(*  ,2004)
              go to 101
            end if
          else  ! default values
            a(2,i) = 1.e+8
            a(4,i) = 1.
            ik(i)  = 1
          end if
                      write(iow,2001) i,ik(i),(a(m,i),m=1,6),iexp(i)
          if(ior.lt.0)write(*  ,2001) i,ik(i),(a(m,i),m=1,6),iexp(i)
100     continue
        nprop = j
        propld = prop
      else
c....   compute value at time t
        propl = 0.0
        do 200 i = 1,nprop

          if(t.ge.a(1,i).and.t.le.a(2,i)) then

            if(ik(i).eq.1) then
c....         functional proportional loading
              tt    = t - a(1,i)
              l     = max(iexp(i),1)
              propl = propl + a(3,i)
     1              + a(4,i)*tt + a(5,i)*((sin(a(6,i)*tt))**l)

            else if(ik(i).eq.2) then
c....         sawtooth proportional loading
              l = int((t+a(4,i))/(a(4,i)+a(4,i)))
              m = max(1,iexp(i))
              m = int((t + (m-1)*a(4,i))/(m*a(4,i))) + 1
              propl = propl +
     1               (1 - 2*mod(m,2))*a(3,i)*(l+l-t/a(4,i))

            else if(ik(i).eq.2) then
c....         polynomial proportional loading
              tt    = t - a(1,i)
              propl = propl + a(3,i)
     1              + a(4,i)*tt + a(5,i)*tt**2 + a(6,i)*tt**3

            end if

          end if
200     continue
        propld = propl
      end if
      return
c.... formats
2000  format(30x,'Proportional Load Table')
2001  format(/,' number    type      tmin',10x,'tmax',/i3,i10,7x,g10.4,
     1 4x,g10.4,/6x,'a(1)',10x,'a(2)',10x,'a(3)',10x,'a(4)',10x,
     2 'exp',/4x,4(g10.4,4x),i5/)
2003  format(' Input: type, exponent, tmin, tmax, a(i),i=1,4')
2004  format( ' Zero a(2) in sawtooth function - reinput prop')
      end
c
      subroutine prtdii1(ieq,isw)
c----------------------------------------------------------------------
c
c      Purpose: print node number/dof for nodes with negative
c               diagonal elements
c
c      Inputs:
c         ieq       -  equation number
c         isw       -  =1:  sign of diagonal changed
c                   -  =2:  reduced diagonal zero
c                   -  =3:  loss of at least 7 digits
c
c      Outputs:
c         None
c
c.... set field for plot via 'NDII'
c     w.wagner 5/94
c----------------------------------------------------------------------
      USE mdata
      call prtdii2(psid,econ,ieq,isw)
      return
      end
c
      subroutine prtdii2(id,ix,ieq,isw)
c----------------------------------------------------------------------
c
c      Purpose: print node number/dof for nodes with negative
c               diagonal elements
c
c      Inputs:
c         id(ndf,*)   - Equation numbers for each active dof
c         ix(nen1,*)  - Element nodal connections of mesh
c         ieq         - equation number
c         isw         -  =1:  sign of diagonal changed
c                     -  =2:  reduced diagonal zero
c                     -  =3:  loss of at least 7 digits
c
c      Outputs:
c         None
c
c     only 50 nodes can be plotted, see common ndii
c----------------------------------------------------------------------
      USE cdata
      USE dii
      USE fdata
      USE iofile
      USE sdata
      USE soltyp
      integer id(ndf,numnp),ix(nen1,numel)
c
c.... find node/dof
      do i = 1,numnp
        do k = 1,ndf
          j = (id(k,i))
          if(j.eq.ieq) go to 200
        end do
      end do
      return
200   continue
c
c...  print/plot node number/dof
      if(isw.eq.1) then
                             write(iow,2001) i,k,ieq
        if(pfr.and.ior.lt.0) write(  *,2001) i,k,ieq
c...    plot
        if(ii(1).lt.50) then
          ii(1) = ii(1)+1
          ndii(1,ii(1))=i
        end if
c
      else if(isw.eq.2) then
cww     iprin=1
        if (istyp.eq.1.or.istyp.eq.2) then
c         in case of smsolv, cheque if node belongs to an element
          if(inodeused(ix,i).eq.0) iprin=0
        end if
cww     if(iprin.eq.1) then
                             write(iow,2002) i,k,ieq
        if(pfr.and.ior.lt.0) write(  *,2002) i,k,ieq
c..     plot
        if(ii(2).lt.50) then
          ii(2) = ii(2)+1
          ndii(2,ii(2))=i
        end if
cww   end if
c
      else if(isw.eq.3) then
                             write(iow,2003) i,k,ieq
        if(pfr.and.ior.lt.0) write(  *,2003) i,k,ieq
c...    plot
        if(ii(3).lt.50) then
          ii(3) = ii(3)+1
          ndii(3,ii(3))=i
        end if
      end if
2001  format('**datri warning: sign of diagonal changed  ',
     1 ' node',i5,' dof ',i2,' eq.',i6)
2002  format('**datri warning: red. diagonal zero ',
     1 ' node',i5,' dof ',i2,' eq.',i6)
2003  format('**datri warning: loss of at least 7 digits ',
     1 ' node',i5,' dof ',i2,' eq.',i6)
      return
      end
c
      function inodeused(ix,n)
c-------------------------------------------------------------
c
c      Purpose:  test if node n is connected to an element
c
c      Input:
c         ix(nen1,*)  - Element nodal connections of mesh
c         n           - Node number
c
c      Output
c         inodeused   - 0 not connected
c                       1     connected
c
c-------------------------------------------------------------
      USE cdata
      USE sdata
      integer ix(nen1,numel)
c.... find element to node n
      inodeused = 0
      do j=1,numel   ! loop over elements
        do k = 1,nen ! loop over element-nodes
          if( ix(k,j) .eq. n) go to 100
        end do
      end do
      return
100   inodeused=1
      return
      end
c
      subroutine prtdis(x,b,angl,ttim,prop,numnp,ndm,ndf,n1,n2,n3,ii,
     + ipola,iline)
c----------------------------------------------------------------------
c
c      Purpose: Output nodal values for real solutions
c
c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh
c         b(*)      - Current value of solution
c         angl      - rotation angles for each node
c         ttim      - Value of solution time
c         prop      - Value of total proportional load
c         ndm       - Spatial dimension of mesh
c         ndf       - Number dof/node
c         n1        - First node to output
c         n2        - Last noed to output
c         n3        - Increment to n1
c         ii        - Type of output: 1 = displacement; 2 = velocity;
c                                     3 = acceleration
c         ipola     - base system for output
c         iline     - print only for nodes on line
c
c      Outputs:
c         None      - Outputs to file
c
c      Comments:
c      ipola = 0   x-y   displacements mixed in case of angl
c      ipola = 12  r-phi displacements  in 1-2 plane
c      ipola = 13  r-phi displacements  in 1 3 plane
c      ipola = 23  r-phi displacements  in 2-3 plane
c
c      w.wagner IBNM UH 17.5.93
c-----------------------------------------------------------------------
      USE bdata
      USE fdata
      USE iofile
      USE mdat2
      USE pnodn
      implicit double precision (a-h,o-z)
      character*6 cd,di(3)
      character*25 di1(3)
      dimension x(ndm,*),b(ndf,*),bl(6),uu(3),uur(3),angl(*)
      data cd/' coord'/,di/' displ',' veloc',' accel'/
      data di1/'d i s p l a c e m e n t s',
     +         'v e l o c i t i e s      ','a c c e l e r a t i o n s'/
c
      pi=datan(1.d0)*4.d0
c...  dofs for angl
      ij1 = ia(1)
      ij2 = ia(2)

c.... for printing nodes on line
      if(iline.eq.1) call node_line(x,numnp,ndm,n,ipos,1)
c
      kount = 0
      do 100 ni = n1,n2,n3
c....   look for tied nodes
        n = iprttie(gtie,ni)

c....   for printing nodes on line
        if(iline.eq.1) then
          call node_line(x,numnp,ndm,n,ipos,2)
          if(ipos.eq.0) go to 100
        end if

        kount = kount - 1
        if(kount.le.0) then
cww       write(iow,2000)  o,head,ttim,ipola,prop,(i,cd,i=1,ndm),
          write(iow,2000) di1(ii),ttim,ipola,prop,(i,cd,i=1,ndm),
     1                    (i,di(ii),i=1,ndf)
          if(ior.lt.0.and.pfr) then
cww         write(*  ,2000)  o,head,ttim,ipola,prop,(i,cd,i=1,ndm),
cww         write(*  ,2000) di1(ii),ttim,ipola,prop,(i,cd,i=1,ndm),
cww   1                           (i,di(ii),i=1,ndf)
cww         write(*  ,2002)  o,head,ttim,ipola,prop,(i,di(ii),i=1,ndf)
            write(*  ,2002) di1(ii),ttim,ipola,prop,(i,di(ii),i=1,ndf)
          end if
          kount = 48
        end if
        if(ipola.eq.0) then
c.....    print global(if angl mixed) displacements
          write(iow,2001) ni,(x(j,n),j=1,ndm),(b(j,n),j=1,ndf)
          if(ior.lt.0.and.pfr) then
c           write(*  ,2001) ni,(x(j,n),j=1,ndm),(b(j,n),j=1,ndf)
            write(*  ,2003) ni,                 (b(j,n),j=1,ndf)
          end if
        else if(ipola.ne.0) then
c.....    print 'polar' displacements
c.....    if(angl) transform mixed displacements to global displ.
          if(angl(n).ne.0.0d0) then
            ang = angl(n)*pi/180.d0
            cn  = cos(ang)
            sn  = sin(ang)
          else
            cn  = 1.0
            sn  = 0.0
          end if

c....     displacements 1-3
          mdf = 2
          if(ndf.gt.2) mdf=3
          call pzero(uu,3)
          do i = 1,mdf
            uu(i) = b(i,n)
          end do
          ut       = uu(ij1)*cn - uu(ij2)*sn
          uu(ij2)  = uu(ij1)*sn + uu(ij2)*cn
          uu(ij1)  = ut
          if(ndf.ge.6) then
           call pzero(uur,3)
c....       displacements 4-6
            do i = 1,3
              uur(i) = b(i+3,n)
            end do
            if(itrot.eq.0) then
              urt      = uur(ij1)*cn - uur(ij2)*sn
              uur(ij2) = uur(ij1)*sn + uur(ij2)*cn
              uur(ij1) = urt
            end if
          end if
c
          if(ipola.eq.12) then
            i=1
            k=2
            j=3
          else if(ipola.eq.13) then
            i=1
            k=3
            j=2
          else if(ipola.eq.23) then
            i=2
            k=3
            j=1
          end if
          xx = x(i,n)
          xy = x(k,n)
          phi = datan2(xy,xx)
          sn = dsin(phi)
          cs = dcos(phi)
c....     displacements 1-3
          call pzero(bl,6)
          ur    = cs*uu(i) + sn*uu(k)
          ut    =-sn*uu(i) + cs*uu(k)
          bl(i) = ur
          bl(k) = ut
          bl(j) = uu(j)
c....     displacements 4-6 (transform always to polar coordinates)
          if(ndf.ge.6) then
            urr    = cs*uur(i) + sn*uur(k)
            utr    =-sn*uur(i) + cs*uur(k)
            bl(i+3) = urr
            bl(k+3) = utr
            bl(j+3) = uur(j)
          end if
c
          if(ndf.le.6) then
            write(iow,2001) ni,(x(l,n),l=1,ndm),(bl(l),l=1,ndf)
            if(ior.lt.0.and.pfr) then
c             write(*  ,2001) ni,(x(l,n),l=1,ndm),(bl(l),l=1,ndf)
              write(*  ,2003) ni,                 (bl(l),l=1,ndf)
            end if
          else
            write(iow,2001) ni,(x(l,n),l=1,ndm),bl,(b(l,n),l=7,ndf)
            if(ior.lt.0.and.pfr) then
c             write(*  ,2001) ni,(x(l,n),l=1,ndm),bl,(b(l,n),l=7,ndf)
              write(*  ,2003) ni,                 bl,(b(l,n),l=7,ndf)
            end if
          end if
        end if
100   continue
      return
cww2000  format(a1,19a4,a3/'  n o d a l   d i s p l a c e m e n t s',5x,
2000  format(/'  n o d a l   ',a,5x,
     1  'time',e18.5/'  polar = ',i3,18x,'prop. ld. (eigenvalue)',
     2  1pe13.5/ '  node',9(i6,a6):/(6x,9(i6,a6)))
2001  format(i6,1p9e12.5:/(6x,1p9e12.5))
cww2002  format(a1,19a4,a3/'  n o d a l   d i s p l a c e m e n t s',5x,
2002  format(/'  n o d a l   ',a,5x,
     1  'time',e18.5/'  polar = ',i3,18x,'prop. ld. (eigenvalue)',
     2  1pe13.5/ '  node',9(i6,a6):/(6x,9(i6,a6)))
2003  format(i6,1p9e12.5:/(6x,1p9e12.5))
      end
c
      subroutine prterr(numel)
c----------------------------------------------------------------------
c
c      Purpose: print global error values
c
c      Input:
c
c      Output: u_om(i),e_om(i)
c
c----------------------------------------------------------------------
      USE errin1
      USE errin2
      USE errnam
      USE iofile
      implicit double precision (a-h,o-z)
cww   dimension eta_rel(2)
c.... output the error indicator values
cww   do i= 1,numerr
cww     eta_rel(i) = sqrt(e_om(i)/(u_om(i)+e_om(i)))
cww   end do

      write(*,2005)'* * E R R O R   A N A L Y S I S  GLOBAL VALUES * *'
      write(*,2005)'Norm           ','System-Energy  ','System-Error   '
     +      ,'accepted value/element'

      dnumel=numel 
      do i = 1,numerr
        v1 = sqrt(u_om(i))
        v2 = sqrt(e_om(i))
        v3 = e_bar(i)
        if(i.eq.3) v1=v1/sqrt(dnumel)  ! then Y0 
        if(i.eq.3) v2=v2/sqrt(dnumel)  ! then Sigma_v middle
        if(ior.lt.0) write(*,2010) e_name(i),v1,v2,v3  
      end do
c
cww   write(*,2005)'* relative error measures'
cww   do i = 1,numerr
cww     if(ior.lt.0) write(*,2010)  e_name(i),eta_rel(i)
cww   end do
2005  format(1x,A,A,A,A,A)
2010  format(1x,A13,3(e15.7))
      end
c
      subroutine prtom(omega,kount)
c----------------------------------------------------------------------
c
c      Purpose: print eigenvalues
c
c      Input:
c        omega   - Eigenvalue
c        kount   - No. of eigenvalue
c
c      Output:
c
c-----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      if(kount.eq.1) then
                     write(iow,2000)
        if(ior.lt.0) write(*  ,2000)
      end if
      write(*,2001) kount, omega
      return
2000  format(/'  No.   Eigenvalue')
2001  format(i4,3x,e15.8)
      end
c
      subroutine prtjint(r,ndf,n1,n2,n3,flg,irdf,ieps,reps)
c----------------------------------------------------------------------
c
c      Purpose: Output material forces for current solution
c
c      Inputs:
c         r(*)      - Current value of reactions
c         ndf       - Number dof/node
c         n1        - First node to output
c         n2        - Last  node to output
c         n3        - Increment to n1
c         flg       - true   r=-r
c         irdf      - dof to print  macro: REAC,eps,irdf,reps
c         ieps      - 0
c                     1 plot r only if  r > reps
c         reps      - tolerance for printing r
c
c      Outputs:
c         None      - Outputs to file
c
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE fdata
      USE iofile
      USE pnodn
      implicit double precision (a-h,o-z)
      logical flg
      dimension r(ndf,*)
      kount = 0
      if(ieps.eq.0) then
        do 100 ni = n1,n2,n3
            kount = kount - 1
            if(kount.le.0) then
              if(pfr)              write(iow,2000)        (k,k=1,ndf)
              if(ior.lt.0.and.pfr) write(*  ,2000)        (k,k=1,ndf)
              kount = 50
            end if
c....     look for tied nodes
          n = iprttie(gtie,ni)
            if(ior.lt.0.and.pfr) write(*  ,2001) ni,(r(k,n),k=1,ndf)
            if(pfr)              write(iow,2001) ni,(r(k,n),k=1,ndf)
100     continue
      else if(ieps.eq.1) then
        do 101 ni = n1,n2,n3
            if(dabs(r(irdf,ni)).lt.reps) go to 101
            kount = kount - 1
            if(kount.le.0) then
              if(pfr)              write(iow,2000)        (k,k=1,ndf)
              if(ior.lt.0.and.pfr) write(*  ,2000)        (k,k=1,ndf)
              kount = 50
            end if
c....     look for tied nodes
          n = iprttie(gtie,ni)
            if(ior.lt.0.and.pfr) write(*  ,2001) ni,(r(k,n),k=1,ndf)
            if(pfr)              write(iow,2001) ni,(r(k,n),k=1,ndf)
101     continue
      end if
      return
2000  format(/'  n o d a l    m a t e r i a l  f o r c e s'/5x,
     1  'node',6(i7,' dof'):/(9x,6(i7,' dof')))
2001  format(i9,1p6e11.3:/(9x,1p6e11.3))
      end
c
      subroutine prtrea(r,ndf,n1,n2,n3,flg,irdf,ieps,reps,mnip)
c----------------------------------------------------------------------
c
c      Purpose: Output nodal reactions for current solution
c
c      Inputs:
c         r(*)      - Current value of reactions
c         ndf       - Number dof/node
c         n1        - First node to output
c         n2        - Last  node to output
c         n3        - Increment to n1
c         flg       - true   r=-r
c         irdf      - dof to print  macro: REAC,eps,irdf,reps
c         ieps      - 0
c                     1 plot r only if  r > reps
c         reps      - tolerance for printing r
c         mnip      - array of node numbers with respect to tie
c
c      Outputs:
c         None      - Outputs to file
c
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE fdata
      USE iofile
      implicit double precision (a-h,o-z)
      logical flg
      dimension r(ndf,*),rsum(100),asum(100),psum(100),mnip(*)
      if(ndf.gt.100) then
        call drawmess('prtrea: rsum',1,0)
        return
      end if
      do 50 k = 1,ndf
        psum(k) = 0.
        rsum(k) = 0.
50      asum(k) = 0.
      do 75 i = 1,numnp
        do 75 k = 1,ndf
          if(flg) r(k,i) = -r(k,i)
          rsum(k) = rsum(k) + r(k,i)
75        asum(k) = asum(k) + dabs(r(k,i))
      kount = 0
      if(ieps.eq.0) then
        do 100 ni = n1,n2,n3
            kount = kount - 1
            do 90 k = 1,ndf
              psum(k) = psum(k) + r(k,ni)
90          continue
            if(kount.le.0) then
cww                                write(iow,2000) o,head,(k,k=1,ndf)
cww           if(ior.lt.0.and.pfr) write(*  ,2000) o,head,(k,k=1,ndf)
cww                                write(iow,2000)        (k,k=1,ndf)
              if(pfr)              write(iow,2000)        (k,k=1,ndf)
              if(ior.lt.0.and.pfr) write(*  ,2000)        (k,k=1,ndf)
              kount = 50
            end if
c....     look for tied nodes
          n = iprttie(mnip,ni)
            if(ior.lt.0.and.pfr) write(*  ,2001) ni,(r(k,n),k=1,ndf)
            if(pfr)              write(iow,2001) ni,(r(k,n),k=1,ndf)
cww                              write(iow,2001) ni,(r(k,n),k=1,ndf)
100     continue
      else if(ieps.eq.1) then
        do 101 ni = n1,n2,n3
            if(dabs(r(irdf,ni)).lt.reps) go to 101
            kount = kount - 1
            do 91 k = 1,ndf
              psum(k) = psum(k) + r(k,ni)
91      continue
            if(kount.le.0) then
cww                                write(iow,2000) o,head,(k,k=1,ndf)
cww           if(ior.lt.0.and.pfr) write(*  ,2000) o,head,(k,k=1,ndf)
cww                                write(iow,2000)        (k,k=1,ndf)
              if(pfr)              write(iow,2000)        (k,k=1,ndf)
              if(ior.lt.0.and.pfr) write(*  ,2000)        (k,k=1,ndf)
              kount = 50
            end if
c....     look for tied nodes
            n = iprttie(mnip,ni)
            if(ior.lt.0.and.pfr) write(*  ,2001) ni,(r(k,n),k=1,ndf)
            if(pfr)              write(iow,2001) ni,(r(k,n),k=1,ndf)
cww                              write(iow,2001) ni,(r(k,n),k=1,ndf)
101     continue
      end if
c.... print statics check
cww   write(iow,2003) (psum(k),k=1,ndf)
cww   write(iow,2004) (rsum(k),k=1,ndf)
      if(pfr) then
        write(iow,2003) (psum(k),k=1,ndf)
        write(iow,2004) (rsum(k),k=1,ndf)
      end if
cww   write(iow,2005) (asum(k),k=1,ndf)
      if(ior.lt.0.and.pfr) then
        write(*,2003) (psum(k),k=1,ndf)
        write(*,2004) (rsum(k),k=1,ndf)
cww   write(*,2005) (asum(k),k=1,ndf)
      end if
      return
cww2000  format(a1,19a4,a3/'  n o d a l    r e a c t i o n s'/5x,
2000  format(/'  n o d a l    r e a c t i o n s'/5x,
     1  'node',6(i7,' dof'):/(9x,6(i7,' dof')))
2001  format(i9,1p6e11.3:/(9x,1p6e11.3))
2003  format( 2x,'prt sum',1p6e11.3:/(9x,1p6e11.3))
2004  format( 6x,'sum',1p6e11.3:/(9x,1p6e11.3))
cww2005  format( 2x,'abs sum',1p6e11.3:/(9x,1p6e11.3))
      end
c
      subroutine prtreas(x,r,ndm,ndf,n1,n2,n3,flg)
c----------------------------------------------------------------------
c
c      Purpose: Output stresses for SIGQ=1/V*A^T R  
c               only 3D, only for equilibrium
c      Inputs:
c         x(ndm,*)  - coordinate array 
c         r(*)      - Current value of reactions
c         ndf       - Number dof/node
c         n1        - First node to output
c         n2        - Last  node to output
c         n3        - Increment to n1
c         flg       - true   r=-r
c
c      Outputs:
c         SIGQ      - Outputs to file
c
c      Comments     
c         tied nodes have r=0 -> no influence  
c         valid only 3D
c         needs to be documented in Manual 
c         only for displacement based BCs: A=A(x) and not A=A(delta x)   
c
c         only in case of convergence correct:
c          # F_a = 0 = equilibrium
c          # -> loop over all boundaries  not necessary 
c          # -> loop over all nodes sufficient
c
c     WW IBS KIT 03/2015
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE fdata
      USE iofile
      implicit double precision (a-h,o-z)
      logical flg
      dimension x(ndm,*),r(ndf,*),sigq(6),dxyz(3),xmin(3),xmax(3)
      data bl/-999.d0/

      if(ndm.eq.3.and.ndf.eq.3) then      

        d5=0.5e0
        do i = 1,numnp
          do  k = 1,ndf
            if(flg) r(k,i) = -r(k,i)
          end do
        end do  
        
c....   calculate volume of RVE =lx*ly*lz
        xmin = 0.d0
        xmax = 0.d0
        dxyz = 0.d0
      
c....   find start node without tie
        do n = 1,numnp 
          if(x(1,n).ne. bl)goto 12  
        end do
      
12      do i = 1,ndm
          xmin(i) = x(i,n)
          xmax(i) = x(i,n)
        end do
      
c....   extreme values
        do n = 1,numnp
          if(x(1,n).ne. bl) then ! tie
            do i = 1,ndm
              xmin(i) = min(xmin(i),x(i,n))
              xmax(i) = max(xmax(i),x(i,n))
            end do ! i
          end if 
        end do ! n
      
c....   length l_x,l_y,l_z,vol
        do i = 1,ndm
          dxyz(i) = xmax(i)-xmin(i) 
        end do
        vol=dxyz(1)*dxyz(2)*dxyz(3)
         
c....   SIGQ=1/V*A^T R    
        sigq = 0 
        do ni = n1,n2,n3
          sigq(1)=sigq(1) + x(1,ni)*r(1,ni) 
          sigq(2)=sigq(2) + x(2,ni)*r(2,ni) 
          sigq(3)=sigq(3) + x(3,ni)*r(3,ni) 
                      
          sigq(4)=sigq(4) + d5*x(2,ni)*r(1,ni) + d5*x(1,ni)*r(2,ni) 
          sigq(5)=sigq(5) + d5*x(3,ni)*r(1,ni) + d5*x(1,ni)*r(3,ni) 
          sigq(6)=sigq(6) + d5*x(3,ni)*r(2,ni) + d5*x(2,ni)*r(3,ni) 
        end do 
        
        sigq=sigq/vol 
        
        if(pfr)              write(iow,2000)  (sigq(k),k=1,6)
        if(ior.lt.0.and.pfr) write(*  ,2000)  (sigq(k),k=1,6)
      else
                     write(iow,'(a)')  'SIGQ only for ndm=ndf=3'
        if(ior.lt.0) write(*  ,'(a)')  'SIGQ only for ndm=ndf=3'
      end if
      return
2000  format(/'  SIGQ [S_11, S_22, S_33, T_12, T_13, T_23]',
     +        ' (only for converged solution!)',/2x,1p6e14.6)
      end
c
      subroutine prtrsum(react,nodesrf,dr,ndf,ndfrs,nfs1,pfr,iow)
c-----------------------------------------------------------------------
c      Purpose print and plot(tplo) rsum(nfs1) of reactions for dof(ndfrs)
c              if macro rsum is used in mesh
c
c      Inputs:
c         nodesrf   - array of node numbers for rsum
c         dr        - array of reaction forces
c         ndf       - Number dof/node
c         ndfrs     - dof for rsum
c         nfs1      - Number of nodes for rsum
c         pfr       - Print generated data if true
c
c      Outputs:
c         react     - Sum of reaction forces
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical pfr
      dimension nodesrf(*),dr(*)
      react = 0.d0
      do irss = 1,nfs1
        irs = nodesrf(irss)
        jreacs = ndf*(irs-1)+ndfrs
        react = react +dr(jreacs)
      end do
      if(pfr) write(*  ,*) ' reaction force from RSUM', react
              write(iow,*) ' reaction force from RSUM', react
      return
      end
c
c
      subroutine prtstr(x,dt,ds,numnp,ndm,n1,n2,n3,na,ne,mnip,iline)
c----------------------------------------------------------------------
c
c      Purpose: Output nodal projected stresses
c
c      Inputs:
c         x(ndm,*)    - nodal coordinates
c         dt(numnp,*) - Weights at nodes
c         ds(numnp,*) - Stress values at nodes
c         numnp       - Number of nodes in mesh
c         ndm         - Spatial dimension of mesh
c         n1          - Number of first node to output
c         n2          - Number of last  node to output
c         n3          - Increment to node from n1
c         na          - first stress number
c         ne          - last stress number
c         mnip        - array of node numbers with respect to tie
c         iline       - print only for nodes on line
c
c      Outputs:
c         None        - Outputs to file/screen
c
c----------------------------------------------------------------------
      USE bdata
      USE fdata
      USE iofile
      USE strnam
      implicit double precision (a-h,o-z)
      dimension dt(numnp),ds(numnp,*),x(*),mnip(*)
      kount = 0

c.... for printing nodes on line
      if(iline.eq.1) call node_line(x,numnp,ndm,n,ipos,1)

      do 200 ni = n1,n2,n3
c....   for printing nodes on line
        if(iline.eq.1) then
          call node_line(x,numnp,ndm,ni,ipos,2)
          if(ipos.eq.0) go to 200
        end if

        kount = kount - 1
        if(kount.le.0) then
          if(strsus(1).eq.'               ') then
c...        use numbers
                               write(iow,2000) (i,i=na,ne)
            if(ior.lt.0.and.pfr) write(*,2000) (i,i=na,ne)
          else
c...        use names
                               write(iow,2002)(strsus(i),i=na,ne)
            if(ior.lt.0.and.pfr) write(*,2002)(strsus(i),i=na,ne)
          end if
          kount = 50
        end if

c....     look for tied nodes
          n = iprttie(mnip,ni)

                             write(iow,2001) ni,(ds(n,i),i=na,ne)
        if(ior.lt.0.and.pfr) write(*,  2001) ni,(ds(n,i),i=na,ne)
200   continue
c

2000  format(/'   n o d a l   s t r e s s e s'/
     1 ' node ',26(3x,i7,'-value'))
2002  format(/'   n o d a l   s t r e s s e s'/
     1 ' node ',26(1x,a15))
2001  format(i5,1p26e16.5)
      end
c
      subroutine node_line(x,numnp,ndm,n,ipos,isw)
c----------------------------------------------------------------------
c
c      Purpose: check if node is on line
c
c      Inputs:
c         x(ndm,*)    - Nodal coordinates of mesh
c         numnp       - Number of nodes in mesh
c         ndm         - Spatial dimension of mesh
c
c      Outputs:
c         ipos        - 1/0  on line/ not on line
c
c      Comments:
c         line p=p1+(p2-p1)*lambda  0<lambda<1
c
c----------------------------------------------------------------------
      USE dspos
      implicit double precision (a-h,o-z)
      dimension x(ndm,numnp)
      dimension p(3),dp(3),dx(3),t(3)
      data tol/1000.d0/

      if(dsdcor2.lt.1.e-8) then ! no line set
        ipos=0
        return
      end if

      goto (1,2) isw

c...  set initial values
1     call pzero(dx,3)
      do i = 1,ndm
        dx(i) = pdiff(x(i,1),ndm,numnp)
      end do
c.... circle tolerance
      tolc = dsqrt(dot(dx,dx,3))/tol
      return

c.... check node
2     call pzero(p,3)
      call pzero(dp,3)
      call pzero(t,3)

c.... check node
      do i = 1,ndm
        p(i) = x(i,n)
        dp(i) = p(i) - dscor(i,1)
      end do

      dlamb = dot(dp,dsdcor,3)/dsdcor2
      if(dlamb.lt.-1.d0/tol.or.dlamb.gt.(1.d0+1.d0/tol)) goto 20

      do i=1,3
        dp(i) = dp(i) - dlamb*dsdcor(i)
      end do
      dd = dsqrt(dot(dp,dp,3))

      if(dd.gt.tolc) goto 20
c.... on line
      ipos=1
      return
c.... not on line
20    ipos=0
      return
      end
c
      integer function iprttie(ip,ni)
c-----------------------------------------------------------------------
c
c      Purpose: look for tied nodes
c
c      Inputs:
c       ip(*)     - Node number list tied nodes
c       ni        - node number input, eventually tied
c
c      Output
c        n        - real node number
c                                              |
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer*4 ni,n
      dimension ip(ni)
      n = ip(ni)
      iprttie = n
      return
      end
c
      subroutine prxtie(ip,ni,x,x1,ndm)
c-----------------------------------------------------------------------
c      Purpose: reset cordinate x1 for tied nodes
c
c      Inputs:
c       ip(*)        - Node number list tied nodes
c       ni           - node number input, eventually tied
c       x(ndm,numnp) - nodal coordinates
c
c      Output
c        1           - coordinate x1 of node ni
c                                              |
c-----------------------------------------------------------------------
      USE cdata
      implicit double precision(a-h,o-z)
      integer*4 ni,n
      dimension ip(numnp),x(*)
      n = ip(ni)
      x1 = x((n-1)*ndm+1)
      return
      end
c
      subroutine pseta(na,ni,ip,afl,name)
c----------------------------------------------------------------------
c
c      Purpose: set pointer for arrays stored in blank common
c
c      Input:
c        ni      - lenght of array to be added
c        ip      - iprecision
c        afl     - Flag
c
c      Output
c        na      - adress of first entry of this array in blank common
c        afl     - Flag set to .false.
c        kmax    - final entry in blank common
c        kpset   - number of used arrays
c        ipset   - list of starting adresses of used arrays
c        lpset   - list of length            of used arrays
c        npset   - list of names             of used arrays
c
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
      USE cdata
      USE debugs
      USE iofile
      USE plong
      USE psethis
      USE psize
      implicit double precision (a-h,o-z)
      logical afl
      character yyy*200
      character*(*) name
c.... set values
      na   = kmax
      kmax = na + ni*ip + mod(ni*ip,ipr)
      afl = .false.

c.... save values for history
      kpset=kpset+1
      ipset(kpset)=na
      lpset(kpset)=kmax-na
      npset(kpset)=name

      if(debug.eq.1) then
        write(iow,1002) dbgtxt
        write(*  ,1002) dbgtxt
        write(iow,1001) na,kmax-na,kmax,name
        write(*  ,1001) na,kmax-na,kmax,name
        dbgtxt = ' '
      end if
      if(kmax.le.maxm) return
      write(yyy,1000) name,kmax,maxm
      call drawmess(yyy,1,6)
      afl = .true.

      stop 'stop in SR PSETA'

1000  format('insufficient storage for array ',a,'(required:',i12,
     +       ' available:',i12,')')
      return
1001  format(' Debug Pseta: Begin: ',i10,' Length: ',i10,' Free: ',i10
     +, 3x,a)
1002  format(a80)
      end
c
      subroutine pseta_out(isw)
c----------------------------------------------------------------------
c
c      Purpose: print pointer for arrays stored in blank common
c
c      Input:
c       isw       -  1 from MACRO   2 at end
c
c      Output:
c        na       - adress of first entry of this array in blank common
c        kmax-na  - length of array
c        name     - name of array
c        kmax    - final entry in blank common
c
c      Comment:
c        300 sets possible
c
c     W. Wagner BS UKA 04/09
c----------------------------------------------------------------------
      USE iofile
      USE plong
      USE psethis
      USE psize
      implicit double precision (a-h,o-z)
      fac  = 4./1024.
      lmax = maxm-kmax
      if(isw.eq.1) then
        write(*,2000)
        do i = 1,kpset
          write(*,2001) i,npset(i),ipset(i),lpset(i),lpset(i)*fac
        end do
        write(*,2002) 'Used  space    ',0,kmax,kmax*fac
        write(*,2002) 'Free  space    ',0,lmax,lmax*fac
        write(*,2002) 'Total space    ',0,maxm,maxm*fac
      else
        write(iow,2000)
        do i = 1,kpset
          write(iow,2001) i,npset(i),ipset(i),lpset(i),lpset(i)*fac
        end do
        write(iow,2002) 'Used  space  ',0,kmax,kmax*fac
        write(iow,2002) 'Free  space  ',0,lmax,lmax*fac
        write(iow,2002) 'Total space  ',0,maxm,maxm*fac
      end if

      return
c.... formats
2000  format('  List of used arrays in Common M',/
     + 2x,'Array No',2x,'Name of Array',8x,'Beginn at',
     + 1x,' Length of M ',3x,' Length in [KB]')
2001  format(2x,i5,5x,a15,i15,i15,f15.3)
2002  format(12x,a,i15,i15,f15.3)
      end
c
      subroutine pmemo_out(isw)
c-----------------------------------------------------------------------
c
c      Purpose: print memory for allocated arrays
c
c      Input:
c       isw       -  1 from MACRO   2 at end
c
c      Output:
c        na       - adress of first entry of this array in blank common
c        kmax-na  - length of array
c        name     - name of array
c        kmax    - final entry in blank common
c
c      Comment:
c        500 sets possible
c
c     S. Klarmann TUD 11/14
c-----------------------------------------------------------------------
      USE allocdata
      USE iofile
      implicit none
      integer i,n,isw
      real*8  absize, sum,sum1
      character(len=2) nsize(5)
      character(len=30) oname
      nsize(1) = 'B '
      nsize(2) = 'KB'
      nsize(3) = 'MB'
      nsize(4) = 'GB'
      sum = 0.d0
      oname = 'Name of Array'
      oname = adjustl(oname)
      if(isw.eq.1) then  
        write(  *,2003)
        write(  *,2001) oname,'|','Size','|','No. of. Elements'
        write(  *,2002)
      else
        write(iow,2003)
        write(iow,2001) oname,'|','Size','|','No. of. Elements'
        write(iow,2002)
      end if
      do i=1,size(asizes)
        if(isalloc(i)) then
          n=1
          sum = sum + asizes(i)
          absize = asizes(i)
          do while(absize.ge.1024.d0)
            absize = absize/1024.d0
            n=n+1
          end do
          if(isw.eq.1) then  
            write(  *,2000) 
     +      allocname(i),'|',absize,nsize(n),'|',asizes(i)/8
          else
            write(iow,2000) 
     +      allocname(i),'|',absize,nsize(n),'|',asizes(i)/8
          end if
        end if
      end do

      sum1 = sum 
      n = 1
      do while(sum.ge.1024.d0)
        sum = sum/1024.d0
        n=n+1
      end do
      oname = 'Sum of all arrays:'
      oname = adjustl(oname)
      if(isw.eq.1) then  
        write(  *,2002)
        write(  *,2000) oname,'|',sum, nsize(n),'|',sum1/8
      else 
        write(iow,2002)
        write(iow,2000) oname,'|',sum, nsize(n),'|',sum1/8
      end if  
 2000 format(3x,A30,A1,F10.3,A4,A1,F16.0)
 2001 format(3x,A30,A1,A14,A1,A16)
 2002 format(3x,62('-'))
 2003 format('   List of allocated arrays ')

      return
      end
c
      subroutine pstores(ns,nv,ixl,valu,ixm,valm)
c----------------------------------------------------------------------
c
c      Purpose: moves int. + real vector to int. + real vector
c
c      Input:
c        ns    - length of integer array
c        nv    - length of    real array
c        ixl   - integer array
c        valu  -    real array
c
c      Ouput:
c        ixm   - integer array
c        valm  -    real array
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ixl(ns),ixm(ns),valu(nv),valm(nv)
      do 100 i = 1,ns
100   ixm(i) = ixl(i)
      do 110 i = 1,nv
110   valm(i) = valu(i)
      return
      end
c
      subroutine pstr3d(bb,bpr)
c----------------------------------------------------------------------
c
c      Purpose:  principal value routine for symmetric tensor
c
c      Input:
c        bb(6)   - tensor components
c
c            | bb(1) bb(4) bb(6) |
c       3d = | bb(4) bb(2) bb(5) |
c            | bb(6) bb(5) bb(3) |
c
c      Output:
c        bpr(3)  - principal values
c
c----------------------------------------------------------------------
c
c..... Declare variable types
      integer i
      real*8 pi23, tol, al, b1, b2, b3, c1, c2,c3
c     real*8 pi4
c..... Declare array types
      real*8  bd(6),bb(6),bpr(3)
c..... Intrinsics
      intrinsic atan2, cos, sqrt
c     data pi4  /0.7853981633974483d0/
      data pi23 /2.0943951023931955d0/
      data tol /1.d-12/
c....  compute mean and deviatoric (upper trianglular part) tensors
      b1  = (bb(1) + bb(2) + bb(3))/3.
      do 100 i = 1,6
        bd(i) = bb(i)
100   continue
      do 110 i = 1,3
        bd(i) = bd(i) - b1
110   continue
c....  compute 2nd and 3rd invariants of deviator
      c1 = bd(4)*bd(4)
      c2 = bd(5)*bd(5)
      c3 = bd(6)*bd(6)
      b2 = 0.5d0*(bd(1)*bd(1)+bd(2)*bd(2)+bd(3)*bd(3))
     1   + c1 + c2 + c3
      if(b2.le.tol*b1*b1) then
        bpr(1) = b1
        bpr(2) = b1
        bpr(3) = b1
      else
        b3 = bd(1)*bd(2)*bd(3)+(bd(4)+bd(4))*bd(5)*bd(6)
     1     + bd(1)*(c1-c2) + bd(2)*(c1-c3)
c....  set constants
        c1 = 2.d0*sqrt(b2/3.d0)
        c2 = 4.d0*b3
        c3 = c1*c1*c1
        al = atan2(sqrt(c3*c3-c2*c2),c2)/3.d0
c....  set principal values
        bpr(1) = b1 + c1*cos(al)
        bpr(2) = b1 + c1*cos(al-pi23)
        bpr(3) = b1 + c1*cos(al+pi23)
      end if
c
      end
c
      subroutine pstres(sig,p1,p2,p3)
c----------------------------------------------------------------------
c
c      Purpose: Compute principal stresses for 2-d problems.
c
c      Input:
c         sig(3) - Stresses in order: s_x, s_xy, s_y
c
c      Output:
c         p1     - Principal stress s_1
c         p2     - Principal stress s_2
c         p3     - angle (degrees): s_1 to s_x
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension sig(3)
      smax=sqrt(dot(sig,sig,3))
      eps=1.e-8*smax

      xi1 = (sig(1) + sig(3))/2.d0
      xi2 = (sig(1) - sig(3))/2.d0

c.... prevent numerical garbage
      if(dabs(xi1).lt.eps)    xi1   = 0.d0
      if(dabs(xi2).lt.eps)    xi2   = 0.d0
      if(dabs(sig(2)).lt.eps) sig(2)= 0.d0

c.... phi = arctan(s_xy/0.5(s_x-s_y)
c.... prevent 0/0
      if(xi2.eq.0.d0.and.sig(2).eq.0.d0) then
c        special cases Mohr-circle->point (p3 is arbitrary!!)
        if(xi1.gt.0.d0) p3 = -2.d0*datan(1.d0)  ! -45
        if(xi1.lt.0.d0) p3 =  2.d0*datan(1.d0)  !  45
        if(xi1.eq.0.d0) p3 =  0.d0              !   0
      else
        p3 = datan2(sig(2),xi2)
      end if

c.... formulation from transformation
      rho = xi2*dcos(p3) + sig(2)*dsin(p3)
      p1 = xi1 + rho
      p2 = xi1 - rho
      p3 = 22.5d0/datan(1.d0)*p3

      return
      end
c
      subroutine pstres1(sig,sig0,z)
c----------------------------------------------------------------------+
c
c      Purpose: Compute principal stresses for 3-d problems.
c
c      Input:
c         sig(6)   - Stresses Sx,Sy,Sz,Sxy,Sxz,Syz
c
c      Scratch:
c         sig3(3,3)- Stress tensor
c
c      Output:
c         sig0(3)  - Principal stress s_1,s_2,s_3
c         z(3,3)   - associated directions
c
c      Comments:
c       solv  [a -   lambda 1] z = 0    3*3
c                              |
c----------------------------------------------------------------------+
c
      implicit double precision (a-h,o-z)
      dimension sig(6),sig3(3,3),b(3,3),sig0(3),z(3,3),d(3,3),
     +          fv1(3),fv2(3),zh(3)
      data ev /0.d0/
c...  set stress tensor
      sig3(1,1) = sig(1)
      sig3(2,2) = sig(2)
      sig3(3,3) = sig(3)
      sig3(1,2) = sig(4)
      sig3(1,3) = sig(5)
      sig3(2,3) = sig(6)
      sig3(2,1) = sig3(1,2)
      sig3(3,1) = sig3(1,3)
      sig3(3,2) = sig3(2,3)
      call pzero(sig0,3)
      call pzero(b,9)
      b(1,1) = 1.d0
      b(2,2) = 1.d0
      b(3,3) = 1.d0
      call pzero(z,9)
      call pzero(d,9)
      call pzero(fv1,3)
      call pzero(fv2,3)
c.....solve EV-Problem from EISPACK
      matz = 1
      ierr = 0
      nm = 3
      call rsg(nm,3,sig3,b,sig0,matz,z,fv1,fv2,ierr)
      if(matz.ne.0) then
c....   scale to length 1
        do k = 1,3
          ev = 0.d0
          do i = 1,3
            ev = ev + z(i,k)*z(i,k)
          end do
          ev = dsqrt(ev)
          do i = 1,3
            z(i,k) = z(i,k)/ev
          end do
        end do
      end if
c.... sort stresses    2-3
10    if(sig0(2).lt.sig0(3)) then
         shelp   = sig0(3)
         do i = 1,3
           zh(i) = z(i,3)
         end do
         sig0(3) = sig0(2)
         sig0(2) = shelp
         do i = 1,3
           z(i,3) = z (i,2)
           z(i,2) = zh(i)
         end do
      end if
c.... sort stresses    1-2
      if(sig0(1).lt.sig0(2)) then
         shelp   = sig0(2)
         do i = 1,3
           zh(i) = z(i,2)
         end do
         sig0(2) = sig0(1)
         sig0(1) = shelp
         do i = 1,3
           z(i,2) = z (i,1)
           z(i,1) = zh(i)
         end do
      end if
c.... sort stresses    2-3
      if(sig0(2).lt.sig0(3)) go to 10
c
      return
      end
c
      subroutine ptrans(ia,itrot,angl,ul,p,s,nel,nen,ndf,nst,isw)
c----------------------------------------------------------------------
c
c      Purpose: transform mixed displacements(if angl) to pure global
c               displacements and other transformations
c               with sloping boundary conditions
c               for angl,eang,vang
c               if angl: global displacements are mixed
c               (angl. dir+ global dirs)
c
c      Inputs:
c         ia(*)     - Degrees of freedom to rotate
c         itrot     - Rotation for dof 4-6
c         angl(*)   - Array of element nodal angles
c         ul(6,6)   - Element solution variables
c         p(*)      - Element vector
c         s(*,*)    - Element matrix
c         nel       - Number of nodes on element
c         nen       - Maximum number nodes/element
c         ndf       - Number dof/node
c         nst       - Dimension of element arrays
c         isw       - Switch:
c                    =1 transform the mixed displ. quantities to global directions
c                    =2 transform the element arrays s+p to mixed values
c                    =3 transform the element array    p to mixed values
c
c      Outputs:
c         ul(*)     - Element solution variables
c         p(*)      - Element vector
c         s(*,*)    - Element matrix
c
c
c      Comments:
c       transformation acts on variables ia(1) and ia(2)
c
c----------------------------------------------------------------------
c.....type declaration for variables
      integer nel,nen,ndf,nst,isw, i1,ij1,ij2,ij3,ij4, i, j,itrot
      real*8  ang,cs,sn,tm,pi
c.....type declaration for arrays
      integer ia(2)
      real*8  angl(*),ul(nst,6),p(ndf,*),s(nst,nst)
      pi=datan(1.d0)*4.d0
      if(ndf.eq.1) return ! added ww
      ij1 = ia(1)
      ij2 = ia(2)
c.... for ndf .ge.6 (transform for rotations)
      ij3 = ij1 + 3
      ij4 = ij2 + 3
c
      if(isw.eq.1) then
c.... transform the mixed to global displacement quantities
        do 100 i = 1,nel
          if(angl(i).ne.0.0d0) then
            ang = angl(i)*pi/180.d0
            cs  = cos(ang)
            sn  = sin(ang)
            do 110 j = 1,6
                tm        = cs*ul(ij1,j) - sn*ul(ij2,j)
                ul(ij2,j) = sn*ul(ij1,j) + cs*ul(ij2,j)
                ul(ij1,j) = tm
              if(ndf.ge.6.and.itrot.eq.0) then
                tm        = cs*ul(ij3,j) - sn*ul(ij4,j)
                ul(ij4,j) = sn*ul(ij3,j) + cs*ul(ij4,j)
                ul(ij3,j) = tm
              end if
110         continue
          end if
          ij1 = ij1 + ndf
          ij2 = ij2 + ndf
          ij3 = ij3 + ndf
          ij4 = ij4 + ndf
100     continue

      else if(isw.eq.2 .or. isw.eq.3) then

c....   transform the element arrays to mixed directions
        i1 = 0
        do 220 i = 1,nel
          if(angl(i).ne.0.0d0) then
            ang = angl(i)*pi/180.d0
            cs  = cos(ang)
            sn  = sin(ang)
c....       transform load vector
              tm       = cs*p(ij1,i) + sn*p(ij2,i)
              p(ij2,i) =-sn*p(ij1,i) + cs*p(ij2,i)
              p(ij1,i) = tm
            if(ndf.ge.6.and.itrot.eq.0) then
              tm       = cs*p(ij3,i) + sn*p(ij4,i)
              p(ij4,i) =-sn*p(ij3,i) + cs*p(ij4,i)
              p(ij3,i) = tm
            end if
            if(isw.eq.2) then
c....         postmultiply s by the transformation
              do 210 j = 1,nst
                  tm         = s(j,i1+ij1)*cs + s(j,i1+ij2)*sn
                  s(j,i1+ij2)=-s(j,i1+ij1)*sn + s(j,i1+ij2)*cs
                  s(j,i1+ij1)= tm
210           continue
              if(ndf.ge.6.and.itrot.eq.0) then
                do 211 j = 1,nst
                  tm         = s(j,i1+ij3)*cs + s(j,i1+ij4)*sn
                  s(j,i1+ij4)=-s(j,i1+ij3)*sn + s(j,i1+ij4)*cs
                  s(j,i1+ij3)= tm
211             continue
              end if
c....         premultiply s by the transformation
              do 215 j = 1,nst
                  tm         = cs*s(i1+ij1,j) + sn*s(i1+ij2,j)
                  s(i1+ij2,j)=-sn*s(i1+ij1,j) + cs*s(i1+ij2,j)
                  s(i1+ij1,j)= tm
215           continue
              if(ndf.ge.6.and.itrot.eq.0) then
                do 216 j = 1,nst
                  tm         = cs*s(i1+ij3,j) + sn*s(i1+ij4,j)
                  s(i1+ij4,j)=-sn*s(i1+ij3,j) + cs*s(i1+ij4,j)
                  s(i1+ij3,j)= tm
216             continue
              end if
            end if
          end if
          i1 = i1 + ndf
220     continue
      end if
c
      end
c
      subroutine ptrans1(ia,itrot,angl,ul,nel,ndf,nst,isw)
c----------------------------------------------------------------------
c
c      Purpose: transform  displacements(if angl)
c               for angl,eang,vang
c
c      Inputs:
c         ia(*)   - Degrees of freedom to rotate
c         itrot   - Rotation for dof 4-6
c         angl(*) - Array of element nodal angles
c         ul(6)   - Element solution variables
c         nel     - Number of nodes on element
c         ndf     - Number dof/node
c         nst     - Dimension of element arrays
c         isw     - Switch:
c                 =1 transform the       mixed displ. to pure global displ.
c                 =2 transform the pure global displ. to       mixed displ.
c
c      Outputs:
c         ul(*)     - Element solution variables
c
c
c      Comments:
c       transformation acts on variables ia(1) and ia(2)
c
c----------------------------------------------------------------------
c.....type declaration for variables
      integer nel,ndf,nst,isw,ij1,ij2,ij3,ij4, i,itrot
      real*8  ang,cs,sn,tm,pi
c.....type declaration for arrays
      integer ia(2)
      real*8  angl(*),ul(nst)
      pi=datan(1.d0)*4.d0
      if(ndf.eq.1) return !
      ij1 = ia(1)
      ij2 = ia(2)
c.... for ndf .ge.6 (transform for rotations)
      ij3 = ij1 + 3
      ij4 = ij2 + 3
c
      do i = 1,nel
        if(angl(i).ne.0.0d0) then
          ang = angl(i)*pi/180.d0
          cs  = cos(ang)
          sn  = sin(ang)
          if(isw.eq.2) sn = -sn
          tm      = cs*ul(ij1) - sn*ul(ij2)
          ul(ij2) = sn*ul(ij1) + cs*ul(ij2)
          ul(ij1) = tm
          if(ndf.ge.6.and.itrot.eq.0) then
            tm      = cs*ul(ij3) - sn*ul(ij4)
            ul(ij4) = sn*ul(ij3) + cs*ul(ij4)
            ul(ij3) = tm
          end if
        end if
        ij1 = ij1 + ndf
        ij2 = ij2 + ndf
        ij3 = ij3 + ndf
        ij4 = ij4 + ndf
      end do
      end
c
      subroutine pzero(v,nn)
c----------------------------------------------------------------------
c
c      Purpose: Zero real array of data
c
c      Inputs:
c         nn     - Length of array
c
c      Outputs:
c         v(*)   - Array with zero values
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension v(nn)
      do 100 n = 1,nn
100   v(n) = 0.d0
      return
      end
c
      subroutine pzeroi(ii,nn)
c----------------------------------------------------------------------
c
c      Purpose: Zero integer array of data
c
c      Inputs:
c         nn     - Length of array
c
c      Outputs:
c         ii(*)  - Array with zero values
c
c----------------------------------------------------------------------
      integer ii(nn)
      do 100 n = 1,nn
        ii(n) = 0
100   continue
      return
      end
c
      subroutine pconsi(ii,nn)
c----------------------------------------------------------------------
c
c      Purpose: Set integer array to itself (e.g. for TIE)
c
c      Inputs:
c         nn     - Length of array
c
c      Outputs:
c         ii(n)   - n
c
c----------------------------------------------------------------------
      integer ii(nn)
      do 100 n = 1,nn
        ii(n) = n
100   continue
      return
      end
c
      subroutine pconsi2(ii,nn,k)
c----------------------------------------------------------------------
c
c      Purpose: Set integer array to value k
c
c      Inputs:
c         nn     - Length of array
c
c      Outputs:
c         ii(n)   - k
c
c----------------------------------------------------------------------
      integer ii(nn)
      do 100 n = 1,nn
        ii(n) = k
100   continue
      return
      end
c
      subroutine prfrst(id,idl,ndf,mm,nad)
c----------------------------------------------------------------------
c
c      Purpose: compute the column height for this equation for SOLV 0
c
c      Input:
c        id(ndf)
c        ndf

c      Output:
c        idl(*)
c
c
c----------------------------------------------------------------------
      integer id(ndf),idl(*)
c
      do 10 j = 1,ndf
        jj = id(j)
        if(jj.gt.0) then
          if(mm.eq.0) mm = jj
          mm = min(mm,jj)
          nad = nad + 1
          idl(nad) = jj
        end if
10    continue
c
      return
      end
c
      subroutine nwprof(jd,neq)
c----------------------------------------------------------------------
c
c      Purpose: set initial column pointers
c
c      Input:
c        neq          - No of unknowns in problem
c
c      Output:
c         jd(ndf,*)   - Equation numbers for each active dof
c
c      Comments:
c        only for standard solver
c
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension jd(*)
c
      jd(1) = 0
      if (neq .le. 1) return
      do 10 n = 2, neq
        jd(n) = jd(n) + jd(n-1)
   10 continue
c
      return
      end
c
      subroutine rstprf(jd,idl,id,ix,ien,inn,ieb,nde,
     1                  ndf,nen1,nen,neq,numnp,numel,edfl)
c----------------------------------------------------------------------
c
c      Purpose: reset profile of element/edge arrays
c
c
c      Inputs:
c       idl            -
c       id(ndf,*)      - Equation numbers for each active dof
c       ix(nen1,numel) - Element conectivity array
c       ien            -
c       inn            -
c       ieb            -
c       nde            -
c       ndf            - number dofs of freedom
c       nen1           - nen+4
c       nen            - max. number of nodes/element
c       neq            - number of equations
c       numnp          - number of nodes
c       numel          - number of elements
c       edfl           - true=edge profile
c
c      Output:
c         jd(*)        - Pointer to row/column ends in 'al' and 'au'.
c                        for solver 0
c
c----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      logical edfl
      dimension jd(*),idl(*),id(ndf,*),ix(nen1,*),ien(*),inn(*),
     +          ieb(nde,*)
c.... zero pointer array
      call pzeroi(jd,neq)
      call perform(0,numel,1,30)
c.... compute column heights
      do 50 n = 1,numel
        mm  = 0
        nad = 0
        do 30 i = 1,nen
          ii = ix(i,n)
c.... set element profile
          if(ii.gt.0) then
            call prfrst(id(1,ii),idl,ndf,mm,nad)
c.... set edge profile
            if(edfl) then
              i1 = ien(ii) +1
              i2 = min(numnp,ii+1)
              i2 = ien(i2)
              do 20 ie = i1,i2
                do 10 j = 1,nen
                  if(ix(j,n).eq.inn(ie)) then
                  call prfrst(ieb(1,ie),idl,nde,mm,nad)
                  end if
10              continue
20            continue
            end if
          end if
30      continue
        if(nad.gt.0) then
          do 40 i = 1,nad
            ii = idl(i)
            jd(ii) = max(jd(ii),ii-mm)
40        continue
        end if
      call perform(n,numel,2,30)
50    continue
      call perform(numel,numel,3,30)
c
      return
      end
c
      subroutine reader(ct,b,nneq)
c----------------------------------------------------------------------
c
c      Purpose: Read nodal displacement and stress values from disk
c
c      Inputs:
c         ct        - Name of array, file rewind, or file name
c         b(*)      - displacement array to read
c         nneq      - Length of array b
c
c      Outputs:
c         none
c
c----------------------------------------------------------------------
      USE cdata
      USE fdata
      USE iodata
      USE iofile
      USE pdata3
cww   USE psize
      USE strnam
      USE doalloc
      logical lflg
      character*4 cc,ct,fname1,yyy*80
      double precision b(*)
      save lflg
      data lflg/.false./
c.... go to processor
      if(ct.eq.'disp') go to 100
      if(ct.eq.'stre') go to 200
      if(ct.eq.'wind') go to 300
      if(ct.eq.'clos') go to 400
c.... open file set with name 'ct'
      fname1 = ct
      inquire(file=fname1,exist=lflg)
      if(lflg) then
                       write(iow,2001) fname1
          if(ior.lt.0) write(*  ,2001) fname1
          open(ird,file=fname1,status='old',form='unformatted',err=990)
      else
          write(yyy,2060) fname1
          call drawmess(yyy,1,0)
      end if
      return
c.... read displacement state from disk
100   if(lflg) then
          read(ird,end=980,err=990) cc
          if(cc.ne.ct) go to 970
          read(ird,end=980,err=990) (b(i),i=1,nneq)
      else
          go to 999
      end if
      return
c.... read nodal stress state from disk
200   if(lflg) then
          read(ird,end=980,err=990) cc
          if(cc.ne.ct) go to 970
          ne = -1
          if(plfl) then
            call ralloc(strea,numnp*npstr,'READ-STRE',plfl)
          end if
          npp = numnp*npstr
          read(ird,end=980,err=990) istv,(strea(i),i=1,npp)
          fl(11) = .true.
      else
          go to 999
      end if
      return
c.... rewind file
300   if(lflg) rewind ird
      return
c.... close file
400   if(lflg) close(ird)
      lflg = .false.
      return
970   write(yyy,2070) ct,cc
        call drawmess(yyy,1,0)
      return
980   write(yyy,2080) ct
        call drawmess(yyy,1,0)
      return
990   write(yyy,2090) ct
        call drawmess(yyy,1,0)
      return
999   write(iow,2099)
        call drawmess(yyy,1,0)
      return
2001  format('   File name for a read has been set to ',a4)
2060  format(' ** ERROR ** File ',a4,' does not exist')
2070  format(' ** ERROR ** READ requested ',a4,' but found ',a4)
2080  format(' ** ERROR ** End of file on a read command for ',a4)
2090  format(' ** ERROR ** on a read command for ',a4)
2099  format(' ** ERROR ** No read file is open.')
      end
c
c
      subroutine readerb(x,ix,ndm,numnp,nen,nen1,numel)
c----------------------------------------------------------------------
c
c      Purpose: Read coordinates and elements from binary file Bname
c               on ios=12
c
c      Inputs:
c          x(ndm,numnp)    - coordinate array
c         ix(nen1,numel)   - element array
c
c      Outputs:
c         none
c      Comments
c      # File ITEST  
c        * coor, elem, bloc einfügen 
c
c      # FEAP(ITEST) 
c        * auf Macroebene: ‘BSYS’ -> schreibt Daten nach BTEST 
c        * NUMNP,NUMEL ablesen
c        * beenden
c      # File ITEST
c        * coor, elem, bloc entfernen, stattdessen ‘BSYS’ einfügen 
c        * NUMNP,NUMEL unter FEAP einfügen! 0,0 geht nicht!!
c      # FEAP(ITEST) 
c        * starten, Daten werden automatisch gelesen   
c
c----------------------------------------------------------------------
      USE comfil
      USE iodata
      USE iofile
      logical lflg
      character*4 ct,yyy*80
      character*229   fbout
      double precision x(ndm,numnp)
      integer ix(nen1,numel)
      data lflg/.true./
c.... open file bname = iname+1=B
      call dochar2(finp,ipos)
      fbout = finp
      call dochar1(fbout,'B',ipos)
      inquire(file=fbout,exist=lflg)
      if(lflg) then
        open(ios,file=fbout,status='old',form='unformatted',err=990)
      else
        write(iow,2001) fbout
        write(*  ,2001) fbout
        stop
      end if
c
c.... read coordinates
      call pzero(x,ndm*numnp)
      read(ios,end=980,err=990) ct
      if(ct.ne.'coor') then
        write(yyy,2070) 'coor',ct
        call drawmess(yyy,1,0)
        return
      end if
      read(ios,end=980,err=990) ((x(i,j),j=1,numnp),i=1,ndm)
c
c.... read elements
      call pzeroi(ix,nen1*numel)
      read(ios,end=980,err=990) ct
      if(ct.ne.'elem') then
        write(yyy,2070) 'elem',ct
        call drawmess(yyy,1,0)
        return
      end if
      read(ios,end=980,err=990) ((ix(   i,j),j=1,numel),i=1,nen) ! nodes
      read(ios,end=980,err=990)  (ix(nen1,j),j=1,numel)          ! mate
c
c.... close file
      close(ios)
      return
980   write(yyy,2080) ct
      call drawmess(yyy,1,0)
      return
990   write(yyy,2090) ct
      call drawmess(yyy,1,0)
      return
2001  format('   File does not exist: ',a229)
2070  format(' ** ERROR ** READ requested ',a4,' but found ',a4)
2080  format(' ** ERROR ** End of file on a read command for ',a4)
2090  format(' ** ERROR ** on a read command for ',a4)
      end
c
!      subroutine reshis(ix,nen1,numel,n1,n2)
!c----------------------------------------------------------------------
!c
!c      Purpose: Initialize t_n+1 history variables from final value
!c               of variables at t_n
!c               this is: copy h2 to h1
!c
!c      Inputs:
!c         ix(nen1,*)  - Element connection/history pointer array
!c                       at position nen+1
!c         nen1        - Dimension of ix array
!c         numel       - Number of elements in mesh
!c         n1          - Pointer in ix to t_n   data
!c         n2          - Pointer in ix to t_n+1 data
!
!c      Outputs:
!c         none        - Output is retained in blank common
!c
!c      Comments:      - nothing to do for nh3, thus back does not work!!
!c
!c----------------------------------------------------------------------
!      USE debugs
!      USE hdata
!      USE iofile
!      implicit double precision (a-h,o-z)
!      dimension ix(nen1,numel)
!
!      if(debug.eq.1) then
!       write(dbgtxt,1000) n1,n2
!1000   format(
!     + ' Debug: move history variables (time,back) from ',i3,' to ',i3)
!       write(*  ,1001) dbgtxt
!       write(iow,1001) dbgtxt
!1001      format(a80)
!          dbgtxt = ' '
!       end if
!      if(n1.eq.1) then
!          gh2=gh1
!      else
!          gh1=gh2
!      endif
!
!      !do n = 1,numel
!      !  nt1 = (n-1)*nhmax+1
!      !  nt1 = ix(n1,n) -1
!      !  nt2 = ix(n2,n) -1
!      !  if(nt2.ne.nt1) then
!      !    nhd = abs(nt2 - nt1)
!      !    do nh = 1,nhd
!      !      m(nt2+nh) = m(nt1+nh)
!      !    end do ! nh
!      !  end if
!      !end do ! n
!      return
!      end
c
      subroutine restrt(fres,b,ix,ndm,ndf,nen1,nneq,isw,asc,asc2)
c----------------------------------------------------------------------
c
c      Purpose: Read/save restart files for resolutions
c
c      Inputs:
c         fres        - Name of restart file to read/save
c         b(*)        - Solution state to save
c         ix(nen1,*)  - Element nodal connections of mesh
c         ndm         - Spatial dimension of mesh
c         ndf         - Number dof/node
c         nen1        - Dimension for ix array
c         nneq        - Total number of parameters in solutions
c         isw         - Switch: = 1 for read; =2 for save.
c         asc         - input/output  in ascii for asc=1
c         asc2        - shared memory FE^2
c
c      Outputs:
c         b(*)        - Solution state read
c         none        - from/to blank common
c
c
c      Comments:
c      Restart data are
c        # general values
c        # displacements
c        # arc length values
c        # <transient fields v,a>
c        # history fields
c        # edge data
c        # <crack values>
c        # <director values>
c
c     calculate length of Restart-File (binary)
c
c----------------------------------------------------------------------
#ifdef __INTEL_
      USE KERNEL32
#endif
      USE arcl
      USE cdata
      USE ddata
      USE dirdat
      USE edgdat
      USE fdata
      USE hdata
      USE iodata
      USE iofile
      USE mdata  ! cww??
      USE ndata
      USE pcrit
      USE prlod
      USE psize
      USE subdt
      USE tdata
      USE doalloc
      implicit double precision (a-h,o-z)
      logical exst,sfl,ldummy
      character fres*229,y*80

#ifdef __INTEL_
      CHARACTER(len=25) :: mapname  ! cww?? 24->25
      INTEGER*8 :: rsmap, rsview
      INTEGER :: lo ! local offset
#endif

cww   character yorn*1
      dimension b(*),ix(nen1,*)
      common m(maxm)
      irecl=0
      iasc = asc
      iasc2 = asc2

      irtyp = 0
      if(iasc2.eq.2) then
        irtyp = 2
      end if

C.... set number of displacement vectors to write or read
      iu = 1      ! short u vector ( only u )
c     iu = 3      ! long  u vector ( u, du, ddu)
c.... check file status
1     inquire(file=fres,exist=exst)

#ifdef __INTEL_

c...  mapping name and record length for shared memory restart file
        CALL feap_shmem_ParseRestartMappingName(fres,mapname)


        if (mapname(8:11).eq.'hrve') then
cDBG            WRITE (*,*) "Skipping hrve restart mapping"
cDBG            READ (*,*)
            RETURN
        end if
      if(irtyp.eq.2)then
c...  open restart file mapping and view (created by macro feap earlier)
        call feap_shmem_OpenMapping(rsmap,mapname)

cDH
cDH   This is a workaround!
cDH
cDH   Micro feap will want to write to other restart files than the
cDH   normal rrve-files when using commands such as 'geom'.
cDH   These restart file names are unknown to macro feap, which must
cDH   open the restart file mappings (they would be closed when exiting
cDH   micro FEAP otherwise).
cDH
cDH   As the data written in hrve-files by 'geom' is overwritten by each
cDH   element, is does not seem to be necessary. Thus, we will skip
cDH   these files altogether, ignoring failed results from OpenMapping.
cDH
      end if
#endif

      if(irtyp.ne.2) then ! only check for file if not using shared memory
      if(.not.exst.and.isw.eq.1) then
        if(ior.gt.0) then
        call drawmess('Restart-File could not be opened.',icolo,1)
        return
        end if
cww     ncol = ipos(fres,229)
cww       write(iow,3002) fres(1:ncol)
cww       if(ior.lt.0) then
cww         write(*,3002) fres(1:ncol)
cww         write(*,3003)
cww10       read (*,1000,err=11,end=12) yorn
cww         go to 13
cww11       call errclr ('RESTRT')
cww         go to 10
cww12       call endclr ('RESTRT',yorn)
cww13       if(yorn.eq.'y' .or. yorn.eq.'Y') then
cww           write(*,3004)
cww20         read(*,1001,err=21,end=22) fres
cww           go to 1
cww21         call errclr ('RESTRT')
cww           go to 20
cww22         call endclr ('RESTRT',fres)
cww           go to 1
cww         end if
cww       end if
        call file_rest(fres)
        ipath = ipos1(fres,229)
        ir    = ipath+1
        if(fres(ir:ir).eq.'r'.or.fres(ir:ir).eq.'R') go to 1
        call drawmess('Restart-File could not be opened.',icolo,1)
        return
      end if
      end if ! irtyp.ne.2

c.... open file
c.DH... only if not using shared memory

      if(irtyp.ne.2) then !not SHMEM
        if(iasc.eq.0) then
          if(exst) then
            open(ios,file=fres,form='unformatted',status='old',err=23)
          else
            open(ios,file=fres,form='unformatted',status='new')
          end if
        else if(iasc.ne.0) then
          if(exst) then
            open(ios,file=fres,form='formatted'  ,status='old',err=23)
          else
            open(ios,file=fres,form='formatted'  ,status='new')
          end if
        end if
        go to 24
23      write(*,25) fres
        return
25      format('Error, the file is not un-/formatted!  ',a229)
24      rewind ios
#ifdef __INTEL_
      else if (irtyp.eq.2) then !SHMEM

c     open full restart file mapping view
        call feap_shmem_OpenView(rsmap,0,rsview)
        irecl = 0 ! starting offset

        if(isw.eq.1) then
c     read from restart mapping
c....     general values,displacements
          call feap_shmem_readGenValsDispl(rsview,irecl,lo,nnpo,nnlo,
     &          nnmo,ndmo,ndfo,nrt,fl(9),ttim,b,nneq,iu)
          irecl = irecl + lo

          if(.not. ((nnpo.eq.numnp).and.(nnlo.eq.numel).and.
     &      (nnmo.eq.nummat).and.(ndmo.eq.ndm).and.(ndfo.eq.ndf)) ) then
                call drawmess(
     &      ' **ERROR** Incorrect information in a restart',1,0)
          else
c....       arc length values
            call feap_shmem_readArcLenVals(rsview,irecl,lo,prop,rlnew,
     &          c0,cs1,cs2,ds0,r,det0,xn)
            irecl = irecl + lo

c....       transient fields v,a,...
            if(fl(9)) then
              call ralloc(trans,nrt*nneq,'READ-TRANS',sfl)  
              call feap_shmem_readTransFields(rsview,irecl,lo,
     &              trans,nv,nrt,nneq,ipr) ! cww?? nv=?
              irecl = irecl + lo
            end if

c....       history fields
            call feap_shmem_readHistFields(rsview,irecl,lo,ttim)
            irecl = irecl + lo

c....       edge data
            if(nde*ned.gt.0) then
              ii = ne5 + 5*m(ne1+numnp-1)*nde*ipr -1
              call feap_shmem_readEdgeData(rsview,irecl,lo,ne5,ii,m)
              irecl = irecl + lo
            end if
            refl = .true.
c....       crack values
            ! read clfl
            call feap_shmem_readSingleInt(rsview,irecl,lo,clfl)
            irecl = irecl + lo

            if(.not.clfl) then
              ! read ncs
              call feap_shmem_readSingleInt(rsview,irecl,lo,ncs)
              irecl = irecl + lo
              clfl = .true.
              nc1 = 1
              nc2 = nc1 + numnp
              nc3 = nc2 + numnp*ncs
              nc4 = nc3 + numnp*2
              ncmax = nc1+numnp*(ncs+5)
              allocate(crit(ncmax))  ! READ-CRIT-nc1 .. nc4
              clfl=.false.
              read(ios) (crit,i=nc1,ncmax)
              call feap_shmem_readCrackValues(rsview,irecl,lo,nc1,ncmax,
     &                  crit)
              irecl = irecl + lo
            end if !.not.clfl

c.....      director values (only 1), all values defined in input-file!
            if(ldir.eq.1) then
              mdirmax = mdir+10*knode*ipr
              call feap_shmem_readDirectorVals(rsview,irecl,lo,mdir,
     &              mdirmax,m)
              irecl = irecl + lo
              call pmove(basea,basea(1+knode),knode)
            end if

          end if !correct restart information
          end if !isw.eq.1

c     write to restart mapping
        if(isw.eq.2) then
c....     save information for restart during iteration

c....     general values,displacements
          call feap_shmem_writeGenValsDispl(rsview,irecl,lo,numnp,numel,
     &          nummat,ndm,ndf,nrt,fl(9),ttim,b,nneq,iu)
          irecl = irecl + lo

c....     arc length values
          call feap_shmem_writeArcLenVals(rsview,irecl,lo,prop,rlnew,c0,
     &          cs1,cs2,ds0,r,det0,xn)
          irecl = irecl + lo

c....     transient fields v,a,...
          if(fl(9)) then
            call feap_shmem_writeTransFields(rsview,irecl,lo,trans,nv
     &              ,nrt,nneq,ipr)  ! cww?? nv=?
            irecl = irecl + lo
          end if
c....     history fields:  numel*2*nh1*ipr
          call feap_shmem_writeHistFields(rsview,irecl,lo,ttim)
          irecl = irecl + lo

c....     edge data 
          if(nde*ned.gt.0) then
            ii = ne5 + 5*m(ne1+numnp-1)*nde*ipr -1
            call feap_shmem_writeEdgeData(rsview,irecl,lo,ne5,ii,m)
            irecl = irecl + lo
          end if

c.....    crack values
          ! write clfl
          call feap_shmem_writeSingleInt(rsview,irecl,lo,clfl)
          irecl = irecl + lo

          if(.not.clfl) then
            ! write ncs separately (needs to be read first later)
            call feap_shmem_writeSingleInt(rsview,irecl,lo,ncs)
            irecl = irecl + lo

            ncmax = nc1+numnp*(ncs+5)
            call feap_shmem_writeCrackValues(rsview,irecl,lo,nc1,
     &              ncmax,crit)
            irecl = irecl + lo
          end if

c.....    director values (only 1)
          if(ldir.eq.1) then
            mdirmax = mdir+10*knode*ipr
            call feap_shmem_writeDirectorVals(rsview,irecl,lo,mdir,
     &              mdirmax,m)
            irecl = irecl + lo
          end if
        end if !isw.eq.2


        call feap_shmem_CloseHandle(rsview)
        call feap_shmem_CloseHandle(rsmap)

        return
#endif
      end if !end SHMEM fork


      if(iasc.eq.0) then
c....   unformatted files
c....   read restart files
        if(isw.eq.1) then
c....     general values,displacements
          read(ios) nnpo,nnlo,nnmo,ndmo,ndfo,nrt,fl(9)
          if((nnpo.eq.numnp).and.(nnlo.eq.numel).and.(nnmo.eq.nummat)
     1        .and.(ndmo.eq.ndm).and.(ndfo.eq.ndf)) then
            read(ios) ttim,(b(i),i=1,nneq*iu)
c....       eigenvalues and eigenvectors    
cww         haben wir nicht gespeichert, neue Version fehlt noch
cww            read(ios) mf,mq,mfmax
cww            if(md.ne.0) then
cww              call pseta(md,mq,         ipr,sfl,'EIGVAL')
cww              call pseta(mv,mq*neq,ipr,sfl,'EIGVEC')
cww              read(ios) (m(i),i=md,md+mq*ipr)
cww              read(ios) (m(i),i=mv,mv+mq*neq*ipr)
cww            end if

c....       arc length values
            read(ios) prop,rlnew,c0,cs1,cs2,ds0,r,det0,xn
c....       transient fields v,a,...
            if(fl(9)) then
              call ralloc(trans,nrt*nneq,'READ-TRANS',sfl)
              read(ios) trans
            end if
c....       history fields
            read(ios) ttim,isgh1,isgh3
            if(size(gh1).ne.isgh1.or.size(gh3).ne.isgh3) then
              write(*,*) 'History data error'
              stop
            end if
            read(ios) gh1,gh2,gh3
c....       edge data
            if(nde*ned.gt.0) then
              ii = ne5 + 5*m(ne1+numnp-1)*nde*ipr -1
              read(ios) edge1,edge2,edge3,edge4
            end if
            refl = .true.
c....       crack values
            read(ios) clfl
            if(.not.clfl) then
              read(ios) ncs
              nc1 = 1
              nc2 = nc1 + numnp
              nc3 = nc2 + numnp*ncs
              nc4 = nc3 + numnp*2
              ncmax = nc1+numnp*(ncs+5)
              allocate(crit(ncmax))  ! READ-CRIT-nc1 .. nc4
              clfl=.false.
              read(ios) crit
            end if
c.....      director values (only 1), all values defined in input-file!
            if(ldir.eq.1) then
              mdirmax = 10*knode
              read(ios) (basea(i),i=1,mdirmax)
              call pmove(basea,basea(1+knode),knode)
            end if
          else
            call drawmess(
     +      ' **ERROR** Incorrect information in a restart',1,0)
          end if
        else if(isw.eq.2) then
c....     save information for restart during iteration
c....     general values,displacements
          write(ios) numnp,numel,nummat,ndm,ndf,nrt,fl(9)
          irecl=irecl+28
          write(ios) ttim,(b(i),i=1,nneq*iu)
          irecl=irecl+8+(nneq*iu)*8
c....     eigenvalues and eigenvectors
cww          write(ios) mf,mq,mfmax
cww          mq = min0(mf+mf,mf+8,neq)
cww          write(ios) mf,mq,mfmax
c            irecl
cww          if(md.ne.0) then
cww            write(ios) (m(i),i=md,md+mq*ipr)
c            irecl
cww            write(ios) (m(i),i=mv,mv+mq*neq*ipr)
c            irecl
cww          end if
c....     arc length values
          write(ios) prop,rlnew,c0,cs1,cs2,ds0,r,det0,xn
          irecl=irecl+72
c....     transient fields v,a,...
          if(fl(9)) then
            write(ios) trans
            irecl=irecl+nrt*nneq*8
          end if
c....     history fields:  numel*(2*nh1+nh3) 
          write(ios) ttim,size(gh1),size(gh3)
          write(ios) gh1,gh2,gh3
          irecl=irecl+12+(size(gh1)*2+size(gh3))*8
c....     edge data
          if(nde*ned.gt.0) then
            ii = ne5 + 5*m(ne1+numnp-1)*nde*ipr -1  
            write(ios) edge1,edge2,edge3,edge4
            irecl=irecl+(ii-ne5)*4
          end if
c.....    crack values
          write(ios) clfl
          irecl=irecl+4
          if(.not.clfl) then
            write(ios) ncs
            irecl=irecl+4
            ncmax = nc1+numnp*(ncs+5)
            write(ios) crit
            irecl=irecl+(ncmax-nc1)*4
          end if
c.....    director values (only 1)
          if(ldir.eq.1) then
            mdirmax = 10*knode
            write(ios) (basea(i),i=1,mdirmax)
            irecl=irecl+mdirmax*8 
          end if
        end if
        write(iow,*) 'IRECL',irecl
c
      else if(iasc.ne.0) then
c....   formatted files
c....   read restart files
        if(isw.eq.1) then
c....     general values,displacements
          read(ios,4020) y
          read(ios,4002) nnpo,nnlo,nnmo,ndmo,ndfo,nrt,fl(9),ttim
          if((nnpo.eq.numnp).and.(nnlo.eq.numel).and.(nnmo.eq.nummat)
     1       .and.(ndmo.eq.ndm).and.(ndfo.eq.ndf)) then
            read(ios,4020) y
            call readf(b,ndf,iu*numnp,ios)
c....       eigenvalues and eigenvectors
cww            read(ios,4020) y
cww            read(ios,4005) mf,mq,mfmax
cww            if(md.ne.0) then
cww              call pseta(md,mq,    ipr,sfl,'READ-EIGVAL')
cww              call pseta(mv,mq*neq,ipr,sfl,'READ-EIGVEC')
cww              read(ios,4020) y
cww              call readf(eigd,1,mq,ios)
cww              read(ios,4020) y
cww              call readf(eigv,mq,neq,ios)
cww            end if
c....       arc length values
            read(ios,4020) y
            read(ios,4009) prop,rlnew,c0,cs1,cs2,ds0,r,det0,xn
            if(fl(9)) then
c....         transient fields v,a,...
              call ralloc(trans,nrt*nneq,'READ-TRANS',sfl)
              read(ios,4020) y
              call readf(trans,ndf,nrt*numnp,ios)
            end if
c....       history fields
            read(ios,4020) y
            read(ios,4013) ttim
            if(allocated(gh1)) then
              call readf(gh1,1,size(gh1),ios)
              call readf(gh2,1,size(gh2),ios)
            end if
            if(allocated(gh3)) then
              call readf(gh3,1,size(gh3),ios)
            end if
c....       edge data
            if(nde*ned.gt.0) then
              ii = ne5 + 5*m(ne1+numnp-1)*nde*ipr -1
              read(ios,4020) y
              call readf(edge4(ne4),1,ii-ne5,ios)
            end if
            refl = .true.
c....       crack values
            read(ios,4020) y
            read(ios,4016) clfl
            if(.not.clfl) then
              read(ios,4005) ncs
              nc1 = 1
              nc2 = nc1 + numnp
              nc3 = nc2 + numnp*ncs
              nc4 = nc3 + numnp*2
              allocate(crit(nc4+2*numnp))  ! READ-CRIT-nc1 .. nc4  cww?? allocate?
              clfl = .false.
              call readf (crit,numnp,ncs+5,ios)
            end if
c.....      director values (only 1)
            if(ldir.eq.1) then
              read(ios,4020) y
              call readf(basea,10,knode,ios)
              call pmove(basea,basea(1+knode),knode)
            end if
          else
            call drawmess(
     +      ' **ERROR** Incorrect information in a restart',1,0)
          end if
        else if(isw.eq.2) then
c....     save information for restart during iteration
c....     general values,displacements
          write(ios,4001)
          write(ios,4002) numnp,numel,nummat,ndm,ndf,nrt,fl(9),ttim
          write(ios,4003)
          call writef(b,ndf,iu*numnp,ios)
c....     eigenvalues and eigenvectors
cww          mq = min0(mf+mf,mf+8,neq)
cww          write(ios,4004)
cww          write(ios,4005) mf,mq,mfmax
cww          if(md.ne.0) then
cww            write(ios,4006)
cww            call writef(eigd,1,mq,ios)
cww            write(ios,4007)
cww            call writef(eigv,mq,neq,ios)
cww          end if
c....     arc length values
          write(ios,4008)
          write(ios,4009) prop,rlnew,c0,cs1,cs2,ds0,r,det0,xn
c....     transient fields v,a,...
          if(fl(9)) then
            write(ios,4010)
            call writef(trans,ndf,nrt*numnp,ios)
          end if
c....     history fields  numel*(2*nhmax+nh3max)   
          write(ios,4012)
          write(ios,4013) ttim
          if(allocated(gh1)) then
            call writef(gh1,1,size(gh1),ios)
            call writef(gh2,1,size(gh2),ios)
          end if
          if(allocated(gh3)) then
            call writef(gh3,1,size(gh3),ios)
          end if
c....     edge data
          if(nde*ned.gt.0) then
            ii = ne5 + 5*m(ne1+numnp-1)*nde*ipr -1
            write(ios,4014)
            call writef(edge4(ne4),1,ii-ne5,ios)
          end if
c....     crack values
          write(ios,4015)
          write(ios,4016) clfl
          if(.not.clfl) then
            write(ios,4005) ncs
            call writef(crit,numnp,ncs+5,ios)
          end if
c.....    director values (only 1)
          if(ldir.eq.1) then
            write(ios,4017)
            call writef(basea,10,knode,ios)
          end if
        end if
      end if
c.... close the file
      close(ios)
      return
cww1000  format(a1)
cww1001  format(a229)
cww3002  format(' **ERROR** Restart file ',a,' does not exist')
cww3003  format(11x,'Specify new name for restart file? (y or n) >',$)
cww3004  format(11x,'New Restart File Name >',$)
4001  format(' numnp,numel,nummat,ndm,ndf,nrt,fl(9),ttim')
4002  format(6i10,l5,e18.9)
4003  format(' displacements (nneq)')
4004  format(' mf,mq,mfmax')
4005  format(5i7)
4006  format(' eigenvalues  (neq)')
4007  format(' eigenvectors (neq)')
4008  format(' prop,rlnew,c0,cs1,cs2,ds0,r,det0,xn')
4009  format(9e18.9)
4010  format(' transient values v,a,... stored with neq, length=nneq')
4012  format(' ttim')
4013  format(e18.9,3i12)
4014  format(' edge data')
4015  format(' crack data')
4016  format(l10)
4017  format(' director data')
4020  format(a80)
      end
c
      subroutine readf(f,n1,n2,ios)
c----------------------------------------------------------------------
c
c      Purpose: read an array from restart file with formats
c
c      Input:
c        f(n1,n2)  - array to read
c        n1        - first  index of array f
c        n2        - second index of array f
c        ios       - Filenumber of file to read
c
c      Output:
c        none
c
c----------------------------------------------------------------------
      real*8 f(n1,n2)
      do i = 1,n2
        read(ios,1000) (f(k,i),k=1,n1)
      end do
1000  format(1x,8e12.5,/,1x,8e12.5)
      return
      end
c
      subroutine writef(f,n1,n2,ios)
c----------------------------------------------------------------------
c
c      Purpose: write an array to restart file with formats
c
c      Input:
c        f(n1,n2)  - array to write
c        n1        - first  index of array f
c        n2        - second index of array f
c        ios       - Filenumber of file to write
c
c      Output:
c        none
c
c----------------------------------------------------------------------
      real*8 f(n1,n2)
      do i = 1,n2
        write(ios,1000) (f(k,i),k=1,n1)
      end do
1000  format(1x,8e12.5,/,1x,8e12.5)
      return
      end
c
c      subroutine readf1(f,nz,numel,ios)
c----------------------------------------------------------------------
c
c      Purpose: read history array from restart file with formats
c
c      Input:
c        f(nz,numel) - history array to read
c        nz          - first  index of array f, nz = (nh*ngp)*2
c        numel       - second index of array f
c        ios         - Filenumber of file to read
c
c      Output:
c        none
c
c----------------------------------------------------------------------
c      real*8 f(nz,numel)
c      nre = int(numel/7)+1
c      do i = 1,nre
c         nsp1 = (i-1)*7+1
c         nsp = 7
c         if(i.eq.nre) nsp=numel-7*(nre-1)
c         do j = 1,nz
c           read(ios,1000) (f(j,k),k=nsp1,nsp1+nsp-1)
c         end do
c      end do
c1000  format(1x,7e12.5)
c      return
c      end
c
c      subroutine writef1(f,nz,numel,ios)
c----------------------------------------------------------------------
c
c      Purpose: write history array to restart file with formats
c
c      Input:
c        f(nz,numel) - history array to write
c        nz          - first  index of array f, nz = (nh*ngp)*2
c        numel       - second index of array f
c        ios         - Filenumber of file to write
c
c      Output:
c        none
c
c----------------------------------------------------------------------
c      real*8 f(nz,numel)
c      nre = int(numel/7)+1
c      do i = 1,nre
c         nsp1 = (i-1)*7+1
c         nsp = 7
c         if(i.eq.nre) nsp=numel-7*(nre-1)
c         do j = 1,nz
c           write(ios,1000) (f(j,k),k=nsp1,nsp1+nsp-1)
c         end do
c      end do
c1000  format(1x,7e12.5)
c      return
c      end
c
      subroutine rot(x,ndm,numnp,prt)
c----------------------------------------------------------------------
c
c      Purpose: rotate cartesian coordinates around axis
c               Xnew = T*Xold
c
c      Input:
c         x(ndm,*)    - Nodal coordinates of mesh
c         ndm         - Spatial dimension of mesh
c         numnp       - Number of nodes in mesh
c      Output:
c         x(ndm,*)    - Nodal coordinates of mesh

c
c----------------------------------------------------------------------
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      dimension x(ndm,*),td(6),tr(3,3),r(3)
      character yyy*80
      mct = 0
      pi  = 4.d0*datan(1.0d0)
      call pzero(tr,9)
      do i = 1,3
        tr(i,i) = 1.d0
      end do
100   if(ior.lt.0) write(*,3001)
      call dinput(td,6)
      if(errck) go to 100
      ni  = td(1)
      ne  = td(2)
      inc = td(3)
      iax = td(4)
      ang = td(5)
      angd= ang*pi/180.d0
      if(ni.le.0) return
      if(ni.gt.numnp.or.ne.gt.numnp) go to 300
      inc = isign(max(iabs(inc),1),ne-ni)
      if(ne.eq.0) ne = ni
      n = ni
c.... rotation matrix
      cn  = cos(angd)
      sn  = sin(angd)
      if(ndm.eq.2.and.(iax.ne.3)) then
          write(yyy,3002) iax
          call drawmess(yyy,1,0)
          return
      end if
      j   = mod(iax,3) + 1
      k   = mod(  j,3) + 1
      do i = 1,3
          te      =  cn*tr(j,i) + sn*tr(k,i)
          tr(k,i) = -sn*tr(j,i) + cn*tr(k,i)
          tr(j,i) =  te
      end do
c...  calculate new values
200   continue
        if(x(1,n).ne. -999.d0) then
        if(ndm.eq.2) then
          do i = 1,2
              r(i) = tr(1,i)*x(1,n)+tr(2,i)*x(2,n)
          end do
        else
            do i = 1,3
              r(i) = tr(1,i)*x(1,n)+tr(2,i)*x(2,n)+tr(3,i)*x(3,n)
          end do
        end if
          x(1,n) = r(1)
          x(2,n) = r(2)
          if(ndm.eq.3) x(3,n) = r(3)
      end if
c
      if(mct.gt.0) go to 250
c
      if(prt)              write(iow,2000) iax,ang,(i,i=1,ndm)
      if(prt.and.ior.lt.0) write(*,  2000) iax,ang,(i,i=1,ndm)
      mct = 50
250   if(prt)              write(iow,2001) n,(x(i,n),i=1,ndm)
      if(prt.and.ior.lt.0) write(*,  2001) n,(x(i,n),i=1,ndm)
      mct = mct - 1
      n = n + inc
      if((ne-n)*inc.ge.0) go to 200
      if(mod(ne-ni,inc).eq.0) go to 100
      ni = ne
      n  = ne
      go to 200
c.... error
300   write(yyy,3000) ni,ne
      call drawmess(yyy,1,0)
      return
c.... formats
2000  format(/,
     1 '  New cartesian coordinates computed from ROT',/,
     2 '  with axis = ',i1,' angl = ',f12.4,/,
     3 4x,'node',3(i6,'-coord'))
2001  format(i8,3f12.4)
3000  format(
     +' Try to use nodes in ROT from',i4,' to ',i4,'(> max node)')
3001  format(' Input: node-1,node-2,inc, axis, angl'/'   >',$)
3002  format(' Rotate around axis ',i4,'not allowed!')
      end
c
      subroutine sblk(nr,ns,xl,ixl,shp,x,ix,dr,ds,ni,ne,ndm,nel1,
     1   nodinc,ntyp,nm,ma,prt)
c----------------------------------------------------------------------
c
c      Purpose: Generate nodes and elements for 2-d problems
c
c      generate a block of elements
c      ntyp = 0   4-node quadrilaterals
c      ntyp = 1   3-node triangles - diags ll to ur
c      ntyp = 2   3-node triangles - diags ul to lr
c      ntyp = 3   3-node triangles - diags to right over 2 el
c      ntyp = 4   3-node triangles - diags to left  over 2 el
c      ntyp = 5   3-node triangles - diags union jack
c      ntyp = 6   3-node triangles - diags union jack invers
c      ntyp = 7   6-node triangles - diags ll to ur
c      ntyp = 8   8-node quadrilaterals
c      ntyp = 9   9-node quadrilaterals
c      ntyp = 16 16-node quadrilaterals
c
c      Inputs:
c         nr        - Number nodes in 1-local coordinate dir.
c         ns        - Number nodes in 2-local coordinate dir.
c         xl(ndm,*) - Block nodal coordinate array
c         ixl(*)    - Block nodal connection list
c         shp(3,*)  - Shape functions for block
c         x
c         ix
c         dr        - 1-local coordinate increment
c         ds        - 2-local coordinate increment
c         ni        - Initial node number for block
c         ne        - Initial element number for block
c         ndm       - Spatial dimension of mesh
c         nel1      - Dimension for ix array
c         nodinc    - Increment to node numbers at end of each line
c         ntyp      - Block type
c         nm        - Number of nodes on block
c         ma        - Material number
c         ctype     - Type of block coordinates
c         prt       - Output generated data if true

c      Outputs:
c         ne        - Final node number for block
c         x(ndm,*)  - Nodal coordinates for block
c         ix(nel1,*)- Element nodal connections of mesh
c
c      Comments:
c         15.05.15 WW IBS KIT to do: set b.c. for unused node for ntyp=8 see vblk
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt,ityp
      character*6, parameter :: xh  = ' coord'
      dimension xl(3,*),ixl(1),x(ndm,*),ix(nel1,*),shp(3,*)

c.... generate nodes

      n = ni
      mct = 0
      s = -1.0d0
      do 200 j = 1,ns
        r = -1.0d0
        do 100 i = 1,nr
          call shape(r,s,xl,shp,xsj,3,nm,ixl,.true.)
          call cfunc(shp,xl,ixl,ndm,x(1,n))
          if(prt) then
            mct = mct + 1
cww         if(mod(mct,50).eq.1) write(iow,2003) o,head,(k,xh,k=1,ndm)
            if(mod(mct,50).eq.1) write(iow,2003) (k,xh,k=1,ndm)
            write(iow,2004) n,(x(k,n),k=1,ndm)
            if(ior.lt.0) then
cww           if(mod(mct,50).eq.1) write(*,2003) o,head,(k,xh,k=1,ndm)
              if(mod(mct,50).eq.1) write(*,2003) (k,xh,k=1,ndm)
              write(*,2004) n,(x(k,n),k=1,ndm)
            end if
          end if
          n = n + 1
        if(n.gt.numnp) go to 201
          r = r + dr
100     continue
        n = n + nodinc
      if(n.gt.numnp) go to 201
        s = s + ds
200   continue
201   continue

c.... generate elements

      if(ne.gt.0) then
        me = ne - 1
        inc = 1
        if(ntyp.ge.7.and.ntyp.le.9) inc = 2
        if(ntyp.eq.16) inc = 3
        do 400 j = 1,ns-1,inc
          n = (nr+nodinc)*(j-1) + ni
          do 300 i = 1,nr-1,inc
            n = n + 1
            me = me + 1
          if(me.gt.numel) go to 401
            ix(nel1,me) = ma
            if(ntyp.eq.0) then
              ix(1,me)    = n - 1
              ix(2,me)    = n
              if(ndm.ne.1) then
                    ix(3,me)    = n + nr + nodinc
                    ix(4,me)    = n + nr - 1 + nodinc
              end if
            else if(ntyp.eq.7) then
              ix(1,me) = n-1
              ix(4,me) = n
              ix(2,me) = n+1
              ix(6,me) = nr+nodinc + n
              ix(5,me) = nr+nodinc + n+1
              ix(3,me) = 2*(nr+nodinc) + n+1
              me = me + 1
            if(me.gt.numel) go to 401
              ix(1,me) = n-1
              ix(6,me) = nr+nodinc + n-1
              ix(4,me) = nr+nodinc + n
              ix(3,me) = 2*(nr+nodinc) + n-1
              ix(5,me) = 2*(nr+nodinc) + n
              ix(2,me) = 2*(nr+nodinc) + n+1
              ix(nel1,me) = ma
              n = n+1
            else if(ntyp.eq.8.or.ntyp.eq.9) then
              ix(1,me) = n-1
              ix(5,me) = n
              ix(2,me) = n+1
              ix(8,me) = nr+nodinc + n-1
              if(ntyp.gt.8) ix(9,me) = nr+nodinc + n
              ix(6,me) = nr+nodinc + n+1
              ix(4,me) = 2*(nr+nodinc) + n-1
              ix(7,me) = 2*(nr+nodinc) + n
              ix(3,me) = 2*(nr+nodinc) + n+1
              n = n+1
            else if(ntyp.eq.16) then

              ix( 1,me) = n - 1
              ix( 5,me) = n
              ix( 6,me) = n + 1
              ix( 2,me) = n + 2

              ix(12,me) =    nr+nodinc  + n - 1
              ix(13,me) =    nr+nodinc  + n
              ix(14,me) =    nr+nodinc  + n + 1
              ix( 7,me) =    nr+nodinc  + n + 2

              ix(11,me) = 2*(nr+nodinc) + n - 1
              ix(16,me) = 2*(nr+nodinc) + n
              ix(15,me) = 2*(nr+nodinc) + n + 1
              ix( 8,me) = 2*(nr+nodinc) + n + 2

              ix( 4,me) = 3*(nr+nodinc) + n - 1
              ix(10,me) = 3*(nr+nodinc) + n
              ix( 9,me) = 3*(nr+nodinc) + n + 1
              ix( 3,me) = 3*(nr+nodinc) + n + 2

crlt neu
c             ix( 4,me) = n - 1
c             ix( 3,me) = n + 2
c             ix( 2,me) = 3*(nr+nodinc) + n + 2
c             ix( 1,me) = 3*(nr+nodinc) + n - 1
c             ix(10,me) = n
c             ix( 9,me) = n + 1
c             ix( 8,me) =    nr+nodinc  + n + 2
c             ix( 7,me) = 2*(nr+nodinc) + n + 2
c             ix( 6,me) = 3*(nr+nodinc) + n + 1
c             ix( 5,me) = 3*(nr+nodinc) + n
c             ix(12,me) = 2*(nr+nodinc) + n - 1
c             ix(11,me) =    nr+nodinc  + n - 1
c             ix(16,me) =    nr+nodinc  + n
c             ix(15,me) =    nr+nodinc  + n + 1
c             ix(14,me) = 2*(nr+nodinc) + n + 1
c             ix(13,me) = 2*(nr+nodinc) + n

              n = n + 2

            else
              ityp = (ntyp.eq.1) .or. (ntyp.eq.3.and.mod(j,2).eq.1)
     1          .or. (ntyp.eq.4.and.mod(j,2).eq.0) .or. (ntyp.eq.5.and.
     2          mod(i+j,2).eq.0) .or. (ntyp.eq.6.and.mod(i+j,2).eq.1)
              if(ityp) then
                    ix(1,me) = n-1
                    ix(2,me) = n + nr + nodinc
                    ix(3,me) = n + nr + nodinc - 1
                    me = me + 1
              if(me.gt.numel) go to 401
                    ix(1,me) = n-1
                    ix(2,me) = n
                    ix(3,me) = n + nr + nodinc
                    ix(nel1,me) = ma
              else
                    ix(1,me) = n-1
                    ix(2,me) = n
                    ix(3,me) = n + nr + nodinc - 1
                    me = me + 1
              if(me.gt.numel) go to 401
                    ix(1,me) = n
                    ix(2,me) = n + nr + nodinc
                    ix(3,me) = n + nr + nodinc - 1
                    ix(nel1,me) = ma
              end if
            end if
300       continue
400     continue
401   continue
      end if
      return
cww2003  format(a1,19a4,a3/'  n o d a l   c o o r d i n a t e s'/
cww     1    6x,'node',5(i7,a6))
2003  format(/'  n o d a l   c o o r d i n a t e s'/
     1    6x,'node',5(i7,a6))
2004  format(i10,5f13.4)
      end
c
c
      subroutine sblkco(nr,ns,xl,ixl,shp,x,ix,dr,ds,ni,ne,ndm,nel1,
     1   nm,ma1,ma2,ma3,mtyp,prt)
c----------------------------------------------------------------------
c
c      Purpose: Generate nodes and elements for 2-d problems cohesive zone
c
c      Inputs:
c         nr        - Number nodes in 1-local coordinate dir.
c         ns        - Number nodes in 2-local coordinate dir.
c         xl(ndm,*) - Block nodal coordinate array
c         ixl(*)    - Block nodal connection list
c         shp(3,*)  - Shape functions for block
c         x
c         ix
c         dr        - 1-local coordinate increment
c         ds        - 2-local coordinate increment
c         ni        - Initial node number for block
c         ne        - Initial element number for block
c         ndm       - Spatial dimension of mesh
c         nel1      - Dimension for ix array
c         nm        - Number of nodes on block
c         ma1       - Material number
c         ma2       - Material number
c         ma3       - Material number
c         mtyp      - 1: 4-node shell elmts and  8-node interface elmt
c                    -2: 4-node shell elmts and 18-node interface elmt
c                    -3: 4-node shell elmts and 16-node interface elmt
c         ctype     - Type of block coordinates
c         prt       - Output generated data if true

c      Outputs:
c         ne        - Final node number for block
c         x(ndm,*)  - Nodal coordinates for block
c         ix(nel1,*)- Element nodal connections of mesh
c
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      character*6, parameter :: xh  = ' coord'
      dimension xl(3,*),ixl(1),x(ndm,*),ix(nel1,*),shp(3,*)

      n   = ni                ! nodes bottom
      n2  = n + nr*ns         ! nodes top

      mct = 0

c.... coordinates

      s = -1.0d0
      do j = 1,ns
        r = -1.0d0
        do i = 1,nr
c....     set coordinates of bottom node
          call shape(r,s,xl,shp,xsj,3,nm,ixl,.true.)
          call cfunc(shp,xl,ixl,ndm,x(1,n))

c....     set coordinates of top node
          do idm = 1,ndm
            x(idm,n2) = x(idm,n)
          end do

          if(prt) then
            mct = mct + 1
            if(mod(mct,50).eq.1) write(iow,2003) (k,xh,k=1,ndm)
            write(iow,2004) n , n2, (x(k,n),k=1,ndm)
            if(ior.lt.0) then
              if(mod(mct,50).eq.1) write(*,2003) (k,xh,k=1,ndm)
              write(*,2004) n, n2, (x(k,n),k=1,ndm)
            end if
          end if

          n   = n + 1
          n2  = n + nr*ns
          r = r + dr
        end do ! i
        s = s + ds
      end do ! j

c.... elements

      nd =  nr   * ns    ! difference node    numbers
      md = (nr-1)*(ns-1) ! difference element numbers

      if(ne.gt.0) then
        me  = ne - 1

c....   generate shell elements
        inc = 1
        do j = 1,ns-1,inc
          n  = nr*(j-1) + ni

          do i = 1,nr-1,inc
            n  = n + 1
            n2 = n + nd

            me  = me + 1
            me2 = me + md

c....       element bottom
            ix(nel1,me) = ma1
            ix(1,me)    = n - 1
            ix(2,me)    = n
            ix(3,me)    = n + nr
            ix(4,me)    = n + nr - 1

c....       element top
            ix(nel1,me2) = ma2
            ix(1,me2)    = n2 - 1
            ix(2,me2)    = n2
            ix(3,me2)    = n2 + nr
            ix(4,me2)    = n2 + nr - 1

          end do ! i
        end do ! j


c....   generate interface elements
        me3 = me2
        if(mtyp.eq.1) then ! 8-node
          inc = 1
          do j = 1,ns-1,inc
            n  = nr*(j-1) + ni

            do i = 1,nr-1,inc
              n  = n + 1
              n2 = n + nd
              me3 = me3 + 1
c....         element interface
              ix(nel1,me3) = ma3
c....         bottom edge
              ix(1,me3)    = n - 1
              ix(2,me3)    = n
              ix(3,me3)    = n + nr
              ix(4,me3)    = n + nr - 1
c....         top edge
              ix(5,me3)    = n2 - 1
              ix(6,me3)    = n2
              ix(7,me3)    = n2 + nr
              ix(8,me3)    = n2 + nr - 1
            end do ! i
          end do ! j
        else if(mtyp.eq.2.or.mtyp.eq.3) then ! 18/16-node
          inc = 2
          do j = 1,ns-1,inc
            n0 = nr*(j-1) + ni

            do i = 1,nr-1,inc
              n  = n0 + i
              n2 = n + nd
              me3 = me3 + 1
c....         element interface
              ix(nel1,me3) = ma3
c....         bottom edge
              ix( 1,me3)   =  n - 1
              ix( 2,me3)   =  n + 1
              ix( 3,me3)   =  n + 1 + nr*2
              ix( 4,me3)   =  n - 1 + nr*2
c....         top edge
              ix( 5,me3)   = n2 - 1
              ix( 6,me3)   = n2 + 1
              ix( 7,me3)   = n2 + 1 + nr*2
              ix( 8,me3)   = n2 - 1 + nr*2
c....         bottom interior
              ix( 9,me3)   =  n
              ix(10,me3)   =  n + 1 + nr
              ix(11,me3)   =  n     + nr*2
              ix(12,me3)   =  n - 1 + nr
              if(mtyp.eq.2) then ! 18-node
c....           bottom central
                ix(13,me3)   =  n     + nr
c....           top interior
                ix(14,me3)   = n2
                ix(15,me3)   = n2 + 1 + nr
                ix(16,me3)   = n2     + nr*2
                ix(17,me3)   = n2 - 1 + nr
c....           top central
                ix(18,me3)   = n2     + nr
              else if(mtyp.eq.3) then ! 16-node
c....           top interior
                ix(13,me3)   = n2
                ix(14,me3)   = n2 + 1 + nr
                ix(15,me3)   = n2     + nr*2
                ix(16,me3)   = n2 - 1 + nr
              end if
            end do ! i
          end do ! j
        end if
      end if
      return
2003  format(/'  n o d a l   c o o r d i n a t e s'/
     1    4x,'node B','    node T',5(i7,a6))
2004  format(i10,i10,5f13.4)
      end
c
c
      subroutine sblkdx(nr,ns,xl,ixl,shp,x,ix,dr,ds,ni,ne,ndm,nel1,
     1   nodinc,ntyp,nm,ma,prt,dr1,ds1,dmar,dmas)
c----------------------------------------------------------------------
c
c      Purpose: Generate nodes and elements for 2-d problems for BLOX
c
c      generate a block of elements
c      ntyp = 0   4-node quadrilaterals
c
c      Inputs:
c         nr        - Number nodes in 1-local coordinate dir.
c         ns        - Number nodes in 2-local coordinate dir.
c         xl(ndm,*) - Block nodal coordinate array
c         ixl(*)    - Block nodal connection list
c         shp(3,*)  - Shape functions for block
c         x
c         ix
c         dr        - 1-local coordinate increment
c         ds        - 2-local coordinate increment
c         ni        - Initial node number for block
c         ne        - Initial element number for block
c         ndm       - Spatial dimension of mesh
c         nel1      - Dimension for ix array
c         nodinc    - Increment to node numbers at end of each line
c         ntyp      - Block type
c         nm        - Number of nodes on block
c         ma        - Material number
c         prt       - Output generated data if true
c         dr1(100)  - increments in r-direction
c         ds1(100)  - increments in s-direction
c         dmar(100) - material   in r-direction
c         dmas(100) - material   in s-direction
c
c      Outputs:
c         ne        - Final node number for block
c         x(ndm,*)  - Nodal coordinates for block
c         ix(nel1,*)- Element nodal connections of mesh
c
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      character*6 xh
      dimension xl(3,*),ixl(1),x(ndm,*),ix(nel1,*),shp(3,*)
      dimension dr1(100),ds1(100),dmar(100),dmas(100)
      data xh/' coord'/
      n = ni
      mct = 0
      s = -1.0d0
      do 200 j = 1,ns
        r = -1.0d0
        do 100 i = 1,nr
          call shape(r,s,xl,shp,xsj,3,nm,ixl,.true.)
          call cfunc(shp,xl,ixl,ndm,x(1,n))
          if(prt) then
            mct = mct + 1
            if(mod(mct,50).eq.1) write(iow,2003) (k,xh,k=1,ndm)
            write(iow,2004) n,(x(k,n),k=1,ndm)
            if(ior.lt.0) then
              if(mod(mct,50).eq.1) write(*,2003) (k,xh,k=1,ndm)
              write(*,2004) n,(x(k,n),k=1,ndm)
            end if
          end if
          n = n + 1
          if(n.gt.numnp) go to 201
cww          r = r + dr1(i)
          if(i.lt.nr) r = r + dr1(i)
100     continue
        n = n + nodinc
        if(n.gt.numnp) go to 201
cww        s = s + ds1(j)
        if(j.lt.ns) s = s + ds1(j)
200   continue
201   continue
      if(ne.gt.0) then
        me = ne - 1
        inc = 1
        do 400 j = 1,ns-1,inc
          n = (nr+nodinc)*(j-1) + ni
          mat2 = dmas(j)
          do 300 i = 1,nr-1,inc
            n = n + 1
            me = me + 1
            if(me.gt.numel) go to 401
            mat1 = dmar(i)
            mat=max(mat1,mat2)
            ix(nel1,me) = mat
            ix(1,me)    = n - 1
            ix(2,me)    = n
            if(ndm.ne.1) then
              ix(3,me)  = n + nr + nodinc
              ix(4,me)  = n + nr - 1 + nodinc
            end if
300       continue
400     continue
401   continue
      end if
      return
2003  format(/'  n o d a l   c o o r d i n a t e s'/
     1    6x,'node',5(i7,a6))
2004  format(i10,5f13.4)
      end
c
      subroutine scalev(v,nn)
c----------------------------------------------------------------------
c
c      Purpose: Scale vector to have maximum element of 1.0
c
c      Inputs:
c         v(*)   - Vector of values
c         nn     - Space dimension of v
c
c      Outputs:
c         v(*)   - Vector of values
c
c----------------------------------------------------------------------
      implicit none
      real*8 v(*), vmax
      integer n, nn
      vmax = dabs(v(1))
      do 100 n = 1,nn
100   vmax = dmax1(vmax,dabs(v(n)))
      do 110 n = 1,nn
110   v(n) = v(n)/vmax
      return
      end
c
      subroutine serchl(g0,f0,f,id,rsd,u,d,stol,t,neq,nneq, step,pfr)
c----------------------------------------------------------------------
c
c      Purpose: Line search driver program, routine makes a line search
c      in direction 'd' and returns the steplength size in 'step'
c
c      Inputs:
c         g0        - Line search reference value
c         f0(*)     - Nodal initial force values
c         f(*)      - Nodal load vector f
c         id(*)     - Equation numbers for each dof
c         rsd       - Current residual
c         u         - Current solution
c         stol      - Line search convergence tolerance
c         t(*)      - Working solution storage
c         neq       - Number of active equations
c         nneq      - Total number of equations in problem
c         pfr       - Print flag

c      Outputs:
c         d(*)      - Line search vector * step
c         step      - Step size for line search
c
c----------------------------------------------------------------------
      USE iofile
      implicit none
      real*8 f0(*),f(*),rsd(*),u(*),d(*),t(*)
      real*8 g,g0,ga,gb,gamma1,stol, step, sb,sa
      logical pfr
      integer j,id(*),nneq,neq
      integer, parameter :: linmax=10
      sa=1.0d0
      g = gamma1(f0,f,id,u,rsd,d,t,sa,neq,nneq)
c.... find bracket on zero
      if(g*g0.gt.0.0d0) then
                                write(iow,3000)
         if(pfr .and. ior.lt.0) write(*  ,3000)
      else
c.... perform 'linmax' steps of line-search
        j = 0
        ga=g
        gb=g0
        sb = 0.0d0
10      j = j + 1
        if (j.le.linmax) then
          step = sa - ga*(sa-sb)/(ga-gb)
          g = gamma1(f0,f,id,u,rsd,d,t,step,neq,nneq)
c.... output line-search parameters
                                 write(iow,3001) j,step,g,sa,sb
          if(pfr .and. ior.lt.0) write(*  ,3001) j,step,g,sa,sb
c.... update postions for next iteration
          gb = 0.5d0*gb
          if (g*ga.lt.0.0d0) then
            sb = sa
            gb = ga
          end if
          sa = step
          ga = g
c.... check convergence
          if(dabs(g).gt.stol*dabs(g0)) go to 10
c.... THIS CHECK CAUSED PROBLEMS AT TIMES WHEN AN ADEQUATE SOLN FOUND
c.... CAN NOT REMEMBER WHY ADDED IN FIRST PLACE!
c         if(dabs(sb-sa).gt.stol*0.5*(sa+sb)) go to 10
        end if
c.... multiply solution increment by step size computed
        do 20 j = 1,neq
          d(j) = step*d(j)
20      continue
      end if
      return
c.... format statements
3000  format(' -> No line search - Energy positive at full step.')
3001  format(' -> Iteration',i3,' Step Size =',e12.5,' Energy =',e12.5,
     1       ' sa=',f5.3,' sb=',f5.3)
      end
c
      subroutine setdibc(ix,u,numnp,numel,ndf,nen,nen1,n1,n2)
c-------------------------------------------------------------
c
c.... Purpose: set values of displacement array
c
c     Input:
c         ix(nen1,numel) - Element nodal connections of mesh
c         u(ndf,3*numnp) - displacement vector
c         numnp          - Number of nodes in mesh
c         numel          - Number of elements in mesh
c         ndf            - Number dof/node
c         nen1           - Dimension for ix array  nen1 = nen+4
c         n1,n2          - Material numbers for change
c
c     Output:
c         u(ndf,3*numnp) - displacement vector modified
c
c     W. Wagner BS KIT 05/11
c-------------------------------------------------------------
      implicit none
      real *8 u(ndf,3*numnp)
      integer ix(nen1,numel)
      integer numnp, numel, ndf, nen, nen1, n1, n2
      integer ma, kk, jj, ii
c
cww      include 'phasefield.h' ! phiA,phiB
c
      do ii=1,numel          ! element loop
        ma=ix(nen1,ii)       ! mat.no.
        if(ma.eq.n1) then
          do jj=1,nen
            kk = ix(jj,ii)   ! element node
cww            u(7,kk) = phiA   ! set phiA to displacement 7
          end do ! jj
        elseif(ma.eq.n2) then
          do jj=1,nen
            kk = ix(jj,ii)   ! element node
cww            u(7,kk) = phiB   ! set phiB to displacement 7
          end do ! jj
        end if ! ma
      end do ! numel
      return
      end
c
      subroutine sethis(ie,ix,nie,nen,nen1,numel,nummat,prt)
c----------------------------------------------------------------------
c
c      Purpose: Set up history addresses in ix array
c
c      Inputs:
c         ie(nie,*) - Material set assembly information
c         nie       - Dimension of ie array
c         nen       - Number of nodes/element
c         nen1      - Dimension of ix array
c         numel     - Number of elements in mesh
c         nummat    - Number of material sets in mesh
c         prt       - Flag, output results if true

c      Outputs:
c         ix(nen1,*)- History data pointers added to positions nen+1
c                     and nen+2
c      Comments
c         nie = ndf+3
c         ie(1-ndof,*) - local element dofs 
c         ie(nie-1,ma) - iel
c         ie(nie-2,ma) - nh3
c         ie(nie,  ma) - nh1
c
c         nen1 = nen+4
c         ix(1-nen,n)    - Element nodes
c         ix(nen+1,n)    - nt1  global adress history terms h1m
c         ix(nen+2,n)    - nt2  global adress history terms h2m
c         ix(nen+3,n)    - nt3  global adress history terms h3m
c         ix(nen+4,n)    - ma   nen+4=nen1
c
c----------------------------------------------------------------------
      USE doalloc
      USE hdata
      USE hdatam
      USE iofile
      implicit none
      logical flg,prt
cww   integer ie(nie,*),ix(nen1,*)
      integer ie(nie,nummat),ix(nen1,numel)
      integer iel,ma,n,nummat,numel,nie,nen1,nen
      integer ndh,nh0,nh30,nh1e,nh3e

c.... output the amount of memory used in each element
      if(prt) then
        write(  *,2000)
        write(iow,2000)
        do n = 1,nummat
          nh1e = ie(nie,n)
          iel  = ie(nie-1,n)
          nh3e = ie(nie-2,n)
          write(  *,2001) n,iel,nh1e,nh3e
          write(iow,2001) n,iel,nh1e,nh3e
        end do
      end if

c.... compute maximum length necessary to store history variables
      nhmax  = 0
      nh3max = 0
      do n = 1,nummat
        nh1e = ie(nie,n)
        nh3e = ie(nie-2,n)
        nhmax  = max(nhmax,nh1e)
        nh3max = max(nh3max,nh3e)
      end do

c.... set the pointers for history variables
      nh0  = 1
      nh30 = 1

c.... absolute position of nt1,nt2,nt3 for each element in GH1/GH2/GH3
cww      do n = 1,numel
cww        ma   = ix(nen1,n)
cww        nh1e = ie(nie,ma)
cww        nh3e = ie(nie-2,ma)
cww        ix(nen+1,n) = nh0      ! nt1
cww        ix(nen+2,n) = nh0      ! nt2
cww        nh0  = nh0 + nh1e      ! compressed storage, else .. + nhmax   
cww        ix(nen+3,n) = nh30     ! nt3
cww        nh30 = nh30 + nh3e     ! compressed storage, else .. + nh3max
cww      end do

      do n = 1,numel
        ma   = ix(nen1,n)
        ix(nen+1,n) = nh0      ! nt1
        ix(nen+2,n) = nh0      ! nt2
        nh0  = nh0 + nhmax  
        ix(nen+3,n) = nh30     ! nt3
        nh30 = nh30 + nh3max
      end do

c.... pointers + space for history terms
      call ralloc(gh1, nhmax*numel,'ELement History Data NH1',flg)
      call ralloc(gh2, nhmax*numel,'ELement History Data NH2',flg)
      call ralloc(gh3,nh3max*numel,'ELement History Data NH3',flg)

      return

2000  format(10x,'Material  Element Type  Hist.terms H1  Hist.terms H3')
2001  format(4i15)
      end
c
      subroutine setval(xi,num, val)
c----------------------------------------------------------------------
c
c      Purpose: Represent character constants by free inputs in strings
c      Inputs:
c        xi - input string - num-characters long
c        xs - string to search
c        xt - temporary string
c
c      Outputs:
c         val       - Value of string
c
c----------------------------------------------------------------------
      USE codat
      USE conval
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical errco
      character*1 xt(75),xs(75),xi(num)
      dimension v(25)
c.... check length of string
      if(num.gt.75) then
         write(*,*) 'Warning: string too long'
         write(*,*) xi
         stop
      end if
c.... read the value if no constants have been set
      if(.not.coflg) then
        do 40 i = 1,15
          xs(i) = ' '
40      continue
        nex = 15 - num
        do 50 i = 1,num
          xs(i+nex) = xi(i)
50      continue
        errco = .false.
        call pinval(xs,val,errco)
        if(errco) go to 60
        return
      end if
c.... find the letter number for this parameter
60    do 100 i = 1,75
        xs(i) = ' '
100   continue
      do 110 i = 1,num
        xs(i) = xi(i)
110   continue
c.... evaluate expression
1     nex = 0
      call parexp(xs,xt,v,nex,errck)
      if(errck) go to 150
      call pfuncs(xs,v,val,nex,errck)
      if(errck) go to 150
      return
c.... an error has been detected in the statement respecify
150   if(ior.lt.0) then
        write(*,2001) (xi(i),i=1,num)
151     read (*,1000,err=152,end=153) xt
        go to  154
152     call  errclr ('SETVAL')
        go to  151
153     call  endclr ('SETVAL',xt)
154     write(iow,2002) xt
        call pcheck(1,xt,errck)
        do 160 i = 1,74
          xs(i) = xt(i)
160     continue
        xs(75) = ' '
        errck = .false.
        go to 1
c.... error on a read
      else
        call  errclr ('SETVAL')
      end if
      return
c.... formats
 1000 format(75a1)
 2001 format(2x,a1,' = ',74a1/'   >',$)
 2002 format('  Correction:>',75a1)
      end
c
c
      subroutine setval2(xi,num, val)
c----------------------------------------------------------------------
c
c      Purpose: Represent character constants by free inputs in strings
c               Version without CONS
c
c      Inputs:
c        xi - input string - num-characters long
c        xs - string to search
c        xt - temporary string
c
c      Outputs:
c         val       - Value of string
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character*1 xs(75),xi(num)
c.... check length of string
      if(num.gt.75) then
         write(16,*) 'Warning: string too long'
         write(16,*) xi
         stop
      end if
c.... read the value if no constants have been set
      do 40 i = 1,15
        xs(i) = ' '
40    continue
      nex = 15 - num
      do 50 i = 1,num
        xs(i+nex) = xi(i)
50    continue
      call pinval2(xs,val)
      return
      end
c
      subroutine pinval(xss,val,errco)
c----------------------------------------------------------------------
c
c      Purpose:  Moves character string into real value
c
c      Inputs:
c         xss(*)  - Character string
c
c      Outputs:
c         val     - Value extracted from character string
c         error   - Flag, true if error occurs
c
c----------------------------------------------------------------------
      logical errco
      character*15 xs
      character*1  xss(15)
      real*8 val
      do 10 i=1,15
10    xs(i:i) = xss(i)
      read(xs,1000,err=100) val
      return
100   errco = .true.
      return
1000  format(f15.0)
      end
c
      subroutine pinval2(xss,val)
c----------------------------------------------------------------------
c
c      Purpose:  Moves character string into real value
c                Version without error flag
c      Inputs:
c         xss(*)  - Character string
c
c      Outputs:
c         val     - Value extracted from character string
c
c----------------------------------------------------------------------
      character*15 xs
      character*1  xss(15)
      real*8 val
      do 10 i=1,15
10    xs(i:i) = xss(i)
      read(xs,1000) val
      return
1000  format(f15.0)
      end
