      subroutine nodew(ndw,msum,ix,nen,nen1,numnp,numel,nstart)
c----------------------------------------------------------------------
c      Purpose: Compute node and element weights for profile
c               minimizations

c      Inputs:
c         ix(nen1,*)     - Element nodal connection list
c         nen            - Maximum number nodes/element
c         nen1           - Dimension of ix  array
c         numnp          - Number of nodes in mesh
c         numel          - Number of elements in mesh

c      Outputs:
c         ndw(*)         - Node weights
c         msum(*)        - Element weights
c         nstart
c----------------------------------------------------------------------
c
      integer k, l, m,minw, n,nel,nen,nen1,numnp,numel,nstart,nsum    
      integer ndw(numnp),msum(numel),ix(nen1,numel)
      intrinsic abs, min

c---- evaluation of node weights -----------------------
      do 100 n=1,numnp
        ndw(n) = 1
  100 continue
      minw = 1

c---- loop five times to evaluate node weights ----------
      do 600 k=1,5
c---- (a) normalize node weights ------------------------
        do 150 n=1,numnp
          ndw(n) = ndw(n) / minw
  150   continue
c---- (b) compute element weights -----------------------
        do 200 m=1,numel
          if(msum(m).gt.0) then
c
            nsum = 0
            nel  = 0
            do 180 l=1,nen
              n = abs(ix(l,m))
              if(n.gt.0) then
                nel  = nel + 1
                nsum = nsum + ndw(n)
              end if
  180       continue
c
            nsum     = min(nsum,numnp)
            msum(m)  = nsum / nel
          end if
c
  200   continue
c---- (c) compute node weights --------------------------
        call pzeroi(ndw, numnp)
c
        do 400 m=1,numel
          if (msum(m).gt.0) then
            nsum = 0
c
            do 380 l=1,nen
              n = abs(ix(l,m))
              if(n.gt.0) then
                ndw(n) = ndw(n) + msum(m)
              end if
  380       continue
          end if
c
  400   continue
c---- (d) find minimum value ---------------------------
        minw = 32000
        do 500 n=1,numnp
          if (ndw(n).gt.0 .and. ndw(n).le.minw) then
            minw   = ndw(n)
            nstart = n
          end if
  500   continue
c
  600 continue
c
      end
c
      subroutine nodel(ix,nd,ln,ne,numnp,numel,nsum,nen,nen1,nbn)
c-----------------------------------------------------------------------
c      Purpose: Formation of "nd" and "ln" arrays for profile
c               optimization.
c
c      Inputs:
c         ix(nen1,*)     - Element connection array
c         numnp          - Number of nodes in mesh
c         numel          - Number of elements in mesh
c         nen            - Maximum number nodes/element
c         nen1           - Dimension of ix  array

c      Outputs:
c         nd(*)          - Nodal-element array
c         ln(*)          - Location array
c         ne(*)          - Element order to minimize front
c         nsum           - Maximum value in location array
c
c-----------------------------------------------------------------------
c
      integer l, m, n,nn,nen,nen1,numnp,numel,nsum
      integer nd(numnp),ix(nen1,numel),ln(numnp),ne(*)
      intrinsic abs

c---- formation of "nd" and "ln" arrays -------------
      call pzeroi(nd,numnp)
c---- count elements or constraints attached to nodes
      do 200 m=1,numel
          do 150 l=1,nen
            n = abs(ix(l,m))
            if(n.gt.0) then
              nd(n) = nd(n) + 1
            end if
  150     continue
  200 continue
c---- form location array ----------------------------
      nsum = 0
      do 300 n=1,numnp
          nsum  = nsum + nd(n)
          ln(n) = nsum
  300 continue
c---- form nodel-element array -----------------------
      do 400 m=1,numel
          do 350 l=1,nen
            n  = abs(ix(l,m))
            if(n.gt.0) then
              nd(n)  = nd(n) - 1
              nn     = ln(n) - nd(n)
 
            if(nn.gt.numnp*nbn) 
     +      call drawmess('Too few storage for OPTI in NODEL',1,0) 

              ne(nn) = m
            end if
  350     continue
  400 continue
      end

c
      function padd(val)
c----------------------------------------------------------------------
c     for pfuncs: calculate increase of value      
c----------------------------------------------------------------------
c
      implicit  none
      real*8    padd, val, xval
      data      xval /0.0d0/
c     Look at parameter
      if(val.eq.0.0d0) then
        xval = 0.0d0
      else
        xval = xval + val
      endif
      padd = xval
      end
c
c
      function psub(val)
c----------------------------------------------------------------------
c     for pfuncs: calculate decrease of value      
c----------------------------------------------------------------------
c
      implicit  none
      real*8    psub, val, xval
      data      xval / 0.0d0 /
c     Look at parameter
      if(val.eq.0.0d0) then
        xval = 0.0d0
      else
        xval = xval - val
      endif
      psub = xval
      end
c
c
      subroutine pangl(ix,nen,angl,angg,nrot)
c----------------------------------------------------------------------
c      Purpose: Set value of angle for each element node

c      Inputs:
c         ix(nen) - Element nodal connection list
c         nen     - Number of nodes connected to element
c         angg(*) - Nodal angle array

c      Outputs:
c         angl(*) - Element angle array
c         nrot    - Number of nodes needing modification
c
c----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      dimension ix(nen), angl(nen),angg(*)
      nrot = 0
      do 100 n = 1,nen
         angl(n) = 0.0
         ii = ix(n)
         if (ii.gt.0) then
            if (angg(ii).ne.0.0d0) then
               angl(n) = angg(ii)
               nrot = nrot + 1
            end if
         end if
100   continue
      return
      end
c
      subroutine panglee(x,ang,ndm,numnp,prt,iang1,iang2)
c----------------------------------------------------------------------
c      Purpose: Set value of angle at edges

c      Inputs:
c         x(ndm,*)    - Nodal coordinates of mesh
c         ndm         - Spatial dimension of mesh
c         numnp       - Number of nodes in mesh
c         prt         - print flag
c         iang1       - first  axis for angl
c         iang2       - second axis for angl
c 
c      Outputs:
c         ang(*)  - Nodal angle array
c
c----------------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      dimension x(ndm,*),ang(*),td(3)
      data blank/-999.d0/
c
c.... read input of angles at edges 
100   if(ior.lt.0) write(*,3001)
      call dinput(td,3)
      i  = td(1)
      if(i.le.0.or.i.gt.ndm) go to 4
      x0  = td(2)
      val = td(3)
      dx = pdiff(x(i,1),ndm,numnp)/1000.
      do 200 n = 1,numnp
         if(x(i,n).ne.blank.and.abs(x(i,n)-x0).le.dx) then ! original
cww      if(                    abs(x(i,n)-x0).le.dx) then ! skip nodes(-999) via ebou,1,-999,1
            ang(n) = val
         end if
200   continue
      go to 100
4     if(prt)              write(iow,2000) iang1,iang2
      if(prt.and.ior.lt.0) write(*  ,2000) iang1,iang2
      do n = 1,numnp
        if(ang(n).ne.0.d0) then     
          if(prt)              write(iow,2001) n,ang(n)
          if(prt.and.ior.lt.0) write(  *,2001) n,ang(n)
        end if 
      end do        
      return
2000  format(/'  e d g e    n o d a l    a n g l e s ',
     +        '(set for directions ',i1,',',i1,')',/,
     +     4x,'node',5x,'angle')
2001  format(i8,3x,g12.5)
3001  format(' Input: ndir,x(ndir),angle','   >',$)
      end
c
      subroutine panglev(x,ang,ndm,numnp,prt,iang1,iang2,inpc)
c----------------------------------------------------------------------
c      Purpose: Set value of angle for block region
c      x_1a < x_1 < x_1e, x_2a < x_2 < x_2e, x_3a < x_3 < x_3e         
c
c      Inputs:
c         x(ndm,*)- Nodal coordinates of mesh
c         ndm     - Spatial dimension of mesh
c         numnp   - Number of nodes in mesh
c         prt     - print flag
c         iang1   - first  axis for angl
c         iang2   - second axis for angl
c         inpc    - defines type of input
c         x_i in  - cartesian coordinate system: inpc = 1                   
c                   xL1=phi   , xL2=phi+90                                  
c         x_i in  - cartesian coordinate system: inpc = 2                   
c                   xL1=phi-90, xL2=phi                                     
c         x_i in  - polar     coordinate system: inpc = 3                   
c                   xL1=phi   , xL2=phi+90                                  
c         x_i in  - polar     coordinate system: inpc = 4                   
c                   xL1=phi-90, xL2=phi                                     
c 
c      Outputs:
c         ang(*)  - Nodal angle array
c
c----------------------------------------------------------------------
c
      USE bdata
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      dimension x(ndm,*),ang(*),td(6)
      logical prt
      character ctyp(4)*9,cphi(4)*6
      data blank/-999.d0/
      data ctyp/'cartesian','cartesian','polar','polar'/
      data cphi/'phi   ','phi-90','phi   ','phi-90'/
      th = 45.0d0/datan(1.0d0)
c
c.... read input for boundary block 
100   if(ior.lt.0) write(*,3001)
      call dinput(td,2*ndm)
      if(errck) go to 100
      sum = 0.0d0
      do i = 1,2*ndm
         sum=sum+td(i)*td(i)
      end do
      if(sum.eq.0.0d0) goto 110 ! end
      xa  = td(1)
      xe  = td(2)
      ya  = td(3)
      ye  = td(4)
      za  = 0.d0
      ze  = 0.d0
      if(ndm.eq.3) then
        za  = td(5)
        ze  = td(6)
      end if
      if(prt) then  ! control input
                  write(iow,2002)xa,xe,ya,ye,za,ze,ctyp(inpc),cphi(inpc) 
        if(ior.lt.0)write(*,2002)xa,xe,ya,ye,za,ze,ctyp(inpc),cphi(inpc)
      end if  
c
      if(ndm.eq.2) then
        do  n = 1,numnp
          if(x(1,n).ne.blank .and. x(2,n).ne.blank .and.
     +       x(1,n).ge.xa    .and. x(1,n).le.xe    .and.
     +       x(2,n).ge.ya    .and. x(2,n).le.ye) then
c....       set angle
            x1 = x(iang1,n)
            x2 = x(iang2,n)
            if(inpc.eq.1.or.inpc.eq.2) then
c....         cartesian input
              if(x1.eq.0d0 .and. x2.eq.0.d0) then 
                ang(n) = 0.d0
              else
                ang(n) = datan2(x2,x1)*th
                if(inpc.eq.2) ang(n)=ang(n)-90.d0
              end if
            else if(inpc.eq.3.or.inpc.eq.4) then
c....         polar input
              ang(n) = x2
              if(inpc.eq.4) ang(n)=ang(n)-90.d0
            end if
          end if
        end do
      else if(ndm.eq.3) then
        do  n = 1,numnp
          if(x(1,n).ne.blank .and. x(2,n).ne.blank .and. x(3,n).ne.blank
     +      .and. x(1,n).ge.xa .and. x(1,n).le.xe .and.
     +            x(2,n).ge.ya .and. x(2,n).le.ye .and.
     +            x(3,n).ge.za .and. x(3,n).le.ze) then
c....       set angle
            x1 = x(iang1,n)
            x2 = x(iang2,n)
            if(inpc.eq.1.or.inpc.eq.2) then
c....         cartesian input
              if(x1.eq.0d0 .and. x2.eq.0.d0) then 
                ang(n) = 0.d0
              else
                ang(n) = datan2(x2,x1)*th
                if(inpc.eq.2) ang(n)=ang(n)-90.d0
              end if
            else if(inpc.eq.3.or.inpc.eq.4) then
c....         polar input
              ang(n) = x2
              if(inpc.eq.4) ang(n)=ang(n)-90.d0
            end if
          end if
        end do
      end if
      goto 100
c...  print
110   if(prt)              write(iow,2000) iang1,iang2
      if(prt.and.ior.lt.0) write(*  ,2000) iang1,iang2
      do n = 1,numnp
        if(ang(n).ne.0.d0) then     
          if(prt)              write(iow,2001) n,ang(n)
          if(prt.and.ior.lt.0) write(  *,2001) n,ang(n)
        end if 
      end do        
      return
2000  format(/'  b l o c k     n o d a l    a n g l e s ',
     +        '(set for directions ',i1,',',i1,')',/,
     +     4x,'node',5x,'angle')
2001  format(i8,3x,g12.5)
2002  format(/'  ANGL Generation with VANG:',/
     +  '  x_1a ',e12.5,'  x_1e ',e12.5,/
     +  '  x_2a ',e12.5,'  x_2e ',e12.5,/
     +  '  x_3a ',e12.5,'  x_3e ',e12.5,/,
     +     a9,' input, calculate ',a6)     
     
3001  format(' Input: X1a,X1e,X2a,X2e,X3a,X3e',' >',$)
      end
c
      subroutine panglb(b,angl,numnp,ndf)
c----------------------------------------------------------------------
c      Purpose: transform the parts defined by angl of displacement 
c               vector b to global directions, see pdefm b is modified!!
c
c      Inputs:
c         b(ndf,*) - displacement vector
c         angl(*)  - coordinates
c         numnp    - Number of nodes in mesh
c         ndf      - Number dof/node
c
c      Ouputs:
c         b(ndf,*) - displacement vector
c
c----------------------------------------------------------------------
      USE mdat2
      implicit double precision (a-h,o-z)
      dimension b(ndf,*),angl(*)
      pi=datan(1.d0)*4.d0
      if(ndf.eq.1) return  ! added ww
c...  dofs for angl 
      ij1 = ia(1)
      ij2 = ia(2)
c.... for ndf .ge.6 (transform for rotations)
      ij3 = ij1 + 3
      ij4 = ij2 + 3
c
      do 120 n = 1,numnp
        if(angl(n).ne.0.0d0) then
          ang = angl(n)*pi/180.d0
          cn  = cos(ang)
          sn  = sin(ang)
        else
          cn  = 1.0
          sn  = 0.0
        end if
c...    dof 1-3
          tm       = cn*b(ij1,n) - sn*b(ij2,n)
          b(ij2,n) = sn*b(ij1,n) + cn*b(ij2,n)
          b(ij1,n) = tm
c...    dof 4-6
        if(ndf.ge.6.and.itrot.eq.0) then
          tm       = cn*b(ij3,n) - sn*b(ij4,n)
          b(ij4,n) = sn*b(ij3,n) + cn*b(ij4,n)
          b(ij3,n) = tm
        end if
120   continue
      return
      end
c
      subroutine panglb1(xl,bl,angl,ndm,ndf,ipola)
c----------------------------------------------------------------------
c      Purpose: transformation of a nodal displacement vector b
c               to polar coordinates 
c
c      Inputs:
c         x(ndm,*) - Nodal coordinates of mesh
c         bl(ndf)  - Nodal displacement vector
c         angl(*)  - Nodal vector of angles
c         ndf      - Number dof/node
c         ndm      - Spatial dimension of mesh
c         ipola    - polar acts to axis ipola=12,13,23 
c
c      Ouputs:
c         b(ndf,*) - Nodal displacement vector
c
c----------------------------------------------------------------------
      USE mdat2
      implicit double precision (a-h,o-z)
      dimension xl(ndm),bl(ndf)
      if(ndf.eq.1) return  ! added ww
c
c.... back transformation due to angl 
      if(angl.ne.0.0d0) then
c...    dofs for angl 
        ij1 = ia(1)
        ij2 = ia(2)
c....   for ndf .ge.6 (transform for rotations)
        ij3 = ij1 + 3
        ij4 = ij2 + 3
        pi  = datan(1.d0)*4.d0
        ang = angl*pi/180.d0
        cn  = dcos(ang)
        sn  = dsin(ang)
      else
        cn  = 1.0d0
        sn  = 0.0d0
      end if
c...  dof 1-3
        tm      = cn*bl(ij1) - sn*bl(ij2)
        bl(ij2) = sn*bl(ij1) + cn*bl(ij2)
        bl(ij1) = tm
c...  dof 4-6
      if(ndf.ge.6.and.itrot.eq.0) then
        tm      = cn*bl(ij3) - sn*bl(ij4)
        bl(ij4) = sn*bl(ij3) + cn*bl(ij4)
        bl(ij3) = tm
      end if
c
c.... forward transformation due to pola 
c
      if(ipola.eq.12) then
        i=1
        k=2
      else if(ipola.eq.13) then
        i=1
        k=3
      else if(ipola.eq.23) then
        i=2
        k=3
      end if
      xk = xl(k)
      xi = xl(i)
      if(xi.eq.0d0 .and. xk.eq.0.d0) then ! special case sphere at top
        sn = 0.d0
        cn = 0.d0
      else
        phi = datan2(xk,xi)
        sn = dsin(phi)
        cn = dcos(phi)
      end if
      ur = cn*bl(i) + sn*bl(k)
      ut =-sn*bl(i) + cn*bl(k)
      bl(i) = ur
      bl(k) = ut
      if(ndf.ge.6) then ! always for rotations
        ur = cn*bl(i+3) + sn*bl(k+3)
        ut =-sn*bl(i+3) + cn*bl(k+3)
        bl(i+3) = ur
        bl(k+3) = ut
      end if 
      return
      end
c
      subroutine parexp(x,xs,v,nex,error)
c----------------------------------------------------------------------
c      Purpose: Identify parenthetical expressions and evaluate
c
c      Inputs:
c         x(*)     - String containing expression to evaluate

c      Scratch:
c         xs(*)    - Array used to temporarily store expression
c         v(*)     - Array to hold values

c      Outputs:
c         x(*)     - Expression replaced by upper case letter
c         nex      - Number of upper case letters used
c         error    - Flag, true if error occurs
c      Common returns:
c         www(*)   - Upper case letters with values assigned
c----------------------------------------------------------------------

      USE conval
      implicit double precision (a-h,o-z)
      logical error
      character*1 x(*),xs(*)
      dimension v(*)
c.... find parenthetical expressions and remove
      do 130 i = 1,75
        if(x(i).eq.'(') then
          i1 = i + 1
          do 120 j = i1,75
            if(x(j).eq.'(') then
              call errclr('PAREXP')
cww           stop
            return
            else if(x(j).eq.')') then
              do 50 l = 1,j-i+1
                    xs(l) = ' '
50            continue
              i2 = j - 1
              if(i2.lt.i1) then
                call errclr('PAREXP')
cww                 stop
              return
              else
                k = 0
                    do 100 l = i1,i2
                      k = k + 1
                      xs(k) = x(l)
                      x(l)      = ' '
100                 continue
                    x(i2+1)  = ' '
c.... evaluate the expression in the parenthesis
                    call evalex(xs,v,val,k,error)
                    if(error) return
                    nex = nex + 1
                    www(nex) = val
c.... put an upper case letter in the expression and close up remainder
                    x(i) = char(nex +64)
                    i2 = i2 -i1 + 2
                    do 110 l = i1,75
                      x(l) = ' '
                      if(l+i2.le.75) then
                        x(l) = x(l+i2)
                      end if
110                 continue
              end if
              go to 125
            end if
120       continue
125       continue
        end if
130   continue
      return
      end
c

      subroutine parse(wd1,wd2,v)
c----------------------------------------------------------------------
c.... input the values of the macro commands
c----------------------------------------------------------------------
      USE iofile
      USE yydata
      double precision v(3)
      character*4 wd1,wd2,x*75,x1*75
      CHARACTER*1 XX(75),YY(75)
      logical error
      EQUIVALENCE (X,XX)

c.... parse to format a macro command - check for errors
 101  error = .false.
      if(ior.gt.0) read(ior,'(a)',err=401,end=402) x1
      if(ior.lt.0) read(  *,'(a)',err=401,end=402) x1

C..   REMOVE LEADING BLANKS
      DO 150 I=1,75
        IF( X1(I:I).NE.' ') THEN
         II=1
         DO 160 J=I,75
           X(II:II) = X1(J:J)
160        II=II+1
         GOTO 170
        end if
150   CONTINUE
170   CONTINUE
c     Parse xx to find fields separated by commas: output yy   
      call acheck(xx,yy,15,75,75)
c
      DO 200 I=1,75
200   yyy(I:I) = YY(I)
c
c     write FEAP input macro data in asscoiated arrays
      call intio (yyy,wd1,wd2,v,error)
      if(error) go to 401
      return
 401  call  errclr ('PMACIO')
      write(*,*) 'input was: ',yyy 
      go to 101
 402  call  endclr ('PMACIO',x)
      return
      end
c
      subroutine pblock(x,id,idl,ndm,ndf,numnp,prt)
c----------------------------------------------------------------------
c      Purpose: Set boundary constraints for a  block region
c      x_1a < x_1 < x_1e, x_2a < x_2 < x_2e, x_3a < x_3 < x_3e         
c
c      x_i may be defined in any coordinate system
c
c      Inputs:
c         x(ndm,*)    - Nodal coordinates of mesh
c         id(ndf,*)   - Equation numbers for each active dof
c         ndm         - Spatial dimension of mesh
c         ndf         - Number dof/node
c         numnp       - Number of nodes in mesh
c         prt         - print flag
c
c      Scratch:
c         idl(ndf*ndf)- Local Equation numbers input
c 
c      Outputs:
c         id(ndf,*)   - Equation numbers for each active dof
c
c----------------------------------------------------------------------

      USE bdata
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      dimension x(ndm,numnp),id(ndf,numnp),idl(ndf*ndf),td(16)
      data blank/-999.d0/
c.... read input for boundary block - limit is 16 nos. / record
100   if(ior.lt.0) write(*,3001)
      il = min(ndf+2*ndm,16)
      call dinput(td,il)
      if(errck) go to 100
      sum = 0.0d0
      do i = 1,2*ndm
         sum=sum+td(i)*td(i)
      end do
      if(sum.eq.0.0d0) goto 110
      xa  = td(1)
      xe  = td(2)
      ya  = td(3)
      ye  = td(4)
      za  = 0.d0
      ze  = 0.d0
      if(ndm.eq.3) then
         za  = td(5)
         ze  = td(6)
      end if
      do 101 j = 1,min(ndf,10)
         idl(j) = td(j+2*ndm)
101   continue
      if(prt) then  ! control input
         write(iow,2002) xa,xe,ya,ye,za,ze, (idl(j),j=1,ndf)
         if(ior.lt.0) then 
         write(  *,2002) xa,xe,ya,ye,za,ze, (idl(j),j=1,ndf)
         end if
       end if  
      if(ndf.gt.10) then ! control input must be added
         do 105 ii = 1,(ndf+2*ndm)/16
          is = il+1
          il = min(is+15,ndf+2*ndm)
103       call dinput(td,il-is+1)
          if(errck) go to 103
          do 104 k = 1,il-is+1
           idl(k+is-1-2*ndm) = td(k)
104       continue
105      continue
      end if
      if(ndm.eq.2) then
        do  n = 1,numnp
          if(x(1,n).ne.blank.and.x(2,n).ne.blank.and.
     3      x(1,n).ge.xa .and. x(1,n).le.xe .and.
     4      x(2,n).ge.ya .and. x(2,n).le.ye) then
           do j = 1,ndf
            id(j,n) = max(abs(id(j,n)),abs(idl(j)))
           end do
          end if
        end do
      else if(ndm.eq.3) then
        do  n = 1,numnp
          if(x(1,n).ne.blank.and.x(2,n).ne.blank.and.x(3,n).ne.blank
     3    .and. x(1,n).ge.xa .and. x(1,n).le.xe .and.
     4          x(2,n).ge.ya .and. x(2,n).le.ye .and.
     5          x(3,n).ge.za .and. x(3,n).le.ze) then
           do j = 1,ndf
            id(j,n) = max(abs(id(j,n)),abs(idl(j)))
           end do
          end if
        end do
      end if
      goto 100
cww110   if(prt)              write(iow,2000) o,head,(i,i=1,ndf)
cww      if(prt.and.ior.lt.0) write(*  ,2000) o,head,(i,i=1,ndf)
110   if(prt)              write(iow,2000) (i,i=1,ndf)
      if(prt.and.ior.lt.0) write(*  ,2000) (i,i=1,ndf)
      do 410 n = 1,numnp
         do 408 i = 1,ndf
408      if(id(i,n).ne.0) go to 409
         go to 410
409      if(prt)              write(iow,2001) n,(id(i,n),i=1,ndf)
         if(prt.and.ior.lt.0) write(*  ,2001) n,(id(i,n),i=1,ndf)
410   continue
      return
cww2000  format(a1,19a4,a3//'  b l o c k   n o d a l    b . c .'/
2000  format(/'  b l o c k   n o d a l    b . c .'/
     1       /4x,'node',9(i3,'-b.c.')/(8x,9(i3,'-b.c.')))
2001  format(10i8/(8x,9i8))
2002  format(/'  BOUN Generation with VBOU:',/
     +  '  x_1a ',e12.5,'  x_1e ',e12.5,/
     +  '  x_2a ',e12.5,'  x_2e ',e12.5,/
     +  '  x_3a ',e12.5,'  x_3e ',e12.5,/,'  Code ',10(i1,1x))
3001  format(' Input: X1a,X1e,X2a,X2e,X3a,X3e,(idl(i),i=1,ndf)',' >',$)
      end
c
      subroutine pblockr(x,id,idl,ndm,ndf,numnp,prt)
c----------------------------------------------------------------------
c      Purpose: Set boundary constraints for a  block region in cyl. coordinates
c      r_a < r < r_e, phi_a < phi < phi_e, z_a < z < z_e
c      phi in degree, always >0
c
c      Inputs:
c         x(ndm,*)    - Nodal coordinates of mesh
c         id(ndf,*)   - Equation numbers for each active dof
c         ndm         - Spatial dimension of mesh
c         ndf         - Number dof/node
c         numnp       - Number of nodes in mesh
c         prt         - print flag
c
c      Scratch:
c         idl(ndf*ndf)- Local Equation numbers input
c 
c      Outputs:
c         id(ndf,*)   - Equation numbers for each active dof
c
c----------------------------------------------------------------------

      USE bdata
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      dimension x(ndm,numnp),id(ndf,numnp),idl(ndf*ndf),td(16)
      data blank/-999.d0/
      fac=45.d0/datan(1.0d0)
c.... read input for boundary block - limit is 16 nos. / record
100   if(ior.lt.0) write(*,3001)
      il = min(ndf+2*ndm,16)
      call dinput(td,il)
      if(errck) go to 100
      sum = 0.0d0
      do i = 1,2*ndm
         sum=sum+td(i)*td(i)
      end do
      if(sum.eq.0.0d0) goto 110
      ra  = td(1)
      re  = td(2)
      pa  = td(3)
      pe  = td(4)
      za  = 0.d0
      ze  = 0.d0
      if(ndm.eq.3) then
        za  = td(5)
        ze  = td(6)
      end if
      do 101 j = 1,min(ndf,10)
        idl(j) = td(j+2*ndm)
101   continue
      if(ndf.gt.10) then
        do 105 ii = 1,(ndf+2*ndm)/16
        is = il+1
        il = min(is+15,ndf+2*ndm)
103     call dinput(td,il-is+1)
        if(errck) go to 103
        do 104 k = 1,il-is+1
          idl(k+is-1-2*ndm) = td(k)
104     continue
105     continue
      end if
      if(prt) then  ! control input
          write(iow,2002) ra,re,pa,pe,za,ze, (idl(j),j=1,ndf)
        if(ior.lt.0) then 
          write(  *,2002) ra,re,pa,pe,za,ze, (idl(j),j=1,ndf)
        end if
      end if  
cww   if(ndm.eq.3) then
      do 200 n = 1,numnp
        r = dsqrt(x(1,n)*x(1,n)+x(2,n)*x(2,n))
        if(x(1,n).eq.0.0d0.and.x(2,n).eq.0.0d0) then
          p = 0.0d0
        else
          p = fac*datan2(x(2,n),x(1,n)) ! angles = +-...  
          if(p.lt.0.d0) p = p + 360.d0
        end if
        if(ndm.eq.2) then 
          if(x(1,n).ne.blank.and.x(2,n).ne.blank.and.
     1      r.ge.ra .and. r.le.re .and. p.ge.pa .and. p.le.pe ) then
            do j = 1,ndf
               id(j,n) = max(abs(id(j,n)),abs(idl(j)))
            end do
          end if
        else if(ndm.eq.3) then   
          if(x(1,n).ne.blank.and.x(2,n).ne.blank.and.x(3,n).ne.blank
     1      .and. r.ge.ra .and.      r.le.re .and.
     2            p.ge.pa .and.      p.le.pe .and.
     3       x(3,n).ge.za .and. x(3,n).le.ze) then
            do j = 1,ndf
               id(j,n) = max(abs(id(j,n)),abs(idl(j)))
            end do
          end if
        end if
200   continue
      goto 100
cww110   if(prt)                   write(iow,2000) o,head,(i,i=1,ndf)
cww      if(prt.and.ior.lt.0) write(*  ,2000) o,head,(i,i=1,ndf)
110   if(prt)              write(iow,2000) (i,i=1,ndf)
      if(prt.and.ior.lt.0) write(*  ,2000) (i,i=1,ndf)
      do 410 n = 1,numnp
        do 408 i = 1,ndf
408       if(id(i,n).ne.0) go to 409
        go to 410
409     if(prt) write(iow,2001) n,(id(i,n),i=1,ndf)
        if(prt.and.ior.lt.0) write(*,2001) n,(id(i,n),i=1,ndf)
410   continue
      return
cww2000  format(a1,19a4,a3//'radial block  n o d a l    b . c .'/
2000  format(/'radial block  n o d a l    b . c .'/
     1       /4x,'node',9(i3,'-b.c.')/(8x,9(i3,'-b.c.')))
2001  format(10i8/(8x,9i8))
2002  format(/'  BOUN Generation with RBOU:',/
     +  '  r_a  ',e12.5,'  r_e  ',e12.5,/
     +  '  p_a  ',e12.5,'  p_e  ',e12.5,/
     +  '  z_a  ',e12.5,'  z_e  ',e12.5,/,'  Code ',10(i1,1x))
3001  format(' Input: Ra,Re,Pa,Pe,Za,Ze,(idl(i),i=1,ndf)',' >',$)
      end
c
c----------------------------------------------------------------------
c
      subroutine pcfour(u,du,aa,ab,ct,ndf,numfs,numnp,isw)
c----------------------------------------------------------------------
c
c.... purpose:  fourier series solution of linear problems
c
c     nf - fourier harmonic
c     displacement
c     u-th    :  u-nf           ; v-nf         ; w-nf
c     multipliers
c     nf > 0  :  cos nf-theta ; sin nf-theta ; cos nf-theta
c     nf = 0  :       1       ;      1       ;      1
c     nf < 0  :  sin nf-theta ; cos nf-theta ; sin nf-theta
c
c 
c     isw  - 1  fourier solution
c          - 2  fourier solution step - displacements
c          - 3  fourier solution step - stresses
c          - 4  sum displacements 
c          - 5  sum stresses
c
c     Comment: Macro needs to be tested 
c
c----------------------------------------------------------------------
      USE fodata
      USE iodata
      USE iofile
      implicit double precision (a-h,o-z)
      logical lfis
      character ffile(2)*8
      integer fleng(2)
      real*8  ct(3),fact
      real*8  aa(8,numfs),ab(8,numfs)
      real*8  u(ndf,numnp),du(ndf,numnp)
      save factor
      data factor /1.0d0/, pi4 /0.0d0/, nfterm /0/
      if(isw.gt.1) then
        ii = 1 + mod(isw,2)
        inquire(file=ffile(ii),exist=lfis)
        if(lfis) then
          open(ios, file=ffile(ii), access='direct', recl=fleng(ii),
     1         status='old')
        else
          open(ios, file=ffile(ii), access='direct', recl=fleng(ii),
     1         status='new')
        end if
      end if
c.... fourier solution
      if(isw.eq.1) then
        foout  = .false.
        nf     = ct(1)
        factor = ct(2)
        if(factor.eq.0.0d0) factor = 1.0
        write(iow,2000) nf,factor
        if(ior.lt.0) write(*,2000) nf,factor
        if(.not.foflg) then
          foflg  = .true.
          nfterm = 0
          pi4    = atan(1.d0)/45.d0
c.... set file names and lengths
          ffile(1) = 'four_dis'
          fleng(1) = ndf*numnp*8 + 16
          ffile(2) = 'four_str'
          fleng(2) =   8*numfs*nfs*8 + 16
        end if
c.... fourier solution step - displacements
      else if(isw.eq.2. and. foflg) then
        nfterm = nfterm + 1
        foout  = .false.
        fact   = factor
        if(ct(1).ne.0.0d0) fact = ct(1)
        write(iow,2001) nfterm,nf,fact
        if(ior.lt.0) write(*,2001) nfterm,nf,fact
        write(ios,rec=nfterm) nf,nf1,fact,u
        close(ios)
c.... fourier solution step - stresses
      else if(isw.eq.3. and. foflg) then
        write(ios,rec=nfterm) nf,nf1,fact,aa
        close(ios)
c.... sum displacements
      else if(isw.eq.4) then
        if(lfis) then
          theta = ct(1)*pi4
          write(iow,2002) nfterm,ct(1),theta
          if(ior.lt.0) write(*,2002) nfterm,ct(1),theta
          ii = ndf*numnp
          call pzero (u,ii)
          do 110 i = 1,nfterm
            read(ios,rec=i,err=400) nf,nf1,fact,du
            if(ior.lt.0) write(*,2003) nf,fact
            if(nf.gt.0) then
              t1 = fact*cos(nf*theta)
              t2 = fact*sin(nf*theta)
            else if(nf.eq.0) then
              t1 = fact
              t2 = fact
            else if(nf.lt.0) then
              t1 = fact*sin(-nf*theta)
              t2 = fact*cos( nf*theta)
            end if
            do 100 n = 1,numnp
              u(1,n) = u(1,n) + t1*du(1,n)
              u(2,n) = u(2,n) + t2*du(2,n)
              if(ndf.ge.3) u(3,n) = u(3,n) + t1*du(3,n)
100         continue
110       continue
        end if
        close(ios)
c
c.... sum stresses
      else if(isw.eq.5) then
        if(lfis) then
          theta = ct(1)*pi4
          ii = 8*numfs
          call pzero (ab,ii)
          do 210 i = 1,nfterm
            read(ios,rec=i,err=400) nf,nf1,fact,aa
            if(nf.gt.0) then
              t1 = fact*cos(nf*theta)
              t2 = fact*sin(nf*theta)
            else if(nf.eq.0) then
              t1 = fact
              t2 = fact
            else if(nf.lt.0) then
              t1 = fact*sin(-nf*theta)
              t2 = fact*cos( nf*theta)
            end if
            do 200 n = 1,nfs
              ab(1,n) = ab(1,n) + t1*aa(1,n)
              ab(2,n) = ab(2,n) + t1*aa(2,n)
              ab(3,n) = ab(3,n) + t1*aa(3,n)
              ab(4,n) = ab(4,n) + t1*aa(4,n)
              ab(5,n) = ab(5,n) + t2*aa(5,n)
              ab(6,n) = ab(6,n) + t2*aa(6,n)
              ab(7,n) = aa(7,n)
              ab(8,n) = aa(8,n)
200         continue
210       continue
          foout  = .true.
        end if
        close(ios)
      end if
      return
400   call errclr('PCFOUR')
      return
2000  format(/'   Fourier harmonic = ',i3,'  Factor =',e12.5/1x)
2001  format(/'   Fourier solution:',i3,' terms.'/
     1        '     Current harmonic = ',i3,'  Factor =',e12.5/1x)
2002  format(/'   Fourier sum for ',i3,' terms at ',f7.2,' degrees;'/
     1        '                             or ',f7.2,' radians.'/1x)
2003  format('   Fourier summation: Harmonic',i3,', factor =',e12.5)
      end
c
c

      subroutine pcharr(y,ivd,n0,nt)
c-----------------------------------------------------------------
c      Purpose: Identify 'character' variables in the string
c               N.B. Lower case variables have been input,
c                    upper case computed

c      Inputs:
c         y(*)     - String to search
c         n0       - Field width
c         nt       - Number of characters in y-array

c      Outputs:
c         ivd(2,*) - Number of characters found
c-----------------------------------------------------------------
      character*1 y(*)
      integer   n0,nt,n1,k,i,j,n,ivd(2,*),kk

      n1 = n0 - 1
      k = 0
      do i = 1,nt,n0
        k = k + 1
        ivd(1,k) = 0
        ivd(2,k) = 0
        n = ichar( y(i) ) - 64
        if(n.gt.0) then
          if(n.gt.58) go to 200
          if(n.gt.26) then
            ivd(1,k) = n - 32
            j = ichar(y(i+1))
            if(j.eq.32) then
              ivd(2,k) = 0
            elseif(j.ge.ichar('a') .and. j.le.ichar('z')) then
              ivd(2,k) = j - ichar('a') + 1
            elseif(j.ge.ichar('0') .and. j.le.ichar('9')) then
              ivd(2,k) = j - ichar('0') + 27
            endif
          else
            ivd(1,k) = - n
          end if


          do j = i,i+n1
            y(j) = ' '
c           1 arbitrary statement necessary for SALFORD without DEBUG Option!!?? 
            kk = 42           
          end do  

c....  set value
          y(i+n1) = '0'

        end if

      end do   


      return
c.... error
 200  call errclr('PCHARR')
      write(*,*) 'error character no: ',n,'field ',i, 'value',y(i) 
 
      return
      end
c
      subroutine pcheck(nc,xs,error)
c-----------------------------------------------------------------
c      Purpose: Check that input string contains admissible data
c               and parentheses match.  Convert all input letters
c               to lower case for further processing

c      Inputs:
c         nc     - Number of characters to check
c         xs(*)  - Character array

c      Outputs:
c         error  - Flag, true if error occurs
c-----------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical error
      character*1 x(75),xs(75)
c.... make sure that all of x is in lower case and blanks are removed
      i = 0
      do 100 j = 1,75
        x(j) = ' '
        if(xs(j).ne.' ' .and. xs(j).ne.'=' .and. xs(j).ne.',') then
          i = i + 1
          x(i)  = xs(j)
          xs(j) = ' '
          n = ichar( x(i) )
          if(n.ge.65 .and. n.le.90) x(i) = char(n + 32)
        end if
100   continue
c.... move back and check the characters for incorrect parenthesis
      error = .false.
      n = 0
      do 200 j = 1,i
        xs(j) = x(j)
        if(xs(j).eq.'(') n = n+1
        if(xs(j).eq.')') n = n-1
        if(n.lt.0 .or. n.gt.1 ) error = .true.
200   continue
      if(n.ne.0) error = .true.
      n = ichar(xs(1))
      if(n.lt.97 .or. n.gt.122) error = .true.
c.... check the characters for incorrect parameters
      if(.not.error) then
        do 210 j = 2,i
          n = ichar(xs(j))
          if(.not.(n.ge.97 .and. n.le.122) .and.
     1       .not.(n.ge.40 .and. n.le.57) ) then
            error = .true.
          end if
210     continue
      end if
      if(error) then
        write(*,2000)
      else
        write(*,2001) nc,(xs(j),j=1,i)
      end if
      return
 2000 format(' Incorrect statement - reinput ')
 2001 format('   No.',i3,'>',a1,' = ',74a1)
      end
c
      subroutine pcktie(id,dt,st,ndf,numnp)
c-----------------------------------------------------------------
c      Purpose: search for tied nodes to average stress projections
c
c      Inputs:
c         id(ndf,*)   - Equation numbers for each active dof
c         dt(*)       - Nodal Sum of weights 
c         st(numnp,*) - Nodal stress values
c         ndf         - Number dof/node
c         numnp       - Number of nodes in mesh
c
c      Outputs:
c         dt(*)       - Nodal Sum of weights 
c         st(numnp,*) - Nodal stress values
c
c      WW no longer used 09.04.2009: no projection necessary,
c                                    tied nodes are not used  
c-----------------------------------------------------------------
      USE strnam
      implicit double precision (a-h,o-z)
      dimension id(ndf,*),dt(*),st(numnp,*)
      ii = iabs(istv)
      do 140 n = 1,numnp-1
         nc = n
         do 110 m = n+1,numnp
            do 50 i = 1,ndf
               if(id(i,nc).ne.id(i,m)) go to 110
50          continue
c.... accumulate values

            write(*,*) 'accumulate at node - testww', nc

            dt(m) = dt(m) + dt(nc)
            do 100 i = 1,ii
100         st(m,i) = st(m,i) + st(nc,i)
c.... reset nc to node 'm'
            nc = m
110      continue
c.... reset other values if ties exist
         if (nc.gt.n) then
            do 130 m = 1,nc-1
               if(id(1,m).eq.id(1,nc)) then
                  dt(m) = dt(nc)
                  do 120 i = 1,ii
120               st(m,i) = st(nc,i)
               end if
130         continue
         end if
140   continue
      return
      end
c

      logical function pcomp(a,b,n)
c-----------------------------------------------------------------
c      Purpose: Compare character strings for match
c               Ignores upper/lower case differences.
c
c      Inputs:
c         a(*)   - Character string 1
c         b(*)   - Character string 2
c         n      - Number of characters to compare
c
c      Outputs:
c         pcomp  - Flag, true if a = b
c-----------------------------------------------------------------

      character*(*) a,b
c      character*4 a,b
      pcomp = .false.
c.... compute the increment between an upper and a lower case letter
      inc = ichar('A') - ichar('a')
c.... compare for a match
      do 100 i = 1,n
         ia = ichar(a(i:i))
         ib = ichar(b(i:i))
c.... test all permutations of characters for a match
         if(ia.ne.ib .and. ia+inc.ne.ib .and. ia.ne.ib+inc ) return
100   continue
      pcomp = .true.
      return
      end
C
c-----------------------------------------------------------------
c
      subroutine pconst(prt)
c-----------------------------------------------------------------
c
c     Purpose: Input parameter expressions:  let = expression
c              arithmetic free input routine
c              start from PMESH 
c
c      Inputs:
c         prt    - Print input values if true
c
c      Outputs:
c        Values of parameters a-z are stored in array vvv(i,j)
c        vvv(i, 1-26) Character a-z 
c        vvv(i,27-36) Character 0-9 
c        vvv(i,    0) Character ' ' 
c        vvv(1-26,k)  Character a-z
c        character position:
c        0=048 ... 9=057   ... i-47    =0
c        A=065 ... Z=090   ... i-64    =A 
c        a=097 ... z=122   ... i-64-32 =a
        
c        modified WW 11/05 to two characters
c-----------------------------------------------------------------
      USE conval
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      character*1 let*2,eql,x(75),xx*4
      logical pcomp,prt
cww   integer n,i,j,ial,izl,iau,izu,id,iq, i0,i9

c.... character positions
      i0  = ichar('0') ! 48 
      i9  = ichar('9') ! 57 
      iq  = ichar('=') ! 61
      ial = ichar('a') ! 97
      izl = ichar('z') !122
      iau = ichar('A') ! 65
      izu = ichar('Z') ! 90
      id  = ial - iau  ! 32
c
      if (prt) then
                     write(iow,2001)
        if(ior.lt.0) write(*  ,2001)
      end if
c.... input a record from file ior or the keyboard *
1     if(ior.gt.0) then
          read (ior,1000,err=901,end=902) let,eql,x
      else
          write(  *,3000)
          read (  *,1000,err=901,end=902) let,eql,x
      end if

c.... check length of string
      num=ipos(x,75) 
      if(num.eq.75) then 
         write(*,*) 'Warning: string too long'
         write(*,*) x      
         stop   
      end if

c     check for macro with one character 
      if(let(2:2).eq.'=') then
        do i = 75,2,-1
          x(i) = x((i-1))
        end do
        x(1)     = eql   
        eql      = '='
        let(2:2) = ' '
      end if

              
c...  check for a blank character or the null character = blank line
10    if(let(1:1).eq.' '.or.ichar(let(1:1)).eq.0) then
          x(1) = ' '
          let  = ' '
          return
      end if
c...  compare 'xx' for a match to 'list' = list values to screen
      xx(1:1) = let(1:1)
      xx(2:2) = let(2:2)
      xx(3:3) = eql
      xx(4:4) = x(1)
      if(pcomp(xx,'list',4)) then
        if(ior.lt.0) then
          write(*,2001)
          do i = 1,26
            if(vvv(i,0).ne.0.0d0) then  ! character,-
              write(*,2000) char(i+96),' ',vvv(i,0)
            end if
            do j = 1,26                 ! character,character
              if(vvv(i,j).ne.0.0d0) then
                write(*,2000) char(i+96),char(j+96),vvv(i,j)
              end if
            end do
            do j = 27,36
              if(vvv(i,j).ne.0.0d0) then ! character,digit
                write(*,2000) char(i+96),char(j+i0-27),vvv(i,j)
              end if
            end do
          end do
        end if
        go to 1
      end if
c.... Locate correct location for the addition in lower case letters
      i = ichar(let(1:1))
      if(ial.le.i .and. i.le.izl) i = i-96
      if(iau.le.i .and. i.le.izu) i = i-64
      if(i.gt.26) go to 901

      if(let(2:2).eq.' ') then
        j        = 0
      else
        j = ichar(let(2:2))
        if(ial.le.j .and. j.le.izl) then
          j = j-96
          if(j.gt.26) go to 901
        end if
        if(iau.le.j .and. j.le.izu) then
          j = j-64
          if(j.gt.26) go to 901
        end if
        if( i0.le.j .and. j.le.i9 ) then 
          j = j-21  !-47+26
          if(j.lt.27.or.j.gt.36) go to 901
        end if
      end if

      errck = .false.
      call setval(x,75, val)


c.... store value at calculated position
      vvv(i,j) = val

      if (prt) then
                     write(iow,2000) let(1:1),let(2:2),vvv(i,j)
        if(ior.lt.0) write(  *,2000) let(1:1),let(2:2),vvv(i,j)
      end if
      go to 1
c.... error on a read
901   call  errclr ('PCONST')
      if (ior.lt.0)  goto 1
      return
c.... eof encountered
902   call  endclr ('PCONST',let)
      goto  10
c.... formats
 1000 format(a2,76a1)
 2000 format(5x,'Parameter ',a1,a1,' = ',e15.8)
 2001 format(/'  p a r a m e t e r   v a l u e s')
 3000 format(' Use "list" to give current values - <CR> to exit'/
     1       ' Input: letter=expression (no blanks)'/'  -->',$)
      end
c
      subroutine pconstm(lct,prt)
c-----------------------------------------------------------------
c
c     Purpose: Input parameter expressions:  let = expression
c              arithmetic free input routine
c              start from PMACR 
c
c               Set parameter 'par' to 'value'.

c      Inputs:
c        lct(1:2) - Parameter name to set
c        lct(2/3) - = 
c        lct (3/4 - value 
c        prt      - Echo value set if true
c
c      Outputs:
c        Values of parameters a-z are stored in array vvv(i,j)
c        vvv(i, 1-26) Character a-z 
c        vvv(i,27-36) Character 0-9 
c        vvv(i,    0) Character ' ' 
c        vvv(1-26,k)  Character a-z
c        character position:
c        0=048 ... 9=057   ... i-47    =0
c        A=065 ... Z=090   ... i-64    =A 
c        a=097 ... z=122   ... i-64-32 =a
c        
c       Comments:
c       similar to pconst, only INPUT is different
c        modified WW 11/05 to two characters
c-----------------------------------------------------------------
      USE conval
      USE iofile
      implicit double precision (a-h,o-z)

      character*1  lct(4),eql,x(75)
      character*2  let 
      character*4  xx
      logical pcomp,prt

c.... character positions
      i0  = ichar('0') ! 48 
      i9  = ichar('9') ! 57 
      iq  = ichar('=') ! 61
      ial = ichar('a') ! 97
      izl = ichar('z') !122
      iau = ichar('A') ! 65
      izu = ichar('Z') ! 90
      id  = ial - iau  ! 32

c...  compare 'xx' for a match to 'list' = list values to screen
      xx(1:1) = lct(1)
      xx(2:2) = lct(2)
      xx(3:3) = lct(3)
      xx(4:4) = lct(4)
      if(pcomp(xx,'list',4)) then
        write(*,2001)
        do i = 1,26
          if(vvv(i,0).ne.0.0d0) then  ! character,-
            write(*,2000) char(i+96),' ',vvv(i,0)
          end if
          do j = 1,26                 ! character,character
            if(vvv(i,j).ne.0.0d0) then
              write(*,2000) char(i+96),char(j+96),vvv(i,j)
            end if
          end do
          do j = 27,36
            if(vvv(i,j).ne.0.0d0) then ! character,digit
              write(*,2000) char(i+96),char(j+i0-27),vvv(i,j)
            end if
          end do
        end do
  
      else
c...    set parameter
c       Put data into 'let','eql','x'
c       macro with two characters 
        let = ' '
        eql = '='
        x   = ' '

        let(1:1)=lct(1)

        if(lct(2).ne.'=') then 
          let(2:2) = lct(2)
          x(1)     = lct(4)
        else
          x(1) = lct(3) 
          x(2) = lct(4)
        end if

c....   Locate correct location for the addition in lower case letters
        i = ichar(let(1:1))
        if(ial.le.i .and. i.le.izl) i = i-96
        if(iau.le.i .and. i.le.izu) i = i-64
        if(i.gt.26) stop 'i.gt.26 pconstm'
        
        if(let(2:2).eq.' ') then
          j        = 0
        else
          j = ichar(let(2:2))
          if(ial.le.j .and. j.le.izl) then
            j = j-96
            if(j.gt.26) stop 'j.gt.26 pconstm'
          end if
          if(iau.le.j .and. j.le.izu) then
            j = j-64
            if(j.gt.26) stop 'j.gt.26 pconstm'
          end if
          if( i0.le.j .and. j.le.i9 ) then 
            j = j-21  !-47+26
            if(j.lt.27.or.j.gt.36) stop 'j.gt.26 pconstm'
          end if
        end if
        
        call setval(x,75, val)
        
c....   store value at calculated position
        vvv(i,j) = val
       
c....   print value 
        if (prt) then
                       write(iow,2001)
                       write(iow,2000) let(1:1),let(2:2),vvv(i,j)
          if(ior.lt.0) write(*  ,2001)
          if(ior.lt.0) write(  *,2000) let(1:1),let(2:2),vvv(i,j)
        end if

      end if

      return
c.... formats
 2000 format(5x,'Parameter ',a1,a1,' = ',e15.8)
 2001 format(/'  p a r a m e t e r   v a l u e s')
      end
c

      subroutine pcontr
c-----------------------------------------------------------------------
c      Purpose: Control program for FEAP problem input and solution.

c      Inputs:
c        none

c      Outputs:
c        none
c
c-----------------------------------------------------------------------
      USE bdata
      USE bisw
      USE cdat1
      USE cdata
      USE dirdat
      USE edgdat
      USE errchk
      USE hdatam
      USE iofile
      USE iosave
      USE maclg
      USE mdata
      USE mdat2
      USE mxasz
      USE ndata
      USE nolink
      USE pdata6
      USE pdata7
      USE pdata8
      USE plodfb
      USE plong
      USE pnodn
      USE printh
      USE psize
      USE qload
      USE rfeap
      USE sdata
      USE slid4
      USE soltyp
      USE strnam
      USE aunr
      USE vdata
      USE working
      USE doalloc
      implicit double precision (a-h,o-z)
      logical tfl,pcomp,errs, ldummy
      character*4 titl(20),titll,yyy*80
      dimension td(7), tl(9)
      real*4  tary,tary1
      integer, allocatable, dimension(:):: opth1,opth2,opth3,opth4,opth5
c
c.... set files back to default values
1     ior = iabs(ior)
c.... read a card and compare first 4 columns with macro list
      read(ior,1000,err=601,end=700) (titl(i),i=1,20)
      if(pcomp(titl(1),'feap',4)) go to 100
      if(pcomp(titl(1),'edge',4)) go to 150
      if(pcomp(titl(1),'inte',4)) go to 200
      if(pcomp(titl(1),'macr',4)) go to 300
c     if(pcomp(titl(1),'batc',4)) go to 300
cww----------------------------------------            
      if(pcomp(titl(1),'batc',4)) then
c...    for batc,TANG and batc,UPDH 
c          execute if FEAP is started in BATCH-mode      
c          dummy if not   
c...    for batc,MICR
c          execute if FEAP is started interactive      
c          dummy if not   
c...    for batc,-  normal

        call acheck(titl,yyy,15,80,80)
        read(yyy,1002,err=600,end=700) titl(1),titl(2)

        if(pcomp(titl(2),'tang',4)) then  
          if(iplot.eq.0) then
            if(iswb.eq.3.or.iswb.eq.4.or.iswb.eq.6.or.iswb.eq.8)then  
              go to 300 ! execute
            else
              goto 10   ! dummy
            end if ! iswb
          else if(iplot.eq.1) then
            goto 10     ! dummy
          end if ! iplot
        else if(pcomp(titl(2),'updh',4)) then  
          if(iplot.eq.0) then
            if(iswb.eq.15) then 
              go to 300 ! execute
            else
              goto 10   ! dummy
            end if
          else if(iplot.eq.1) then
            goto 10     ! dummy 
          end if ! iplot
        else if(pcomp(titl(2),'micr',4)) then   
          if(iplot.eq.0) then
            goto 10     ! dummy
          else if(iplot.eq.1) then
            go to 300   ! execute
          end if ! iplot
        else ! standard case: execute always 
          go to 300
        end if ! titl(2)
c....   read batch macros as dummys
10      read(ior,1000) (titl(i),i=1,20)
        if(pcomp(titl(1),'stop',4)) go to 1
        go to 10
      end if ! batc
cww----------------------------------------            
      if(pcomp(titl(1),'macr',4)) go to 300
      if(pcomp(titl(1),'link',4)) go to 400
      if(pcomp(titl(1),'stop',4)) return
      if(pcomp(titl(1),'tie ',3)) go to 500
      if(pcomp(titl(1),'opti',4)) go to 550
      if(pcomp(titl(1),'optn',4)) go to 550

      go to 1

c.... read and print control information
c.....[feap]
100   do 101 i = 1,20
101   head(i) = titl(i)
c
c.... set initial parameters
      call ini
c
      call dinput(tl,9)
      numnp  =tl(1)
      numel  =tl(2)
      nummat =tl(3)
      ndm    =tl(4)
      ndf    =tl(5)
      nen    =tl(6)
      nad    =tl(7)
      ndd    =tl(8)
      nthread=tl(9)
      nde = 0
      ned = 0
      ndd = max(ndd,50)
c      if(nthread.eq.0) nthread=omp_get_max_threads ()
c      call omp_set_num_threads(nthread)

c     If number of nodes, or elements is zero compute number from data
      if(numnp.eq.0 .or. numel.eq.0  .or. nummat.eq.0) then
cww     call drawmess('calculate number of nodes/elements/materials',1,4)
        call pnums(numnp,numel,ndm,nen,nummat,prt)
cww     call drawmess(' ',1,5)
      end if 
      write(iow,2000) head,versn,
     1                numnp,numel,nummat,ndm,ndf,nen,nad,ndd

c.... set parameters for page eject, and edge/rotation dof
      o   = '    '
      errck = .false.
      flnk  = .false.
      lsave = .false.
c...  edgdat
      ne1   = 1
      ne2   = 1
      ne3   = 1
      ne4   = 1
      ne5   = 1

      ia(1) = 1
      ia(2) = 2
      itrot = 0

c.... plot data for 100 elements with max 40 nodes to plot 
      do i = 1,100
        inord(i) = 0
        nedg(i)  = 0
        do j = 1,12
          iedg(j,i)  = 0
        end do  
        do j = 1,40
          ipord(j,i) = 0
        end do  
      end do  
c.... set pointers for allocation of data arrays
      nen1 = nen + 4
      nie  = ndf + 3
      nst  = nen*ndf + ned*nde + nad
cww   kmax = noff + (numnp*max(ndf,ndm)+ned*nde*numel)*ipr
      kmax=noff
      nmaxdr=(numnp*max(ndf,ndm)+ned*nde*numel)

c     allocation of main arrays   
      call ralloc(   dr,nmaxdr,'DR-ndr',tfl)
      call ralloc( edis,      6*nst,    'Elmt-DISP-nn',tfl)
      call ralloc( ecor,    nen*ndm,    'Elmt-COOR-n0',tfl)
      call ralloc( etem,        nen,    'Elmt-TEMP-n1',tfl)
      call ialloc( eeqn,        nst,    'Elmt-EQNO-n2',tfl)
      call ralloc( epve,        nst,    'Elmt-Pvec-n3',tfl)
      call ralloc( ekma,    nst*nst,    'Elmt-Kmat-n4',tfl)
      call ialloc( nmat, nummat*nie,     'MAT-List-n5',tfl)
      call ralloc( edma, nummat*ndd,         'DMAT-n6',tfl)
      call ialloc( psid,  ndf*numnp,     'ID-Array-n7',tfl)
      call ralloc( coor,  ndm*numnp,         'COOR-n8',tfl)
      call ialloc( econ, nen1*numel,         'ELEM-n9',tfl)
      call ralloc( gloa,  numnp*ndf,        'LOAD-n10',tfl)
      call ralloc( gtem,      numnp,        'TEMP-n11',tfl)
      call ralloc( aang,        nen, 'Elmt-ANGLE-n11a',tfl)
      call ralloc( bang,      numnp, 'Elmt-ANGLE-n11b',tfl)
      call ralloc( glo0,  numnp*ndf,     'LOAD F0-n13',tfl)
      call ralloc(   gu,3*numnp*ndf,        'DISP-n14',tfl)
      call ialloc( gtie,      numnp,   'TIE-Nodes-nip',tfl)
      call ralloc(  plo,    10*nplo,       'TPLO-mplo',tfl)
      call ialloc(aipma,     nummat,           'IPMA',tfl)

c.... set initial values ipma
      call pconsi(aipma,nummat)
c.... set initial values gtie
      call pconsi(gtie,numnp)
            
c.... dummy allocation of needed arrays for check pointer alloc intel
      call ialloc( edge1,1, 'Edge1-dummy',tfl)
      call ialloc( edge2,1, 'Edge2-dummy',tfl)
      call ialloc( edge3,1, 'Edge3-dummy',tfl)
      call ialloc( edge4,1, 'Edge4-dummy',tfl)
      call ralloc( trans,1, 'Trans-dummy',tfl)
      call ralloc( strea,1, 'Strea-dummy',tfl)
      call ralloc( basea,1, 'Base-dummy' ,tfl)
      call ralloc( dampm,1, 'Mass-dummy' ,tfl)
      call ralloc(   eh1,1,   'EH1-dummy',tfl)
      call ralloc(   eh2,1,   'EH2-dummy',tfl)
      call ralloc(   eh3,1,   'EH3-dummy',tfl)
      call ralloc(gstiff,1,'gstiff-dummy',tfl)
      call ralloc( aqloa,1,        'QLOA',tfl)

      nal=1
      nau=1

c.... call mesh input subroutine to read and print all mesh data
      iii = 0
      call pmesh1(numnp,numel,numat,ndm,ndf,nen)
      call pmesh(eeqn,nmat,edma,psid,coor,econ,gloa,gtem,glo0,
     1     ndd,nie,ndf,ndm,nen1,iii,prt)

c.... perform simple check on mesh to ensure basic data was input
      call meshck(coor,econ,nen,nen1,ndm,numnp,numel,nummat,errs)
      if(errs) 
     + call drawmess(' Error in Mesh data, modify input file!!',1,-2)
      call pmesh2(n1,n2,n3,n4,n5,n6,n7,n8,n9,nmat,numnp,numel,nummat,
     +        nie)
      tfl = .true.
      go to 1

c.... determine edge profile for compact storage
c.....[edge]
150   if(nde*ned.gt.0) then
        call ialloc(edge1,numnp,'Edge1',tfl)
        ne2 = numnp + mod(numnp,ipr)
        ne3 = ned*numel
        call ialloc(edge2,ne2,'Edge2',tfl)
        call ialloc(edge3,ne3,'Edge3',tfl)
        
        call esiden(econ,nmat,edge1,edge2,edge3,
     1              nie,nen1,numel,numnp)
        if(prt) call eprint(edge1,edge2,numnp)
        
        iii = edge1(numnp-1)
        ne4 = iii*nde+1
        ! EDGE-ne4 + EDGE-ne5
        call ialloc(edge2,iii,'Edge2',tfl)
        call ialloc(edge3,nde*iii,'Edge3',tfl)
        call ialloc(edge4,6*nde*iii,'Edge4',tfl)

c.... set edge boundary conditions and loads
        call esideb(coor,edge1,edge2,edge3,edge4,
     1              ndm,nde,numnp,prt)
        tfl = .true.
      end if
      go to 1

c.... set files for interactive macro execution (ior negative)
c.....[inte]
200   close(ior)  ! close input file to have access by other programs
      ior = -iabs(ior)

c.... solution of problem using macros 
c.....[macr]
300   if(tfl) then

c....   determine the maximum storage that can occur with no b.c. fixed
        neq = ndf*numnp + nde*ned*numel
        call ialloc(jdt12,neq+1,'JD-TEST-n12',tfl) ! for JD (eq.pointer)

c....   determine current profile
        call perform(0,1,1,30)
        call etimef(tary)
c....   set up equation numbers of dof's
        call profil(jdt12,eeqn,psid,econ,edge1,edge2,edge3,1,prt)
        call etimef(tary1)
c....   change profile with respect to LINK
        if(flnk) call plink(eeqn,psid,coor,ndm,ndf,numnp,neq,
     +           link1,link2,prt)
c....   compute the column/row lengths and profile
        call profil(jdt12,eeqn,psid,econ,edge1,edge2,edge3,2,prt)
        call etimef(tary1)
        if(prt) write(iow,2018) tary1-tary  ! time for Profile
        if(prt) write(*  ,2018) tary1-tary
2018    format('     Time for Profile',37x,'t=',0pf9.4)
        call perform(1,1,3,30)

c....   store up to ncmds macro commands (batch mode or interactive loops)
        ncmds=200
        call ralloc(macl,3*ncmds,'MACRO-List-nct',tfl)

c....   set up stress history addresses
        call sethis(nmat,econ,nie,nen,nen1,numel,nummat,prt)
      end if

c.... zero the solution vectors
cww   call pzero(m(n13),  ndf*numnp) ! F0
      gu = 0.d0        !DISP

c.... macro module for establishing solution algorithm
      call pmacr(edis,ecor,etem,eeqn,epve,ekma,nmat,edma,psid,
     1     coor,econ,gloa,gtem,jdt12,glo0,gu,dr,macl,
     2     ncmds,ndf,ndm,nen1,nst,prt,basea,plo)

      if(irfeap.gt.0) return ! RINP RINP,NEW REME,NEW REME,OLD
      go to 1

c.... reset id list to link dof's on different nodes
c.....[link]
400   call perform(0,1,1,31)
      call ialloc(link1,numnp,'LINK-nli1',tfl)      ! array to plot linked nodes
      call ialloc(link2,numnp*ndf,'LINK-nli2',tfl)  ! array to plot linked nodes
      call ialloc(link3,nen*ndf,'LINK-nli3',tfl)    ! local array for periodic B.C.s
      call plinka ! input of link data and save to file   
      flnk = .true.
      call perform(1,1,3,31)
      tfl = .true.
      go to 1

c.... tie nodes within tolerance of one another, no parameter input possible!!
c.... [tie,,<tol>]
c.... [tie,all,<tol>]
c.... [tie,node,<tol>,n2,n3] 
c.... [tie,line,<tol>]  cont. x1,y1,<z1>,x2,y2,<z2> 
c.... [tie,mlin,<tol>],n2,n3  cont. x1,y1,<z1>,x2,y2,<z2> 
c.... [tie,dir,<tol>,j,xj] 
c.... [tie,mate,<tol>,n2,n3] 
c.... [tie,prin] 
c.... [tie,nopr] 

500   continue
      if(pcomp(titl(2),'prin',4)) then
        prt=.true.
        goto 1
      end if   
      if(pcomp(titl(2),'nopr',4)) then
        prt=.false.
        goto 1
      end if   
c.... make copy of original array for use in paraview
      call ialloc(tecon,nen1*numel,'ELEM with TIE',tfl)
      tecon=econ
      ixtie = 1 ! = tie exist for paraview 
      call pzero(td,7)
      call perform(0,1,1,28)
      call acheck(titl,yyy,15,80,80)
      read(yyy,1002,err=600,end=700) titl(1),titl(2),(td(i),i=1,3)

      if(pcomp(titl(2),'node',4)) then
c       tie node between numbers l1 and l2
        j  = 1
        l1 = max(    1,int(td(2)))
        l2 = min(numnp,int(td(3)))
        l3 = 0
        l4 = 0
        write(iow,2011) l1,l2
2011    format(/5x,'Tie nodes from',i8,' to ',i8/1x)

      else if(pcomp(titl(2),'line',4)) then
c       tie node along line between x1,y1,<z1> and x2,y2,<z2>
        j  = 2
        l1 = 1
        l2 = numnp
        l3 = 0
        l4 = 0
        call dinput(td(2),6)
        if(ndm.eq.2) then
          write(iow,2012) (td(i),i=2,5)
        else if(ndm.eq.3) then
          write(iow,2016) (td(i),i=2,7)
      end if    
2012    format(/5x,'Tie nodes on line from node 1 to node 2',/,
     +  3x,'  x1=',e15.5,'  y1=',e15.5,/,
     +  3x,'  x2=',e15.5,'  y2=',e15.5)
2016    format(/5x,'Tie nodes on line from node 1 to node 2',/,
     +  3x,'  x1=',e15.5,'  y1=',e15.5,'  z1=',e15.5,/,
     +  3x,'  x2=',e15.5,'  y2=',e15.5,'  z2=',e15.5)

      else if(pcomp(titl(2),'mlin',4)) then
c       tie node along line between x1,y1,<z1> and x2,y2,<z2> only for same material n2=mat
        j  = 3
        l1 = 1
        l2 = numnp
        l3 = int(td(2))
        l4 = int(td(3))
        call dinput(td(2),6)
        write(iow,2014) l3,l4
        if(ndm.eq.2) then
          write(iow,2012) (td(i),i=2,5)
        else if(ndm.eq.3) then
          write(iow,2016) (td(i),i=2,7)
      end if    

      else if(pcomp(titl(2),'dir',3)) then
c       tie all nodes in dir td(2) with coordinate td(3)
        j  = 4
        l1 = 1
        l2 = numnp
        l3 = 0
        l4 = 0
        write(iow,2013) td(2),td(3)
2013    format(/5x,'Tie: direction =',f3.0,' X =',1p,1e12.5/1x)

      else if(pcomp(titl(2),'mate',4)) then
c       tie nodes for mate n3 and n4   
        j  = 5
        l1 = 1
        l2 = numnp
        l3 = int(td(2))
        l4 = int(td(3))
        write(iow,2014) l3,l4
2014    format(/5x,'Tie from material',i4,' to material',i4/1x)

      else
c       tie all nodes
        j  = 0
        l1 = 1
        l2 = numnp
        l3 = 0
        l4 = 0
        write(iow,2015)
2015    format(/5x,'Tie all nodes with common coordinates'/1x)
      endif

      call tienod(psid,econ,coor,gloa,gtie,
     1         ndm,ndf,nen,nen1,numnp,numel,ntied,l1,l2,l3,l4,j,td,prt)
      tfl = .true.
      call perform(1,1,3,28)
      go to 1

c.... optimize profile
c.....[opti,optn]                     optn must be checked !!!!!
550   continue
      if(istyp.gt.0) goto 1   ! do not use for any other solver
      call perform(0,1,1,29)
      call ialloc(optin,numnp,'OPTI-nren',tfl)
      nren=1
c.... adresses for scratch files
      allocate(opth1(numnp))
      allocate(opth2(numnp))
      allocate(opth3(numnp))
      allocate(opth4(numel))
      allocate(opth5(numnp*numnp))
      nbn = numnp
c.... front optimization
      if(pcomp(titl(1),'opti',4))
     1 call opnum(econ,opth1,opth2,opth5,opth3,opth4,optin,
     2            numnp,numel,nen,nen1,nbn,.false.)
c.... node numbering optimization (Cuthill McKee)
      if(pcomp(titl(1),'optn',4))
     1  call optn(econ,optin,opth1,opth2,opth3,opth5,numnp,numel,
     2            nen,nen1,nren,nbn,.true.)
      tfl = .true.
      deallocate(opth1)
      deallocate(opth2)
      deallocate(opth3)
      deallocate(opth4)
      deallocate(opth5)
      call perform(1,1,3,29)
      go to 1
600   call  errclr ('PCONTR')
      return
c...  provisorial ww
601   return
c
700   call  endclr ('PCONTR',titl(1))
      return
c.... input/output formats
1000  format(20a4)
1001  format(10i8)
1002  format(2(a4,11x),3f15.0)
2000  format(1h ,19a4,a3//14x,a16/14x,a16/14x,a16
     1               //5x,'Number of nodal points        =',i6
     2                /5x,'Number of elements            =',i6
     3                /5x,'Number of material sets       =',i6
     4                /5x,'Dimension of coordinate space =',i6
     5                /5x,'Degrees-of-freedom/node       =',i6
     6                /5x,'Nodes per element (maximum)   =',i6
     7                /5x,'Added degrees-of-freedom/elmt =',i6
     8                /5x,'Maximum matl. properties/elmt =',i6)
cww     9                     /5x,'Degrees-of-freedom/edge       =',i6/
cww     x                     /5x,'Edges per element (maximum)   =',i6)
c
cww     7                     /5x,'Degrees-of-freedom/edge      =',i6
cww     8                     /5x,'Edges per element (maximum)  =',i6/
cww     9                     /5x,'Maximum matl. properties/elmt=',i6/
cww     x                     /5x,'Added degrees-of-freedom/elmt=',i6)
c.kne2001  format(10x,'Maximum Array Storage =',i8/)
      end
c
      double precision function pdiff(x,ndm,numnp)
c-----------------------------------------------------------------
c      Purpose: Compute the difference between maximum and minimum
c               nodal coordinates.
c
c      Inputs:
c         x(ndm,* ) - Nodal coordinates for mesh
c         ndm       - Spatial dimension of mesh
c         numnp     - Number of nodes in mesh
c
c      Outputs:
c         pdiff     - Difference between maximum and minimum
c-----------------------------------------------------------------
      USE iofile
      implicit double precision (a-h,o-z)
      dimension x(ndm,numnp)
      data blank /-999./
      do 100 n = 1,numnp
         if(x(1,n).ne.blank) go to 110
100   continue
      call drawmess(' error all coodinates are still blank',1,0)
      pdiff = 0.0d0
      return
110   xmx = x(1,n)
      xmn = x(1,n)
      do 200 n = 1,numnp
         if(x(1,n).ne.blank) then
            xmx = dmax1(xmx,x(1,n))
            xmn = dmin1(xmn,x(1,n))
         end if
200   continue
      pdiff = xmx - xmn
      return
      end
c
      subroutine pdirec(x,dr,dir,numnp,numel,ndm,prt,u)
c-----------------------------------------------------------------
c
c      Purpose: calculate nodal triads in field: T=[t_1,t_2,t_3]
c
c      Inputs:
c         x(ndm,* )       - Nodal coordinates for mesh
c         dir(10,knode,2) - triad field 1-9 = dir, 10 = no.of el.per node
c         dir( , ,1)      - save   
c         dir( , ,2)      - test, update at TIME
c         knode           - numnp(idtyp=>0,3-17)   or nen*numel(idtyp=1)
c         numnp           - Number of nodes in mesh
c         nen             - Max. Number of nodes per element
c         numel           - Number of elements in mesh
c         ndm             - Spatial dimension of mesh
c         prt             - Print flag
c         u

c      idtyp =  type of triad
c      --------------------------------
c      idtyp =  0  triad calculated on element-level(isw=18) for each
c                  global node, + averageing 
c      idtyp =  1  triad calculated on element-level(isw=18) for each
c                  element node,  no averageing on global level
c      idtyp =  2  triad calculated on element-level(isw=18) for each
c                  element node,  no averageing on global level
c                  e_3 from triad, input t_1, calculate e_1=t_1/|t_1|, e_2=e_3xe_1  
c      idtyp =  3-17 analytical solutions for global nodes                                                 
c      --------------------------------
c      * field dir defined in pmesh
c      * director calculated in pmacr isw=18 
c      * director plotted in plot 
c
c      Scratch:
c         dr(*)
c
c      Outputs:
c
c         dir(10,knode,2) - triad field 1-9 = dir, 10 = no.of el.per node
c
c.... ww bs uka 12/96+2/04+10/05
c-----------------------------------------------------------------
      USE bdata
      USE dirdat
      USE hdatam
      USE iofile
      USE pnodn
      implicit double precision (a-h,o-z)
      logical prt,fa
      character head1*72,yyy*80
      dimension dr(*),x(ndm,*),dir(10,knode,2),r(3,3),u(*),xl(3),add(3)
      fa=.false.
      call pzero(dir,knode*10*2)
      idtyp =xdir(1)
      if(idtyp.lt.0.or.idtyp.gt.20) then
        write(yyy,2003) idtyp
        call drawmess(yyy,1,0)
        return
      end if
      if(idtyp.eq.0) then
c....   triad calculated at nodes -avering all elements values
        head1='  Triad at nodes from averaging element values'
        hflgu  = .false.
        h3flgu = .false.
        call formfe(u,dr,dr,dr,fa,fa,fa,fa,18,1,numel,1)
c....   average director field - lamina base Hughes
        do i = 1,numnp
          node = dir(10,i,1)
          if(node.gt.0) then
            do k = 1,9
              dir(k,i,1) = dir(k,i,1) / node
            end do
            call plamina_base(dir(1,i,1))
          end if
        end do
      end if
      if(idtyp.eq.1) then
c....   triad calculated at all element nodes
        head1='  Triad at all element nodes, no averaging'
        hflgu  = .false.
        h3flgu = .false.
        call formfe(u,dr,dr,dr,fa,fa,fa,fa,18,1,numel,1)
c....   norm if necessary
        do i = 1, knode
          call vnorm (dir(1,i,1),dummy)
          call vnorm (dir(4,i,1),dummy)
          call vnorm (dir(7,i,1),dummy)
        end do
        goto 100
      end if 
      if(idtyp.eq.2) then
c....   triad calculated at all element nodes
        head1='  Triad at all element nodes, no averaging, given t_1'
        hflgu  = .false.
        h3flgu = .false.
        call formfe(u,dr,dr,dr,fa,fa,fa,fa,18,1,numel,1)
c....   input t_1
        add(1) = xdir(2)
        add(2) = xdir(3)
        add(3) = xdir(4)
        call vnorm (add,dummy)
c....   modify triad
        do i = 1, knode
c....     norm on e_3
          call vnorm (dir(7,i,1),dummy) 
c....     e_1
          dir(1,i,1) = add(1)  
          dir(2,i,1) = add(2) 
          dir(2,i,1) = add(3) 
c....     e_2
          call vecp (dir(3,i,1),dir(1,i,1),dir(2,i,1))
        end do
        goto 100
      end if 
      If(idtyp.gt.2.and.idtyp.lt.18) then
c....   analytical solutions
        if(idtyp.eq. 3) head1=' for plate (x-y plane)'
        if(idtyp.eq. 4) head1=' for plate (y-z plane)'
        if(idtyp.eq. 5) head1=' for plate (x-z plane)'
        if(idtyp.eq. 6) head1=' for cylindrical shell(x-axis)'
        if(idtyp.eq. 7) head1=' not defined' 
        if(idtyp.eq. 8) head1=' for cylindrical shell(y-axis)'
        if(idtyp.eq. 9) head1=' not defined' 
        if(idtyp.eq.10) head1=' for cylindrical shell(z-axis) outside'
        if(idtyp.eq.11) head1=' for cylindrical shell(z-axis)  inside'
        if(idtyp.eq.12) head1=' for sphere'
        if(idtyp.eq.13) head1=' for hypar'
        if(idtyp.eq.14) head1=' for rotations hyperboloid'
        if(idtyp.eq.15) head1=' for rotations paraboloid'
        if(idtyp.eq.16) head1=' for parabola'
        if(idtyp.eq.17) head1=' for twisted beam'

c
c.....  loop over nodes
        do i = 1,numnp
          call pzero(r,9)
          do k=1,3
            xl(k)=x(k,i)
          end do   
c....     triad
          call pdirec3(r,xdir(2),xl,idtyp)
c
c.....    store nodal values in field dir (1+2)  
          do k = 1,3
            dir(k  ,i,1) = r(k,1)
            dir(k+3,i,1) = r(k,2)
            dir(k+6,i,1) = r(k,3)

            dir(k  ,i,2) = r(k,1)
            dir(k+3,i,2) = r(k,2)
            dir(k+6,i,2) = r(k,3)

          end do
        end do
      end if
c>>wd   Director for NURBS shells
      if (idtyp.ge.18.and.idtyp.le.20) then
        if(idtyp.eq.18) head1=' NURBS control point net lin. int.'
        if(idtyp.eq.19) head1=' NURBS control point net lin. int2.'
        if(idtyp.eq.20) head1=' NURBS derivatives'
        if (idtyp.eq.18.or.idtyp.eq.19)
     +    call pdireciga1(idtyp,x,dir,numnp,ndm,1)
        if (idtyp.eq.20) call pdireciga2(idtyp,x,dir,numnp,ndm,1,gtie)
        do i = 1,numnp
          node = dir(10,i,1)
          if(node.gt.0) then
            do k = 1,9
              dir(k,i,1) = dir(k,i,1) / node
            end do
            call plamina_base(dir(1,i,1))
          end if
        end do
      end if
c<<wd
c.... print director
cww100   if(prt)                write(iow,2001) o,head,head1
cww      if(prt .and. ior.lt.0) write(*  ,2001) o,head,head1
100   if(prt)                write(iow,2001) head1
      if(prt .and. ior.lt.0) write(*  ,2001) head1
c     if(prt) write(iow,2000) o,head,head1,(i,i=1,9)
c             write(*  ,2000) o,head,head1,(i,i=1,9)
c     do n = 1,knode
c         if(prt) write(iow,2002) n,(dir(i,n,1),i=1,10)
c                     write(*  ,2002) n,(dir(i,n,1),i=1,10)
c     end do
      return
cww2001  format(a1,19a4,a3,// ' n o d a l  d i r e c t o r  ',a72,/)
2001  format(/ ' n o d a l  d i r e c t o r  ',a72,/)
c2000  format(a1,19a4,a3,//
c     1 ' n o d a l  d i r e c t o r  ',a72,//,
c     2  4x,'node',9(2x,i4,' -dir   '),(2x,'el. per node'))
c2002  format(i8,10(2x,e12.5))
2003  format(' Error: invalid type ',i3,' of director under base')
      end
c
      subroutine plamina_base(diri)
c-----------------------------------------------------------------------
c
c      Purpose: calculate 'best' nodal triad from Lamina Basis (Hughes)
c.....          see Book Belytschko Nonlinear FEM p.541
c
c      Inputs:
c         diri(i,3 ) - director field for node T=[t_1|t_2|t_3]
c         i          - component of base vector i=1..3
c
c      Outputs:
c         diri(i,3 ) - director field for node
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension diri(3,3),a(3),b(3),g1(3),g2(3),ex(3),ey(3),ez(3)      
c.... g_alpha
      do i = 1,3
        g1(i)=diri(i,1)  
        g2(i)=diri(i,2)  
      end do
c.... e_z      
      call vecp (g1,g2,ez)
      call vnorm(ez,ezb)     
c.... a
      call vnorm(g1,g1b)     
      call vnorm(g2,g2b)     
      do i = 1,3
        a(i) = g1(i) + g2(i)
      end do
c.... b
      call vecp (ez,a,b)
c.... e_x,e_y
      do i = 1,3
        ex(i) = a(i) - b(i)
        ey(i) = a(i) + b(i)
      end do
      call vnorm(ex,exb)
      call vnorm(ey,eyb)
c.... averaged triad       
      do i = 1,3  
        diri(i,1) = ex(i) 
        diri(i,2) = ey(i) 
        diri(i,3) = ez(i) 
      end do
      return
      end  
c
      subroutine pdirec1(ix,dir,diri,i,isw,ifdir)
c-----------------------------------------------------------------------
c
c      Purpose: read/save nodal triad for global(!) node i on actual 
c               element for averaging of global nodal values
c
c      Inputs:
c
c         ix(nen1,*)      - Element nodal connections of mesh
c         dir(10,knode,2) - triad field 1-9 = dir, 10 = no.of el.per node
c         dir( , ,ifdir)  - save   
c         dir( , ,ifdir)  - test, update at TIME
c         diri(3,3)       - triad for node i
c         i               - local  node number
c         isw             - 1=write 2=read
c         ifdir           - 1=save 2=test, update at TIME
c
c      Outputs:
c
c         dir(10,knode,2) - triad field 1-9 = dir, 10 = no.of el.per node
c
c.... ww bs uka 12/96+2/04+10/05
c-----------------------------------------------------------------
      USE dirdat
      implicit double precision (a-h,o-z)
      dimension dir(10,knode,2),diri(3,3),ix(*)
      if(isw.eq.1) then
c.....  write
        ii = abs(ix(i))
        if(ii.eq.0) return
        do k = 1,3
          dir(k  ,ii,ifdir) = dir(k  ,ii,ifdir) + diri(k,1)
          dir(k+3,ii,ifdir) = dir(k+3,ii,ifdir) + diri(k,2)
          dir(k+6,ii,ifdir) = dir(k+6,ii,ifdir) + diri(k,3)
        end do
        dir(10,ii,ifdir) = dir(10,ii,ifdir) + 1
      else if(isw.eq.2) then
c.....  read
        call pzero(diri,9)
        ii = abs(ix(i))
        if(ii.eq.0) return
        do k = 1,3
          diri(k,1) = dir(k  ,ii,ifdir)
          diri(k,2) = dir(k+3,ii,ifdir)
          diri(k,3) = dir(k+6,ii,ifdir)
        end do
      end if
      return
      end
c
      subroutine pdirec2(dir,diri,i,n,nen,isw,ifdir)
c-----------------------------------------------------------------------
c
c      Purpose: read/save nodal triad for element node i on actual 
c               element n
c
c      Inputs:
c
c         dir(10,knode,2) - triad field 1-9 = dir, 10 = no.of el.per node
c         dir( , ,ifdir)  - save   
c         dir( , ,ifdir)  - test, update at TIME
c         diri(3,3)       - triad for node i
c         i               - local  node number
c         n               - actual element number
c         isw             - 1=write 2=read
c         ifdir           - 1=save 2=test, update at TIME
c
c      Outputs:
c
c         dir(10,knode,2) - triad field 1-9 = dir, 10 = no.of el.per node
c
c.... ww bs uka 12/96+2/04+10/05
c-----------------------------------------------------------------
      USE dirdat
      implicit double precision (a-h,o-z)
      dimension dir(10,knode,2),diri(3,3)
      if(isw.eq.1) then
c.....  write
        ii = (n-1) * nen + i
        do k = 1,3
          dir(k  ,ii,ifdir) = diri(k,1)
          dir(k+3,ii,ifdir) = diri(k,2)
          dir(k+6,ii,ifdir) = diri(k,3)
        end do
      else if(isw.eq.2) then
c.....  read
        call pzero(diri,9)
        ii = (n-1) * nen + i
        do k = 1,3
          diri(k,1) = dir(k  ,ii,ifdir)
          diri(k,2) = dir(k+3,ii,ifdir)
          diri(k,3) = dir(k+6,ii,ifdir)
        end do
      end if
      return
      end
c
      subroutine pdirec3(r,a,x,idtyp)
c-----------------------------------------------------------------------
c
c      Purpose: calculate nodal triad for an analytical description 
c               of a surface, 
c               used from macro BASE and from elements
c
c      Inputs:
c
c         r(3,3)  - Nodal triad [t_1|t|2|t_3] (in columns) 
c         a(3)    - additional values to describe surface 
c         x(3)    - Coordinate of sampling point
c         idtyp   - Typ of triad
c
c         idtyp   -  3 plate x-y plane e_1(1,0,0),e_2(0,1,0),e_3(0,0,1)
c         idtyp   -  4 plate y-z plane e_1(0,1,0),e_2(0,0,1),e_3(1,0,0)
c         idtyp   -  5 plate x-z plane e_1(1,0,0),e_2(0,0,1),e_3(0,-1,0)
c         idtyp   -  6 plate x-y plane polar

c         idtyp   -  8 cylinder(x=axis)      y0,z0
c         idtyp   -  9 cylinder(y=axis)      x0,z0
c         idtyp   - 10 cylinder(z=axis)      x0,y0 outside
c         idtyp   - 11 cylinder(z=axis)      x0,y0  inside
c         idtyp   - 12 sphere                x0,y0,z0
c         idtyp   - 13 hypar                 lx,ly,lz
c         idtyp   - 14 rotations hyperboloid r0,c
c         idtyp   - 15 rotations paraboloid  r0
c         idtyp   - 16 parabola              a
c         idtyp   - 17 twisted beam          lx

c
c      Outputs:
c
c         r       - Nodal triad [t_1|t|2|t_3] (in columns) 
c
c.... ww bs uka 2/04+10/05
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension r(3,3),a(3),x(3) 

      call pzero(r,9)  
      if(idtyp.eq.3) then
c....   plate in x-y plane
c       t_1 = (1,0,0), t_2 = (0,1,0), t_3 = (0,0,1)
c     
        r(1,1) = 1.0d0

        r(2,2) = 1.0d0

        r(3,3) = 1.0d0
c
      else if(idtyp.eq.4) then
c....   plate in y-z plane
c       t_1 = (0,1,0), t_2 = (0,0,1), t_3 = (1,0,0)
c     
        r(2,1) = 1.0d0

        r(3,2) = 1.0d0

        r(1,3) = 1.0d0
c
      else if(idtyp.eq.5) then
c....   plate in x-z plane
c       t_1 = (1,0,0), t_2 = (0,0,1), t_3 = (0,-1,0)
c     
        r(1,1) = 1.0d0

        r(3,2) = 1.0d0

        r(2,3) =-1.0d0
c
      else if(idtyp.eq.6) then
c.....  plate x-y plane polar
c       t_1=(cos a, sin a,0), t_3=(-sin a,cos a,0), t_3=(0,0,1) 
c
        xi = x(1)
        yi = x(2)
        dphi = datan2(yi,xi)
        sn = dsin(dphi)
        cs = dcos(dphi)

        r(1,1) =  cs
        r(2,1) =  sn

        r(1,2) = -sn
        r(2,2) =  cs

        r(3,3) =  1.d0
c
      else if(idtyp.eq.8) then
c....   cylindrical shell (axis in x-direction)
c....   x=x, y=y_0+rcos a, z=z_0+rsin a
c       t_1=(1,0,0), t_2=(0,-sin a,cos a), t_3=(0,cos a,sin a) 
c
        yi = x(2) + a(1)
        zi = x(3) + a(2)
        dphi = datan2(zi,yi)
        sn = dsin(dphi)
        cs = dcos(dphi)

        r(1,1) =  1.d0

        r(2,2) =  sn
        r(3,2) = -cs

        r(2,3) =  cs
        r(3,3) =  sn
c 
      else if(idtyp.eq.9) then
c....   cylindrical shell (axis in y-direction)
c....   x=x_0+rcos a, y=y,z=z_0+rsin a
c       t_1=(sin a,0,-cos a), t_2=(0,1,0), t_3=(cos a,0,sin a) 
c
        xi = x(1) + a(1)
        zi = x(3) + a(2)
        dphi = datan2(zi,xi)
        sn = dsin(dphi)
        cs = dcos(dphi)

        r(1,1) =  sn
        r(3,1) = -cs

        r(2,2) =  1.d0

        r(1,3) =  cs
        r(3,3) =  sn
c 


c
      else if(idtyp.eq.10) then
c.....  cylindrical shell (axis in z-direction)  outside view
c....   x=x_0+rcos a, y=y,z=z_0+rsin a
c       t_1=(-sin a,cos a,0), t_2=(0,0,1), t_3=(cos a,sin a,0) 
c
        xi = x(1) + a(1)
        yi = x(2) + a(2)
        dphi = datan2(yi,xi)
        sn = dsin(dphi)
        cs = dcos(dphi)

        r(1,1) = -sn
        r(2,1) =  cs

        r(3,2) =  1.d0

        r(1,3) =  cs
        r(2,3) =  sn

c
      else if(idtyp.eq.11) then
c.....  cylindrical shell (axis in z-direction)   inside view
c....   x=x_0+rcos a, y=y,z=z_0+rsin a
c       t_1=(sin a,-cos a,0), t_2=(0,0,1), t_3=(-cos a,-sin a,0) 
c
        xi = x(1) + a(1)
        yi = x(2) + a(2)
        dphi = datan2(yi,xi)
        sn = dsin(dphi)
        cs = dcos(dphi)

        r(1,1) =  sn
        r(2,1) = -cs

        r(3,2) =  1.d0

        r(1,3) = -cs
        r(2,3) = -sn

c
      else if(idtyp.eq.12) then
c.....  spherical shell  x0x+x0,y=y+y0, vertical = z=z+z0
c       sphere (z-axis): 1=in theta d. , 2=in phi-d.(circum.), 3=rad. 
c       t_1=[-snp,      csp,         0], 
c       t_2=[-csth*csp,-csth*snp, snth], 
c       t_3=[ snth*csp, snth*snp, csth]
c
        xi = x(1) + a(1)
        yi = x(2) + a(2)
        zi = x(3) + a(3)
        rs = dsqrt(xi*xi+yi*yi+zi*zi)
        if(yi.eq.0.0.and.xi.eq.0.0) then
          phi = 0.0
        else
          phi  = datan2(yi,xi)
        end if
        snp  = dsin(phi)
        csp  = dcos(phi)
        csth = x(3)/rs
        th   = dacos(csth)
        snth = dsin(th)
c     
        r(1,1) =-snp
        r(2,1) = csp
        r(3,1) = 0.d0

        r(1,2) =-csth*csp
        r(2,2) =-csth*snp
        r(3,2) = snth
      
        r(1,3) = snth*csp
        r(2,3) = snth*snp
        r(3,3) = csth
c     
      else if(idtyp.eq.13) then
c....   hypar shell: plane in x,y KOS at 0,0 = midpoint of hypar, z=axy
c       high: (-1,-1)  (l,l)  low:  ( l,-1)  (-1,l)      
c       x=x,y=y,z=4*lz*(x/lx)*(y/ly)=a*x*y
c
c       t_1=[1,0,a*y]/t1 t_2=[0,1,a*x]/t2 
c       t_3=[-ay,-ax,1)/t3
c       t1=(1+a^2x^2)^1/2
c       t2=(1+a^2y^2)^1/2
c       t3=(1+a^2(x^2+y^2))^1/2
c----------------------------------------------------------------------
c       these base vectors are not othogonal!! 
c         a) solution t_2 = t_3 x t_1
c         b) solution t_1 = t_2 x t_3
c       >>c) lamina base Belytschko p. 541
c----------------------------------------------------------------------
c     
        ax = a(1)
        ay = a(2)
        az = a(3) 
        aq = 4.d0*az/ax/ay
        a2 = aq*aq
        xi = x(1)       
        x2 = xi*xi         
        yi = x(2)       
        y2 = yi*yi         
c       length of vectors
        t1 = 1.d0/sqrt(1.d0 + a2*y2)
        t2 = 1.d0/sqrt(1.d0 + a2*x2)
        t3 = 1.d0/sqrt(1.d0 + a2*(x2+y2))

        r(1,1) =  1.d0*t1 
        r(3,1) =  aq*yi*t1

        r(2,2) =  1.d0*t2 
        r(3,2) =  aq*xi*t2

        r(1,3) = -aq*yi*t3
        r(2,3) = -aq*xi*t3
        r(3,3) =  1.d0*t3

c....   a) orthogonal vector t_2     
cww     call vecp (r(1,3),r(1,1),r(1,2))
c....   b) orthogonal vector t_1     
cww     call vecp (r(1,2),r(1,3),r(1,1))
c....   c) lamina base     
        call plamina_base(r)
c     
      else if(idtyp.eq.14.or.idtyp.eq.15) then
c....   11  rotations hyperboloid r = r0/c * sqrt(c**2 + z**2)
c....   12  rotations paraboloid  r = r0 * sqrt(2z)
c....   t_1=[-sn;cs;0], t_2=[r,z*cs;r,z*sn;1]/|t|, t_3=[cs;sn;-r,z]/|t|
c....                   |t|=(1+r,z**2)^0.5
        xi = x(1)
        yi = x(2)
        zi = x(3)
        dphi = datan2(yi,xi)
        sn   = dsin(dphi)
        cs   = dcos(dphi)
        r0   = a(1)
        if(idtyp.eq.14) then
          cc   = a(2)
          rkz  = r0/cc*zi/dsqrt(cc*cc+zi*zi)
        else if(idtyp.eq.15) then
          rkz  = r0/dsqrt(2.d0*zi)
        end if
        rn   = dsqrt(1.0d0+rkz*rkz)
c
        r(1,1) = -sn
        r(2,1) =  cs

        r(1,2) =  rkz*cs/rn
        r(2,2) =  rkz*sn/rn
        r(3,2) =  1.d0/rn

        r(1,3) =  cs/rn
        r(2,3) =  sn*rn
        r(3,3) = -rkz/rn
c     
      else if(idtyp.eq.16) then
c....   parabola  z = a*x**2
c....   t_1=[1;0;2ax]/|t|, t_2=[0;1;0], t_3=[-2ax;0;1]/|t|
c....   |t|=(1+(2ax)**2)^0.5
        aq = a(1)
        dz = 2.d0*aq*x(1)
        rn = dsqrt(1.d0+dz*dz)
        sn = dz/rn
        cs = 1.d0/rn
      
        r(1,1) =  cs
        r(3,1) =  sn

        r(2,2) =  1.d0

        r(1,3) = -sn
        r(3,3) =  cs
c
      else if(idtyp.eq.17) then
c....   twisted beam 0 < (lx=a(1) <  lx                           
c....   angle phi = pi/(2*lx)*x
c....   f(x,r) = [x, r cos(phi), r sin(phi)]
c....   f,x    = [1, -r*sin(phi)*pi/(2*lx), r*cos(phi)pi/(2*lx)]->t_1 + norm             
c....        r = sign(y)*(y*y+z*z)^0.5
c....   f,r    = [0, cos(phi), sin(phi)]                        ->t_2 
c....   t_3 = t_1 x t_2             
c
        pi   = 4.d0*datan(1.d0)
        pia  = pi/2.d0/a(1)
        dphi = pia*x(1) 
        sn   = dsin(dphi)
        cs   = dcos(dphi)
        yi   = x(2)    
        zi   = x(3)    
        if(yi.eq.0.d0) then
          r(1,1) =  1.d0
        else
          rq2  = yi*yi+zi*zi
          r2   = dsqrt(rq2)
          t1   = dsqrt(1.d0+pia*pia*rq2) 
          t2   = x(2)/dsqrt(yi*yi) 
          r(1,1) =  1.d0/t1
          r(2,1) = -sn*pia/t1*t2*r2
          r(3,1) =  cs*pia/t1*t2*r2
        end if
          
        r(2,2) =  cs
        r(3,2) =  sn

        call vecp (r(1,1),r(1,2),r(1,3))
c 
      else
        write(*,*) 'wrong director type', idtyp
      end if
c.... write directors on screen
c     write(*,*) 'director'                
c     write(*,*) 't1',r(1,1),r(2,1),r(3,1)
c     write(*,*) 't2',r(1,2),r(2,2),r(3,2)
c     write(*,*) 't3',r(1,3),r(2,3),r(3,3)

      return
      end
c
      subroutine pedgeb(x,id,idl,ndm,ndf,numnp,prt)
c-----------------------------------------------------------------------
c
c      Purpose: Set boundary constraints along lines        
c               ltyp=1: line
c               ltyp=2: parabola
c               ltyp=3: circle
c
c      Inputs:
c         x(ndm,*)    - Nodal coordinates of mesh
c         id(ndf,*)   - Equation numbers for each active dof
c         ndm         - Spatial dimension of mesh
c         ndf         - Number dof/node
c         numnp       - Number of nodes in mesh
c         prt         - print flag
c
c      Scratch:
c         idl(ndf*ndf)- Local Equation numbers input
c 
c      Outputs:
c         id(ndf,*)   - Equation numbers for each active dof
c
c.... ww,fg bs uka 11/95                                               |
c-----------------------------------------------------------------------
      USE bdata
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      character*8  ltyp(3)
      character*20 ltyp1(3)
      dimension x(ndm,numnp),id(ndf,numnp),idl(ndf*ndf),td(16)
      dimension dp(3),p(3),p1(3),p2(3),p3(3),a(3),dx(3),dn(3),dh1(3),
     +          dh2(3),dh3(3)
      dimension idt(6)
      data blank/-999.d0/,tol/1000.d0/
      data ltyp/'line    ','parabola','circle  '/
      data ltyp1/' intermediate point ',' intermediate point ',
     +           ' origin of circle   '/
      data eps /0.001d0/
      pi2 = datan(1.0d0)*8.0d0
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
      call pzero(td ,16)
      call pzeroi(idt,6)
      ltd = 3*ndm + 3 
c.... read input of boundary edge data
100   if(ior.lt.0) write(*,3001)
      call dinput(td,ltd)
      if(errck) go to 100
      dii = dsqrt(dot(td,td,ndm*3))
      if(dii.le.0) go to 4
      toluser = td(ltd-2)
      iadd    = td(ltd-1)
      if(dabs(toluser).gt.1.e-10) tol = toluser
      radius  = td(ltd)
      do i = 1,ndm
        p1(i) = td(i      )
        p2(i) = td(i+  ndm)
        p3(i) = td(i+2*ndm)
      end do

c.... read input of boundary restraints - limit is 16 nos. / record
106   if(ior.lt.0) write(*,3002)
      il = min(ndf,16)
      call dinput(td,il)
      if(errck) go to 106
      do j = 1,min(ndf,16)
        idl(j) = td(j)
      end do     
      if(ndf.gt.16) then
        do 105 ii = 1,(ndf)/16
          is = il+1
          il = min(is+15,ndf)
103       call dinput(td,il-is+1)
          if(errck) go to 103
          do 104 k = 1,il-is+1
            idl(k+is-1) = td(k)
104       continue
105     continue
      end if
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
cww     if(ndm.eq.3) stop 'edge for circle not in 3d in SR pedgeb'
c....   calculate boundary point 1 of circle line
        dxp1 = p1(1)-p3(1)
        dyp1 = p1(2)-p3(2)
        fr1 = dabs(dxp1*dxp1 + dyp1*dyp1 - radius**2)
        fr1 = dsqrt(fr1)
        if(fr1.ge.tolc) then
          write(iow,2003)
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
          write(iow,2003)
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
c....   line
        ityp = 1
      end if
c.... print base data
      if(prt) write(iow,2002)              ltyp(ityp),p1,p2,ltyp1(ityp),
     +        p3,iadd,toluser,radius,(idl(i),i=1,ndf)
      if(prt.and.ior.lt.0) write(*  ,2002) ltyp(ityp),p1,p2,ltyp1(ityp),
     +        p3,iadd,toluser,radius,(idl(i),i=1,ndf)
c
c.... find nodes 
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
c.....    special case
          if(phin.eq.0.d0.and.phi1.gt.0.d0) phin = phin+pi2
c
          if(phin.lt.phi1-eps) goto 200
          if(phin.gt.phi2+eps) goto 200
        end if
        if(dd.lt.tolc) then
c......   set boundary conditions for point
          do j = 1,ndf
            if(iadd.eq.0) id(j,n) = max(abs(id(j,n)),abs(idl(j))) ! add
            if(iadd.eq.1) id(j,n) = idl(j)                        ! set
          end do
        end if
200   continue
      go to 100
cww4  if(prt)              write(iow,2000) o,head,(i,i=1,ndf)
cww   if(prt.and.ior.lt.0) write(*  ,2000) o,head,(i,i=1,ndf)
4     if(prt)              write(iow,2000) (i,i=1,ndf)
      if(prt.and.ior.lt.0) write(*  ,2000) (i,i=1,ndf)
      do 410 n = 1,numnp
        do 408 i = 1,ndf
408       if(id(i,n).ne.0) go to 409
        go to 410
409     continue
        do i=1,ndf
          idt(i) =idt(i)+id(i,n) ! count no of bc/dof
        end do
        if(prt)              write(iow,2001) n,(id(i,n),i=1,ndf)
        if(prt.and.ior.lt.0) write(*  ,2001) n,(id(i,n),i=1,ndf)
410   continue

      if(prt)              write(iow,2004) (idt(i),i=1,ndf)
      if(prt.and.ior.lt.0) write(*  ,2004) (idt(i),i=1,ndf)

      return
cww2000  format(a1,19a4,a3//'  n o d a l  b. c. along line    '/
2000  format(/'  n o d a l  b. c. along line    '/
     1       /4x,'node',9(i3,'-b.c.')/(8x,9(i3,'-b.c.')))
2001  format(10i8/(8x,9i8))
2002  format(/,' BOUN generation with EDGE on ',a8,/,
     +         ' From         Point ',3(1x,e12.5),/,
     +         ' To           Point ',3(1x,e12.5),/,
     +           a20                 ,3(1x,e12.5),/,
     +         ' add bc(0)set bc(1):',   i5,/,
     +         ' User defined Tol   ',e12.5,/,
     +         ' Radius(only ityp=3)',e12.5,/,
     +         ' B.C. Code:         ',6(1x,i2))
2003  format(' EDGE: Point is not on defined circle')
2004  format(' No. of bcs/dof',/,2x,6i8)
3001  format(' Input: P_1(x,y,(z)) P_2(x,y,(z)),<P_3(x,y,(z))>',
     +        'tol,add,r  >',$)
3002  format(' Input: (idl(i),i=1,ndf) >',$)
      end
c

      subroutine pedgeg(x,id,idl,ndm,ndf,numnp,prt)
c----------------------------------------------------------------------
c
c      Purpose: Set boundary constraints at edges + generation of higher       
c               values 
c               set b.c. at edges    1 < idf < 3  like pedges
c               additional: generate values for idf > 3  (repetition!)  
c
c      Inputs:
c         x(ndm,*)    - Nodal coordinates of mesh
c         id(ndf,*)   - Equation numbers for each active dof
c         ndm         - Spatial dimension of mesh
c         ndf         - Number dof/node
c         numnp       - Number of nodes in mesh
c         prt         - print flag
c
c      Scratch:
c         idl(ndf*ndf)- Local Equation numbers input
c 
c      Outputs:
c         id(ndf,*)   - Equation numbers for each active dof
c
c----------------------------------------------------------------------
      USE bdata
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      dimension x(ndm,numnp),id(ndf,numnp),idl(ndf*ndf),td(16)
      data blank/-999.d0/
c.... read input of boundary edge for restraints
100   if(ior.lt.0) write(*,3001)
      il = 3 + 2
      call dinput(td,il)
      if(errck) go to 100
      i  = td(1)
      if(i.le.0.or.i.gt.ndm) go to 4
      x0 = td(2)
      idl(1) = td(3)
      idl(2) = td(4)
      idl(3) = td(5)
c.... generate boundary conditions by repetition
      do j = 4,ndf,3
          idl(j)         = idl(1)
          idl(j+1) = idl(2)
          idl(j+2) = idl(3)
      end do
c
      dx = pdiff(x(i,1),ndm,numnp)/1000.
      do 200 n = 1,numnp
         if(x(i,n).ne.blank.and.abs(x(i,n)-x0).le.dx) then
            do 210 j = 1,ndf
               id(j,n) = max(abs(id(j,n)),abs(idl(j)))
210         continue
         end if
200   continue
      go to 100
cww4     if(prt)              write(iow,2000) o,head,(i,i=1,ndf)
cww      if(prt.and.ior.lt.0) write(*  ,2000) o,head,(i,i=1,ndf)
4     if(prt)              write(iow,2000) (i,i=1,ndf)
      if(prt.and.ior.lt.0) write(*  ,2000) (i,i=1,ndf)
      do 410 n = 1,numnp
        do 408 i = 1,ndf
408       if(id(i,n).ne.0) go to 409
        go to 410
409     if(prt) write(iow,2001) n,(id(i,n),i=1,ndf)
        if(prt.and.ior.lt.0) write(*,2001) n,(id(i,n),i=1,ndf)
410   continue
      return
cww2000  format(a1,19a4,a3//'  e d g e(g) n o d a l    b . c .'/
2000  format(/'  e d g e(g) n o d a l    b . c .'/
     1       /4x,'node',9(i3,'-b.c.')/(8x,9(i3,'-b.c.')))
2001  format(10i8/(8x,9i8))
3001  format(' Input: ndir,x(ndir),(idl(i),i=1,ndf)','   >',$)
      end
c
      subroutine pedges(x,id,idl,ndm,ndf,numnp,etyp,prt)
c----------------------------------------------------------------------
c
c      Purpose: Set boundary constraints at edges        
c
c      Inputs:
c         x(ndm,*)    - Nodal coordinates of mesh
c         id(ndf,*)   - Equation numbers for each active dof
c         ndm         - Spatial dimension of mesh
c         ndf         - Number dof/node
c         numnp       - Number of nodes in mesh
c         etyp        - option for userdefined GAP
c         prt         - print flag
c
c      Scratch:
c         idl(ndf*ndf)- Local Equation numbers input
c 
c      Outputs:
c         id(ndf,*)   - Equation numbers for each active dof
c
c----------------------------------------------------------------------
      USE bdata
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt,gapfl,pcomp
      character etyp*4
      dimension x(ndm,numnp),id(ndf,numnp),idl(ndf*ndf),td(16)
      data blank/-999.d0/
      gap   = 1.0d-3/sqrt(dble(max(1,numnp))) 
      gapfl = .false.
      if(pcomp(etyp,'gap',3)) then
        call dinput(td,1)
        gap   = td(1)
        gapfl = .true.
      end if  
c.... read input of boundary edge for restraints - limit is 16 nos. / record
100   if(ior.lt.0) write(*,3001)
      il = min(ndf+2,16)
      call dinput(td,il)
      if(errck) go to 100
      i  = td(1)
      if(i.le.0.or.i.gt.ndm) go to 4
      x0 = td(2)
      do 101 j = 1,min(ndf,14)
          idl(j) = td(j+2)
101   continue
      if(ndf.gt.14) then
        do 105 ii = 1,(ndf+2)/16
          is = il+1
          il = min(is+15,ndf+2)
103       call dinput(td,il-is+1)
          if(errck) go to 103
          do 104 k = 1,il-is+1
            idl(k+is-3) = td(k)
104       continue
105     continue
      end if
c.... search circle        
      if(gapfl) then
        dx = gap
      else
        dx = pdiff(x(i,1),ndm,numnp)*gap
      endif
c...  control input
      if(prt) then  
        write(iow,2002) i,x0,dx,(ii,ii=1,ndf)
        write(iow,2003) (idl(ii),ii=1,ndf)     
        if(ior.lt.0) then 
          write(*,2002) i,x0,dx,(ii,ii=1,ndf)
          write(*,2003) (idl(ii),ii=1,ndf)     
        end if
      end if  
c.... look for nodes
      do 200 n = 1,numnp
        if(x(i,n).ne.blank.and.dabs(x(i,n)-x0).le.dx) then 
          do 210 j = 1,ndf
            id(j,n) = max(abs(id(j,n)),abs(idl(j)))
210       continue
        end if
200   continue
      go to 100
cww4  if(prt)              write(iow,2000) o,head,(i,i=1,ndf)
cww   if(prt.and.ior.lt.0) write(*  ,2000) o,head,(i,i=1,ndf)
4     if(prt)              write(iow,2000)        (i,i=1,ndf)
      if(prt.and.ior.lt.0) write(*  ,2000)        (i,i=1,ndf)
c407   do 410 n = 1,numnp
      do 410 n = 1,numnp
        do 408 i = 1,ndf
408       if(id(i,n).ne.0) go to 409
        go to 410
409     if(prt) write(iow,2001) n,(id(i,n),i=1,ndf)
        if(prt.and.ior.lt.0) write(*,2001) n,(id(i,n),i=1,ndf)
410   continue
      return
cww2000  format(a1,19a4,a3//'  e d g e    n o d a l    b . c .'//
2000  format(/'  e d g e    n o d a l    b . c .'/
     1       4x,'node',9(i3,'-b.c.')/(8x,9(i3,'-b.c.')))
2001  format(10i8/(8x,9i8))
2002  format(/'  BOUN Generation with EBOUN on axis ',i3,
     1        ' at value = ',e12.5,' with gap = ',e12.5/
     2       4x,9(i3,'-b.c.')/(4x,9(i3,'-b.c.')))
2003  format(4x,9i8/(4x,9i8))
3001  format(' Input: ndir,x(ndir),(idl(i),i=1,ndf)','   >',$)
      end
c
      subroutine pedgesd(x,f,ndm,ndf,numnp,prt)
c----------------------------------------------------------------------
c
c      Purpose: Set prescribed values for boundary constraints at edges        
c
c      Inputs:
c         x(ndm,*)    - Nodal coordinates of mesh
c         f(ndf,*)    - Load-vector
c         ndm         - Spatial dimension of mesh
c         ndf         - Number dof/node
c         numnp       - Number of nodes in mesh
c         prt         - print flag
c
c      Outputs:
c         id(ndf,*)   - Equation numbers for each active dof
c
c      Open: 
c      possible only for 14 dofs
c 
c----------------------------------------------------------------------
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt
      dimension x(ndm,numnp),f(ndf,numnp),td(16)
      data blank/-999.d0/
      td = 0
c.... read input of boundary edge for restraints - limit is 16 nos. / record
100   if(ior.lt.0) write(*,3001)
      il = min(ndf+2,16)
      call dinput(td,il)
      if(errck) go to 100
c...  axis
      i  = td(1)
      if(i.le.0.or.i.gt.ndm) go to 4
c...  value
      x0 = td(2)
c...  code      
      do j = 1,14
        td(j) = td(j+2) 
      end do        
cww      do 101 j = 1,min(ndf,14)
cww          idl(j) = td(j+2)
101   continue
cww      if(ndf.gt.14) then
cww          do 105 ii = 1,(ndf+2)/16
cww            is = il+1
cww            il = min(is+15,ndf+2)
cww103         call dinput(td,il-is+1)
cww            if(errck) go to 103
cww            do 104 k = 1,il-is+1
cww              idl(k+is-3) = td(k)
cww104         continue
cww105       continue
cww      end if
c...  control input
      if(prt) then  
         write(iow,2002) i,x0,(ii,ii=1,ndf)
         write(iow,2003)  (td(ii),ii=1,ndf)     
         if(ior.lt.0) then 
           write(*,2002) i,x0,(ii,ii=1,ndf)
           write(*,2003)  (td(ii),ii=1,ndf)     
         end if
       end if  
c.... search circle        
      dx = pdiff(x(i,1),ndm,numnp)/1000.
c.... look for nodes
      do 200 n = 1,numnp
      if(x(i,n).ne.blank.and.abs(x(i,n)-x0).le.dx) then ! original
cww   if(                    abs(x(i,n)-x0).le.dx) then ! skip nodes(-999) via ebou,1,-999,1
            do 210 j = 1,ndf
               f(j,n)= td(j)
210         continue
         end if
200   continue
      go to 100
cww4  if(prt)              write(iow,2000) o,head,(i,i=1,ndf)
cww   if(prt.and.ior.lt.0) write(*  ,2000) o,head,(i,i=1,ndf)
4     if(prt)              write(iow,2000)        (i,i=1,ndf)
      if(prt.and.ior.lt.0) write(*  ,2000)        (i,i=1,ndf)
c407   do 410 n = 1,numnp
      do 410 n = 1,numnp
        do 408 i = 1,ndf
408       if(f(i,n).ne.0.d0) go to 409
        go to 410
409     if(prt) write(iow,2001) n,(f(i,n),i=1,ndf)
        if(prt.and.ior.lt.0) write(*,2001) n,(f(i,n),i=1,ndf)
410   continue
      return
cww2000  format(a1,19a4,a3//'  e d g e  n o d a l  prescribed values'//
2000  format(/'  e d g e    n o d a l    prescribed b.c. values'/
     1       4x,'node',9(i3,'-b.c.')/(8x,9(i3,'-b.c. val.')))
2001  format(i8,9f8.4/(8x,9f8.4))
2002  format(/'  DISP Generation with EDIS on axis ',i3,
     1        ' at value = ',e12.5/
     2       4x,9(i3,'-b.c.')/(4x,9(i3,'-b.c.')))
2003  format(4x,9f8.4/(4x,9f8.4))
3001  format(' Input: ndir,x(ndir),(f(i),i=1,ndf)','   >',$)
      end
c
      subroutine ppload(x,f,id,ndm,ndf,numnp,prt)
c----------------------------------------------------------------------
c
c      Purpose: set loads and boundary conditions for a point P(x,y,z),                  

c      Comment: - loads and b.c.s are set for all (!) nodes within a 
c               - certain circle. 
c               - values are added 
c               - use MACRO POIN only after nodes are defined 
c               - factor for radius of circle is input, def=1  
c
c      Inputs:
c         x(ndm,*)    - Nodal coordinates of mesh
c         f(ndf,*)    - Nodal force values
c         id(ndf,*)   - Equation numbers for each active dof
c         ndm         - Spatial dimension of mesh
c         ndf         - Number dof/node
c         numnp       - Number of nodes in mesh
c         prt         - print flag
c
c      Outputs:
c         f(ndf,*)    - Nodal force values
c         id(ndf,*)   - Equation numbers for each active dof
c
c      Comments:  
c         return      - if xl=0    
c                     - input of a point with xl=0: fact must be set!!
c
c      Open:
c         ndf <16!    - but no problem to change|
c 
c----------------------------------------------------------------------
      USE bdata
      USE iofile
      implicit double precision (a-h,o-z)
      logical prt,poinf
      dimension x(ndm,numnp),f(ndf,numnp),id(ndf,numnp),
     +          xl(4),fl(16),idl(16),td(16),dx(3),dxp(3)
c.... write title
      if(prt) write(iow,2000)(i,i=1,ndf)
c.... read input 
100   if(ior.lt.0) write(*,3001)
c.... coordinates 
      call dinput(xl,ndm+1)
      fact = xl(ndm+1)
      r1 = dot(xl,xl,ndm)
      if(r1.eq.0.d0.and.fact.eq.0.d0) return
      if(fact.eq.0.d0) fact=1.d0
c.... loads 
      call dinput(fl,ndf)
c.... b.c.
      call dinput(td,ndf)
      call pzeroi(idl,ndf)
      do i=1,ndf
        idl(i) = td(i)
      end do  
c.... search circle
      call pzero(dx ,3)
      do i = 1,ndm
        dx(i) = pdiff(x(i,1),ndm,numnp)/1000.d0/fact
      end do
      tol = dsqrt(dot(dx,dx,3))
c.... find node
      poinf=.false.
      nsp = 0
      do  n = 1,numnp
        call pzero(dxp,3)
        do i = 1,ndm
          dxp(i) = (x(i,n)-xl(i)) * (x(i,n)-xl(i)) 
        end do
        tolp = dsqrt(dot(dxp,dxp,3))
        if(tolp.lt.tol) then
c....     node found 
c....     check, if more than one node is found
          if(poinf) then
            if(prt) write(  *,2003) nsp,n
                    write(iow,2003) nsp,n
          end if ! poinf 
          poinf = .true. ! set node as used, TIE not available here!
          nsp = n
c....     set values for this node
          do i=1,ndf
c....       from input file id=0/1            
            if(ior.ge.0.and.id(i,n).eq.0) f(i,n) = f(i,n) + fl(i)
c....       interactive  id=-i/k via macro>mesh         
            if(ior.lt.0.and.id(i,n).gt.0) f(i,n) = f(i,n) + fl(i) 
c....       for bcs  only from input file       
            if(ior.ge.0) id(i,n) = max( abs(id(i,n)),idl(i) )
          end do
c...      print values
          if(prt) write(iow,2001) n, ( fl(i),i=1,ndf)
          if(prt) write(iow,2002) n, (idl(i),i=1,ndf)
        end if ! tolp
      end do ! n
      go to 100
      return
2000  format(/'  p o i n t    l o a d   and    b . c .'/
     1       2x,'node',9(i3,'-dof')/(8x,9(i3,'-dof')))
2001  format(i6,2x,8(e12.5,1x)/(8x,8(e12.5,1x)))
2002  format(i6,2x,8(i12,1x)  /(8x,8(i12,1x)))
2003  format(2x,'#### more than 1 node found in Macro POIN',/,
     +       2x,' org node =',i5,' second node found =',i5,/,
     +       2x,' this is ok, if second node is skipped by macro TIE!')  

3001  format(' Input: x(1,ndm)//F(1,ndf)//bc(1,ndf)  >',$)
      end
c
c----------------------------------------------------------------------
c
      subroutine pform(ul,xl,tl,ld,p,s,h1,h2,h3,ie,d,id,x,ix,f,t,jp,f0,
     1  u,ud,b,aq,c,ang,angl,ndd,nie,ndf,ndm,nen1,nst,afl,bfl,cfl,dfl,
     2    isw,nn1,nn2,nn3)
c----------------------------------------------------------------------
c
c      Purpose: Compute element arrays and assemble global arrays
c
c      Inputs:
c         ie(nie,numat)  - Assembly information for material set
c         nie =ndf+3
c         ie(1-ndof,ma)  - local element dofs 
c         ie(nie-2, ma)  - nh3=no of time ind. history terms for mat ma   
c         ie(nie-1, ma)  - iel - Element type number
c         ie(nie  , ma)  - nh1=no of history terms for material ma     
c         d(ndd,*)       - Material set parameters
c         id(ndf,*)      - Equation numbers for each active dof
c         x(ndm,numnp)   - Nodal coordinates of mesh
c         ix(nen1,numel) - Element nodal connections of mesh
c         nen1 = nen+4   
c         ix(1-nen,*)    - Element nodes
c         ix(nen+1,*)    - nt1  global adress history terms h1m
c         ix(nen+2,*)    - nt2  global adress history terms h2m
c         ix(nen+3,*)    - nt3  global adress history terms h3m
c         ix(nen+4,*)    - ma   nen+4=nen1
c         f(ndf,*,2)     - Nodal force and displacement values
c         t(numnp)       - Nodal temperature values
c         jp(*)          - Pointer array for row/columns of tangent
c         f0(ndf)        - Nodal initial force values
c         u(ndf,numnp*3) - Nodal solution values
c         ud(*)          - Nodal rate values
c         ang(numnp)     - Rotation angles at global nodes   
c         ndd            - Dimension for d array
c         nie            - Dimension for ie array
c         ndf            - Number dof/node
c         ndm            - Spatial dimension of mesh
c         nen1           - Dimension for ix array
c         nst            - Dimension for element array
c         afl            - Flag, assemble coefficient array if true
c         bfl            - Flag, assemble vector if true
c         cfl            - Flag, coefficient array unsymmetric if true ??
c         dfl            - Flag, assemble reactions if true
c         isw            - Switch to control quantity computed
c         nn1            - First element number to process
c         nn2            - Last element number to process
c         nn3            - Increment to nn1
c         m              - m-array
c    
c      Scratch:
c         ul(nst,6)      - Element solution and rate values
c         xl(ndm,nen)    - Element nodal coordinates
c         tl(nen)        - Element nodal temperatures
c         angl(nen)      - Rotation angles at element nodes   
c         ld(nst)        - Element local/global equation numbers
c         p(nst)         - Element vector
c         s(nst,nst)     - Element array
c         h1(nhmax)      - History data element
c         h2(nhmax)      - History data element
c         h3(nh3max)     - History data element
c
c      Outputs:
c         b(*)           - Global vector
c         aq(*)          - Global matrix, diagonal and upper part
c         al(*)          - Global matrix, lower part
c
c      Comments:
c         ul(nst,1)      - Element displacements                   =   u
c         ul(nst,2)      - Element displacements since last TIME   =  du
c         ul(nst,3)      - Element displacements of last iteration = ddu
c
c         extended System k>0
c         ul(nst,4)      - Element displacements                 c = h2 
c         ul(nst,5)      - Element displacements                 d = h1 
c         ul(nst,6)      - Element displacements                 e = h3 
c
c         dynamic         k>0
c         ul(nst,1)      - Element displacements                   = u 
c         ul(nst,4)      - Element displacements                   = v  
c         ul(nst,5)      - Element displacements                   = a  
c
c         k<0(B.C.)
c         ul(nst,6)      - initial displacements at b.c.s
c
c
c      Open:
c
               
c
c----------------------------------------------------------------------
      USE arcext
      USE cdata
      USE ddata
      USE debugs
      USE edgdat
      USE eldata
      USE ext1
      USE ext2
      USE fdata
      USE feapmpi
      USE hdata
      USE hdatam
      USE iofile
      USE iscsr
      USE mdat2
      USE prlod
      USE smpak
      USE soltyp
      USE tdata
      USE uneig
      implicit real*8 (a-h,o-z)
c.... Declare variable types
      logical afl,bfl,cfl,dfl,efl,gfl,edf
      logical parforms
c.... Declare array types
      integer ld(nst), ie(nie,*), id(ndf,*), ix(nen1,*), jp(*)
      real*8  x(ndm,*),b(*),aq(*),c(*),t(*),f(ndf,*),f0(ndf,*),ang(*)
      real*8  xl(ndm,nen),p(nst),s(nst,nst),ul(nst,6),tl(nen)
      real*8  h1(nhmax),h2(nhmax),h3(nh3max),angl(nen)
      real*8  d(ndd,*),u(ndf,*),ud(*)
      real*8  xl1(ndm,nen),p1(nst),s1(nst,nst),temp(nhmax),temp3(nh3max)
c

      !dimension m(*) 
c
c.... set up local arrays before calling element library
      iel = 0
      edf = nde*ned.gt.0
      efl = .false.
      if(.not.dfl.and.isw.eq.6) efl = .true.  ! assemble p + not REAC
      if(bfl.and.isw.eq.3)      efl = .true.  ! assemble p + TANG
      gfl = (efl.and.fl(1)).and.fl(9)
      nl1 = ndf*nen + 1
      numnp2 = numnp + numnp
      nneq = numnp*ndf
      nrkn = nrk*nneq - nneq
      nrcn = nrc*nneq - nneq
      nrmn = nrm*nneq - nneq
c.... show state ! kann parallel nicht laufen, stoert aber nicht! ww
      call perform(0,nn2,1,isw) 
      npos = (nn2-nn1)/nn3+1

c.... Parallelization of element loop 
c     Exceptions:
c.... 02=check: 
c.... 04=prin: Reihenfolge geht verloren 


c     Exceptions if not ordered:
c     08=stre: Problem dt=dt+xsji auf El.ebene, dto für st, 
c     09=erro: wie 08
c     13=forc: wie 08
c     14=str1: wie 08

      if (isw.eq. 2 .or. isw.eq. 4 .or.isw.eq. 8 .or.isw.eq. 9
     +   .or.isw.eq.13 .or. isw.eq.14 .or.nn1.eq.nn2) then
        parforms = parform ! save
        parform=.false.
      end if     

!$OMP  PARALLEL DEFAULT(NONE) 
!$OMP& IF(parform) 
!$OMP& PRIVATE(
!$OMP&    nt1,nt3,ndh1,ndh3,nrot,dbgtxt,un,dun,dan,
!$OMP&    h1,h2,h3,i,ii,iid,ild,j,jj,k,ul,tl,angl,ld,s,p,xl
!$OMP&    ,eh1,eh2,eh3,temp,temp3)
!$OMP&    SHARED(nn1,nn2,nn3,ix,nen1,nen,ie,nie,debug,iow,
!$OMP&    nst,ndm,ndf,t,x,id,u,extflg,numnp,
!$OMP&    numnp2,nmode,b,ang,aq,kex,c,fl,nrk,ud,nrkn,nrc,nrcn,nrmn,
!$OMP&    neq,f0,f,prop,dfl,edf,nedg,nl1,ne4,nde,ia,
!$OMP&    itrot,d,isw,efl,gfl,npos,jp,
!$OMP&    afl,bfl,istyp,japt,lodrv,cfl,hflgu,h3flgu,
!$OMP&    parform,s1,p1,xl1,
!$OMP&    nhmax,nh3max,gh1,gh2,gh3,csrja,edge1,edge2,edge3,edge4,
!$OMP&    extkc,extkd,extke,eigmma,smsperm1)      

c.... ordered schedule ohne Critical  
cww!$OMP DO ORDERED SCHEDULE(DYNAMIC)  
c.... ordered mit Critical  
!$OMP DO SCHEDULE(DYNAMIC) 

c.... main loop over elements  
      do 110 n = nn1,nn2,nn3

c....   local material and local history adresses
        ma  = ix(nen1,n)
cww     nt1 = ix(nen+1,n)  ! same as below  
cww     nt2 = ix(nen+2,n) 
cww     nt3 = ix(nen+3,n) 
        nt1 = (n-1)*nhmax+1 
        nt3 = (n-1)*nh3max+1
c
c....   in case of element does not exist(node(1)=0) ww        
        if(ix(1,n).eq.0) goto 110
c        
c....   history variables: move el. values into working arrays h1,h2,h3
        if(ie(nie,ma).gt.0) then
          if(debug.eq.1.and.n.eq.1) then
            dbgtxt = ' Debug: move history variables h1,h2'
            write(*  ,1000) dbgtxt
            write(iow,1000) dbgtxt
            dbgtxt = ' '
          end if
          h1=gh1(nt1:nt1+nhmax-1)
          h2=gh2(nt1:nt1+nhmax-1)
        end if
        if(ie(nie-2,ma).gt.0) then
          if(debug.eq.1.and.n.eq.1) then
            dbgtxt = ' Debug: move history variables h3'
            write(*  ,1000) dbgtxt
            write(iow,1000) dbgtxt
            dbgtxt = ' '
          end if
          h3=gh3(nt3:nt3+nh3max-1)
        end if
        if(ie(nie-1,ma).ne.iel) mct = 0
c....   Element type number
        iel = ie(nie-1,ma)
        
c....   zero local arrays 
        call pzero (ul,6*nst)
        call pzero (xl,nen*ndm)
        call pzero (tl,nen)
        call pzero (angl,nen)
        call pzeroi(ld,nst)

        un = 0.d0
        dun= 0.d0
        dan= 0.d0

c....   local angles on element nodes
        call pangl(ix(1,n),nen,angl,ang,nrot)

c....  local temp, coord,displ,veloc,accel on element nodes
        do 108 i = 1,nen
          ii = ix(i,n)
          if(ii.gt.0) then
            iid = ii*ndf - ndf
            ild =  i*ndf - ndf
            nel = i
            tl(i) = t(ii)
            do j = 1,ndm
              xl(j,i) = x(j,ii)
            end do  
            do 107 j = 1,ndf
              jj = ie(j,ma)
              if(jj.gt.0) then
                k = id(jj,ii)
                ul(j+ild,1) = u(jj,ii)
                if(.not.extflg) then
                  ul(j+ild,2) = u(jj,ii+numnp)
                  ul(j+ild,3) = u(jj,ii+numnp2)
                end if
                if(extflg.and.nmode.eq.1.and.k.gt.0) then
                  ul(j+ild,4) = b(k)
                  ul(j+ild,5) = aq(k)
                  if(kex.ne.0) ul(j+ild,6) = c(k)
                end if
                if(k.gt.0.and.fl(9)) then
                  if(nrk.gt.0) ul(j+ild,1) = ud(nrkn+k)
                  if(nrc.gt.0) ul(j+ild,4) = ud(nrcn+k)
                               ul(j+ild,5) = ud(nrmn+k)
                else if(k.lt.0) then
                  if(fl(9)) then
                    if(nrk.gt.0) ul(j+ild,1) = ud(nrkn+neq-k)
                    if(nrc.gt.0) ul(j+ild,4) = ud(nrcn+neq-k)
                                 ul(j+ild,5) = ud(nrmn+neq-k)
                    dan  = max(dan,abs(ul(j+ild,5)))
                  end if
                  ul(j+ild,6) = f0(jj,ii) + f(jj,ii)*prop - u(jj,ii)
                  dun         = max(dun,abs(ul(j+ild,6)))
                end if
                un  = max( un,abs(ul(j+ild,1)))
                if(dfl) k = iid + jj
                ld(j+ild) = k
              end if
107         continue
          end if
108     continue

c....   set up edge values if necessary
        if(edf.and.nedg(iel).gt.0) then
          call pforme(iel,ix(1,n),ld(nl1),ul(nl1,1),
     1         edge1,edge2,edge3,edge4,
     2         edge4(ne4),un,dun,nde,nst)
        end if

c....   modify local element solution array for transient algorithms
        if(fl(9).and.(dun.gt.1.0d-7*un)) call modu(ul,nst)

c....   transform all(!) mixed displ. to pure global dir.
        if(nrot.gt.0)
     1  call ptrans(ia,itrot,angl,ul,p,s,nel,nen,ndf,nst,1)

        dm = prop ! for what ww??

c....   calculate data on element level
        call elmlib(d(1,ma),ul,xl,ix(1,n),tl,s,p,h1,h2,h3,
     +              ndf,ndm,nst,iel,isw)

        if(isw.eq.19) goto 110

c....   for extended system: compute derivatives
cww!$OMP ORDERED
csk!$OMP CRITICAL
        if(extflg .and. nmode .eq. 1) then
c....     vector c = h2
          call pmult(s,p,nst,ul(1,4),2)
          if(nrot.gt.0)
     1    call ptrans(ia,itrot,angl,ul,p,s,nel,nen,ndf,nst,3)
          if(istyp.ne.0 ) stop 'extended system only with solver 0' 
          call dasbly(s,p,ld,jp,nst,.false.,.false.,.true.,
     1                extkc,c,aq,aq)

c....     vector d = h1
          call pmult(s,p,nst,ul(1,5),2)
          if(nrot.gt.0)
     1    call ptrans(ia,itrot,angl,ul,p,s,nel,nen,ndf,nst,3)
          call dasbly(s,p,ld,jp,nst,.false.,.false.,.true.,
     1                extkd,c,aq,aq)

c....     vector e = h3 (in case of singularity)
          if(kex.ne.0) then
            call pmult(s,p,nst,ul(1,6),2)
            if(nrot.gt.0)
     1      call ptrans(ia,itrot,angl,ul,p,s,nel,nen,ndf,nst,3)
            call dasbly(s,p,ld,jp,nst,.false.,.false.,.true.,
     1                  extke,c,aq,aq)
          end if

        else ! transform s,p to mixed global displacements back

          if(nrot.gt.0)
     1    call ptrans(ia,itrot,angl,ul,p,s,nel,nen,ndf,nst,2)

c....     add to total array with respect to solver type
          if(isw.eq.21) then ! ueig
            call assemb(s,ix,id,ie,ma,neq,nen1,ndf,nie,nen,n,nst,
     1                  eigmma)

            go to 111 
          else
c...        assemble coefficient array (afl) or vector (bfl) if true
            if(afl.or.bfl) then
c...          assemble with respect to solver
                     
              if(istyp.eq.0 ) then  ! Standard
                call dasbly(s,p,ld,jp,nst,cfl,afl,bfl,b,c,aq(neq+1),aq)
              else if(istyp.eq.1.or.istyp.eq.2)then  ! SM
                call dasbl2(s,p,ld,jp,nst,afl,bfl,b,aq,csrja,
     1                      smsperm1,lodrv)
              else if(istyp.ge.3.and.istyp.le.8)then ! all other CSR
                call dasbl_csr(s,p,ld,jp,nst,afl,bfl,b,aq,csrja)
              end if ! istyp
                   
            end if   ! afl,bfl
          end if     ! ueig
        end if       ! extflg 


c....   history variables: move working arrays nh1,nh2,nh3 
c               back into element arrays nt1,nt2,nt3 if update(hflgu)
c               also save in case need to recompute tangent.
        if(hflgu .and. ie(nie,ma).gt.0) then
          if(debug.eq.1.and.n.eq.1) then
            dbgtxt = ' Debug: save history variables h1,h2'
            write(*  ,1000) dbgtxt
            write(iow,1000) dbgtxt
            dbgtxt = ' '
          end if
cww        call copyhis(h1,h2,m(nt1),m(nt2),ndh1,3)
          !call copyhis(h1,h2,eh1(nt1),eh2(nt1),ndh1,3)
          temp=gh1(nt1:nt1+nhmax-1)
          gh1(nt1:nt1+nhmax-1)=h1
          h1=temp
          temp=gh2(nt1:nt1+nhmax-1)
          gh2(nt1:nt1+nhmax-1)=h2
          h2=temp
        end if
        if(h3flgu .and. ie(nie-2,ma).gt.0) then  
          if(debug.eq.1.and.n.eq.1) then
            dbgtxt = ' Debug: save history variables h3'
            write(*  ,1000) dbgtxt
            write(iow,1000) dbgtxt
            dbgtxt = ' '
          end if
cww          call copyhis(h3,h3,m(nt3),m(nt3),ndh3,4)
          temp3=gh3(nt3:nt3+nh3max-1)
          gh3(nt3:nt3+nh3max-1)=h3
          h3=temp3
        end if
c....   modify for non-zero displacement boundary conditions

        if(efl.and.dun.gt.1.0d-7*un) then
c....     get current element tangent matrix
          if (.not.afl) then
            dm = prop ! for what ww??
            call elmlib(d(1,ma),ul,xl,ix(1,n),tl,s,p,h1,h2,h3,
     1                  ndf,ndm,nst,iel,3)
            if(nrot.gt.0)
     1      call ptrans(ia,itrot,angl,ul,p,s,nel,nen,ndf,nst,2)
          end if
c....     perform modify for displacements
c....     ul(1,6)  must be  mixed(!) displacements, added ww      
          if(nrot.gt.0)
     1      call ptrans1(ia,itrot,angl,ul(1,6),nel,ndf,nst,2)
                   
          call modify(b,ld,s,ul(1,6),nst)
                   
        end if

c....   modify for non-zero acceleration boundary conditions (cmas only)
        if(gfl.and.dan.gt.0.0d0) then
c....     get current element mass matrix
          dm = prop   ! for what ww??
          call elmlib(d(1,ma),ul,xl,ix(1,n),tl,s,p,h1,h2,h3,
     1                ndf,ndm,nst,iel,5)
          if(nrot.gt.0)
     1    call ptrans1(ia,itrot,angl,ul,p,s,nel,nen,ndf,nst,2)
c....     perform modify for accelerations
c....     ul(1,5)  must be  mixed(!) accelerations, added ww      
          if(nrot.gt.0)
     1      call ptrans(ia,itrot,angl,ul(1,5),nel,ndf,nst,2)
                     
          call modify(b,ld,s,ul(1,5),nst)
                    
        end if

 111    continue
cww!$OMP END ORDERED
csk!$OMP END CRITICAL
c....   show actual state
        call perform(n,npos,2,isw)

       if (nn1.eq.nn2) then       
        s1  = s
        p1  = p
        xl1 = xl
       endif 
110   continue ! end loop elements

!$OMP END DO
!$OMP END PARALLEL
      if (nn1.eq.nn2) then       
       s  = s1
       p  = p1
       xl = xl1
      endif 

      if (isw.eq. 2 .or. isw.eq. 4 .or.isw.eq. 8 .or.isw.eq. 9
     +.or.isw.eq.13 .or. isw.eq.14 .or.nn1.eq.nn2) parform = parforms ! reset

      call perform(npos,npos,3,isw)
      return 
c
c     formats
1000  format(a80)
      end
c
      subroutine modu(ul,nst)
c----------------------------------------------------------------------
c
c      Purpose: modify the localized solution vector for each dynamic 
c               algorithm
c
c      Inputs:
c        ul(nst,*) - element vectors
c        ul(nst,1) - solution vector or ur-nrk
c        ul(nst,2) - current total increment for t-n to t-n+1
c        ul(nst,3) - last incremental correction in step t-n+1
c        ul(nst,4) - rate vector ur-nrc
c        ul(nst,5) - rate vector ur-nrm
c        ul(nst,6) - current increment to boundary solutions.
c         nst      - Dimension of element arrays

c      Outputs:
c        ul(nst,*) - element vectors
c----------------------------------------------------------------------
c
c..... Declare variable types
      USE ddata
      USE tdata
      integer i, nst

c..... Declare array types
      real*8  ul(nst,6)
c
c.... modify transient algorithms for essential boundary conditions
      if(nop.eq.1.or.nop.eq.3.or.nop.eq.4) then  ! Newmark, HHT???,Genalpha?? ww
        do 100 i = 1,nst
          ul(i,4) = ul(i,4) + c2*ul(i,6)
          ul(i,5) = ul(i,5) + c1*ul(i,6)
100     continue
      end if
c
      end
c
cww      subroutine pglobx(igi,xg,x,ngn,ndf1,ndm)
cwwc----------------------------------------------------------------------
cwwc
cwwc      Purpose: set the global coordinates
cwwc
cwwc      Inputs:
cwwc
cwwc      Outputs:
cwwc----------------------------------------------------------------------
cwwc..... Declare variable types
cww      integer i, in, n, ngn, ndf1, ndm
cwwc..... Declare array types
cww      integer igi(ndf1,*)
cww      real*8  xg(ndm,*), x(ndm,*)
cwwc.... set the global coordinates
cww      do 100 n = 1,ngn
cww        in = igi(1,n)
cww        do 50 i = 1,ndm
cww          xg(i,n) = x(i,in)
cww50      continue
cww100   continue
cwwc
cww      end
c
      subroutine pforme(iel,ix,ld,ul,ien,inn,ieb,ef,u,un,dun,ndg,nst)
c----------------------------------------------------------------------
c
c      Purpose: PFORM for edges
c
c      Inputs:
c
c      Outputs:
c----------------------------------------------------------------------
c..... Declare variable types
      USE cdata
      USE edgdat
      USE prlod
      integer iel, ndg, nst

      integer i, j, k, i1,i2, ie, jj, n1, n2, nl2, nl3, nl4, nl5
      real*8  un, dun
c..... Declare array types
      integer ix(*), ld(*), ien(*), inn(*), ieb(ndg,*)

      real*8  ul(nst,6), ef(ndg,*), u(ndg,*)
c.... set up pointers
      nl2 = ien(numnp)
      nl3 = nl2 + nl2
      nl4 = nl3 + nl2
      nl5 = nl4 + nl2
      k  = 0
      do 130 i = 1,nedg(iel)
c.... get nodal pairs for the edge
        i1 = iedg(i,iel)
        i2 = mod(i,nedg(iel)) + 1
        i2 = iedg(i2,iel)
        n1 = min(ix(i1),ix(i2))
        n2 = max(ix(i1),ix(i2))
        if(n1.gt.0) then
          i1 = ien(n1) +1
          i2 = ien(n1+1)
          do 120 ie = i1,i2
            if(inn(ie).eq.n2) then
              do 110 j = 1,nde
                ul(k+j,1) = u(j,ie)
                ul(k+j,2) = u(j,ie+nl2)
                ul(k+j,3) = u(j,ie+nl3)
                ul(k+j,4) = u(j,ie+nl4)
                ul(k+j,5) = u(j,ie+nl5)
                jj = ieb(j,ie)
                ld(k+j)   = jj
                if(jj.le.0) then
                  ul(k+j,6) = prop*ef(j,ie) - u(j,ie)
                  dun         = max(dun,abs(ul(k+j,6)))
                end if
                un  = max( un,abs(ul(k+j,1)))
110           continue
            end if
120       continue
        end if
        k = k + nde
130   continue
c
      end
c
      subroutine pfuncs(x,v,val,nex,error)
c----------------------------------------------------------------------
c
c      Purpose: Evaluate expressions with functions in input records for
c               inc dec  int  abs   sqrt
c               sin sind asin asind sinh
c               cos cosd acos acosd cosh
c               tan tand atan atand tanh
c               log exp
c      extension d: input in degrees 
c
c      Inputs:
c         x(*)   - String of input
c         v(*)

c      Outputs:
c         nex    - Number of www(*) used to hold function value
c         error  - Flag, true if error occurs
c         val    - Expression value
c----------------------------------------------------------------------
      USE conval
      implicit double precision (a-h,o-z)
      character*1 x(*),yy(75),xx*5
      logical pcomp,error
      dimension v(*)

c.... rad->deg
      rtod = datan(1.d0)/45.d0 

c.... evaluate functions in expressions

      yy = ' '
      k  = 0
      i  = 1
140   continue
       k = k + 1
        xx(1:1) = x(i)
        xx(2:2) = x(i+1)
        xx(3:3) = x(i+2)
        xx(4:4) = x(i+3)
        xx(5:5) = x(i+4)
        if(     pcomp(xx,'atand',5))then
          kk = 1
        else if(pcomp(xx,'asind',5))then
          kk = 2
        else if(pcomp(xx,'acosd',5))then
          kk = 3
        else if(pcomp(xx,'atan',4))then
          kk = 4
        else if(pcomp(xx,'asin',4))then
          kk = 5
        else if(pcomp(xx,'acos',4))then
          kk = 6
        else if(pcomp(xx,'cosh',4))then
          kk = 7
        else if(pcomp(xx,'sinh',4))then
          kk = 8
        else if(pcomp(xx,'tanh',4))then
          kk = 9
        elseif(pcomp(xx,'cosd',4))then
          kk = 10
        else if(pcomp(xx,'sind',4))then
          kk = 11
        else if(pcomp(xx,'tand',4))then
          kk = 12
        else if(pcomp(xx,'sqrt',4))then
          kk = 13
        else if(pcomp(xx,'exp',3)) then
          kk = 14
        else if(pcomp(xx,'sin',3)) then
          kk = 15
        else if(pcomp(xx,'cos',3)) then
          kk = 16
        else if(pcomp(xx,'tan',3)) then
          kk = 17
        else if(pcomp(xx,'abs',3)) then
          kk = 18
        else if(pcomp(xx,'int',3)) then
          kk = 19
        else if(pcomp(xx,'log',3)) then
          kk = 20
        else if(pcomp(xx,'inc',3)) then
          kk = 21
        else if(pcomp(xx,'dec',3)) then
          kk = 22
        else
          kk = 0
        endif

c       Evaluate functions
        if(kk.ne.0 ) then
c         Functions 1 to 3
          if(kk.le.3) then
            j = 5
c         Functions 4 to 13
          elseif(kk.le.13) then
            j = 4
c         Functions 14 to 22
          else
            j = 3
          endif

          nn = ichar(x(i+j)) - 64
          if(nn.gt.26) then
            jj = nn + ichar(x(i+j+1))
            if    (jj.ge.ichar('a') .and. jj.le.ichar('z')) then
              ii = jj - ichar('a') + 1
            elseif(jj.ge.ichar('0') .and. jj.le.ichar('9')) then
              ii = jj - ichar('0') + 27
            else
              ii = 0
            endif

            val = vvv(nn-32,ii)
          else
            val = www(nn)
          endif

          nex = nex + 1
          if(     kk.eq.1) then
            www(nex) = atan(val)/rtod
          else if(kk.eq.2) then
            www(nex) = asin(val)/rtod
          else if(kk.eq.3) then
            www(nex) = acos(val)/rtod
          else if(kk.eq.4) then
            www(nex) = atan(val)
          else if(kk.eq.5) then
            www(nex) = asin(val)
          else if(kk.eq.6) then
            www(nex) = acos(val)
          else if(kk.eq.7) then
            www(nex) = cosh(val)
          else if(kk.eq.8) then
            www(nex) = sinh(val)
          else if(kk.eq.9) then
            www(nex) = tanh(val)
          else if(kk.eq.10) then
            www(nex) = cos(val*rtod)
          else if(kk.eq.11) then
            www(nex) = sin(val*rtod)
          else if(kk.eq.12) then
            www(nex) = tan(val*rtod)
          else if(kk.eq.13) then
            www(nex) = sqrt(val)
          else if(kk.eq.14) then
            www(nex) = exp(val)
          else if(kk.eq.15) then
            www(nex) = sin(val)
          else if(kk.eq.16) then
            www(nex) = cos(val)
          else if(kk.eq.17) then
            www(nex) = tan(val)
          else if(kk.eq.18) then
            www(nex) = abs(val)
          else if(kk.eq.19) then
            www(nex) = int(val)
          else if(kk.eq.20) then
            www(nex) = log(val)
          else if(kk.eq.21) then
            www(nex) = padd(val)
          else if(kk.eq.22) then
            www(nex) = psub(val)
          endif
          yy(k) = char(nex+64)
          i = i + j
        else
          yy(k) = x(i)
        endif
        i = i + 1
        if(i.lt.75) go to 140   ! fehler mit fullcheck 
c       if(i.lt.73) go to 140   

c.... final evaluation of expression
      call evalex(yy,v,val,k,error)
      return
      end
c
      subroutine pgauss(l,lint,r,z,w)
c----------------------------------------------------------------------
c
c      Purpose: Form Gauss points and weights for two dimensions
c               1-4 points
c-----------------------------------------------------------------------  
c               5-20 points                         TH KIT 06/12                                      
c-----------------------------------------------------------------------
c
c      Inputs:
c         l       - Number of points/direction
c
c      Outputs:
c         lint    - Total number of points
c         r(*)    - 1-direction Gauss point
c         z(*)    - 2-direction Gauss point
c         w(*)    - Gauss weight
c
c      Position
c         1   2      3       4
c         1   4 3    4 7 3    4  3  2  1    
c             1 2    8 9 6    8  7  6  5
c                    1 5 2   12 11 10  9
c                            16 15 14 13
c----------------------------------------------------------------------
      USE eldata
      implicit double precision (a-h,o-z)
      dimension lr(9),lz(9),lw(9),r(*),z(*),w(*),g4(4),h4(4)
      dimension g5(5),h5(5),g6(6),h6(6),g7(7),h7(7),g8(8),h8(8)
      dimension g9(9),h9(9),g10(10),h10(10),g11(11),h11(11)
      dimension g12(12),h12(12),g13(13),h13(13),g20(20),h20(20)
      dimension g14(14),h14(14),g15(15),h15(15),g16(16),h16(16)
      dimension g17(17),h17(17),g18(18),h18(18),g19(19),h19(19)
      data lr/-1,1,1,-1,0,1,0,-1,0/,lz/-1,-1,1,1,-1,0,1,0,0/
      data lw/4*25,4*40,64/
      lint = l*l
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),l
c.... 1x1 integration
1     r(1) = 0.d0
      z(1) = 0.d0
      if(nel.eq.3) z(1) = -1.d0/3.0d0
      w(1) = 4.d0
      return
c.... 2x2 integration
2     g = 1.d0/dsqrt(3.d0)
      do 21 i = 1,4
      r(i) = g*lr(i)
      z(i) = g*lz(i)
21    w(i) = 1.d0
      return
c.... 3x3 integration
3     g = dsqrt(0.6d0)
      h = 1.d0/81.d0
      do 31 i = 1,9
      r(i) = g*lr(i)
      z(i) = g*lz(i)
31    w(i) = h*lw(i)
      return
c.... 4x4 integration
4     g = dsqrt(4.8d0)
      h = dsqrt(30.0d0)/36.d0
      g4(1) = dsqrt((3.d0+g)/7.d0)
      g4(4) = - g4(1)
      g4(2) = dsqrt((3.d0-g)/7.d0)
      g4(3) = -g4(2)
      h4(1) = 0.5d0 - h
      h4(2) = 0.5d0 + h
      h4(3) = 0.5d0 + h
      h4(4) = 0.5d0 - h
      i = 0
      do 41 j = 1,4
      do 41 k = 1,4
      i = i + 1
      r(i) = g4(k)
      z(i) = g4(j)
      w(i) = h4(j)*h4(k)
41    continue
      return
      
c... 5x5 integration
5     g5(1) = 0.9061798459386639927976269d0
      g5(5) = - g5(1)
      g5(2) = 0.5384693101056830910363144d0
      g5(4) = -g5(2)
      g5(3) = 0.d0

      h5(1) = 0.2369268850561890875142640d0
      h5(5) = h5(1)
      h5(2) = 0.4786286704993664680412915d0
      h5(4) = h5(2)
      h5(3) = 0.5688888888888888888888889d0

      i = 0
      do 51 j = 1,5
      do 51 k = 1,5
      i = i + 1
      r(i) = g5(k)
      z(i) = g5(j)
      w(i) = h5(j)*h5(k)
51    continue
      return

c... 6x6 integration
6     g6(1) = 0.9324695142031520278123016d0
      g6(6) = - g6(1)
      g6(2) = 0.6612093864662645136613996d0
      g6(5) = - g6(2)
      g6(3) = 0.2386191860831969086305017d0
      g6(4) = - g6(3)

      h6(1) = 0.1713244923791703450402961d0
      h6(6) = h6(1)
      h6(2) = 0.3607615730481386075698335d0
      h6(5) = h6(2)
      h6(3) = 0.4679139345726910473898703d0
      h6(4) = h6(3)

      i = 0
      do 61 j = 1,6
      do 61 k = 1,6
      i = i + 1
      r(i) = g6(k)
      z(i) = g6(j)
      w(i) = h6(j)*h6(k)
61    continue
      return

c... 7x7 integration
7     g7(1) = 0.9491079123427585245261897d0
      g7(7) = - g7(1)
      g7(2) = 0.7415311855993944398638648d0
      g7(6) = - g7(2)
      g7(3) = 0.4058451513773971669066064d0
      g7(5) = - g7(3)
      g7(4) = 0d0

      h7(1) = 0.1294849661688696932706114d0
      h7(7) = h7(1)
      h7(2) =  0.2797053914892766679014678d0
      h7(6) = h7(2)
      h7(3) = 0.3818300505051189449503698d0
      h7(5) = h7(3)
      h7(4) = 0.4179591836734693877551020d0

      i = 0
      do 71 j = 1,7
      do 71 k = 1,7
      i = i + 1
      r(i) = g7(k)
      z(i) = g7(j)
      w(i) = h7(j)*h7(k)
71    continue
      return

c... 8x8 integration
8     g8(1) = 0.9602898564975362316835609d0
      g8(8) = - g8(1)
      g8(2) = 0.7966664774136267395915539d0
      g8(7) = - g8(2)
      g8(3) = 0.5255324099163289858177390d0
      g8(6) = - g8(3)
      g8(4) = 0.1834346424956498049394761d0
      g8(5) = -g8(4)

      h8(1) = 0.1012285362903762591525314d0
      h8(8) = h8(1)
      h8(2) = 0.2223810344533744705443560d0
      h8(7) = h8(2)
      h8(3) = 0.3137066458778872873379622d0
      h8(6) = h8(3)
      h8(4) = 0.3626837833783619829651504d0
      h8(5) = h8(4)

      i = 0
      do 81 j = 1,8
      do 81 k = 1,8
      i = i + 1
      r(i) = g8(k)
      z(i) = g8(j)
      w(i) = h8(j)*h8(k)
81    continue
      return

c... 9x9 integration
9     g9(1) = 0.9681602395076260898355762d0
      g9(9) = - g9(1)
      g9(2) = 0.8360311073266357942994298d0
      g9(8) = - g9(2)
      g9(3) = 0.6133714327005903973087020d0
      g9(7) = - g9(3)
      g9(4) = 0.3242534234038089290385380d0
      g9(6) = -g9(4)
      g9(5) = 0d0


      h9(1) = 0.0812743883615744119718922d0
      h9(9) = h9(1)
      h9(2) = 0.1806481606948574040584720d0
      h9(8) = h9(2)
      h9(3) = 0.2606106964029354623187429d0
      h9(7) = h9(3)
      h9(4) =  0.3123470770400028400686304d0
      h9(6) = h9(4)
      h9(5) = 0.3302393550012597631645251d0

      i = 0
      do 91 j = 1,9
      do 91 k = 1,9
      i = i + 1
      r(i) = g9(k)
      z(i) = g9(j)
      w(i) = h9(j)*h9(k)
91    continue
      return

c... 10x10 integration
10    g10(1) = 0.9739065285171717200779640d0
      g10(10) = - g10(1)
      g10(2) = 0.8650633666889845107320967d0
      g10(9) = - g10(2)
      g10(3) = 0.6794095682990244062343274d0
      g10(8) = - g10(3)
      g10(4) = 0.4333953941292471907992659d0
      g10(7) = - g10(4)
      g10(5) = 0.1488743389816312108848260d0
      g10(6) = - g10(5)

      h10(1) = 0.0666713443086881375935688d0
      h10(10) = h10(1)
      h10(2) = 0.1494513491505805931457763d0
      h10(9) = h10(2)
      h10(3) = 0.2190863625159820439955349d0
      h10(8) = h10(3)
      h10(4) = 0.2692667193099963550912269d0
      h10(7) = h10(4)
      h10(5) = 0.2955242247147528701738930d0
      h10(6) = h10(5) 

      i = 0
      do 101 j = 1,10
      do 101 k = 1,10
      i = i + 1
      r(i) = g10(k)
      z(i) = g10(j)
      w(i) = h10(j)*h10(k)
101   continue
      return
       
c... 11x11 integration
11    g11(1) = 0.9782286581460569928039380d0
      g11(11) = - g11(1)
      g11(2) = 0.8870625997680952990751578d0
      g11(10) = - g11(2)
      g11(3) = 0.7301520055740493240934163d0
      g11(9) = - g11(3)
      g11(4) = 0.5190961292068118159257257d0
      g11(8) = - g11(4)
      g11(5) = 0.2695431559523449723315320d0
      g11(7) = - g11(5)
      g11(6) = 0d0

      h11(1)  = 0.0556685671161736664827537d0
      h11(11) = h11(1)
      h11(2)  = 0.1255803694649046246346943d0
      h11(10) = h11(2)
      h11(3)  = 0.1862902109277342514260976d0
      h11(9)  = h11(3)
      h11(4)  = 0.2331937645919904799185237d0
      h11(8)  = h11(4)
      h11(5)  = 0.2628045445102466621806889d0   
      h11(7)  = h11(5) 
      h11(6)  = 0.2729250867779006307144835d0         

      i = 0
      do 111 j = 1,11
      do 111 k = 1,11
      i = i + 1
      r(i) = g11(k)
      z(i) = g11(j)
      w(i) = h11(j)*h11(k)
111   continue
      return 

c... 12x12 integration
12    g12(1) = 0.9815606342467192506905491d0
      g12(12) = - g12(1)
      g12(2) = 0.9041172563704748566784659d0
      g12(11) = - g12(2)
      g12(3) = 0.7699026741943046870368938d0
      g12(10) = - g12(3)
      g12(4) = 0.5873179542866174472967024d0
      g12(9) = - g12(4)
      g12(5) = 0.3678314989981801937526915d0
      g12(8) = - g12(5)
      g12(6) = 0.1252334085114689154724414d0
      g12(7) = - g12(6)

      h12(1)  = 0.0471753363865118271946160d0
      h12(12) = h12(1)
      h12(2)  = 0.1069393259953184309602547d0
      h12(11) = h12(2)
      h12(3)  = 0.1600783285433462263346525d0
      h12(10) = h12(3)
      h12(4)  = 0.2031674267230659217490645d0
      h12(9)  = h12(4)
      h12(5)  = 0.2334925365383548087608499d0
      h12(8)  = h12(5) 
      h12(6)  = 0.2491470458134027850005624d0
      h12(7)  = h12(6)      

      i = 0
      do 121 j = 1,12
      do 121 k = 1,12
      i = i + 1
      r(i) = g12(k)
      z(i) = g12(j)
      w(i) = h12(j)*h12(k)
121   continue
      return

c... 13x13 integration
13    g13(1)  = 0.9841830547185881494728294d0
      g13(13) = - g13(1)
      g13(2)  =0.9175983992229779652065478d0
      g13(12) = - g13(2)
      g13(3)  = 0.8015780907333099127942065d0
      g13(11) = - g13(3)
      g13(4)  = 0.6423493394403402206439846d0
      g13(10) =  - g13(4)
      g13(5)  = 0.4484927510364468528779129d0
      g13(9)  =  - g13(5)
      g13(6)  = 0.2304583159551347940655281d0
      g13(8)  =  - g13(6)
      g13(7)  = 0d0
      
      h13(1)  = 0.0404840047653158795200216d0
      h13(13) = h13(1)
      h13(2)  = 0.0921214998377284479144218d0
      h13(12) = h13(2)
      h13(3)  = 0.1388735102197872384636018d0
      h13(11) = h13(3)
      h13(4)  = 0.1781459807619457382800467d0
      h13(10) = h13(4)
      h13(5)  = 0.2078160475368885023125232d0
      h13(9)  = h13(5) 
      h13(6)  = 0.2262831802628972384120902d0
      h13(8)  = h13(6)  
      h13(7)  =  0.2325515532308739101945895d0
    
      i = 0
      do 131 j = 1,13
      do 131 k = 1,13
      i = i + 1
      r(i) = g13(k)
      z(i) = g13(j)
      w(i) = h13(j)*h13(k)
131   continue
      return

c... 14x14 integration
14    g14(1)  = 0.9862838086968123388415973d0
      g14(14) =  - g14(1)
      g14(2)  = 0.9284348836635735173363911d0
      g14(13) =  - g14(2)
      g14(3)  = 0.8272013150697649931897947d0
      g14(12) =  - g14(3)
      g14(4)  = 0.6872929048116854701480198d0
      g14(11) =  - g14(4)
      g14(5)  = 0.5152486363581540919652907d0
      g14(10) =  - g14(5)
      g14(6)  = 0.3191123689278897604356718d0
      g14(9)  =  - g14(6)
      g14(7)  = 0.1080549487073436620662447d0
      g14(8)  =  -g14(7)

      h14(1)  = 0.0351194603317518630318329d0
      h14(14) = h14(1)
      h14(2)  = 0.0801580871597602098056333d0
      h14(13) = h14(2)
      h14(3)  = 0.1215185706879031846894148d0
      h14(12) = h14(3)
      h14(4)  = 0.1572031671581935345696019d0
      h14(11) = h14(4)
      h14(5)  = 0.1855383974779378137417166d0
      h14(10) = h14(5) 
      h14(6)  = 0.2051984637212956039659241d0
      h14(9)  = h14(6)  
      h14(7)  = 0.2152638534631577901958764d0
      h14(8)  = h14(7)
    
      i = 0
      do 141 j = 1,14
      do 141 k = 1,14
      i = i + 1
      r(i) = g14(k)
      z(i) = g14(j)
      w(i) = h14(j)*h14(k)
141   continue
      return

c... 15x15 integration
15    g15(1)  = 0.9879925180204854284895657d0
      g15(15) =  - g15(1)
      g15(2)  = 0.9372733924007059043077589d0
      g15(14) =  - g15(2)
      g15(3)  = 0.8482065834104272162006483d0
      g15(13) =  - g15(3)
      g15(4)  = 0.7244177313601700474161861d0
      g15(12) =  - g15(4)
      g15(5)  = 0.5709721726085388475372267d0
      g15(11) =  - g15(5)
      g15(6)  = 0.3941513470775633698972074d0
      g15(10) =  - g15(6)
      g15(7)  = 0.2011940939974345223006283d0
      g15(9)  =  - g15(7)
      g15(8)  = 0d0

      h15(1)  = 0.0307532419961172683546284d0
      h15(15) = h15(1)
      h15(2)  = 0.0703660474881081247092674d0
      h15(14) = h15(2)
      h15(3)  = 0.1071592204671719350118695d0
      h15(13) = h15(3)
      h15(4)  = 0.1395706779261543144478048d0
      h15(12) = h15(4)
      h15(5)  = 0.1662692058169939335532009d0
      h15(11) = h15(5) 
      h15(6)  = 0.1861610000155622110268006d0
      h15(10) = h15(6)  
      h15(7)  = 0.1984314853271115764561183d0
      h15(9)  = h15(7)
      h15(8)  = 0.2025782419255612728806202d0
    
      i = 0
      do 151 j = 1,15
      do 151 k = 1,15
      i = i + 1
      r(i) = g15(k)
      z(i) = g15(j)
      w(i) = h15(j)*h15(k)
151   continue
      return

c... 16x16 integration
16    g16(1)  = 0.9894009349916499325961542d0
      g16(16) =  - g16(1)
      g16(2)  = 0.9445750230732325760779884d0
      g16(15) =  - g16(2)
      g16(3)  = 0.8656312023878317438804679d0
      g16(14) =  - g16(3)
      g16(4)  = 0.7554044083550030338951012d0
      g16(13) =  - g16(4)
      g16(5)  = 0.6178762444026437484466718d0
      g16(12) =  - g16(5)
      g16(6)  = 0.4580167776572273863424194d0
      g16(11) =  - g16(6)
      g16(7)  = 0.2816035507792589132304605d0
      g16(10) =  - g16(7)
      g16(8)  = 0.0950125098376374401853193d0
      g16(9)  =  - g16(8)
 
      h16(1)  = 0.0271524594117540948517806d0
      h16(16) = h16(1)
      h16(2)  = 0.0622535239386478928628438d0
      h16(15) = h16(2)
      h16(3)  = 0.0951585116824927848099251d0
      h16(14) = h16(3)
      h16(4)  = 0.1246289712555338720524763d0
      h16(13) = h16(4)
      h16(5)  = 0.1495959888165767320815017d0
      h16(12) = h16(5) 
      h16(6)  = 0.1691565193950025381893121d0
      h16(11) = h16(6)  
      h16(7)  = 0.1826034150449235888667637d0
      h16(10)  = h16(7)
      h16(8)  = 0.1894506104550684962853967d0
      h16(9)  = h16(8)
         
      i = 0
      do 161 j = 1,16
      do 161 k = 1,16
      i = i + 1
      r(i) = g16(k)
      z(i) = g16(j)
      w(i) = h16(j)*h16(k)
161   continue
      return

c... 17x17 integration
17    g17(1)  = 0.9905754753144173356754340d0
      g17(17) =  - g17(1)
      g17(2)  = 0.9506755217687677612227170d0
      g17(16) =  - g17(2)
      g17(3)  = 0.8802391537269859021229557d0
      g17(15) =  - g17(3)
      g17(4)  = 0.7815140038968014069252301d0
      g17(14) =  - g17(4)
      g17(5)  = 0.6576711592166907658503022d0
      g17(13) =  - g17(5)
      g17(6)  = 0.5126905370864769678862466d0
      g17(12) =  - g17(6)
      g17(7)  = 0.3512317634538763152971855d0
      g17(11) =  - g17(7)
      g17(8)  = 0.1784841814958478558506775d0
      g17(10) =  - g17(8)
      g17(9)  = 0d0
 
      h17(1)  = 0.0241483028685479319601100d0
      h17(17) = h17(1)
      h17(2)  = 0.0554595293739872011294402d0
      h17(16) = h17(2)
      h17(3)  = 0.0850361483171791808835354d0
      h17(15) = h17(3)
      h17(4)  = 0.1118838471934039710947884d0
      h17(14) = h17(4)
      h17(5)  = 0.1351363684685254732863200d0
      h17(13) = h17(5) 
      h17(6)  = 0.1540457610768102880814316d0
      h17(12) = h17(6)  
      h17(7)  = 0.1680041021564500445099707d0
      h17(11) = h17(7)
      h17(8)  = 0.1765627053669926463252710d0
      h17(10) = h17(8)
      h17(9)  = 0.1794464703562065254582656d0
         
      i = 0
      do 171 j = 1,17
      do 171 k = 1,17
      i = i + 1
      r(i) = g17(k)
      z(i) = g17(j)
      w(i) = h17(j)*h17(k)
171   continue
      return

c... 18x18 integration
18    g18(1)  = 0.9915651684209309467300160d0
      g18(18) =  - g17(1)
      g18(2)  = 0.9558239495713977551811959d0
      g18(17) =  - g17(2)
      g18(3)  = 0.8926024664975557392060606d0
      g18(16) =  - g17(3)
      g18(4)  = 0.8037049589725231156824175d0
      g18(15) =  - g17(4)
      g18(5)  = 0.6916870430603532078748911d0
      g18(14) =  - g17(5)
      g18(6)  = 0.5597708310739475346078715d0
      g18(13) =  - g17(6)
      g18(7)  = 0.4117511614628426460359318d0
      g18(12) =  - g17(7)
      g18(8)  = 0.2518862256915055095889729d0
      g18(11) =  - g17(8)
      g18(9)  = 0.0847750130417353012422619d0
      g18(10) =  - g18(9)
 
      h18(1)  = 0.0216160135264833103133427d0
      h18(18) = h18(1)
      h18(2)  = 0.0497145488949697964533349d0
      h18(17) = h18(2)
      h18(3)  = 0.0764257302548890565291297d0
      h18(16) = h18(3)
      h18(4)  = 0.1009420441062871655628140d0
      h18(15) = h18(4)
      h18(5)  = 0.1225552067114784601845191d0
      h18(14) = h18(5) 
      h18(6)  = 0.1406429146706506512047313d0
      h18(13) = h18(6)  
      h18(7)  = 0.1546846751262652449254180d0
      h18(12) = h18(7)
      h18(8)  = 0.1642764837458327229860538d0
      h18(11) = h18(8)
      h18(9)  = 0.1691423829631435918406565d0
      h18(10) = h18(9)
         
      i = 0
      do 181 j = 1,18
      do 181 k = 1,18
      i = i + 1
      r(i) = g18(k)
      z(i) = g18(j)
      w(i) = h18(j)*h18(k)
181   continue
      return

c... 19x19 integration
19    g19(1)  = 0.9924068438435844031890177d0
      g19(19) =  - g17(1)
      g19(2)  = 0.9602081521348300308527788d0
      g19(18) =  - g17(2)
      g19(3)  = 0.9031559036148179016426609d0
      g19(17) =  - g17(3)
      g19(4)  = 0.8227146565371428249789225d0
      g19(16) =  - g17(4)
      g19(5)  = 0.7209661773352293786170959d0
      g19(15) =  - g17(5)
      g19(6)  = 0.6005453046616810234696382d0
      g19(14) =  - g17(6)
      g19(7)  = 0.4645707413759609457172671d0
      g19(13) =  - g17(7)
      g19(8)  = 0.3165640999636298319901173d0
      g19(12) =  - g17(8)
      g19(9)  = 0.1603586456402253758680961d0
      g19(11) =  - g18(9)
      g19(10) = 0.d0
 
      h19(1)  = 0.0194617882297264770363120d0
      h19(19) = h19(1)
      h19(2)  = 0.0448142267656996003328382d0
      h19(18) = h19(2)
      h19(3)  = 0.0690445427376412265807083d0
      h19(17) = h19(3)
      h19(4)  = 0.0914900216224499994644621d0
      h19(16) = h19(4)
      h19(5)  = 0.1115666455473339947160239d0
      h19(15) = h19(5) 
      h19(6)  = 0.1287539625393362276755158d0
      h19(14) = h19(6)  
      h19(7)  = 0.1426067021736066117757461d0
      h19(13) = h19(7)
      h19(8)  = 0.1527660420658596667788554d0
      h19(12) = h19(8)
      h19(9)  = 0.1589688433939543476499564d0
      h19(11) = h19(9)
      h19(10) = 0.1610544498487836959791636d0
         
      i = 0
      do 191 j = 1,19
      do 191 k = 1,19
      i = i + 1
      r(i) = g19(k)
      z(i) = g19(j)
      w(i) = h19(j)*h19(k)
191   continue
      return

c... 20x20 integration
20    g20(1)  = 0.9931285991850949247861224d0
      g20(20) =  - g20(1)
      g20(2)  = 0.9639719272779137912676661d0
      g20(19) =  - g20(2)
      g20(3)  = 0.9122344282513259058677524d0
      g20(18) =  - g20(3)
      g20(4)  = 0.8391169718222188233945291d0
      g20(17) =  - g20(4)
      g20(5)  = 0.7463319064601507926143051d0
      g20(16) =  - g20(5)
      g20(6)  = 0.6360536807265150254528367d0
      g20(15) =  - g20(6)
      g20(7)  = 0.5108670019508270980043641d0
      g20(14) =  - g20(7)
      g20(8)  = 0.3737060887154195606725482d0
      g20(13) =  - g20(8)
      g20(9)  = 0.2277858511416450780804962d0
      g20(12) =  - g20(9)
      g20(10) = 0.0765265211334973337546404d0
      g20(11) = - g20(10)

      h20(1)  = 0.0176140071391521183118620d0
      h20(20) = h20(1)
      h20(2)  = 0.0406014298003869413310400d0
      h20(19) = h20(2)
      h20(3)  = 0.0626720483341090635695065d0
      h20(18) = h20(3)
      h20(4)  = 0.0832767415767047487247581d0
      h20(17) = h20(4)
      h20(5)  = 0.1019301198172404350367501d0
      h20(16) = h20(5) 
      h20(6)  = 0.1181945319615184173123774d0
      h20(15) = h20(6)  
      h20(7)  = 0.1316886384491766268984945d0
      h20(14) = h20(7)
      h20(8)  = 0.1420961093183820513292983d0
      h20(13) = h20(8)
      h20(9)  = 0.1491729864726037467878287d0
      h20(12) = h20(9)
      h20(10) = 0.1527533871307258506980843d0 
      h20(11) = h20(10)
         
      i = 0
      do 201 j = 1,20
      do 201 k = 1,20
      i = i + 1
      r(i) = g20(k)
      z(i) = g20(j)
      w(i) = h20(j)*h20(k)
201   continue
      return

      end
c
      subroutine plobatto(ngp,pt,wt)
c----------------------------------------------------------------------
c
c      Purpose: Form Lobatto-points and weights for 1D
c               1-10 points
c
c      Inputs:
c         ngp     - Number of points
c
c      Outputs:
c        pt(*)    - 1-direction Lobatto point
c        wt(*)    - asscoiated  weight
c
c      Comments: Lobatto = Newton-Cotes for ngp<=3
c
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension  pt(*), wt(*)
      if ( ngp .lt. 1 ) stop 'no points in plobatto'
      if ( ngp .eq. 1 )  then
        pt( 1) = 0.00000000000000d+00
        wt( 1) = 2.00000000000000d+00
        return
      elseif ( ngp .eq. 2 )  then
        pt( 1) = 1.00000000000000d+00
        wt( 1) = 1.00000000000000d+00
        pt( 2) = -pt( 1)
        wt( 2) =  wt( 1)
        return
      elseif ( ngp .eq. 3 )  then
        pt( 1) = 1.00000000000000d+00
        wt( 1) = 3.33333333333333d-01
        pt( 2) = 0.00000000000000d+00
        wt( 2) = 1.33333333333333d+00
        pt( 3) = -pt( 1)
        wt( 3) =  wt( 1)
        return
      elseif ( ngp .eq. 4 )  then
        pt( 1) = 1.00000000000000d+00
        wt( 1) = 1.66666666666667d-01
        pt( 2) = 4.47213595499958d-01
        wt( 2) = 8.33333333333333d-01
        pt( 3) = -pt( 1)
        wt( 3) =  wt( 1)
        pt( 4) = -pt( 2)
        wt( 4) =  wt( 2)
        return
      elseif ( ngp .eq. 5 )  then
        pt( 1) = 1.00000000000000d+00
        wt( 1) = 1.00000000000000d-01
        pt( 2) = 6.54653670707977d-01
        wt( 2) = 5.44444444444444d-01
        pt( 3) = 0.00000000000000d+00
        wt( 3) = 7.11111111111111d-01
        pt( 4) = -pt( 1)
        wt( 4) =  wt( 1)
        pt( 5) = -pt( 2)
        wt( 5) =  wt( 2)
        return
      elseif ( ngp .eq. 6 )  then
        pt( 1) = 1.00000000000000d+00
        wt( 1) = 6.66666666666667d-02
        pt( 2) = 7.65055323929465d-01
        wt( 2) = 3.78474956297847d-01
        pt( 3) = 2.85231516480645d-01
        wt( 3) = 5.54858377035486d-01
        pt( 4) = -pt( 1)
        wt( 4) =  wt( 1)
        pt( 5) = -pt( 2)
        wt( 5) =  wt( 2)
        pt( 6) = -pt( 3)
        wt( 6) =  wt( 3)
      elseif ( ngp .eq. 7 )  then
        pt( 1) = 1.00000000000000d+00
        wt( 1) = 4.7619047619048d-02
        pt( 2) = 8.30223896278567d-01
        wt( 2) = 2.76826047361566d-01
        pt( 3) = 4.68848793470714d-01
        wt( 3) = 4.31745381209863d-01
        pt( 4) = 0.00000000000000d+00
        wt( 4) = 4.87619047619048d-01
        pt( 5) = -pt( 1)
        wt( 5) =  wt( 1)
        pt( 6) = -pt( 2)
        wt( 6) =  wt( 2)
        pt( 7) = -pt( 3)
        wt( 7) =  wt( 3)
        return
      elseif ( ngp .eq. 8 )  then
        pt( 1) = 1.00000000000000d+00
        wt( 1) = 3.5714285714286d-02
        pt( 2) = 8.71740148509607d-01
        wt( 2) = 2.10704227143506d-01
        pt( 3) = 5.91700181433142d-01
        wt( 3) = 3.41122692483504d-01
        pt( 4) = 2.09299217902479d-01
        wt( 4) = 4.12458794658704d-01
        pt( 5) = -pt( 1)
        wt( 5) =  wt( 1)
        pt( 6) = -pt( 2)
        wt( 6) =  wt( 2)
        pt( 7) = -pt( 3)
        wt( 7) =  wt( 3)
        pt( 8) = -pt( 4)
        wt( 8) =  wt( 4)
        return
      elseif ( ngp .eq. 9 )  then
        pt( 1) = 1.00000000000000d+00
        wt( 1) = 2.77777777777778d-02
        pt( 2) = 8.99757995411460d-01
        wt( 2) = 1.65495361560805d-01
        pt( 3) = 6.77186279510738d-01
        wt( 3) = 2.74538712500162d-01
        pt( 4) = 3.63117463826178d-01
        wt( 4) = 3.46428510973406d-01
        pt( 5) = 0.00000000000000d+00
        wt( 5) = 3.71519274376417d-01
        pt( 6) = -pt( 1)
        wt( 6) =  wt( 1)
        pt( 7) = -pt( 2)
        wt( 7) =  wt( 2)
        pt( 8) = -pt( 3)
        wt( 8) =  wt( 3)
        pt( 9) = -pt( 4)
        wt( 9) =  wt( 4)
        return
      elseif ( ngp .eq. 10 )  then
        pt( 1) = 1.00000000000000d+00
        wt( 1) = 2.22222222222222d-02
        pt( 2) = 9.19533908166459d-01
        wt( 2) = 1.33305990851070d-01
        pt( 3) = 7.38773865105505d-01
        wt( 3) = 2.24889342063126d-01
        pt( 4) = 4.77924949810444d-01
        wt( 4) = 2.92042683679684d-01
        pt( 5) = 1.65278957666387d-01
        wt( 5) = 3.27539761183897d-01
        pt( 6) = -pt( 1)
        wt( 6) =  wt( 1)
        pt( 7) = -pt( 2)
        wt( 7) =  wt( 2)
        pt( 8) = -pt( 3)
        wt( 8) =  wt( 3)
        pt( 9) = -pt( 4)
        wt( 9) =  wt( 4)
        pt(10) = -pt( 5)
        wt(10) =  wt( 5)
        return
      endif
      return
      end
c
      subroutine pnewcot(l,lint,r,s,w)
c----------------------------------------------------------------------
c
c      Purpose: Form Newton-Cotes points and weights for two dimensions
c               1-8 points
c
c      Inputs:
c         l       - Number of points/direction
c
c      Outputs:
c         lint    - Total number of points
c         r(*)    - 1-direction Newton-Cotes point
c         s(*)    - 2-direction Newton-Cotes point
c         w(*)    - Newton-Cotes weight
c
c        l=1 p=   [0]                    w=     [2]    
c        l=2 p=   [-1  1]                w=     [1  1]                                                      Trapez
c        l=3 p=   [-1  0  1]             w=1/3  [1  4   1]                                                  Simpson
c        l=4 p=1/3[-3 -1  1 3]           w=1/4  [1  3   3 1]                                                Pulcherima
c        l=5 p=1/2[-2 -1  0 1 2]         w=1/45 [7  32  12 32 7]                                            Milne
c        l=6 p=1/5[-5 -3 -1 1 3 5]       w=1/144[19 75  50 50 75 19]                                        NC6
c        l=7 p=1/6[-6 -4 -2 0 2 4 6]     w=1/420[41 216 27 272 27 216 41]                                   Weddle
c        l=8 p=1/7[-7 -5 -3 -1 1 3 5 7]  w=1/320[1502 7154 2646 5978 5978 2646 7154 1502]                   NC8
c
c        For Newton Cotes > 8 we have some negative weights, for that reason the results are inferior. TH KIT 06/12
c
c      Position e.g.l=4
c          4  8 12 16    
c          3  7 11 15
c          2  6 10 14
c          1  5  9 13
c
c      WW KIT 10/11 
c
c----------------------------------------------------------------------
      USE eldata
      USE iofile
      implicit double precision (a-h,o-z)
      dimension lp1(1),lw1(1),
     +          lp2(2),lw2(2),
     +          lp3(3),lw3(3),
     +          lp4(4),lw4(4),
     +          lp5(5),lw5(5),
     +          lp6(6),lw6(6),
     +          lp7(7),lw7(7),
     +          lp8(8),lw8(8),    
     +          r(*),s(*),w(*)

      data lp1/0/,                   lw1/2/,
     +     lp2/-1,1/,                lw2/1,1/,
     +     lp3/-1,0,1/,              lw3/1,4,1/,
     +     lp4/-3,-1,1,3/,           lw4/1,3,3,1/,
     +     lp5/-2,-1,0,1,2/,         lw5/7,32,12,32,7/,
     +     lp6/-5,-3,-1,1,3,5/,      lw6/19,75,50,50,75,19/,
     +     lp7/-6,-4,-2,0,2,4,6/,    lw7/41,216,27,272,27,216,41/
     +     lp8/-7,-5,-3,-1,1,3,5,7/ 
     +     lw8/1502,7154,2646,5978,5978,2646,7154,1502/


      lint = l*l
      go to (1,2,3,4,5,6,7,8),l

c.... 1x1 integration
1     i = 0   
      do 11 j = 1,l
       do 11 k = 1,l
         i=i+1 
         r(i) = lp1(j)
         s(i) = lp1(k)
11       w(i) = lw1(j)*lw1(k)
      goto 10

c.... 2x2 integration
2     i = 0   
      do 21 j = 1,l
       do 21 k = 1,l
         i=i+1 
         r(i) = lp2(j)
         s(i) = lp2(k)
21       w(i) = lw2(j)*lw2(k)
      goto 10

c.... 3x3 integration
3     h = 1.d0/9.d0
      i = 0 
      do 31 j = 1,l
       do 31 k = 1,l
         i=i+1 
         r(i) = lp3(j)
         s(i) = lp3(k)
31       w(i) = lw3(j)*lw3(k)*h
      goto 10

c.... 4x4 integration
4     h = 1.d0/16.d0
      g = 1.d0/3.d0
      i = 0 
      do 41 j = 1,l
       do 41 k = 1,l
         i=i+1 
         r(i) = lp4(j)*g
         s(i) = lp4(k)*g
41       w(i) = lw4(j)*lw4(k)*h
      goto 10

c.... 5x5 integration
5     h = 1.d0/2025.d0
      g = 1.d0/2.d0
      i = 0 
      do 51 j = 1,l
       do 51 k = 1,l
         i=i+1 
         r(i) = lp5(j)*g
         s(i) = lp5(k)*g
51       w(i) = lw5(j)*lw5(k)*h
      goto 10

c.... 6x6 integration
6     h = 1.d0/20736.d0
      g = 1.d0/5.d0
      i = 0 
      do 61 j = 1,l
       do 61 k = 1,l
         i=i+1 
         r(i) = lp6(j)*g
         s(i) = lp6(k)*g
61       w(i) = lw6(j)*lw6(k)*h
      goto 10

c.... 7x7 integration
7     h = 1.d0/176400.d0
      g = 1.d0/6.d0
      i = 0 
      do 71 j = 1,l
       do 71 k = 1,l
         i=i+1 
         r(i) = lp7(j)*g
         s(i) = lp7(k)*g
71       w(i) = lw7(j)*lw7(k)*h
      goto 10

c.... 8x8 integration
8     h = 1.d0/298598400.d0
      g = 1.d0/7.d0
      i = 0 
      do 81 j = 1,l
       do 81 k = 1,l
         i=i+1 
         r(i) = lp8(j)*g
         s(i) = lp8(k)*g
81       w(i) = lw8(j)*lw8(k)*h
      goto 10
c.... print
10    continue

c      write(  *,*) 'L',l
c      write(iow,*) 'L',l
c      ws = 0.d0  
c      do i=1,lint
c        write(  *,'(i5,f12.5,f12.5,f12.5)') i,r(i),s(i),w(i)
c        write(iow,'(i5,f12.5,f12.5,f12.5)') i,r(i),s(i),w(i)
c        ws =ws+w(i)
c      end do
c      write(  *,'(a7,i5,f12.5)') 'Lint,WS',lint,ws
c      write(iow,'(a7,i5,f12.5)') 'Lint,WS',lint,ws
 
      return
      end
c
      subroutine phist(wd,clab2,jct,lct,ct,ll,is)
c----------------------------------------------------------------------
c
c      Purpose: Controls history lists for command language execution
c
c      Inputs:
c         wd(*)    - List of command language options
c         clab2    - option for history request
c         ct(3,*)  - Numerical value of command option request

c      Outputs:
c         jct(*)   - Stores list of commands for history
c         lct(*)   - Stores list of options to commands for history
c         ll       - Number of active comands
c         is       - Output: = 0 for command execution, otherwise = 1
c----------------------------------------------------------------------
      USE idata
      USE iodata
      implicit double precision (a-h,o-z)
      character*4 wd,clab2,lct,y*1
      dimension ct(3,*),jct(*),lct(*),wd(*)
      is = 1
      nl1 = abs(ct(1,ll))
      nl2 = abs(ct(2,ll))
      if  (clab2.eq.'read') then
        open(ios,file='feap.his',status='unknown')
        do 100 n = 1,100
          read(ios,1000,end=110) js(n),ljs(n),(vjs(i,n),i=1,3)
100     continue
110     nn = n - 1
        close(ios)
      else if (clab2.eq.'save') then
        open(ios,file='feap.his',status='unknown')
        rewind ios
        do 150 n = 1,nn
          write(ios,1000) js(n),ljs(n),(vjs(i,n),i=1,3)
150     continue
        close(ios)
      else if (clab2.eq.'list'.or.
     1       (clab2.eq.'    '.and.nl1+nl2.eq.0)) then
        nl1 = min(nn,max(1,nl1))
        if(nl1.gt.0) then
          nl2 = max(nl1,min(nn,nl2))
          if(ct(1,ll).eq.0.d0.and.ct(2,ll).eq.0.d0) nl2 = nn
          write(*,2000)
          do 200 n = nl1,nl2
            write(*,2001) n,wd(js(n)),ljs(n),(vjs(i,n),i=1,3)
200       continue
        else
          write(*,3001) nn
        end if
      else if (clab2.eq.'edit') then
        write(*,2002) nl1,wd(js(nl1)),ljs(nl1),(vjs(i,nl1),i=1,3)
        read(*,1001) y
        if(y.eq.'y' .or. y.eq.'Y') then
          if(nl1.gt.0.and.nl1.le.nn) then
            nl2 = max(nl1,min(nn,nl2))
            nl1 = nl2 - nl1 + 1
            if(nl2.lt.nn) then
              do 250 n = nl2+1,nn
                js(n - nl1) = js(n)
                ljs(n - nl1) = ljs(n)
                vjs(1,n-nl1) = vjs(1,n)
                vjs(2,n-nl1) = vjs(2,n)
                vjs(3,n-nl1) = vjs(3,n)
250           continue
            end if
            nn = nn - nl1
          else
            write(*,2005) nn
          end if
        end if
      else if (clab2.eq.'add ') then
        hadd = .true.
      else if (clab2.eq.'noad') then
        hadd = .false.
      else
        if(nl1.gt.0.and.nl1.le.nn) then
          nl2 = max(nl1,min(nn,nl2))
          do 300 n = nl1,nl2
            jct(ll) = js(n)
            lct(ll) = ljs(n)
            ct(1,ll) = vjs(1,n)
            ct(2,ll) = vjs(2,n)
            ct(3,ll) = vjs(3,n)
            ll = ll + 1
300       continue
          is = 0
        else
          write(*,3002) nn
        end if
      end if
      return
1000  format(i5,1x,a4,3e15.5)
1001  format(a1)
2000  format('  no. macro option   value-1      value-2',
     1    '      value-3')
2001  format(i5,2(2x,a4),3g13.4)
2002  format('  Remove command:',i5,2(2x,a4),3g13.4/3x,'(y or n)? >',$)
2005  format(' ** ERROR ** Not that many items in list, nn =',i4)
3001  format(' ** ERROR ** No items in list, nn = ',i3)
3002  format(5x,'Currently history list contains',i4,' items.'/
     1     7x,'Options:',/
     1    16x,'hist,list,n1,n2 - list items n1 to n2'/
     2    16x,'hist            - list all items'/
     3    16x,'hist,,n1,n2     - execute items n1 to n2'/
     4    16x,'hist,,n1        - execute item n1'/
     5    16x,'hist,edit,n1    - remove item n1'/
     6    16x,'hist,add        - add new macros to list'/
     7    16x,'hist,noad       - do not add macros to list'/
     8    16x,'hist,save       - save current list on   disk'/
     9    16x,'hist,read       - read current list from disk'/
     x    34x,'(file = feap_his)'/1x)
      end
c
      subroutine piden(d,ns,ne)
c----------------------------------------------------------------------
c      Purpose: Sets array to identity (1.0)

c      Inputs:
c         ns      - First entry to set
c         ne      - Last  entry to set

c      Outputs:

c         d(*)    - Array set to unity
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension d(*)
      do 100 n = ns,ne
100   d(n) = 1.0d0
      return
      end
c
      subroutine pintio(y,ifld)
c----------------------------------------------------------------------
c      Purpose: Character string input routine for data
c               N.B. This routine has largely been superceded by
c                    dinput functions.

c      Inputs:
c         lfld   - Field width to separate input data items into.

c      Outputs:
c         y      - Character string input, in field widths of ifld
c
c----------------------------------------------------------------------
      USE iofile
      USE iosave
      CHARACTER*80 X,Y
      CHARACTER*1 XX(80),YY(80)
      EQUIVALENCE (X,XX)
c.... read a record from the input file
100   if (ior.gt.0)  then
        read (ior,1000,err=901,end=902) x
      else
        read (  *,1000,err=901,end=902) x
      end if
      if(lsave) write(lfile,1000) x
c.... adjust and move the record to the y-array
c.kneb       call acheck(x,y,ifld,80,80)
      call acheck(xx,yy,ifld,80,80)
      DO 50 I = 1,80
  50  Y(I:I) = YY(I)
      return
c.... read error encountered
901   call  errclr ('PINTIO')
      stop 'ERROR in PINTIO at data input' !ww 
      goto  100
c.... eof encountered
c.kne 902   call  endclr ('PINTIO',x(1))
902   call  endclr ('PINTIO',x)
      goto  100
c.kne 1000  format(80a1)
1000  format(a)
      end
c
      subroutine plink(idl,id,x,ndm,ndf,numnp,neq,iplk1,iplk2,prt)
c----------------------------------------------------------------------
c
c      Purpose: Link degrees of freedom to have same solution value
c
c      Inputs:
c         id(ndf,numnp)   - equation numbers for each active dof
c         x(ndm,numnp)    - Nodal coordinates of mesh
c         ndm             - Spatial dimension of mesh
c         ndf             - Number dof/node
c         numnp           - Number of nodes in mesh
c         prt             - Output links performed if true
c
c      Scratch:
c         idl(ndf*ndf)    - linking condition  
c
c      Outputs:
c         id(ndf,numnp)   - Equation numbers for each dof after link
c         neq             - Number of equations active after link
c         iplk1(numnp)    - List of linked nodes to plot
c         iplk2(ndf,numnp)- List of bc for linked nodes to plot
c         with
c         iplk1(i)        - = -i  master node 
c                             i=k slave node, master node is k 
c                             0=not linked
c         iplk2(j,i)      - = dof j:  
c                             j=0:  not linked
c                             j=k>0 linked, master node is node k 
c
c      Comments
c      # Tie of linked nodes should be ok
c      # Nodes with B.C.s are not linked
c  
c      Types:
c
c      type 1: 1, na, ne,inc,code................ link nodes na-ne       
c      type 2: 2, na,idm,cor,code,eps............ link nodes on a line
c      type 3: 3,na1,ne1,inc1,na2,ne2,inc2,code.. link node ni1 with ni2 
c      type 4: 4,idm,x1,x2,code,eps... link nodes on x2 with nodes on x1  
c      type 5: 5,idmx,x1,x2,idmy,y1,y2,code,eps... link nodes on x1,x2 with nodes on x2,y2  
c      type 6: 6,idm,x1,x2,code,eps... link nodes on x2 with nodes on x1 with y(x2)=-y(x1)   
c      type6=type4, only coordinate check is different (minus -> plus) 
c      type 7: 7,lx,ly,lz  .... link nodes for z=const on border of area in x-y dir

c
c      type 3 added ww ibs uka  2/07
c      type 4 added ww ibs kit 12/10
c      type 5 added ww ibs kit 01/11
c      type 6 added ww ibs kit 01/12
c      type 7 added ww ibs kit 03/13
c      control output added ww ibs kit 10/12
c----------------------------------------------------------------------
      USE comfil
      USE errchk
      USE iodata
      USE iofile
      USE pnodn
      implicit double precision(a-h,o-z)

      logical prt,lsave
      character*229 flink,yyy*80,yyy1*100
      integer id(ndf,numnp),idl(ndf*ndf),iplk1(numnp),iplk2(ndf,numnp)
      real*8  td(16),x(ndm,*)


c      integer ndm, ndf, numnp, neq, 
c     1        isfile, i, ii, inc, j, jj, n1, n2, ni, nelim, nmax,
c     2        na1,ne1,na2,ne2,inc1,inc2,iadd, ir 
c
c.... set intrinsics
c      intrinsic min, max
c
      save   isfile
      flink = finp
      call addext(flink,'lnk ')
      inquire(file=flink,exist=lsave)
      if(lsave) then
        open(ios,file=flink,status='old')
      else
        write(yyy1,2010) flink(1:55)
        call drawmess(yyy1,1,0)
      end if
      if(lsave) then
        isfile = ior
        ior    = ios

        nlk=0
c
10      continue
c....   input the element records - N.B. limit is 16 nos. / record
        call pzero(td,16)
        call dinput(td,16)
        ityp = td(1)
c....   position of link code
        if(ityp.le.2) iadd=4 
        if(ityp.eq.3) iadd=7 
        if(ityp.eq.4) iadd=4 
        if(ityp.eq.5) iadd=7 
        if(ityp.eq.6) iadd=4 
        if(ityp.eq.7) iadd=4 

        il = min(ndf+iadd,16)
c...... input finished 
        if(ityp.eq.0) then
          close(ios)
          ior = isfile

          if(ior.gt.0) then
            write(iow,2008) 
            do i = 1,numnp
              write(iow,2009) i,(iplk2(j,i),j=1,ndf),(id(j,i),j=1,ndf)         
            end do
          else
            write(  *,2008) 
            do i = 1,numnp
              write(  *,2009) i,(iplk2(j,i),j=1,ndf),(id(j,i),j=1,ndf)         
            end do
          end if
          return
        end if

        nlk=nlk+1
        
        if(ityp.eq.1) then
c......   ityp = 1 nodes from n1 to n2 with inc
          n1  = td(2)
          n2  = td(3)
          n1q = iprttie(gtie,n1)
          n2q = iprttie(gtie,n2)
          if(n1.gt.numnp.or.n1.lt.1) then
             write(iow,3001) n1
c             write(  *,3001) n1
          end if
          if(n2.gt.numnp.or.n2.lt.1) then
             if(prt) write(iow,3001) n2
          end if
          inc  = td(4)             
          if(inc.eq.0) inc = n2 - n1
          ni   = n1 + inc            ! first node to link

        else if(ityp.eq.2) then
c......   ityp = 2 nodes  on coor idm  to n1
          n1   = td(2)              
          idm  = td(3)              
          n1q = iprttie(gtie,n1)
          if(n1.gt.numnp.or.n1.lt.1) then
             write(iow,3001) n1
c             write(  *,3001) n1
          end if
          idm  = min(idm,ndm)
          x0   = td(4)               ! coor
          dx   = pdiff(x(idm,1),ndm,numnp)/1000.d0
          if(ndf.gt.11) stop 'eps in PLINK1 out of TD'
          eps  = td(5+ndf) ! geht nur fuer eine Zeile!
          if(eps.ne.0.d0) dx = eps  

        else if(ityp.eq.3) then
c......   ityp = 3 nodes from na2 to ne2 with inc2  to  na1 to ne1 with inc1
          na1  = td(2)    ! do not change order!!
          ne1  = td(3)   
          na1q = iprttie(gtie,na1)
          ne1q = iprttie(gtie,ne1)
          if(na1.gt.numnp.or.na1.lt.1) then
             write(iow,3001) na1
c             write(  *,3001) na1
          end if
          if(ne1.gt.numnp.or.ne1.lt.1) then
             write(iow,3001) na1
c             write(  *,3001) na1
          end if
          inc1 = td(4)              
          na2  = td(5)    
          ne2  = td(6)    
          na2q = iprttie(gtie,na2)
          ne2q = iprttie(gtie,ne2)
          if(na2.gt.numnp.or.na2.lt.1) then
             write(iow,3001) na1
c             write(  *,3001) na1
         end if
          if(ne2.gt.numnp.or.ne2.lt.1) then
             write(iow,3001) na1
c             write(  *,3001) na1
          end if
          inc2 = td(7)               
          ni   = n1 + inc   

          if(inc1.eq.0) inc1 = ne1 - na1
          if(inc2.eq.0) inc2 = ne2 - na2
          ni    = na1

        else if(ityp.eq.4) then
c......   ityp = 4 nodes on coor(idm)=x2 to nodes on coor(idm)=x1
          idm  = td(2)              
          x1   = td(3)
          x2   = td(4)
          dx   = pdiff(x(idm,1),ndm,numnp)/1000.d0
          eps  = td(5+ndf) ! geht nur fuer eine Zeile!
          if(eps.ne.0.d0) dx = eps  

        else if(ityp.eq.5) then
c......   ityp = 5 nodes on coor(idmx,idmy)=x2,y2 to nodes on coor(x1,y1)
          idmx  = td(2)              
          x1    = td(3)
          x2    = td(4)
          idmy  = td(5)              
          y1    = td(6)
          y2    = td(7)
          dx   = pdiff(x(idmx,1),ndm,numnp)/1000.d0
          dy   = pdiff(x(idmy,1),ndm,numnp)/1000.d0
          dx = min(dx,dy)
          eps  = td(8+ndf) ! geht nur fuer eine Zeile!
          if(eps.ne.0.d0) then  dx = eps  

        else if(ityp.eq.6) then
c......   ityp = 6 nodes on coor(idm)=x2 to nodes on coor(idm)=x1 with y(x2)=-y(x1)
          idm  = td(2)              
          x1   = td(3)
          x2   = td(4)
          dx   = pdiff(x(idm,1),ndm,numnp)/1000.d0
          eps  = td(5+ndf) ! geht nur fuer eine Zeile!
          if(eps.ne.0.d0) dx = eps  

        else if(ityp.eq.7) then
c......   ityp = 7 for x3=c all nodes on boundary in x1,x2 are linked 
c                  to one node [-x1/2,-x2/2,c]
          x1   = td(2)              
          x2   = td(3)
          x3   = td(4)
          dx1  = pdiff(x(1,1),ndm,numnp)/1000.d0
          dx2  = pdiff(x(1,1),ndm,numnp)/1000.d0
          dx3  = pdiff(x(1,1),ndm,numnp)/1000.d0
          x1m=-0.5d0*x1  
          x2m=-0.5d0*x2  
          x3m=-0.5d0*x3  
          x1p= 0.5d0*x1  
          x2p= 0.5d0*x2  
          x3p= 0.5d0*x3  
        end if

c...... input code
        do i = 1,min(ndf,16-iadd)
          idl(i) = td(i+iadd)
        end do    
        if(ndf.gt.16-iadd) then        ! more than one line input
          do 105 ii = 1,(ndf+iadd)/16
            is = il+1
            il = min(is+15,ndf+iadd)
103         call dinput(td,il-is+1)
            if(errck) go to 103
            do 104 k = 1,il-is+1
              idl(k+is-iadd-1) = td(k) ! code
104         continue
105       continue
        end if

c...    protocol
        if(ityp.eq.1) then 
c.....    protocol link node1 - node2 to node1
          if(ior.gt.0) then
            if(prt) write(iow,2001) n1,n1q,n2,n2q,inc,(idl(i),i=1,ndf)
            if(prt) write(iow,4000) n1,n1q
          else
            if(prt) write(*  ,2001) n1,n1q,n2,n2q,inc,(idl(i),i=1,ndf)
            if(prt) write(*  ,4000) n1,n1q
          end if

        else if(ityp.eq.2) then
c.....    protocol link nodes on x1 to node1
          if(ior.gt.0) then
            if(prt) write(iow,2002) idm,x0,n1,n1q,(idl(i),i=1,ndf)
            if(prt) write(iow,4000) n1,n1q
          else
            if(prt) write(*  ,2002) idm,x0,n1,n1q,(idl(i),i=1,ndf)
            if(prt) write(*  ,4000) n1,n1q
          end if

        else if(ityp.eq.3) then
c.....    protocol link (node1-node2)_i to (node1-node2)_j
          if(ior.gt.0) then
            if(prt) write(iow,2003) na2,na2q,ne2,ne2q,inc2,
     +                     na1,na1q,ne1,ne1q,inc1,(idl(i),i=1,ndf)
          else
            if(prt) write(*  ,2003) na2,na2q,ne2,ne2q,inc2,
     +                     na1,na1q,ne1,ne1q,inc1,(idl(i),i=1,ndf)
          end if

        else if(ityp.eq.4) then
c.....    protocol link nodes_x2 to nodes_x1
          if(ior.gt.0) then
            if(prt) write(iow,2004) idm,x2,x1,(idl(i),i=1,ndf)
          else
            if(prt) write(  *,2004) idm,x2,x1,(idl(i),i=1,ndf)
          end if

        else if(ityp.eq.5) then
c.....    protocol link nodes_x2,y2 to nodes_x1,y1
          if(ior.gt.0) then
            if(prt)write(iow,2005)idmx,idmy,x2,y2,x1,y1,(idl(i),i=1,ndf)
          else
            if(prt)write(  *,2005)idmx,idmy,x2,y2,x1,y1,(idl(i),i=1,ndf)
          end if

        else if(ityp.eq.6) then
c.....    protocol link nodes_x2 to nodes_x1 with y(x2)=-y(x1)
          if(ior.gt.0) then
            if(prt) write(iow,2006) idm,x2,x1,(idl(i),i=1,ndf)
          else
            if(prt) write(  *,2006) idm,x2,x1,(idl(i),i=1,ndf)
          end if

        else if(ityp.eq.7) then
c.....    protocol link nodes_ob boundary of plnae x3=c to corner node
          if(ior.gt.0) then
            if(prt) write(iow,2007) x1m,x2m,x3m,x3p,(idl(i),i=1,ndf)
          else
            if(prt) write(  *,2007) x1m,x2m,x3m,x3p,(idl(i),i=1,ndf)
          end if

        end if
        
c....   find nodes

        if(ityp.eq.1) then
c.....    link node1 - node2 to node1
          iplk1(n1q) = -n1q        ! master node

c.....    find all nodes to link
          do 81 i = ni,n2,inc  ! loop over other nodes to link
cww         iq = iprttie(gtie,i) ! wrong if tied node is used!   
            iq = i                 ! modifed ww 30/3/11
c....       loop over dofs   
            i2p = 0
            do 71 j = 1,ndf 
              ir = iplk2(j,iq) 
              if (ir.gt. 0 .and. idl(j).eq.0) then ! dof has been used before 
                if(prt) write(iow,3002) i,n1,j,ir,nlk
                if(prt) write(  *,3002) i,n1,j,ir,nlk
                go to 71 
              end if
              nmax = 0
              if(idl(j).eq.0) nmax=max(nmax,id(j,iq))  ! max. eqn no to link at this node
              if(nmax.gt.0) then
                iplk1(iq)   =  n1q     ! slave node
                iplk2(j,iq) =  n1q     ! link cond. = no of master node
                i2p = 1 
                nelim = 0
c.....          link nodes with defined code
                if(id(j,iq).gt.0) nelim = nelim + 1
                idsave   = id(j,iq)
                id(j,iq) = id(j,n1q)  ! link node if idl=0
c.....          look if this eq. number occurs elsewhere (=slave nodes) 
                do ii = 1,numnp 
                  do jj = 1,ndf
                    if( id(jj,ii).eq.idsave ) then
                      id(jj,ii) = id(j,n1q) 
                    end if                   
                  end do
                end do  
c....           loop through all nodes to reduce the equation numbers
                do 61 ii = 1,numnp
                  do 61 jj  = 1,ndf
                    if(id(jj,ii).gt.nmax) then   ! only if eqn no is higher
                      id(jj,ii) = id(jj,ii) - nelim
                    end if
61              continue
                neq = neq - nelim   ! total number of eqn
              end if
71          continue
c.....      print linked node
            if(i2p.eq.1) then  
              if(ior.gt.0) then
                if(prt) write(iow,4000) i,iq
              else
c               if(prt) write(*  ,4000) i,iq
              end if
            end if
81        continue

        else if(ityp.eq.2) then
c.....    link nodes on x1 to node1
          iplk1(n1q) = -n1q  ! master node

c.....    find all nodes to link
          do 82 i = 1,numnp  ! loop over other nodes to link
            if(i.eq.n1) goto 82 
            if(abs(x(idm,i)-x0).gt.dx) goto 82 
cww         iq = iprttie(gtie,i) ! wrong if tied node is used!   
            iq = i                 ! modifed ww 30/3/11
c....       loop over dofs   
            i2p = 0
            do 72 j = 1,ndf 
              ir = iplk2(j,iq) 
              if (ir.gt. 0 .and. idl(j).eq.0) then ! dof has been used before 
                if(prt) write(iow,3002) i,n1,j,ir,nlk
                if(prt) write(  *,3002) i,n1,j,ir,nlk
                go to 72 
              end if

              nmax = 0
              if(idl(j).eq.0) nmax=max(nmax,id(j,iq)) 
              if(nmax.gt.0) then
                iplk1(iq)   = n1q  ! slave node
                iplk2(j,iq) = n1q  ! link cond. = no of master node 
                i2p = 1 
                nelim = 0
c.....          link nodes
                if(id(j,iq).gt.0) nelim = nelim + 1
                idsave   = id(j,iq)
                id(j,iq) = id(j,n1q)  ! link if idl=0
c.....          look if this eq. number occurs elsewhere (=slave nodes) 
                do ii = 1,numnp 
                  do jj = 1,ndf
                    if( id(jj,ii).eq.idsave ) then
                      id(jj,ii) = id(j,n1q) 
                    end if                   
                  end do
                end do  
c....           loop through all nodes to reduce the equation numbers
                do 62 ii = 1,numnp
                  do 62 jj  = 1,ndf
                    if(id(jj,ii).gt.nmax) then
                      id(jj,ii) = id(jj,ii) - nelim
                    end if
62              continue
                neq = neq - nelim
              end if
72          continue
c.....      print linked node
            if(i2p.eq.1) then 
              if(ior.gt.0) then
                if(prt) write(iow,4000) i,iq
              else
                if(prt) write(*  ,4000) i,iq
              end if
            end if
82        continue

        else if(ityp.eq.3) then
c.....    link (node1-node2)_i to (node1-node2)_j
          iic   = 0
          do 83 i1 = na1,ne1,inc1  ! loop over  1. node to link
            i2   = na2 + iic*inc2   ! associated 2. node 
cww         i1q = iprttie(gtie,i1) ! wrong if tied node is used!   
            i1q = i1                 ! modifed ww 30/3/11
            iplk1(i1q) = -i1q  ! master node
cww         i2q = iprttie(gtie,i2) ! wrong if tied node is used!   
            i2q = i2                 ! modifed ww 30/3/11
c....       loop over dofs   
            i2p = 0
            do 73 j = 1,ndf 
              ir = iplk2(j,i2q) 
              if (ir.gt. 0 .and. idl(j).eq.0) then ! dof has been used before 
                if(prt) write(iow,3002) i2,i1,j,ir,nlk
                if(prt) write(  *,3002) i2,i1,j,ir,nlk
                go to 73 
              end if
              
              ir = iplk2(j,i1q)
              if (ir.eq. i2q .and. idl(j).eq.0) then ! circle application on this dof for these nodal combination 
                if(prt) write(iow,3003) i2,i1,j,nlk
c               if(prt) write(  *,3003) i2,i1,j,nlk
                go to 73 
              end if
              
c....         link if not b.c.
              nmax = 0
              if(idl(j).eq.0) nmax=max(nmax,id(j,i2q))  ! max. eqn no to link at this node
              if(nmax.gt.0) then
                iplk1(i2q)   = i1q ! slave node
                iplk2(j,i2q) = i1q ! link cond. = no of master node
                i2p = 1
                nelim = 0
c.....          link node 2 to node 1
                if(id(j,i2q).gt.0) nelim = nelim + 1
                idsave    = id(j,i2q)
                id(j,i2q) = id(j,i1q)  ! link
c.....          look if this eq. number occur elsewhere (=slave nodes) 
                do ii = 1,numnp 
                  do jj = 1,ndf
                    if( id(jj,ii).eq.idsave ) then
                      id(jj,ii) = id(j,i1q) 
                    end if                   
                  end do
                end do  
c....           loop through all nodes to reduce the equation numbers
                do 63 ij = 1,numnp
                  do 63 jj  = 1,ndf
                    if(id(jj,ij).gt.nmax) then
                      id(jj,ij) = id(jj,ij) - nelim
                    end if
63              continue
                neq = neq - nelim
              end if  
73          continue

c.....      print linked node
            if(i2p.eq.1) then 
              if(ior.gt.0) then
                if(prt) write(iow,4001) i2,i2q,i1,i1q
              else
c               if(prt) write(*  ,4001) i2,i2q,i1,i1q
              end if
            end if
            iic = iic+1

83        continue

        else if(ityp.eq.4) then
c....     link nodes_x2 to nodes_x1
          do 84 i1 = 1,numnp  ! loop over nodes on x1
            if(abs(x(idm,i1)-x1).gt.dx) goto 84 
c....       node on x1
cww         i1q = iprttie(gtie,i1) ! wrong if tied node is used!   
            i1q = i1                 ! modifed ww 30/3/11

            iplk1(i1q) = -i1q        ! master node

            do 74 i2 = 1,numnp  ! loop over nodes on x2
              if(abs(x(idm,i2)-x2).gt.dx) goto 74 

c....         compare other coordinates  
              do jdm = 1,ndm
                if(jdm.ne.idm) then
                  if(abs(x(jdm,i1)-x(jdm,i2)) .gt. dx ) go to 74
                end if
              end do ! jdm
c....         node on x2
cww           i2q = iprttie(gtie,i2) ! wrong if tied node is used!   
              i2q = i2                 ! modifed ww 30/3/11
c....         link dofs of node i2q to node i1q
              i2p = 0
              do 64 j = 1,ndf 
                ir = iplk2(j,i2q) 
                if (ir.gt. 0 .and. idl(j).eq.0) then ! dof has been used before 
                  if(prt) write(iow,3002) i2,i1,j,ir,nlk
c                 if(prt) write(  *,3002) i2,i1,j,ir,nlk
                  go to 64 
                end if
  
                ir = iplk2(j,i1q)
                if (ir.eq. i2q .and. idl(j).eq.0) then ! circle application on this dof for these nodal combination 
                  if(prt) write(iow,3003) i2,i1,j,nlk
c                 if(prt) write(  *,3003) i2,i1,j,nlk
                  go to 64 
                end if
                
c....           link if not b.c.
                nmax = 0
                if(idl(j).eq.0) nmax=max(nmax,id(j,i2q)) 
                if(nmax.gt.0) then
                  iplk1(i2q)   = i1q  ! slave node
                  iplk2(j,i2q) = i1q  ! link cond. = no of master node
                  i2p = 1 
                  nelim = 0
c.....            link nodes
                  if(id(j,i2q).gt.0) nelim = nelim + 1
                  idsave    = id(j,i2q)
                  id(j,i2q) = id(j,i1q)  ! link if idl=0

c.....            look if this eq. number occur elsewhere (=slave nodes) 
                  do ii = 1,numnp 
                    do jj = 1,ndf
                      if( id(jj,ii).eq.idsave ) then
                        id(jj,ii) = id(j,i1q) 
                      end if                   
                    end do
                  end do  
c....             loop through all nodes to reduce the equation numbers
                  do 54 ii = 1,numnp
                    do 54 jj  = 1,ndf
                      if(id(jj,ii).gt.nmax) then
                        id(jj,ii) = id(jj,ii) - nelim
                      end if
54                continue
                  neq = neq - nelim
                end if
64            continue
c.....        print linked node
              if(i2p.eq.1) then
                if(ior.gt.0) then
                  if(prt) write(iow,4001) i2,i2q,i1,i1q
                else
c                 if(prt) write(*  ,4001) i2,i2q,i1,i1q
                end if
              end if
74          continue ! i2

84        continue ! i1  

        else if(ityp.eq.5) then
c....     link nodes_x2,y2 to nodes_x1,y1
          do 85 i1 = 1,numnp  ! loop over nodes on x1,y1
            if(abs(x(idmx,i1)-x1).gt.dx) goto 85 
            if(abs(x(idmy,i1)-y1).gt.dx) goto 85 
c....       node on x1
cww         i1q = iprttie(gtie,i1) ! wrong if tied node is used!   
            i1q = i1                 ! modifed ww 30/3/11

            iplk1(i1q) = -i1q        ! master node

            do 75 i2 = 1,numnp  ! loop over nodes on x2,y2
              if(abs(x(idmx,i2)-x2).gt.dx) goto 75 
              if(abs(x(idmy,i2)-y2).gt.dx) goto 75 

c....         compare other coordinate  

              do 45 jdm = 1,ndm
                if(jdm.eq.idmx) go to 45
                if(jdm.eq.idmy) go to 45
                if(abs(x(jdm,i1)-x(jdm,i2)) .gt. dx ) go to 75
45            continue ! jdm
c....         node on x2
cww           i2q = iprttie(gtie,i2) ! wrong if tied node is used!   
              i2q = i2                 ! modifed ww 30/3/11
c....         link dofs of node i2q to node i1q
              i2p = 0
              do 65 j = 1,ndf 
                ir = iplk2(j,i2q) 
                if (ir.gt. 0 .and. idl(j).eq.0) then ! dof has been used before 
                  if(prt) write(iow,3002) i2,i1,j,ir,nlk
c                 if(prt) write(  *,3002) i2,i1,j,ir,nlk
                  go to 65 
                end if

                ir = iplk2(j,i1q)
                if (ir.eq. i2q .and. idl(j).eq.0) then ! circle application on this dof for these nodal combination 
                  if(prt) write(iow,3003) i2,i1,j,nlk
c                 if(prt) write(  *,3003) i2,i1,j,nlk
                  go to 65 
                end if

c....           link if not b.c.
                nmax = 0
                if(idl(j).eq.0) nmax=max(nmax,id(j,i2q)) 
                if(nmax.gt.0) then
                  iplk1(i2q)   = i1q  ! slave node
                  iplk2(j,i2q) = i1q  ! link cond. = no of master node
                  i2p = 1 
                  nelim = 0
c.....            link nodes
                  if(id(j,i2q).gt.0) nelim = nelim + 1
                  idsave    = id(j,i2q)
                  id(j,i2q) = id(j,i1q)  ! link if idl=0

c.....            look if this eq. number occur elsewhere (=slave nodes) 
                  do ii = 1,numnp 
                    do jj = 1,ndf
                      if( id(jj,ii).eq.idsave ) then
                        id(jj,ii) = id(j,i1q) 
                      end if                   
                    end do
                  end do  
c....             loop through all nodes to reduce the equation numbers
                  do 55 ii = 1,numnp
                    do 55 jj  = 1,ndf
                      if(id(jj,ii).gt.nmax) then
                        id(jj,ii) = id(jj,ii) - nelim
                      end if
55                continue
                  neq = neq - nelim
                end if
65            continue
c.....        print linked node
              if(i2p.eq.1) then
                if(ior.gt.0) then
                  if(prt) write(iow,4001) i2,i2q,i1,i1q
                else
c                 if(prt) write(*  ,4001) i2,i2q,i1,i1q
                end if
              end if
75          continue ! i2

85        continue ! i1  

        else if(ityp.eq.6) then
c....     link nodes_x2 to nodes_x1 with y(x2)=-y(x1)
          do 86 i1 = 1,numnp  ! loop over nodes on x1
            if(idm.eq.3) stop 'SR PLINK: for typ 6 only idm=1,2' 
            if(abs(x(idm,i1)-x1).gt.dx) goto 86 

c            if(abs(x(3,i1)).gt.dx
c     +      .and.abs(x(1,i1)).gt.dx.and.abs(x(2,i1)).gt.dx) goto 86

c....       node on x1
cww         i1q = iprttie(gtie,i1) ! wrong if tied node is used!   
            i1q = i1                 ! modifed ww 30/3/11

            iplk1(i1q) = -i1q        ! master node

            do 76 i2 = 1,numnp  ! loop over nodes on x2
              if(abs(x(idm,i2)-x2).gt.dx) goto 76 


c....         compare other coordinates
              if(abs(x(3,i1)-x(3,i2)) .gt. dx ) go to 76
               if(idm.eq.1) then
                 if(abs(x(2,i1)+x(2,i2)) .gt. dx ) go to 76
               else if(idm.eq.2) then
                 if(abs(x(1,i1)+x(1,i2)) .gt. dx ) go to 76
               end if

c....         node on x2
cww           i2q = iprttie(gtie,i2) ! wrong if tied node is used!   
              i2q = i2                 ! modifed ww 30/3/11
c....         link dofs of node i2q to node i1q
              i2p = 0
              do 66 j = 1,ndf 
                ir = iplk2(j,i2q) 
                if (ir.gt. 0 .and. idl(j).eq.0) then ! dof has been used before 
                  if(prt) write(iow,3002) i2,i1,j,ir,nlk
c                 if(prt) write(  *,3002) i2,i1,j,ir,nlk
                  go to 66 
                end if

                ir = iplk2(j,i1q)
                if (ir.eq. i2q .and. idl(j).eq.0) then ! circle application on this dof for these nodal combination 
                  if(prt) write(iow,3003) i2,i1,j,nlk
c                 if(prt) write(  *,3003) i2,i1,j,nlk
                  go to 66 
                end if

c....           link if not b.c.
                nmax = 0
                if(idl(j).eq.0) nmax=max(nmax,id(j,i2q)) 
                if(nmax.gt.0) then
                  iplk1(i2q)   = i1q  ! slave node
                  iplk2(j,i2q) = i1q  ! link cond. = no of master node
                  i2p = 1 
                  nelim = 0
c.....            link nodes
                  if(id(j,i2q).gt.0) nelim = nelim + 1
                  idsave    = id(j,i2q)
                  id(j,i2q) = id(j,i1q)  ! link if idl=0

c.....            look if this eq. number occur elsewhere (=slave nodes) 
                  do ii = 1,numnp 
                    do jj = 1,ndf
                      if( id(jj,ii).eq.idsave ) then
                        id(jj,ii) = id(j,i1q) 
                      end if                   
                    end do
                  end do  
c....             loop through all nodes to reduce the equation numbers
                  do 56 ii = 1,numnp
                    do 56 jj  = 1,ndf
                      if(id(jj,ii).gt.nmax) then
                        id(jj,ii) = id(jj,ii) - nelim
                      end if
56                continue
                  neq = neq - nelim
                end if
66            continue
c.....        print linked node
              if(i2p.eq.1) then
                if(ior.gt.0) then
                  if(prt) write(iow,4001) i2,i2q,i1,i1q
                else
c                 if(prt) write(*  ,4001) i2,i2q,i1,i1q
                end if
              end if
76          continue ! i2

86        continue ! i1  

        else if(ityp.eq.7) then
c....     link nodes for x3=const on border of area in x1-x2 dir
c....     loop on master nodes between -x3/2 and +x3/2   
          do 87 iz = 1,numnp  ! loop over other nodes to link
            if(abs(x(1,iz)-x1m).gt.dx1) goto 87 !  
            if(abs(x(2,iz)-x2m).gt.dx2) goto 87 ! node=[x1m,x2m,x3] 
            x3q = x(3,iz)
            n1q = iprttie(gtie,iz)
            iplk1(n1q) = -n1q  ! master node
c....       find all nodes to link on plane with x3=x3q
c....       nodes could be on x1m,x2 x1p,x2 or x1,x2m x1,x2p  
            do 77 ixy = 1,numnp  ! loop over other nodes to link
              if(abs(x(3,ixy)-x3q).gt.dx3) goto 77 ! x3-coor    

              if(abs(x(1,ixy)-x1p).le.dx1) goto 771
              if(abs(x(1,ixy)-x1m).le.dx1) goto 771
              if(abs(x(2,ixy)-x2p).le.dx2) goto 771
              if(abs(x(2,ixy)-x2m).le.dx2) goto 771
c....         node not to link
              goto 77
771           continue          
c...          do not link for master node
              if(ixy.eq.iz) goto 77
c....         link node
c....         loop over dofs   
              iq  = ixy
              i2p = 0
              do 67 j = 1,ndf 
              ir = iplk2(j,iq) 
              if (ir.gt. 0 .and. idl(j).eq.0) then ! dof has been used before 
                if(prt) write(iow,3002) iq,n1q,j,ir,nlk
                if(prt) write(  *,3002) iq,n1q,j,ir,nlk
                go to 67 
              end if

              nmax = 0
              if(idl(j).eq.0) nmax=max(nmax,id(j,iq)) 
              if(nmax.gt.0) then
                iplk1(iq)   = n1q  ! slave node
                iplk2(j,iq) = n1q  ! link cond. = no of master node 
                i2p = 1 
                nelim = 0
c.....          link nodes
                if(id(j,iq).gt.0) nelim = nelim + 1
                idsave   = id(j,iq)
                id(j,iq) = id(j,n1q)  ! link if idl=0
c.....          look if this eq. number occurs elsewhere (=slave nodes) 
                do ii = 1,numnp 
                  do jj = 1,ndf
                    if( id(jj,ii).eq.idsave ) then
                      id(jj,ii) = id(j,n1q) 
                    end if                   
                  end do
                end do  
c....           loop through all nodes to reduce the equation numbers
                do 57 ii = 1,numnp
                  do 57 jj  = 1,ndf
                    if(id(jj,ii).gt.nmax) then
                      id(jj,ii) = id(jj,ii) - nelim
                    end if
57              continue
                neq = neq - nelim
              end if
67            continue
c.....        print linked node
              if(i2p.eq.1) then 
                if(ior.gt.0) then
                  if(prt) write(iow,4001) iq,iq,iz,n1q
                  
                else
                  if(prt) write(*  ,4001) iq,iq,iz,n1q
                end if
              end if
77          continue ! ixy
87        continue ! iz

        end if   

c       print id-array
c        write(iow,*) ' final  node i, id(j,i), iplk2(j,i)'
c        write(  *,*) ' final  node i, id(j,i), iplk2(j,i)'
c        do i = 1,numnp
c          write(iow,5000) i,(id(j,i),j=1,ndf),(iplk2(j,i),j=1,ndf)         
c          write(  *,5000) i,(id(j,i),j=1,ndf),(iplk2(j,i),j=1,ndf)         
c        end do
c5000    format(10i4)
        
        go to 10
      else
        write(  *,3000)
        write(iow,3000)
        return
      end if
c.... formats
2001  format(1x,'link nodes(tied to)',i6,'(',i6,')',' to',
     +       i6,'(',i6,')',' at increments ',i4,
     +       ' code ',6i4,/,(63x,6i4))
2002  format(1x,'link nodes on ',i2,' coord.(value ',e12.5,')',
     +       ' to master node(tied to)  ',i6,'(',i6,')',' code ',6i4,/,
     +       (61x,6i4))
2003  format(5x,'link nodes(tied to) ',i6,'(',i6,')',' to',i6,'(',i6,')'
     +       ,' at increments ',i4,' to',/, 
     +       1x, 'master nodes(tied to) ',i6,'(',i6,')',' to',i6,'(',i6,
     +       ')',' at increments ',i4,' code ',6i4,/,(63x,6i4))
2004  format(5x,'link nodes(tied to) in',i3,'-dir',/,  
     +       5x,'Nodes on coordinate x2 ',e12.5,/,
     +       5x,'to nodes on coordinate x1 ',e12.5,/,
     +       5x,' code ',6i4,/,(63x,6i4))

2005  format(5x,'link nodes(tied to) in',i3,i3,'-dir',/,  
     +       5x,'Nodes on coordinate x2,y2 ',e12.5,e12.5,/,
     +       5x,'to nodes on coordinate x1,y1 ',e12.5,e12.5,/,
     +       5x,' code ',6i4,/,(63x,6i4))

2006  format(5x,'link nodes(tied to) in',i3,'-dir',/,  
     +       5x,'Nodes on coordinate x2 ',e12.5,/,
     +       5x,'to nodes on coordinate x1 ',e12.5,'with y(x2)=-y(x1)'/, 
     +       5x,' code ',6i4,/,(63x,6i4))

2007  format(5x,'link nodes(tied to) on b.c. in 1-2-dir',/,  
     +       5x,'to node with x1= ',e12.5,'and x2= ',e12.5,/,
     +       5x,'x3min,x3max ',2e12.5,'with'/, 
     +       5x,' code ',6i4,/,(63x,6i4))
2008  format(/,3x,'F i n a l  l i n k  c o n d i t i o n s',/,
     +       5x,'Node, dofs linked to master node (0=no link)',
     +       ' eq. no. of dofs',/)  
2009  format(5x,13i6)
2010  format(' No file named ',a55,' exists.                      ')
3000  format(5x,'**ERROR** link file does not exist')
3001  format(5x,'**Stop** linked node ',i5,' does not exist.')
3002  format(5x,'Link node ',i5,' to node ',i5,' : dof ',i3,
     +       5x,' is always in use by node',i5,' (input card ',i3,')' )
3003  format(5x,'Link node ',i5,' to node ',i5,' : dof ',i3,
     +       5x,' is used before opposite (input card ',i3,')' )
4000  format(5x,' node(tied to) ',i6,'(',i6,')')
4001  format(5x,' node(tied to) ',i6,'(',i6,')',
     +       ' linked to node(tied to) ',i6,'(',i6,')')
      end
c
      subroutine plinka
c----------------------------------------------------------------------
c
c      Purpose: Read data and save on a file for future processing
c
c      Inputs:
c         fext      - Name of file extender for save file
c
c      Outputs:
c         none      - Data saved to disk
c
c      Comments:    - variable input possible, read only complete line
c                     without action
c
c----------------------------------------------------------------------
      USE comfil
      USE iodata
      USE iofile
      logical pcomp,lsave
      character*229 yyy*80
      character*229 flink
c
      flink = finp
      call addext(flink,'lnk ')
      inquire(file=flink,exist=lsave)
      if(lsave) then
          open(ios,file=flink,status='old')
      else
        open(ios,file=flink,status='new')
      end if
c
  10  read(ior,1000,err=901,end=902) yyy
      write(ios,1000) yyy
c
      if(.not.pcomp(yyy,'    ',4)) then
        go to 10
      else
        close(ios)
      end if
      return
c
 901  call errclr('PLINKA')
      return
c
 902  call endclr('PLINKA',yyy)
c
1000  format (a80)
      end
c
      subroutine pload(id,f,f0,b,nn,p)
c----------------------------------------------------------------------
c      Purpose: Form nodal load vector for current time
c
c        b = b + f*p + f0
c
c      Inputs:
c         id(*)    - Equation numbers for degree of freedom
c         f(nneq)  - Nodal load vector f
c         f0(nneq) - Nodal load vector f0
c         nn       - Number of equations nneq 
c         p        - Total proportional load level

c      Outputs:
c         b(neq)   - Current Total Nodal load vector
c
c----------------------------------------------------------------------
      USE cdata
      USE edgdat
      USE fdata
      implicit double precision (a-h,o-z)
c..... type declaration for variables
      integer j, n,nn,ne
      real*8  p
c..... type declaration for arrays
      integer id(*)
      real*8  f(*),f0(*),b(*)
c
      if(nde.eq.0) then
         ne = nn
      else
         ne = nn + nde*edge1(size(edge1))
      end if
      fl(11) = .false.
cww      call pzero (b,ne)
c.... load nodal forces
      do 100 n = 1,nn
        j = id(n)
        if(j.gt.0) b(j) = f(n)*p + f0(n) + b(j)
100   continue
c.... load edge forces
      if(ne.gt.nn) call ploadf(ne-nn,edge3,edge4,p,b)
c
      end
c
      subroutine ploadf(nn,ieb,ef,p,b)
c----------------------------------------------------------------------
c      Purpose: Form load edge forces
c
c      Inputs:

c      Outputs:
c
c----------------------------------------------------------------------
c..... type declaration for variables
      integer j, n,nn
      real*8  p
c..... type declaration for arrays
      integer ieb(*)
      real*8  ef(*),b(*)
c.... load edge forces
      do 100 n = 1,nn
        j = ieb(n)
        if(j.gt.0) b(j) = ef(n)*p + b(j)
100   continue
c
      end
c
      subroutine ploadi
c----------------------------------------------------------------------
c      Purpose:  Input data for surface loading from elements
c
c      parameter: 
c      ns              - no. nodes on element 
c      nv              - no. of load values
c      nl              - load typ always=1
c      nums            - no. of different elmt.typs  (<(26)
c      iels(1,nums)    - el.typ
c      iels(2,nums)    - kmax = position  m-vec
c      iels(3,nums)    - nl = typ->1
c      iels(4,nums)    - nrec = no. of records for the current loading state
c
c      inods(nums)     - ns = no. nodes
c      ivals(nums)     - nv = no. values
c      inc1            - inc nodes
c      inc2            - inc values
c
c----------------------------------------------------------------------
      USE bdata
      USE cdata
      USE errchk
      USE iofile
      USE plong
      USE psize
      USE sldata
      implicit double precision (a-h,o-z)
      logical error
      dimension ixl(9),valu(27),td(52)
      common m(maxm) 

c.... input of surface loading data for each set 
100   if(ior.lt.0) write(*,3000)
      call dinput(td,4)
      if(errck) go to 100

      ns = td(2)              ! no. nodes 
      nv = td(3)              ! no. values
      nl = td(4)              ! typ->1
      if(nl.eq.0) nl = 1

c.... set parameters
      nums = nums + 1         ! up to 26 diff.elmt.typs allowed

      iels(1,nums) = td(1)    ! element typ
      iels(2,nums) = kmax     ! position  m-vec
      iels(3,nums) = nl       ! typ->1

      inods(nums) = ns        ! no. nodes
      ivals(nums) = nv        ! no. values
      inc1 = ns + mod(ns,ipr) ! inc nodes
      inc2 = inc1 + nv*ipr    ! inc values

      nrec = 0
c.... input records for the current loading state
      if(ior.gt.0) then
        if(ns.lt.10.and.nv.lt.3) then
            write(iow,2004) o,head,(i,i=1,9),(i,i=1,3)
        else
            write(iow,2000) o,head,(i,i=1,ns)
            write(iow,2001) (i,i=1,nv)
        end if
      else
          write(*,2000) o,head,(i,i=1,ns)
          write(*,2001) (i,i=1,nv)
      end if
      error = .false.
      nsv = ns + nv 
      nsvq= nsv + 2
      if(nsvq.gt.52) then
        call drawmess('only 52 nodes+values allowed in ploadi',1,0)
        return
      end if
1     continue
c.....read values, increment of 8!
101   do i = 1,nsvq,8
          k = min(nsvq-i+1,8)
          call dinput(td(i),k)
          if(errck) go to 101
      end do
c.....copy nodes at element
      if(ns.gt.9) then
        call drawmess('only 9 nodes allowed in ploadi',1,0)
        return
      end if
      call pzeroi(ixl,9)
      do i = 1,ns
          ixl(i) = td(i)
      end do
c.....copy p-values for loading
      if(nv.gt.27) then
        call drawmess('only 27 values allowed in ploadi',1,0)
        return
      end if
      call pzero(valu,27)
      do i = 1,nv
          valu(i) = td(i+ns)
      end do
c.....copy gen-increment
      ngen  = td(nsv+1)    
c.....copy number of generations
      lgen  = td(nsv+2)    
c
      if(ixl(1).ne.0) then   ! values occur -> store
c.....  in case of error 
2         do i = 1,ns
            if(ixl(i).gt.numnp) then
              if(ior.ge.0) then
              write(iow,2006) 
                write(iow,2002) (ixl(k),k=1,ns)
cww           stop
              return
              else if(ior.lt.0) then
                write(*  ,2006) 
                write(*  ,2002) (ixl(k),k=1,ns)
                error = .true.
                go to 1
              end if
            end if
        end do
c.....  make output
          if(ior.gt.0) then
            if(ns.lt.10.and.nv.lt.3) then
              write(iow,2005) (ixl(i),i=1,9),(valu(i),i=1,3)
            else
              write(iow,2002) (ixl(i),i=1,ns)
              write(iow,2003) (valu(i),i=1,nv)
            end if
          else
            write(*,2002) (ixl(i),i=1,ns)
            write(*,2003) (valu(i),i=1,nv)
          end if
          call pstores(ns,nv,ixl,valu,m(kmax),m(kmax+inc1))
          kmax = kmax + inc2
          nrec = nrec + 1
c.....  repeat for generation
        if(lgen.gt.0) then
          lgen = lgen - 1 
c.....    new nodes
          do i = 1,ns
              ixl(i) = ixl(i) + ngen
          end do
          goto 2
        end if
          go to 1
      end if
      iels(4,nums) = nrec
      return
c.... formats
2000  format(a1,20a4//'     s u r f a c e   l o a d i n g'//
     1  4x,10(i3,'-node'))
2001  format(5x,6(i5,'-value '))
2002  format(/8i8)
2003  format(4x,1p6e12.4)
2004  format(a1,20a4//'     s u r f a c e   l o a d i n g'//
     1 1x,9(i2,'-no'),1x,3(i5,'-value '))
2005  format(1x,9i5,1x,1p3e12.4)
2006  format(1x,'Wrong Node Number in SLOA')
3000  format(' Input: el-type, #-nodes, #-vals, ld-type'/'   >',$)
      end
c
      subroutine ploadl(id,jd,ld,b,ad,al,au,p,s,x,xl,u,ul,
     1                  propq,ndf,ndm,fl,aufl,alfl)
c----------------------------------------------------------------------
c
c      Purpose: Assemble surface loading from element contributions.
c               Data is input using ploadi routine.

c      Inputs:
c         id(ndf,*)  - Active equation numbers
c         jd(*)      - Pointers for row/columns in tangent array
c         x(ndm,*)   - Nodal coordinates of mesh
c         u(ndf,*)   - Nodal solutions of mesh
c         propq      - Current total proportional load level
c         ndf        - Number dof/node
c         ndm        - Spatial dimension of mesh
c         fl         - Flag, true if uncompressed vector assembled
c         aufl       - Flag, form tangent if true
c         alfl       - Flag, tangent is unsymmetric if true
c
c         nums       -
c         ns              - no. nodes on element 
c         nv              - no. of load values
c         nl              - load typ always=1
c         nums            - no. of different elmt.typs  (<(26)
c         iels(1,nums)    - iel=el.typ
c         iels(2,nums)    - mad = position  m-vec
c         iels(3,nums)    - ma = typ->1
c         iels(4,nums)    - nrec = no. of records for the current loading state
c
c         inods(nums)     - nel = no. nodes
c         ivals(nums)     - mct = no. values
c         inc1            - inc nodes
c         inc2            - inc values
c
c
c      Scratch:
c         ld(ndf,*)  - Element local/global equation numbers
c         p(*)       - Element vector
c         s(*)       - Element matrix
c         xl(ndm,*)  - Element nodal coordinates
c         ul(ndf,*)  - Element solution/rate parameters

c      Outputs:
c         b(*)       - Includes surface load effect in residual/reaction
c         ad(*)      - Diagonal part of tangent array from surface load
c         al(*)      - Lower part of tangent array
c         au(*)      - Lower part of tangent array
c
c----------------------------------------------------------------------
      USE cdata
      USE eldata
      USE iscsr
      USE mdat2
      USE psize
      USE sldata
      USE smpak
      USE soltyp
      implicit double precision (a-h,o-z)
      logical fl,alfl,aufl
      dimension id(ndf,*),jd(*),ld(ndf,*),b(*),p(*),xl(ndm,*)
     1,x(ndm,*),u(ndf,*),ul(ndf,*),valu(27),s(*),ad(*),al(*),au(*)
      common m(maxm)
      if(nums.le.0) return
      ir = 3
      if(aufl) ir = 2
      do 150 nn = 1,nums
        iel = iels(1,nn)
        nad = iels(2,nn)
        ma  = iels(3,nn)
        nrec= iels(4,nn)
        nel = inods(nn)
        mct = ivals(nn)
        inc1 = nel + mod(nel,ipr)
        inc2 = inc1 + mct*ipr
        ns1 = ndf*nel
        do 120 j = 1,nrec
cww       do 110 n = 1,nel   ! n is used in eldata  ->nk
          do 110 nk= 1,nel
            ii = m(nad+nk-1)
            do 100 i = 1,ndm
100         xl(i,nk) = x(i,ii)
            do 110 i = 1,ndf
              ul(i,nk) = u(i,ii)
              ld(i,nk) = id(i,ii)
              if(fl) ld(i,nk) = ndf*ii - ndf + i
110       continue
          call pangl(m(nad),nel,aang,bang,nrot)
          call pzero(p,ns1)
          call pzero(s,ns1*ns1)
          call pzero(valu,mct)
c....     multiply load values by actual value of prop
          call daxpty(mct,propq,m(nad+inc1),valu)         

c....     form element array-transform all(!) mixed displ. to pure global dir.
          if(nrot.gt.0)
     1    call ptrans(ia,itrot,aang,ul,p,s,nel,nen,ndf,ns1,1)

          call elmlib(valu,ul,xl,m(nad),valu,s,p,h1,h2,h3,
     1                ndf,ndm,ns1,iel,7)

c....     transform s,p to mixed global displacements back, see ir
          if(nrot.gt.0)
     1    call ptrans(ia,itrot,aang,ul,p,s,nel,nen,ndf,ns1,ir)

c....     add to total array with respect to solver type
          if(istyp.eq.0 ) then ! Standard Solver
            call dasbly(s,p,ld,jd,ns1,alfl,aufl,.true., b,al,au,ad)
          
          else if(istyp.eq.1.or.istyp.eq.2) then ! SM
            call dasbl2(s,p,ld,jd,ns1,.false.,.true.,
     +                  b,ad,csrja,smsperm1,lodrv)

          else if(istyp.ge.3.and.istyp.le.8) then ! all other CSR
            call dasbl_csr(s,p,ld,jd,ns1,aufl,.true.,b,ad,csrja)

          end if
120     nad = nad + inc2
150   continue
      return
      end
c
      subroutine ploads(u,b,propt,flg,aufl,alfl)
c----------------------------------------------------------------------
c
c      Purpose:  Surface loading routine: Calls ploadl/ploadq
c
c      Inputs:
c         u(*)     - Current solution state
c         propt    - actual proportional load parameter
c         flg      - Flag, form reaction (uncompressed) if true
c         aufl     - Flag, form tangent if true
c         alfl     - Flag, tangent is unsymmetric if true

c      Outputs:
c         b(*)     - Residual/reaction array
c                    Tangent returned through blank common address
c
c     Comments: Stiffness matrix must be at m(na) etc!!  
c
c      modified by ww for macros SLOA and QLOA 2/06
c
c----------------------------------------------------------------------
      USE cdata
      USE hdatam
      USE mdata
      USE ndata
      USE sdata
      !implicit double precision (a-h,o-z)
      implicit none
      logical flg,aufl,alfl
      real*8 u(*),b(*), propt

      call pzero(b,ndf*numnp) ! sequence pload and ploads changed, ww 
      hflgu  = .false.        ! no update of history variables 
      h3flgu = .false.
      
c.... for SLOA (original)      
      call ploadl(psid,jdt12,eeqn,b,gstiff,gstiff(nal),gstiff(nau)
     1      ,epve,ekma,coor,ecor,u,edis,propt,ndf,ndm,flg,aufl,alfl)

c.... for QLOA (added)     
      call ploadq(u,b,gstiff,gstiff(nal),numel,flg,aufl,alfl,propt)

      return
      end
c
      subroutine ploadq(u,b,au,al,numel,dflg,aufl,alfl,propt)
c-----------------------------------------------------------------
c
c      Purpose: calculate element load vectors based on values defined 
c               in SR qloadi (MESH>QLOA) and add to load vector dr 
c               general linear and ninlinear case e.g. follower loads
c
c      Inputs:
c         u          - Nodal solution vector
c         id(*)      - Equation numbers for each active dof
c         numel      - Number of elements in mesh
c         dflg       - Flag, form reaction (uncompressed) if true
c         aufl       - Flag, form tangent if true
c         alfl       - Flag, tangent is unsymmetric if true
c
c      Outputs:
c
c         b(*)       - Nodal load vector   from QLOA for mesh
c         au(*)      - Upper Part K (tang) from QLOA for mesh
c         al(*)      - Lower Part K (utan) from QLOA for mesh
c
c         p          - is part of b 
c         s          - is calculated via common adress
c
c.... ww bs uka 02/06
c-----------------------------------------------------------------
      USE hdatam
      USE qload
      implicit double precision (a-h,o-z)
      logical aufl,alfl,dflg
      dimension u(*),b(*),au(*),al(*)

      hflgu  = .false. ! no update of h-array 
      h3flgu = .false.
      propq = propt

      if(mqloa.eq.1) return
c.... calculate load terms in s and p including load factor      
      call formfe(u,b,au,al,aufl,.true.,alfl,dflg,22,1,numel,1)
      return
      end
c
