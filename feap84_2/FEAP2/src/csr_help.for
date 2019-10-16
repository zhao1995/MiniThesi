c----------------------------------------------------------------------+
c     Subroutines for CSR Storage                                      |
c----------------------------------------------------------------------+
      subroutine mkptr_csr(ix,id,ia,nen1,ndf,is_csr)
c-----------------------------------------------------------------------
c
c      Purpose:    Set storage for compressed arrays  CSR
c
c      Inputs:     diag
c                  is_csr
c
c      Outputs:
c        ia=ir      - ia(neq+1): first entry in each row
c                     (=sum of nonzero elements of previous rows
c        ja=jc      - m(japt(ljacsr)): array of columns with nonzero
c                     entries for each row
c        ljacsr     - length of ja stored in ia(neq+1)
c        ka         = m(kapt(neq): array of diagonal entries
c
c      Comments:
c        Modifed version from FEAP82 SR iters and others
c
c     W. Wagner BS UKA 04/09
c     S. Klarmann  TUD 06/15 parallelization and modifactions 
c
c-----------------------------------------------------------------------
      USE cdata
      USE doalloc
      USE iofile
      USE iscsr
      USE pdata8
      USE plong
      USE soltyp
      implicit  double precision(a-h,o-z)
      logical ldummy
      logical kdiag,kall
      integer ia(*), ips
      integer, allocatable, dimension(:) :: csric, csrielc,csriwr
      integer, allocatable, dimension(:) :: cshja,cshia
      ldummy  = .true.

c.... compute ic:  1)holds element degree of each node
c                  2)becomes pointer for array that contains set of
c                    elements connected to each node.

      kdiag  = .true.
      call setup_csr(ix,id,ia,nen1,ndf,kdiag)

c.... store ja
      call ialloc(csrka,ljacsr,'MKPTR-CSR-ja',ldummy)

c.... set array ka and store adress of diagonal entries
      call setka(ia,csrja,csrka,neq)

      end
c
c-----------------------------------------------------------------------
      subroutine isort(ifeld,n)
c-----------------------------------------------------------------------
c
c     Purpose:  Sorts the entries in array ifeld into increasing order.
c
c     Inputs:
c        ifeld(n) - Array of integers to sort
c        n        - Number of entries to sort
c
c     Output:
c        ifeld(n) - Sorted array of integers
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
c
      implicit  none

      logical   flag
      integer   i,ihilf,n,k
      integer   ifeld(n)

c     Sort integer array using bubble sort algorithm

      flag = .true.
      k = 1
      do while( flag )
        flag = .false.
        do i = 1, n-k
          if (ifeld(i).gt.ifeld(i+1)) then
            flag       = .true.
            ihilf      = ifeld(i)
            ifeld(i)   = ifeld(i+1)
            ifeld(i+1) = ihilf
          end if
        end do ! i
        k = k+1
      end do ! while flag = true

      end
c
c-----------------------------------------------------------------------
      subroutine ioprof(ir,jc,neq)
c-----------------------------------------------------------------------
c
c      Purpose: Print ir,jc - only for test problems!
c
c      Inputs:
c         ir(*)  - Array of integer values to print
c         jc(*)  - Column pointers
c         neq    - Number of row/colums
c
c      Outputs:
c         To file and screen
c
c     W. Wagner BS UKA 04/09
c-----------------------------------------------------------------------
c
      USE iofile
      implicit   none

      integer neq,ir(*),jc(*)
      integer n,i

c     Print profile stored array by equation number

      write(iow,2000)
      if(ior.lt.0) write(*,2000)

      do n = 1,neq
        write(iow,2001) n, ir(n),(jc(i),i=ir(n),ir(n+1)-1)
        if(ior.lt.0) then
          write(*,2001) n, ir(n),(jc(i),i=ir(n),ir(n+1)-1)
        end if
      end do ! n

      write(iow,*) 'ia'

      do n=1,neq
        write(iow,2002) ir(n)
      end do

      write(iow,*) 'ja'

      do n = 1,ir(neq+1)-1
        write(iow,2003) jc(n)
      end do

c     Formats

2000  format(/'  S p a r s e    M a t r i x  : '/,
     +       'Row No. ',' ','1. Entry',' Associated columns')
2001  format(i8,':',i8,':',8i8:/(9x,8i8))
2002  format(i7)
2003  format(i7)

      end
c
c-----------------------------------------------------------------------
      subroutine setka(ia,ja,ka,neq)
c-----------------------------------------------------------------------
c
c     Purpose: calculate position of diagonals for CSR-storage
c
c     Inputs:
c         ia(neq+1)  - element number at entry position in each row
c         ja(nja)    - array of columns where nonzero entries occur
c         neq        - number of equations
c
c     Output:
c         ka(neq)    - element number at diagonal in each row
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical col
      dimension ia(neq+1),ja(*),ka(neq)
!$OMP  PARALLEL DEFAULT(NONE) 
!$OMP& PRIVATE(ii,ir,nr,l,col,ipos)
!$OMP& SHARED(neq,ia,ja,ka)   
!$OMP DO SCHEDULE(DYNAMIC) 
      do ii = 1,neq
         ir  = ia(ii)           ! first entry in row
         nr  = ia(ii+1)-ia(ii)  ! no of entries in row
c.....   find column
         l = 1
         col=.false.
10       if (ja(ir).eq.ii) then
           ipos=ir
           col=.true.
         end if
         ir=ir+1
         l=l+1
         if ((.not.col).and.(l.le.nr)) goto 10
c...     store value
         ka(ii) = ipos
      end do
!$OMP END DO
!$OMP END PARALLEL

      return
      end
c
c--------------------------------------------------------------
      subroutine dsolv_csr(n,x,b,a,ka)
c--------------------------------------------------------------
c
c....  Purpose: form x = b/a_ii in CSR-format
c
c      Inputs:
c        x(n)   -  vector
c        b(n)   -  vector
c        a      -  matrix stored in CSR-format
c        ka(n)  -  pointer to diagonal entries
c        n      -  dimenion of arrays
c
c      Outputs:
c        x(n)   -  vector
c
c--------------------------------------------------------------
      integer n,i,ka(*)
      double precision x(*),b(*),a(*)
!$OMP  PARALLEL DEFAULT(NONE) 
!$OMP& PRIVATE(i)
!$OMP& SHARED(x,b,a,ka,n)   
!$OMP DO SCHEDULE(DYNAMIC) 
      do i=1,n
        x(i)=b(i)/a(ka(i))
      end do
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c--------------------------------------------------------------
      subroutine dadd_csr(n,a,b,c,ka)
c--------------------------------------------------------------
c
c....  Purpose: form a_ii = a_ii + c*b_ii
c
c      Inputs:
c        a(*)   -  matrix stored in CSR-format
c        b(*)   -  diagonal matrix
c        ka(n)  -  pointer to diagonal entries
c        c      -  factor
c        n      -  dimension of arrays
c
c      Outputs:
c        a(*)   -  matrix stored in CSR-format
c
c--------------------------------------------------------------
      integer n,i,ka(*)
      double precision a(*),b(*),c
!$OMP  PARALLEL DEFAULT(NONE) 
!$OMP& PRIVATE(i,j)
!$OMP& SHARED(x,b,a,ka,n,c)   
!$OMP DO SCHEDULE(DYNAMIC) 
      do i=1,n
        j = ka(i)
        a(j) = a(j) + c*b(i)
      end do
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c--------------------------------------------------------------
      subroutine dnzero_csr(n,a,ka,nn,ityp,imas)
c--------------------------------------------------------------
c
c....  Purpose: count number of nonzero elements in a_ii in CSR-format
c
c      Inputs:
c        a      -  matrix stored in CSR-format
c        ka(n)  -  pointer to diagonal entries
c        n      -  dimenion of arrays
c        ityp   -  1 .ne. 0, 2 .lt.0
c        imas   -  1 = consist., 2 = diagonal
c
c      Outputs:
c        nn     -  dimenion of arrays
c
c--------------------------------------------------------------
      integer n,nn,i,ka(*)
      double precision a(*)

      nn = 0

      if(ityp.eq.1) then
        do i=1,n
          if(imas.eq.1) then
            if(a(ka(i)).ne.0.0d0) nn = nn+1
          else
            if(a(i)    .ne.0.0d0) nn = nn+1
          end if
        end do
      else if(ityp.eq.2) then
        do i=1,n
          if(imas.eq.1) then
            if(a(ka(i)).lt.0.0d0) nn = nn+1
          else
            if(a(i)    .lt.0.0d0) nn = nn+1
          end if
        end do
      end if
      return
      end
c
c--------------------------------------------------------------
      subroutine pick_csr(v,n,a,ka,imas)
c--------------------------------------------------------------
c
c....  Purpose: pick value of element n of a in CSR-format
c
c      Inputs:
c        a      -  matrix stored in CSR-format
c        ka(n)  -  pointer to diagonal entries
c        n      -  dimenion of arrays
c        imas   -  1 = consist., 2 = diagonal
c
c      Outputs:
c        nn     -  dimenion of arrays
c
c--------------------------------------------------------------
      integer n,ka(*)
      double precision a(*),v

        if(imas.eq.1) v = a(ka(n))
        if(imas.eq.2) v = a(n)

      return
      end
c
c--------------------------------------------------------------
      subroutine prof_csr(pfact,ptone,ptnin,x0,y0,ds,neq,ia,ja,lower)
c-----------------------------------------------------------------------
c
C     Purpose: SR for plotting profile for CSR-based solvers
c
c           neq    - number of equations
c           ix,iy  - number of row/column of eq.
c           x,y    - position of eq. in profil
c           ia,ja  - arrays for CSR storage
c
c
c-----------------------------------------------------------------------
      USE pdata1
      implicit double precision (a-h,o-z)
c
      logical lower
      real*8          x0,y0,x,y,pfact,ptone,ptnin,ds,pix,piy
      integer         ia(*),ja(*)

      do ix  = 1,neq                   ! row
        do k = ia(ix),ia(ix+1)-1      ! column
          iy = ja(k)
          pix = ix-0.5d0
          piy = iy-0.5d0
          x = ptone  + piy*pfact + x0
          y = ptnin  - pix*pfact + y0
c......   plot point
          call ppbox2(x,y,ds)
          if(lower) then
            if(ix.eq.iy) go to 10
            ixq = iy
            iyq = ix
            pix = ixq-0.5d0
            piy = iyq-0.5d0
            x = ptone  + piy*pfact + x0
            y = ptnin  - pix*pfact + y0
c......     plot point
            call ppbox2(x,y,ds)
10        continue
          end if
        end do ! k
      end do ! i

      return
      end
c
c--------------------------------------------------------------
      subroutine promul_csr(a,b,c,ia,ja,neq,add,isymcsr)
c--------------------------------------------------------------
c
c....  Purpose: form c = c +/- a*b
c
c      Inputs:
c        a        -  matrix stored in CSR-format
c        b(neq)   -  vector
c        c(neq)   -  vector
c        ia(neq)  -  pointer to first entries in rows
c        ja(*)    -  pointer to entries in columns
c        neq      -  dimenion of arrays
c        add      -  true c = c+a*b   false c = c-a*b
c
c      Outputs:
c        c(neq)   -  vector
c
c        (c)  ww 5/2006
c--------------------------------------------------------------
c
      logical add
      double precision a(*),b(*),c(*)
      integer neq,ia(*),ja(*),i,j,k
c
!$OMP  PARALLEL DEFAULT(NONE) 
!$OMP& PRIVATE(i,j,k)
!$OMP& SHARED(a,b,c,add,neq,ia,ja,isymcsr)   
      if(isymcsr.eq.1) then  ! symmetry
c....   loop over all rows
!$OMP SINGLE
        do i=1,neq,1
c....     loop over all entries in row, except first=diagonal
          do j=ia(i)+1,ia(i+1)-1
c....       column index
            k=ja(j)
c....       upper + lower part
            if(add) then
              c(i) = c(i) + a(j)*b(k)
              c(k) = c(k) + a(j)*b(i)
            else
              c(i) = c(i) - a(j)*b(k)
              c(k) = c(k) - a(j)*b(i)
            end if
          end do
        end do
!$OMP END SINGLE  
c....   diagonal part
!$OMP DO SCHEDULE(DYNAMIC) 
        do i=1,neq
          if(add) then
            c(i) = c(i) + a(ia(i))*b(i)
          else
            c(i) = c(i) - a(ia(i))*b(i)
          end if
        end do
!$OMP END DO
      else if(isymcsr.eq.2) then ! unsymmetry
!$OMP DO SCHEDULE(DYNAMIC) 
        do i = 1,neq
c....     compute the inner product of row i with vector x
          do j=ia(i), ia(i+1)-1
            k=ja(j)
            if(add) then
              c(i) = c(i) + a(j)*b(k)
            else
              c(i) = c(i) - a(j)*b(k)
            end if
          end do
        end do
!$OMP END DO
      end if
!$OMP END PARALLEL
      return
      end
c
c--------------------------------------------------------------
      subroutine dasbl_csr (s,p,ld,ia,ns,aufl,bfl,b,ad,ja)
c-----------------------------------------------------------------------
c
c      Purpose: Assemble full arrays for CSR storage technique
c
c      Inputs:
c         s(ns,ns) - Element array to assemble
c         p(ns)    - Element vector to assemble
c         ld(ns)   - Local to Global equation numbers for assemble
c         ia(*)    - Pointer to row    CSR sparse technique
c         ja(*)    - Pointer to column CSR
c         ns       - Size of element arrays
c         aufl     - If true, assemble A array
c         bfl      - If true, assemble B vector
c
c      Outputs:
c         b(*)     - Assembled right hand side B vector
c         ad(*)    - Assembled A array
c
c     (c) ww 5/2006
c-----------------------------------------------------------------------
      USE cdata
      USE iscsr
      USE soltyp
      implicit double precision (a-h,o-z)
      logical aufl,bfl,col
      dimension ad(*),b(*),s(ns,ns),p(*),ld(*)
      dimension ia(*),ja(*)
      
c.... loop through the rows to perform the assembly
      do i=1,ns
        ii=ld(i)
        if (ii.gt.0) then
          if (aufl) then
            if(isymcsr.eq.1) jjs = ii ! symmetry
            if(isymcsr.eq.2) jjs = 1  ! full
c....       loop through the columns to perform the assembly
            do j=1,ns
              jj=ld(j)
              if (jj.ge.jjs) then
                do ipos=ia(ii),ia(ii+1)-1
                  if(ja(ipos).eq.jj) then
                    !$OMP ATOMIC
                    ad(ipos)=ad(ipos)+s(i,j) !  OMP ATOMIC
                    exit
                  end if
                end do
              end if
            end do !(j)
          end if
          if (bfl) then
            !$OMP ATOMIC
            b(ii)=b(ii)+p(i) ! OMP ATOMIC
          end if
        end if
      end do !(i)

      return
      end
c
c--------------------------------------------------------------
      subroutine amux (n, x, y, a,ja,ia)
c-----------------------------------------------------------------------
c     y = A * x
c-----------------------------------------------------------------------
c     Purpose: multiplies a full matrix by a vector using
c              the dot product form. Matrix A is stored in CSR format.
c
c     Inputs:
c      n     = row dimension of A
c      x     = real array of length equal to the column dimension of
c              the A matrix.
c      a, ja,ia = input matrix in CSR format.
c
c     Outputs:
c      y     = real array of length n, containing the product y=Ax
c
c     Comments:
c       This SR is similar to promul2 execpt c = c+-A*b
c
c-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      real*8  x(*), y(*), a(*)
      integer n, ja(*), ia(*),i,k

      call pzero(y,n)
!$OMP  PARALLEL DEFAULT(NONE) 
!$OMP& PRIVATE(i,k)
!$OMP& SHARED(y,a,x,ja,n)   
!$OMP DO SCHEDULE(DYNAMIC) 
      do i = 1,n
c....    compute the inner product of row i with vector x
         do k=ia(i), ia(i+1)-1
            y(i) = y(i) + a(k)*x(ja(k))
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c--------------------------------------------------------------
      subroutine atmux (n, x, y, a,ja,ia)
c-----------------------------------------------------------------------
c     y = A^T * x
c-----------------------------------------------------------------------
c     Purpose: multiplies a transpose of a full matrix by a
c              vector using the dot product form. Matrix A is stored in
c             CSR format.
c
c     Inputs:
c      n     = row dimension of A
c      x     = real array of length equal to the column dimension of
c              the A matrix.
c      a, ja,ia = input matrix in CSR format.
c
c     Outputs:
c      y     = real array of length n, containing the product y=Ax
c
c
c     Comments:
c       This SR is similar to promul2 execpt c = c+-A*b
c
c-----------------------------------------------------------------------

      real*8  x(*), y(*), a(*)
      integer n, ja(*), ia(*),i,k

      call pzero(y,n)

      do i = 1,n
c....   compute the inner product of row i with vector x
        do k=ia(i), ia(i+1)-1
          y(ja(k)) = y(ja(k)) + x(i)*a(k)
        end do
      end do
      return
      end
c
c-----------------------------------------------------------------------
      subroutine setup_csr(ix,id,ia,nen1,ndf,wdiag)
c-----------------------------------------------------------------------
c
c      Purpose:    Set storage for compressed arrays  CSR
c
c      Inputs:     diag       - always true
c                  is_csr     - un-/symmetric
c                  numnp      - number of nodes
c                  numel      - number of elements
c                  nen        - max number of nodes per element
c                  nen1       - entries in ix per element
c                  id         - eq. numbers
c                  ix         - element connection list
c
c      Outputs:
c        ia         - ia(neq+1): first entry in each row
c                     (=sum of nonzero elements of previous rows
c        ja         - csrja: array of columns with nonzero
c                     entries for each row
c        ljacsr     - length of ja stored in ia(neq+1)
c
c      Comments:
c        Parallelized version
c        in parallel mode each process processes one equation
c      
c        First two blocks must run sequential
c
c     S. Klarmann TU Darmstadt 6/2015
c
c-----------------------------------------------------------------------    
      USE cdata
      USE doalloc
      USE iofile
      USE iscsr
      USE pdata8
      USE plong
      USE soltyp
      implicit none
      
      integer, allocatable, dimension(:) :: csric, csrielc, csriwr
      integer, allocatable, dimension(:) :: iel_ja, iel_ia
      integer, allocatable, dimension(:) :: csrja_help,csria_help
      integer, allocatable, dimension(:) :: nonzero_count
      integer id(ndf,numnp), ix(nen1,numel), ia(neq+1)
      integer nen1, ndf, is_csr
      integer max1, kp
      logical kall, search,all,wdiag,addeq,addeq2,ldummy
      integer i,j,k,kk,l,m,n, iel_start,iel_end,node_elem,ipos
      integer num_elem,eq_count,jpos
      integer neqi, neqj
      
c.... 1= symmetry(def), 2=un-symmetry
      if(isymcsr.eq.1) all = .false.
      if(isymcsr.eq.2) all = .true.
      
c.... compute row pointer for eq element connection list
      allocate(iel_ia(neq+1))
      iel_ia = 0
      iel_ia(1) = 1
      do i=1,numel
        do j=1,nen
          do k=1,ndf
            if(ix(j,i).gt.0) then
              neqi = id(k,ix(j,i))
              if(neqi.gt.0) iel_ia(neqi+1)=iel_ia(neqi+1)+1
            end if
          end do
        end do
      end do
      do i=2,neq+1
        iel_ia(i)=iel_ia(i)+iel_ia(i-1)
      end do
c.... compute column pointer for node element connection list
      allocate(iel_ja(iel_ia(neq+1)-1))
      iel_ja = 0
      do i=1,numel
        do j=1,nen
          kk=ix(j,i)
          if(kk.gt.0) then
            do k=1,ndf
              neqi=id(k,kk)
              if(neqi.gt.0) then
                l=iel_ia(neqi)
                search=.true.
                do while(search)
                  if(iel_ja(l).eq.0) then
                    iel_ja(l)=i
                    search = .false.
                  end if
                  l=l+1
                end do
              end if
            end do
          end if
        end do
      end do
c.... max space for ja
c.... max connections: numelem*nodes_per_elem*ndf
      allocate(csria_help(neq+1))
      csria_help   =0
      csria_help(1)=1
!$OMP  PARALLEL DEFAULT(NONE) 
!$OMP& PRIVATE(i,j,kk,iel_end,iel_start)
!$OMP& SHARED(neq,iel_ia,id,csria_help,nen,ndf)
!$OMP DO SCHEDULE(DYNAMIC)       
      do i=1,neq
        iel_start=iel_ia(i)
        iel_end  =iel_ia(i+1)
        csria_help(i+1)=(iel_end-iel_start)*nen*ndf
      end do ! i = numnp
!$OMP END DO
!$OMP END PARALLEL
c.... calculate start positions of row in ja
      do i=2,neq+1
        csria_help(i)=csria_help(i)+csria_help(i-1)
      end do
c.... allocate temporary ja field
      allocate(csrja_help(csria_help(neq+1)-1))
      csrja_help=0
      
c.... set up column pointers
      allocate(nonzero_count(neq))
      nonzero_count = 0
!$OMP  PARALLEL DEFAULT(NONE) 
!$OMP& PRIVATE(neqi,iel_start,iel_end,j,num_elem,k,node_elem,l
!$OMP&        ,neqj,kk,addeq,addeq2,search)
!$OMP& SHARED(neq,iel_ia,iel_ja,nen,ix,id,ndf,csria_help,all,wdiag
!$OMP&       ,csrja_help,nonzero_count)
!$OMP DO SCHEDULE(DYNAMIC)           
      do neqi=1,neq
        iel_start=iel_ia(neqi)
        iel_end  =iel_ia(neqi+1)
        do j=iel_start,iel_end-1
          num_elem = iel_ja(j)
          do k=1,nen
            node_elem = ix(k,num_elem)
            if(node_elem.gt.0) then ! element node exists
            do l=1,ndf
              neqj = id(l,node_elem)
              if(neqj.gt.0) then
                if(all) then ! structure of storage
                  addeq = neqj.gt.0
                else
                  if(wdiag) then
                    addeq = neqj.ge.neqi
                  else
                    addeq = neqj.gt.neqi
                  end if
                end if
                if(addeq) then
                  addeq2 = .true.
                  search = .true.
                  kk = csria_help(neqi)
                  do while(search) !.and.kk.le.(csria_help(neqi+1)-1))
                    if(csrja_help(kk).eq.neqj) then
                      search = .false.
                      addeq2 = .false.
                    else if(csrja_help(kk).eq.0) then
                      search = .false.
                    else
                      kk = kk+1
                    end if
                  end do
                  if(addeq2) then
                    csrja_help(kk) = neqj
                    nonzero_count(neqi)=nonzero_count(neqi)+1
                  end if
                end if
              end if !eq exists
            end do !l
            end if ! element node exists
          end do !k
        end do !j
      end do !neqi
!$OMP END DO
!$OMP END PARALLEL      

c...  compressing ja 
c...  set up new row pointer for compressed ja
      ia=0
      ia(1)=1
      do i=2,neq+1
        ia(i)=ia(i-1)+nonzero_count(i-1)
      end do
      
c...  copy non-zero values in new ja field
      ljacsr= ia(neq+1)-1
      call ialloc(csrja,ljacsr,'MKPTR-CSR-jat',ldummy)
!$OMP  PARALLEL DEFAULT(NONE) 
!$OMP& PRIVATE(i,j,k,kk,l)
!$OMP& SHARED(csrja,ia,neq,csria_help,csrja_help)
!$OMP DO SCHEDULE(DYNAMIC)       
      do i=1,neq
        j = csria_help(i)-1
        k = ia(i)-1
        kk=ia(i+1)-ia(i)
        do l=1,kk
          csrja(k+l)=csrja_help(j+l)
        end do
c...  sorting values in ja in ascending order
        call isort(csrja(k+1),kk)
      end do
!$OMP END DO
!$OMP END PARALLEL      

      return
      end
