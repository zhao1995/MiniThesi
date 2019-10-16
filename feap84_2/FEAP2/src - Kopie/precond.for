      subroutine precond(a,ia,ja,ka)
c-----------------------------------------------------------------------
c
c      Purpose: Preconditioner for iterative solvers
c
c      Inputs:
c        a(*)    - Matrix A (non-factored) in CSR storage
c        ia(*)   - pointer to row    of A
c        ja(*)   - pointer to column of A
c        ka(*)   - pointer to diagonals of A
c        neq     - Number of equations
c
c
c      Output:
c        alu(*)  - Matrix ALU in MSR storage
c        jlu(*)  - pointer to non-diagonals of ALU
c        ju(*)   - pointer to column
c
c
c      Comments
c         ip     - 1: Unit matrix
c         ip     - 2: 1/Diagonals of A
c         ip     - 3: Incomplete LU factorization of A
c
c-----------------------------------------------------------------------
      USE cdata
      USE isprec
      USE doalloc
      implicit double precision (a-h,o-z)
      dimension ia(neq+1),ja(*),ka(*)


c.... set parameter for ilux default, reset via macro [prec,,n1,n2,n3]
c     ippc  = 3(in SR ini)      ! preconditioner
c     lfil  = 5(in SR ini)      ! fill-in parameter >0 <=neq,
c                               ! ilut = additional elements in row, recommendation 5-10
c                               ! ilup = total elements in row + diag element, recommendation ?
      lfil  = min(lfil,neq)     ! ilut
c     tolpc = 1.e-3 (in SR ini) ! Thresholding in L and U tolerance, recommendation 1-5 e-3

      ipcwk =  ia(neq+1)        ! maximum length of ALU = ia(neq+1)-1 = nzz-1
                                ! maximum length of JLU = nzz, see description in ILUT
      ipcwk = 10*ipcwk          ! tyipcal


c.....Gauss-Solver as special case of ilut
c      lfil     =  neq
c      tolpc   =  0.0

c.... define necessary arrays
      if (lpreco) then                                        ! conditioner use
        call ralloc( rmpcalu,ipcwk,'PRECO- mpcalu',lpreco)  ! 1,2,3, ,5    ,6
        call ialloc( impcjlu,ipcwk,'PRECO- mpcjlu',lpreco)  !  , ,3, ,5    ,6
        call ialloc(  impcju,  neq,'PRECO-  mpcju',lpreco)  !  , ,3, ,5    ,6
        call ialloc( impcjwl,  neq,'PRECO- mpcjwl',lpreco)  !  , ,3, ,5a=jw,6a=jw
        call ialloc( impcjwu,  neq,'PRECO- mpcjwu',lpreco)  !  , ,3, ,5b   ,6b
        call ialloc(  impcjr,  neq,'PRECO-  mpcjr',lpreco)  !  , ,3, ,5c   ,-
        call ralloc(  rmpcwl,  neq,'PRECO-  mpcwl',lpreco)  !  , ,3, ,5=wl ,6
        call ralloc(  rmpcwu,neq+1,'PRECO-  mpcwu',lpreco)  !  , ,3, ,     ,-
        call ialloc(impclevs,ipcwk,'PRECO-mpclevs',lpreco)  !  , , , ,5    ,6a=iperm(2*neq)
      end if


      call precond1(a,ia,ja,ka,rmpcalu,impcjlu,impcju,neq,ippc,
     +  lfil,tolpc,ipcwk,rmpcwu,impcwl,impcjr,rmpcjwl,impcjwu,
     +  impclevs)

      return
      end


      subroutine precond1(a,ia,ja,ka,alu,jlu,ju,neq,ippc,lfil,tolpc,
     +                    ipcwk,wu,wl,jr,jwl,jwu,levs)
c----------------------------------------------------------------------
c
c     Purpose: Preconditioner for iterative solvers
c
c      Inputs:
c        a(*)    - Matrix A (non-factored) in CSR storage
c        ia(*)   - pointer to row    of A
c        ja(*)   - pointer to column of A
c        ka(*)   - pointer to diagonals of A
c        alu(*)  - Matrix ALU in MSR storage
c        jlu(*)  - pointer to non-diagonals of ALU
c        ju(*)   - pointer to column
c        neq     - Number of equations
c        ippc    - Type of preconditioner
c
c      Outputs:
c        alu(*)  - Matrix ALU in MSR storage
c        jlu(*)  - pointer to non-diagonals of ALU
c        ju(*)   - pointer to column of ALU
c
c         a(*)       - Non factored terms of A
c         neq        - Number of equations
c
c
c      Comments
c         ippc   - 1: Unit matrix                                                          ILU1
c         ippc   - 2: 1/Diagonals of A                                                     ILU2
c         ippc   - 3: Incomplete LU factorization of A with dropout and without pivoting   ILUT
c         ippc   - 4: Incomplete LU factorization of A with dropout and with    pivoting   ILUTP
c         ippc   - 5: Incomplete LU factorization of A with level k-fill in                ILUK
c
c----------------------------------------------------------------------
      USE iscsr
      implicit double precision (a-h,o-z)
      character*5 ip(5)
      dimension a(*),ia(*),ja(*),alu(*),jlu(*),ju(*),wu(*),wl(*),jr(*),
     +          jwl(*),jwu(*),levs(*)
      data ip /'IDEN ','DIAG ','ILUT ','ILUTP','ILUK'/

      if(ippc.eq.1 ) then          ! [prec,,1]
c       Unit matrix preconditioner
        call ilu1(neq,alu)
        ju0 = neq

      else if(ippc.eq.2) then      ! [prec,,2]
c       Diagonal of A preconditioner
        call ilu2(neq,alu,a,ka)
        ju0 = neq

      else if(ippc.eq.3) then       ! [prec,,3,tolpc,lfil] tolpc=1.e-3,lfil = 5-15
c       Incomplete LU factorization of A with dropout and without pivoting

        call ilut(neq,a,ja,ia,lfil,tolpc,alu,jlu,ju,ipcwk,wu,wl,jr,jwl,
     +            jwu,ju0,ierr)

        if(ierr.ne.0 ) go to 10

      else if(ippc.eq.4) then       ! [prec,,4,tolpc,lfil] tolpc=1.e-3,lfil = 5-15
c       Incomplete LU factorization of A with dropout and with    pivoting

c       permtol = 0.01 ! recommendation 0.1 - 0.01
        permtol = tolpc
        mbloc   = neq
c       array iperm   = array levs

        call ilutp(neq,a,ja,ia,lfil,tolpc,permtol,mbloc,alu,jlu,ju,
     +       ipcwk,wl,jwl,levs,ju0,ierr)

c.....differences in iperm
      do k = 1,neq
        idiff = levs(k)-levs(k+neq)
        if(idiff.ne.0) then
          write(*,*) 'difference ',k
          stop
        end if
      end do

c.... set A in original state
c      do k=ia(1), ia(neq+1)-1
c         ja(k) = levs(ja(k))
c      end do

        if(ierr.ne.0 ) go to 10

      else if(ippc.eq.5) then       ! [prec,,5,,level-k] level-k = 2-3
c       Incomplete LU factorization of A with level k-fill in

        call iluk(neq,a,ja,ia,lfil,alu,jlu,ju,levs,ipcwk,wl,jwl,ju0,
     +            ierr)

        if(ierr.ne.0 ) go to 10

      end if

        factju01 = float(ju0)/ia(neq+1)

c        write(*,*) 'work space factor ALU/A for ',ip(ippc), factju01
c        write(*,*) 'tolpc,lfil ',tolpc,lfil

      return

c.... errors
10    continue
c     if(ierr.eq. 0) write(*,*) ip(ippc),':  successful return.'
      if(ierr.gt. 0) write(*,*) ip(ippc),':  zero pivot at step number',
     +                          ierr
      if(ierr.eq.-1) write(*,*) ip(ippc),':  input matrix may be wrong.'
     + ,'(The elimination process has generated a row in L or U whose ',
     +  ' length is >  n.)'
      if(ierr.eq.-2) write(*,*) ip(ippc),':  matrix L overflows the ',
     +                          'array alu.'
      if(ierr.eq.-3) write(*,*) ip(ippc),':  matrix U overflows the ',
     +                          'array alu.'
      if(ierr.eq.-4) write(*,*) ip(ippc),':  Illegal value for lfil.'


      if(ierr.eq.-5.and.ippc.eq.3)
     +           write(*,*) ip(ippc),':  zero pivot encountered.'
      if(ierr.eq.-5.and.ippc.eq.4.or.ippc.eq.5)
     +           write(*,*) ip(ippc),': zero row encountered in A or U:'

      stop 'SR precond1'

      end
c

      subroutine ilu1(neq,ap)
c----------------------------------------------------------------------
c
c      Purpose: Unit matrix as Preconditioner for the PGMRES solver
c
c      Inputs:
c         neq    - Number of equations
c
c      Outputs:
c         ap(*)  - matrix M=ALU stored in MCSR format
c
c----------------------------------------------------------------------

      integer i,neq
      real*8 ap(*)
      do i =1,neq
        ap(i) = 1.d0
      end do
      return
      end

      subroutine ilu2(neq,ap,a,ka)
c----------------------------------------------------------------------
c
c      Purpose: Inverse of Diagonals of A as Preconditioner for the
c               PGMRES solver
c
c      Inputs:
c         a(*)   - matrix A stored in CSR format
c         ka(*)  - pointer to diagonals of matrix A stored in CSR format
c         neq    - Number of equations
c
c      Outputs:
c         ap(*)  - matrix M=ALU stored in MCSR format
c
c----------------------------------------------------------------------
      integer i,neq,ka(*)
      real*8 a(*),ap(*)
      do i =1,neq
        ap(i) = 1.d0/a(ka(i))
      end do
      return
      end


      subroutine ilut ( n, a, ja, ia, lfil, tol, alu, jlu, ju, iwk, wu,
     +                  wl, jr, jwl, jwu, ju0, ierr )
      !*******************************************************************************
      !
      !! ILUT is an ILUT preconditioner.
      !
      !  Discussion:
      !
      !    This routine carries out incomplete LU factorization with dual
      !    truncation mechanism.  Sorting is done for both L and U.
      !
      !    The dual drop-off strategy works as follows:
      !
      !    1) Thresholding in L and U as set by TOL.  Any element whose size
      !       is less than some tolerance (relative to the norm of current
      !       row in u) is dropped.
      !
      !    2) Keeping only the largest lenl0+lfil elements in L and the
      !       largest lenu0+lfil elements in U, where lenl0=initial number
      !       of nonzero elements in a given row of lower part of A
      !       and lenlu0 is similarly defined.
      !
      !    Flexibility: one can use tol=0 to get a strategy based on keeping the
      !    largest elements in each row of L and U. Taking tol /= 0 but lfil=n
      !    will give the usual threshold strategy (however, fill-in is then
      !    unpredictible).
      !
      !    A must have all nonzero diagonal elements.
      !
      !  Modified:
      !
      !    07 January 2004
      !
      !  Author:
      !
      !    Youcef Saad
      !
      !  Parameters:
      !
      !    Input, integer N, the order of the matrix.
      !
      !    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
      !    Compressed Sparse Row format.
      !
      ! lfil    = integer. The fill-in parameter. Each row of L and
      !           each row of U will have a maximum of lfil elements
      !           in addition to the original number of nonzero elements.
      !           Thus storage can be determined beforehand.
      !           lfil must be >= 0.
      !
      ! iwk     = integer. The minimum length of arrays alu and jlu
      !
      ! On return:
      !
      ! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
      !           the L and U factors together. The diagonal (stored in
      !           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
      !           contains the i-th row of L (excluding the diagonal entry=1)
      !           followed by the i-th row of U.
      !
      ! ju      = integer array of length n containing the pointers to
      !           the beginning of each row of U in the matrix alu,jlu.
      !
      ! ierr    = integer. Error message with the following meaning.
      !           ierr  = 0    --> successful return.
      !           ierr  > 0    --> zero pivot encountered at step number ierr.
      !           ierr  = -1   --> Error. input matrix may be wrong.
      !                            (The elimination process has generated a
      !                            row in L or U whose length is >  n.)
      !           ierr  = -2   --> The matrix L overflows the array al.
      !           ierr  = -3   --> The matrix U overflows the array alu.
      !           ierr  = -4   --> Illegal value for lfil.
      !           ierr  = -5   --> zero pivot encountered.
      !
      ! work arrays:
      !
      ! jr,jwu,jwl, integer work arrays of length n.
      ! wu, wl, real work arrays of length n+1, and n resp.
      !
cww---------------------------------------------------------------------
      !  CSR        ->   MSR storage technique
      !  a, ja, ia  ->   ao, jao
      !
      ! ao (1:n)     contains the diagonal of the matrix.
      ! ao( n+2:nnz) contains the nondiagonal elements of the  matrix, stored rowwise.
      ! jao(n+2:nnz) contains their column indices
      ! jao(1:n+1)   contains the pointer array for the nondiagonal
      !              elements in ao(n+1:nnz) and jao(n+2:nnz).
      !              i.e., for i <= n+1 jao(i) points to beginning of row i
      !              in arrays ao, jao.
      ! nnz        = number of nonzero elements+1
cww---------------------------------------------------------------------




        implicit none

        integer n

        real*8 a(*)
        real*8 alu(*)
        real*8 fact
        integer ia(n+1)
        integer idiag
        integer ierr
        integer ii
        integer iwk
        integer j
        integer j1
        integer j2
        integer ja(*)
        integer jj
        integer jlu(*)
        integer jpos
        integer jr(*)
        integer jrow
        integer ju(*)
        integer ju0
        integer jwl(n)
        integer jwu(n)
        integer k
        integer len
        integer lenl
        integer lenl0
        integer lenu
        integer lenu0
        integer lfil
        integer nl
        real*8 s
        real*8 t
        real*8 tnorm
        real*8 tol
        real*8 wl(n)
cww     real*8 wu(n)
        real*8 wu(n+1)


        if ( lfil < 0 ) then
          ierr = -4
          return
        end if
      !
      !  Initialize JU0 (points to next element to be added to ALU, JLU)
      !  and pointer.
      !
        ju0 = n + 2

        jlu(1) = ju0
      !
      !  Integer double pointer array.
      !
        jr(1:n) = 0
      !
      !  The main loop.
      !
        do ii = 1, n

          j1 = ia(ii)
          j2 = ia(ii+1) - 1
          lenu = 0
          lenl = 0

          tnorm = 0.0D+00
          do k = j1, j2
            tnorm = tnorm + abs ( a(k) )
          end do
          tnorm = tnorm / real(j2-j1+1)
      !
      !  Unpack L-part and U-part of row of A in arrays WL, WU.
      !
          do j = j1, j2

            k = ja(j)
            t = a(j)

            if ( tol * tnorm <= abs ( t ) ) then

              if ( k < ii ) then
                lenl = lenl + 1
                jwl(lenl) = k
                wl(lenl) = t
                jr(k) = lenl
              else
                lenu = lenu+1
                jwu(lenu) = k
                wu(lenu) = t
                jr(k) = lenu
              end if

            end if

          end do

          lenl0 = lenl
          lenu0 = lenu
          jj = 0
          nl = 0
      !
      !  Eliminate previous rows.
      !
  150 continue

          jj = jj + 1

          if ( lenl < jj ) then
            go to 160
          end if
      !
      !  In order to do the elimination in the correct order we need to
      !  exchange the current row number with the one that has
      !  smallest column number, among JJ, JJ+1, ..., LENL.
      !
          jrow = jwl(jj)
          k = jj
      !
      !  Determine the smallest column index.
      !
          do j = jj+1, lenl
             if ( jwl(j) < jrow ) then
                jrow = jwl(j)
                k = j
             end if
          end do
      !
      !  Exchange in JWL.
      !
          j = jwl(jj)
          jwl(jj) = jrow
          jwl(k) = j
      !
      !  Exchange in JR.
      !
          jr(jrow) = jj
          jr(j) = k
      !
      !  Exchange in WL.
      !
          s = wl(k)
          wl(k) = wl(jj)
          wl(jj) = s

          if ( ii <= jrow ) then
            go to 160
          end if
      !
      !  Get the multiplier for row to be eliminated: JROW.
      !
          fact = wl(jj) * alu(jrow)
          jr(jrow) = 0

          if ( abs ( fact ) * wu(n+2-jrow) <= tol * tnorm ) then
            go to 150
          end if
      !
      !  Combine current row and row JROW.
      !
          do k = ju(jrow), jlu(jrow+1)-1
             s = fact * alu(k)
             j = jlu(k)
             jpos = jr(j)
      !
      !  If fill-in element and small disregard.
      !
             if ( abs ( s ) < tol * tnorm .and. jpos == 0 ) then
               cycle
             end if

             if ( ii <= j ) then
      !
      !  Dealing with upper part.
      !
                if ( jpos == 0 ) then
      !
      !  This is a fill-in element.
      !
                   lenu = lenu + 1

                   if ( n < lenu ) then
                     go to 995
                   end if

                   jwu(lenu) = j
                   jr(j) = lenu
                   wu(lenu) = - s
                else
      !
      !  No fill-in element.
      !
                   wu(jpos) = wu(jpos) - s
                end if
             else
      !
      !  Dealing with lower part.
      !
                if ( jpos == 0 ) then
      !
      !  This is a fill-in element.
      !
                   lenl = lenl + 1

                   if ( n < lenl ) then
                     go to 995
                   end if

                   jwl(lenl) = j
                   jr(j) = lenl
                   wl(lenl) = -s
                else
      !
      !  No fill-in element.
      !
                   wl(jpos) = wl(jpos) - s
                end if
             end if

        end do

          nl = nl + 1
          wl(nl) = fact
          jwl(nl) = jrow
        go to 150
      !
      !  Update the L matrix.
      !
  160 continue

          len = min ( nl, lenl0 + lfil )

          call bsort2 ( wl, jwl, nl, len )

          do k = 1, len

             if ( iwk < ju0 ) then
               ierr = -2
               write(*,*) 'storage for alu(ierr=-2) exceeded, ju0=',ju0
cww            return
             end if


             alu(ju0) =  wl(k)
             jlu(ju0) =  jwl(k)
             ju0 = ju0 + 1

          end do
      !
      !  Save pointer to beginning of row II of U.
      !
          ju(ii) = ju0

      !
      !  Reset double pointer JR to zero (L-part - except first
      !  JJ-1 elements which have already been reset).
      !
        do k = jj, lenl
          jr(jwl(k)) = 0
        end do
      !
      !  Be sure that the diagonal element is first in W and JW.
      !
          idiag = jr(ii)

          if ( idiag == 0 ) then
            go to 900
          end if

          if ( idiag /= 1 ) then

             s = wu(1)
             wu(j) = wu(idiag)
             wu(idiag) = s

             j = jwu(1)
             jwu(1) = jwu(idiag)
             jwu(idiag) = j

          end if

          len = min ( lenu, lenu0 + lfil )

          call bsort2 ( wu(2), jwu(2), lenu-1, len )
      !
      ! Update the U-matrix.
      !
          t = 0.0D+00

          do k = 2, len

             if ( iwk < ju0 ) then
               ierr = -3
               write(*,*) 'storage for alu(ierr=-3) exceeded, ju0=',ju0
cww            return
             end if

             jlu(ju0) = jwu(k)
             alu(ju0) = wu(k)
             t = t + abs ( wu(k) )
             ju0 = ju0 + 1

          end do
      !
      !  Save norm in WU (backwards). Norm is in fact average absolute value.
      !
          wu(n+2-ii) = t / real (len + 1)
      !
      !  Store inverse of diagonal element of U.
      !
          if ( wu(1) == 0.0D+00 ) then
            ierr = -5
            return
          end if

          alu(ii) = 1.0D+00 / wu(1)
      !
      !  Update pointer to beginning of next row of U.
      !
        jlu(ii+1) = ju0
      !
      !  Reset double pointer JR to zero (U-part).
      !
        do k = 1, lenu
          jr(jwu(k)) = 0
        end do

        end do

        ierr = 0

cww        write(*,*) 'storage used in ILUT', ju0

        return
      !
      !  Zero pivot :
      !
  900    ierr = ii
          return
      !
      !  Incomprehensible error. Matrix must be wrong.
      !
  995    ierr = -1
          return
      end

      subroutine bsort2 ( w, ind, n, ncut )
      !*******************************************************************************
      !
      !! BSORT2 returns the NCUT largest elements of an array, using bubble sort.
      !
      !  Discussion:
      !
      !    This routine carries out a simple bubble sort for getting the NCUT largest
      !    elements in modulus, in array W.  IND is sorted accordingly.
      !    (Ought to be replaced by a more efficient sort especially
      !    if NCUT is not that small).
      !
      !  Modified:
      !
      !    07 January 2004
      !
      !  Author:
      !
      !    Youcef Saad
      !
      !  Parameters:
      !
        implicit none

        integer n

        integer i
        integer ind(*)
        integer iswp
        integer j
        integer ncut
        logical test
        real*8 w(n)
        real*8 wswp

        i = 1

        do

          test = .false.

          do j = n-1, i, -1

            if ( abs ( w(j) ) < abs ( w(j+1) ) ) then
      !
      !  Swap.
      !
              wswp = w(j)
              w(j) = w(j+1)
              w(j+1) = wswp
      !
      !  Reorder the original ind array accordingly.
      !
              iswp = ind(j)
              ind(j) = ind(j+1)
              ind(j+1) = iswp
      !
      !  Set indicator that sequence is still unsorted.
      !
              test = .true.

            end if

          end do

          i = i + 1

          if ( .not. test .or. ncut < i ) then
            exit
          end if

        end do

        return
      end

      subroutine ilutp(n,a,ja,ia,lfil,droptol,permtol,mbloc,alu,jlu,ju,
     + iwk,w,jw,iperm,ju0,ierr)
c----------------------------------------------------------------------*
c       *** ILUTP preconditioner -- ILUT with pivoting  ***            *
c      incomplete LU factorization with dual truncation mechanism      *
c----------------------------------------------------------------------*
c     author Yousef Saad *Sep 8, 1993 -- Latest revision, August 1996. *
c----------------------------------------------------------------------*
c     on entry:
c
c     n       = integer. The dimension of the matrix A.
c
c     a,ja,ia = matrix stored in Compressed Sparse Row format.
c               ON RETURN THE COLUMNS OF A ARE PERMUTED. SEE BELOW FOR
c               DETAILS.
c
c     lfil    = integer. The fill-in parameter. Each row of L and each row
c               of U will have a maximum of lfil elements (excluding the
c               diagonal element). lfil must be .ge. 0.
c               ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
c               EARLIER VERSIONS.
c
c     droptol = real*8. Sets the threshold for dropping small terms in the
c               factorization. See below for details on dropping strategy.
c
c     lfil    = integer. The fill-in parameter. Each row of L and
c               each row of U will have a maximum of lfil elements.
c               WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
c               EARLIER VERSIONS.
c               lfil must be .ge. 0.
c
c     permtol = tolerance ratio used to  determine whether or not to permute
c               two columns.  At step i columns i and j are permuted when
c
c               abs(a(i,j))*permtol .gt. abs(a(i,i))
c
c               [0 --> never permute; good values 0.1 to 0.01]
c
c     mbloc   = if desired, permuting can be done only within the diagonal
c               blocks of size mbloc. Useful for PDE problems with several
c               degrees of freedom.. If feature not wanted take mbloc=n.
c
c
c     iwk     = integer. The lengths of arrays alu and jlu. If the arrays
c               are not big enough to store the ILU factorizations, ilutp
c               will stop with an error message.
c
c     On return:
c
c     alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c              the L and U factors together. The diagonal (stored in
c              alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c              contains the i-th row of L (excluding the diagonal entry=1)
c              followed by the i-th row of U.
c
c     ju      = integer array of length n containing the pointers to
c               the beginning of each row of U in the matrix alu,jlu.
c
c     iperm   = contains the permutation arrays.
c               iperm(1:n) = old numbers of unknowns
c               iperm(n+1:2*n) = reverse permutation = new unknowns.
c
c     ierr    = integer. Error message with the following meaning.
c               ierr  = 0    --> successful return.
c               ierr .gt. 0  --> zero pivot encountered at step number ierr.
c               ierr  = -1   --> Error. input matrix may be wrong.
c                            (The elimination process has generated a
c                            row in L or U whose length is .gt.  n.)
c               ierr  = -2   --> The matrix L overflows the array al.
c               ierr  = -3   --> The matrix U overflows the array alu.
c               ierr  = -4   --> Illegal value for lfil.
c               ierr  = -5   --> zero row encountered.
c
c     work arrays:
c     jw      = integer work array of length 2*n.
c     w       = real work array of length n
c
c     IMPORTANR NOTE:
c
c     TO AVOID PERMUTING THE SOLUTION VECTORS ARRAYS FOR EACH LU-SOLVE,
C     THE MATRIX A IS PERMUTED ON RETURN. [all column indices are
c     changed]. SIMILARLY FOR THE U MATRIX.
c     To permute the matrix back to its original state use the loop:
c
c     do k=ia(1), ia(n+1)-1
c        ja(k) = iperm(ja(k))
c     enddo
c
c-----------------------------------------------------------------------
c     implicit none
      integer n,ja(*),ia(n+1),lfil,jlu(*),ju(n),jw(2*n),iwk,
     *     iperm(2*n),ierr
      real*8 a(*), alu(*), w(n), droptol
c     local variables
c
      integer k,i,j,jrow,ju0,ii,j1,j2,jpos,len,imax,lenu,lenl,jj,mbloc,
     *     icut
      real*8 s, tmp, tnorm,xmax,xmax0, fact, abs, t, permtol
c
      if (lfil .lt. 0) goto 998
c-----------------------------------------------------------------------
c     initialize ju0 (points to next element to be added to alu,jlu)
c     and pointer array.
c-----------------------------------------------------------------------
      ju0 = n+2
      jlu(1) = ju0
c
c.... integer double pointer array.
c
      do 1 j=1, n
         jw(n+j)  = 0
         iperm(j) = j
         iperm(n+j) = j
 1    continue
c-----------------------------------------------------------------------
c.... beginning of main loop.
c-----------------------------------------------------------------------
      do 500 ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
         tnorm = 0.0d0
         do 501 k=j1,j2
            tnorm = tnorm+abs(a(k))
 501     continue
         if (tnorm .eq. 0.0) goto 999
cww      tnorm = tnorm/(j2-j1+1)
         tnorm = tnorm/real(j2-j1+1)
c
c.... unpack L-part and U-part of row of A in arrays  w  --
c
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = 0.0
         jw(n+ii) = ii
c
         do 170  j = j1, j2
            k = iperm(n+ja(j))
            t = a(j)
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n+k) = lenl
            else if (k .eq. ii) then
               w(ii) = t
            else
               lenu = lenu+1
               jpos = ii+lenu-1
               jw(jpos) = k
               w(jpos) = t
               jw(n+k) = jpos
            endif
 170     continue
         jj = 0
         len = 0
c
c.... eliminate previous rows
c
 150     jj = jj+1
         if (jj .gt. lenl) goto 160
c-----------------------------------------------------------------------
c     in order to do the elimination in the correct order we must select
c     the smallest column index among jw(k), k=jj+1, ..., lenl.
c-----------------------------------------------------------------------
         jrow = jw(jj)
         k = jj
c
c.... determine smallest column index
c
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            endif
 151     continue
c
         if (k .ne. jj) then
c.... exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
c.... exchange in jr
            jw(n+jrow) = jj
            jw(n+j) = k
c.... exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
c
c.... zero out element in row by resetting jw(n+jrow) to zero.
c
         jw(n+jrow) = 0
c
c.... get the multiplier for row to be eliminated: jrow
c
         fact = w(jj)*alu(jrow)
c
c.... drop term if small
c
         if (abs(fact) .le. droptol) goto 150
c
c.... combine current row and row jrow
c
         do 203 k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
c.... new column number
            j = iperm(n+jlu(k))
            jpos = jw(n+j)
            if (j .ge. ii) then
c
c.... dealing with upper part.
c
               if (jpos .eq. 0) then
c
c.... this is a fill-in element
c
                  lenu = lenu+1
                  i = ii+lenu-1
                  if (lenu .gt. n) goto 995
                  jw(i) = j
                  jw(n+j) = i
                  w(i) = - s
               else
c.... no fill-in element --
                  w(jpos) = w(jpos) - s
               endif
            else
c
c.... dealing with lower part.
c
               if (jpos .eq. 0) then
c
c.... this is a fill-in element
c
                 lenl = lenl+1
                 if (lenl .gt. n) goto 995
                 jw(lenl) = j
                 jw(n+j) = lenl
                 w(lenl) = - s
              else
c
c.... this is not a fill-in element
c
                 w(jpos) = w(jpos) - s
              endif
           endif
 203  continue
c
c.... store this pivot element -- (from left to right -- no danger of
c.... overlap with the working elements in L (pivots).
c
        len = len+1
        w(len) = fact
        jw(len)  = jrow
       goto 150
 160    continue
c
c.... reset double-pointer to zero (U-part)
c
        do 308 k=1, lenu
           jw(n+jw(ii+k-1)) = 0
 308  continue
c
c.... update L-matrix
c
        lenl = len
        len = min0(lenl,lfil)
c
c.... sort by quick-split
c
        call qsplit (w,jw,lenl,len)
c
c.... store L-part -- in original coordinates ..
c
        do 204 k=1, len
           if (ju0 .gt. iwk) goto 996
           alu(ju0) =  w(k)
           jlu(ju0) = iperm(jw(k))
           ju0 = ju0+1
 204    continue
c
c.... save pointer to beginning of row ii of U
c
        ju(ii) = ju0
c
c.... update U-matrix -- first apply dropping strategy
c
         len = 0
         do k=1, lenu-1
            if (abs(w(ii+k)) .gt. droptol*tnorm) then
               len = len+1
               w(ii+len) = w(ii+k)
               jw(ii+len) = jw(ii+k)
            endif
         enddo
         lenu = len+1
         len = min0(lenu,lfil)
         call qsplit (w(ii+1), jw(ii+1), lenu-1,len)
c
c.... determine next pivot --
c
        imax = ii
        xmax = abs(w(imax))
        xmax0 = xmax
        icut = ii - 1 + mbloc - mod(ii-1,mbloc)
        do k=ii+1,ii+len-1
           t = abs(w(k))
           if (t .gt. xmax .and. t*permtol .gt. xmax0 .and.
     *          jw(k) .le. icut) then
              imax = k
              xmax = t
           endif
        enddo
c
c.... exchange w's
c
        tmp = w(ii)
        w(ii) = w(imax)
        w(imax) = tmp
c
c.... update iperm and reverse iperm
c
        j = jw(imax)
        i = iperm(ii)
        iperm(ii) = iperm(j)
        iperm(j) = i
c
c.... reverse iperm
c
        iperm(n+iperm(ii)) = ii
        iperm(n+iperm(j)) = j
c-----------------------------------------------------------------------
c
        if (len + ju0 .gt. iwk) goto 997
c
c.... copy U-part in original coordinates
c
        do 302 k=ii+1,ii+len-1
           jlu(ju0) = iperm(jw(k))
           alu(ju0) = w(k)
           ju0 = ju0+1
 302  continue
c
c.... store inverse of diagonal element of u
c
        if (w(ii) .eq. 0.0) w(ii) = (1.0D-4 + droptol)*tnorm
        alu(ii) = 1.0d0/ w(ii)
c
c.... update pointer to beginning of next row of U.
c
	jlu(ii+1) = ju0
c-----------------------------------------------------------------------
c     end main loop
c-----------------------------------------------------------------------
 500  continue
c
c.... permute all column indices of LU ...
c
      do k = jlu(1),jlu(n+1)-1
         jlu(k) = iperm(n+jlu(k))
      enddo
c
c.... ...and of A
c
      do k=ia(1), ia(n+1)-1
         ja(k) = iperm(n+ja(k))
      enddo
c
      ierr = 0
      return
c
c.... incomprehensible error. Matrix must be wrong.
c
 995  ierr = -1
      return
c
c.... insufficient storage in L.
c
 996  ierr = -2
      return
c
c.... insufficient storage in U.
c
 997  ierr = -3
      return
c
c.... illegal lfil entered.
c
 998  ierr = -4
      return
c
c.... zero row encountered
c
 999  ierr = -5
      return
      end


        subroutine qsplit(a,ind,n,ncut)
c-----------------------------------------------------------------------
c     does a quick-sort split of a real array.
c     on input a(1:n). is a real array
c     on output a(1:n) is permuted such that its elements satisfy:
c
c     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
c     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
c
c     ind(1:n) is an integer array which permuted in the same way as a(*).
c-----------------------------------------------------------------------
        real*8 a(n)
        integer ind(n), n, ncut
        real*8 tmp, abskey
        integer itmp, first, last
c-----
        first = 1
        last = n
        if (ncut .lt. first .or. ncut .gt. last) return
c
c     outer loop -- while mid .ne. ncut do
c
 1      mid = first
        abskey = abs(a(mid))
        do 2 j=first+1, last
           if (abs(a(j)) .gt. abskey) then
              mid = mid+1
c     interchange
              tmp = a(mid)
              itmp = ind(mid)
              a(mid) = a(j)
              ind(mid) = ind(j)
              a(j)  = tmp
              ind(j) = itmp
           endif
 2      continue
c
c     interchange
c
        tmp = a(mid)
        a(mid) = a(first)
        a(first)  = tmp
c
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
c
c     test for while loop
c
        if (mid .eq. ncut) return
        if (mid .gt. ncut) then
           last = mid-1
        else
           first = mid+1
        endif
        goto 1

        end


      subroutine iluk(n,a,ja,ia,lfil,alu,jlu,ju,levs,iwk,w,jw,ju0,ierr)
c----------------------------------------------------------------------*
c     Purpose:  ILU WITH LEVEL OF FILL-IN OF K (ILU(k))                *
c----------------------------------------------------------------------*
c
c     on entry:
c
c     n       = integer. The row dimension of the matrix A. The matrix
c
c     a,ja,ia = matrix stored in Compressed Sparse Row format.
c
c     lfil    = integer. The fill-in parameter. Each element whose
c               leve-of-fill exceeds lfil during the ILU process is dropped.
c               lfil must be .ge. 0
c
cnot used??ww     tol     = real*8. Sets the threshold for dropping small terms in the
cnot used??ww               factorization. See below for details on dropping strategy.
c
c     iwk     = integer. The minimum length of arrays alu, jlu, and levs.
c
c     On return:
c
c
c     alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c               the L and U factors together. The diagonal (stored in
c               alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c               contains the i-th row of L (excluding the diagonal entry=1)
c               followed by the i-th row of U.
c
c     ju      = integer array of length n containing the pointers to
c               the beginning of each row of U in the matrix alu,jlu.
c
c     levs    = integer (work) array of size iwk -- which contains the
c               levels of each element in alu, jlu.
c
c     ierr    = integer. Error message with the following meaning.
c               ierr  = 0    --> successful return.
c               ierr .gt. 0  --> zero pivot encountered at step number ierr.
c               ierr  = -1   --> Error. input matrix may be wrong.
c                                (The elimination process has generated a
c                                row in L or U whose length is .gt.  n.)
c               ierr  = -2   --> The matrix L overflows the array al.
c               ierr  = -3   --> The matrix U overflows the array alu.
c               ierr  = -4   --> Illegal value for lfil.
c               ierr  = -5   --> zero row encountered in A or U.
c
c     work arrays:
c
c     jw      = integer work array of length 3*n.    !??? ww 3*n or 2*n
c     w       = real work array of length n
c
c     Notes/known bugs: This is not implemented efficiently storage-wise.
c           For example: Only the part of the array levs(*) associated with
c           the U-matrix is needed in the routine.. So some storage can
c           be saved if needed. The levels of fills in the LU matrix are
c           output for information only -- they are not needed by LU-solve.
c
c     ------------------------------------------------------------------
c     w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]
c     jw(n+1:2n)  stores the nonzero indicator.
c
c     Notes:
c     All the diagonal elements of the input matrix must be  nonzero.
c
c----------------------------------------------------------------------*
      implicit none
      real*8 a(*),alu(*),w(n)
      integer ja(*),ia(n+1),jlu(*),ju(n),levs(*),jw(3*n),n,lfil,iwk,ierr

c     locals
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,n2,jlev, min
      real*8 t, s, fact
      if (lfil .lt. 0) goto 998
c-----------------------------------------------------------------------
c     initialize ju0 (points to next element to be added to alu,jlu)
c     and pointer array.
c-----------------------------------------------------------------------
      n2 = n+n
      ju0 = n+2
      jlu(1) = ju0
c
c     initialize nonzero indicator array + levs array --
c
      do 1 j=1,2*n
         jw(j)  = 0
 1    continue
c-----------------------------------------------------------------------
c     beginning of main loop.
c-----------------------------------------------------------------------
      do 500 ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
c
c     unpack L-part and U-part of row of A in arrays w
c
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = 0.0
         jw(n+ii) = ii
c
         do 170  j = j1, j2
            k = ja(j)
            t = a(j)
            if (t .eq. 0.0) goto 170
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n2+lenl) = 0
               jw(n+k) = lenl
            else if (k .eq. ii) then
               w(ii) = t
               jw(n2+ii) = 0
            else
               lenu = lenu+1
               jpos = ii+lenu-1
               jw(jpos) = k
               w(jpos) = t
               jw(n2+jpos) = 0
               jw(n+k) = jpos
            endif
 170     continue
c
         jj = 0
c
c     eliminate previous rows
c
 150     jj = jj+1
         if (jj .gt. lenl) goto 160
c-----------------------------------------------------------------------
c     in order to do the elimination in the correct order we must select
c     the smallest column index among jw(k), k=jj+1, ..., lenl.
c-----------------------------------------------------------------------
         jrow = jw(jj)
         k = jj
c
c     determine smallest column index
c
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            endif
 151     continue
c
         if (k .ne. jj) then
c     exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
c     exchange in jw(n+  (pointers/ nonzero indicator).
            jw(n+jrow) = jj
            jw(n+j) = k
c     exchange in jw(n2+  (levels)
            j = jw(n2+jj)
            jw(n2+jj)  = jw(n2+k)
            jw(n2+k) = j
c     exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
c
c     zero out element in row by resetting jw(n+jrow) to zero.
c
         jw(n+jrow) = 0
c
c     get the multiplier for row to be eliminated (jrow) + its level
c
         fact = w(jj)*alu(jrow)
         jlev = jw(n2+jj)
         if (jlev .gt. lfil) goto 150
c
c     combine current row and row jrow
c
         do 203 k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
            j = jlu(k)
            jpos = jw(n+j)
            if (j .ge. ii) then
c
c     dealing with upper part.
c
               if (jpos .eq. 0) then
c
c     this is a fill-in element
c
                  lenu = lenu+1
                  if (lenu .gt. n) goto 995
                  i = ii+lenu-1
                  jw(i) = j
                  jw(n+j) = i
                  w(i) = - s
                  jw(n2+i) = jlev+levs(k)+1
               else
c
c     this is not a fill-in element
c
                  w(jpos) = w(jpos) - s
                  jw(n2+jpos) = min(jw(n2+jpos),jlev+levs(k)+1)
               endif
            else
c
c     dealing with lower part.
c
               if (jpos .eq. 0) then
c
c     this is a fill-in element
c
                  lenl = lenl+1
                  if (lenl .gt. n) goto 995
                  jw(lenl) = j
                  jw(n+j) = lenl
                  w(lenl) = - s
                  jw(n2+lenl) = jlev+levs(k)+1
               else
c
c     this is not a fill-in element
c
                  w(jpos) = w(jpos) - s
                  jw(n2+jpos) = min(jw(n2+jpos),jlev+levs(k)+1)
               endif
            endif
 203     continue
         w(jj) = fact
         jw(jj)  = jrow
         goto 150
 160     continue
c
c     reset double-pointer to zero (U-part)
c
         do 308 k=1, lenu
            jw(n+jw(ii+k-1)) = 0
 308     continue
c
c     update l-matrix
c
         do 204 k=1, lenl
            if (ju0 .gt. iwk) goto 996
            if (jw(n2+k) .le. lfil) then
               alu(ju0) =  w(k)
               jlu(ju0) =  jw(k)
               ju0 = ju0+1
            endif
 204     continue
c
c     save pointer to beginning of row ii of U
c
         ju(ii) = ju0
c
c     update u-matrix
c
         do 302 k=ii+1,ii+lenu-1
            if (jw(n2+k) .le. lfil) then
               jlu(ju0) = jw(k)
               alu(ju0) = w(k)
               levs(ju0) = jw(n2+k)
               ju0 = ju0+1
            endif
 302     continue

         if (w(ii) .eq. 0.0) goto 999
c
         alu(ii) = 1.0d0/ w(ii)
c
c     update pointer to beginning of next row of U.
c
         jlu(ii+1) = ju0
c-----------------------------------------------------------------------
c     end main loop
c-----------------------------------------------------------------------
 500  continue
      ierr = 0
      return
c
c     incomprehensible error. Matrix must be wrong.
c
 995  ierr = -1
      return
c
c     insufficient storage in L.
c
 996  ierr = -2
      return
c
c     insufficient storage in U.
c
 997  ierr = -3
      return
c
c     illegal lfil entered.
c
 998  ierr = -4
      return
c
c     zero row encountered in A or U.
c
 999  ierr = -5
      return
      end



      subroutine lusol(n, y, x, alu, jlu, ju)
c-----------------------------------------------------------------------
c
c     This routine solves the system (LU) x = y,
c     given an LU decomposition of a matrix stored in (alu, jlu, ju)
c     modified sparse row format
c
c-----------------------------------------------------------------------
c     on entry:
c     n   = dimension of system
c     y   = the right-hand-side vector
c     alu, jlu, ju
c     = the LU matrix as provided from the ILU routines.
c
c     on return
c     x   = solution of LU x = y.
c-----------------------------------------------------------------------
c
c     Note: routine is in place: call lusol (n, x, x, alu, jlu, ju)
c     will solve the system with rhs x and overwrite the result on x .
c
c-----------------------------------------------------------------------
      real*8 x(n), y(n), alu(*)
      integer n, jlu(*), ju(*)

c.... local variables
c
      integer i,k
c
c.... forward solve
c
      do i = 1, n
        x(i) = y(i)
        do k=jlu(i),ju(i)-1
          x(i) = x(i) - alu(k)* x(jlu(k))
        end do
      end do

c
c.... backward solve.
c
      do i = n, 1, -1
        do k=ju(i),jlu(i+1)-1
          x(i) = x(i) - alu(k)*x(jlu(k))
        end do
        x(i) = alu(i)*x(i)
      end do
c
      return
      end


c-----------------------------------------------------------------------
      subroutine lutsol(n, y, x, alu, jlu, ju)
c-----------------------------------------------------------------------
c
c     This routine solves the system  Transp(LU) x = y,
c     given an LU decomposition of a matrix stored in (alu, jlu, ju)
c     modified sparse row format. Transp(M) is the transpose of M.
c-----------------------------------------------------------------------
c     on entry:
c     n   = dimension of system
c     y   = the right-hand-side vector
c     alu, jlu, ju
c     = the LU matrix as provided from the ILU routines.
c
c     on return
c     x   = solution of transp(LU) x = y.
c-----------------------------------------------------------------------
c
c     Note: routine is in place: call lutsol (n, x, x, alu, jlu, ju)
c     will solve the system with rhs x and overwrite the result on x .
c
c-----------------------------------------------------------------------
      real*8 x(n), y(n), alu(*)
      integer n, jlu(*), ju(*)

c.... local variables
c
      integer i,k
c
      do i = 1, n
         x(i) = y(i)
      end do
c
c.... forward solve (with U^T)
c
      do i = 1, n
        x(i) = x(i) * alu(i)
        do k=ju(i),jlu(i+1)-1
          x(jlu(k)) = x(jlu(k)) - alu(k)* x(i)
        end do
      end do
c
c.... backward solve (with L^T)
c
      do i = n, 1, -1
        do k=jlu(i),ju(i)-1
          x(jlu(k)) = x(jlu(k)) - alu(k)*x(i)
        end do
      end do
c
      return
      end
