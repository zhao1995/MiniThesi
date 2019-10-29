subroutine CalcInverse(n,A,invA)!SW: not tested/checked!!!!!!!!!

!---------------- Purpose: calculate inverse of nxn matrix ----------------!

implicit none

!---------------- Transfer parameters ----------------!

integer, intent(in) :: n

double precision, intent(in) :: A(n,n)

double precision, intent(out) :: invA(n,n)

!---------------- Local variables ----------------!

integer :: ipiv(n)

integer :: info

double precision :: work(n)

double precision :: detA

if ( n .EQ. 1 ) then ! 1x1 Matrix

  detA = A(1,1)

  if ( detA .NE. 0.0D+00 ) then

    invA(1,1) = 1.0D+00 / detA
    
  else
  
    write(*,*) 'ERROR: Problems in computing inverse. STOP!'
    STOP
  
  endif

elseif ( n .EQ. 2 ) then ! 2x2 Matrix

  detA = A(1,1) * A(2,2) - A(2,1) * A(1,2)
  
  if ( detA .NE. 0.0D+00 ) then
  
    invA(1,1) =   A(2,2) / detA
    invA(1,2) = - A(1,2) / detA
  
    invA(2,1) = - A(2,1) / detA
    invA(2,2) =   A(1,1) / detA
    
  else
  
    write(*,*) 'ERROR: Problems in computing inverse. STOP!'
    STOP
  
  endif
  
elseif ( n .EQ. 3 ) then ! 3x3 Matrix

  detA = A(1,1) * A(2,2) * A(3,3) + A(1,2) * A(2,3) * A(3,1) + A(1,3) * A(2,1) * A(3,2) - &
         A(3,1) * A(2,2) * A(1,3) - A(3,2) * A(2,3) * A(1,1) - A(3,3) * A(2,1) * A(1,2)
         
  if ( detA .NE. 0.0D+00 ) then
  
    invA(1,1) = ( A(2,2) * A(3,3) - A(2,3) * A(3,2) ) / detA
    invA(1,2) = (-A(1,2) * A(3,3) + A(1,3) * A(3,2) ) / detA
    invA(1,3) = ( A(1,2) * A(2,3) - A(1,3) * A(2,2) ) / detA
  
    invA(2,1) = (-A(2,1) * A(3,3) + A(2,3) * A(3,1) ) / detA
    invA(2,2) = ( A(1,1) * A(3,3) - A(1,3) * A(3,1) ) / detA
    invA(2,3) = (-A(1,1) * A(2,3) + A(1,3) * A(2,1) ) / detA
  
    invA(3,1) = ( A(2,1) * A(3,2) - A(2,2) * A(3,1) ) / detA
    invA(3,2) = (-A(1,1) * A(3,2) + A(1,2) * A(3,1) ) / detA
    invA(3,3) = ( A(1,1) * A(2,2) - A(1,2) * A(2,1) ) / detA
    
  else
  
    write(*,*) 'ERROR: Problems in computing inverse. STOP!'
    STOP
  
  endif

else ! nxn Matrix

  invA(1:n,1:n) = A(1:n,1:n)
  
  ! Call LAPACK routines
     
  call dgetrf(n,n,invA,n,ipiv,info)

  if (info .NE. 0) then

    write(*,*) 'ERROR: Problems in computing inverse. STOP!'
    STOP

  endif

  call dgetri(n,invA,n,ipiv,work,n,info)

endif

end subroutine CalcInverse