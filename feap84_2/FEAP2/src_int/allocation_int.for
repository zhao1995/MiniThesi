      module doalloc
c-----------------------------------------------------------------------
c.... Purpose: allocation/deallocation of arrays
c....          needs to be inserted everwhere an array is allocated
c....
c....     Store name and size of arrays
c....     Each array get its own number.
c....     The arrays are identified by their memory location.
c
c-----------------------------------------------------------------------

      contains
c
      subroutine ralloc(arr,alen,name,flag)
c-----------------------------------------------------------------------
c.... Purpose: allocation of a real array
c-----------------------------------------------------------------------
      USE allocdata
      implicit none
      real*8, allocatable, dimension(:) :: arr
      real*8  testsize
      character(*)name
      logical flag
      integer alen,num,status

c...  setting start values
      status = 0
      num    = 1

c...  check for length overflow
      if(alen.lt.0) then
        write(*,*) 'Something wrong in subroutine ralloc'
        write(*,*) 'Tried to allocate array: ', name, 'with length: '
     +            , alen
        write(*,*) 'Maybe integer length of: 2.147.483.647 exceeded'
        stop
      end if

c...  Array is allocated
      if(allocated(arr)) then
c...    Get current number of the array
        do while(loc(arr).ne.memloc(num))
          num = num + 1
          if(num.ge.maxnum) then
            write(*,*) 'Something wrong in subroutine ralloc'
            write(*,*) 'Array already allocated but has no number'
            write(*,*) 'Cannot allocate array ', name
            write(*,*) 'Maybe array was not allocated with ralloc'
            stop
          end if  
        end do
c...    Check if array has the correct size
        if(size(arr).ne.alen) then
          deallocate(arr)
          allocate(arr(alen),stat=status)
        end if
c...  Array is not allocated
      else
c...    Calculate new number for the array 
        do while(isalloc(num))
          num = num + 1
          if(num.gt.maxnum) then
            write(*,*) 'Too much arrays allocated array list'
            write(*,*) 'Increase maxnum in module allocdata'
            stop
          end if
        end do
        allocate(arr(alen),stat=status)
      end if
      
c...  Update storage informations for [show,memo]
      if(status.eq.0) then
        memloc(num)    = loc(arr)
        arr = 0.d0
        allocname(num) = name            ! stores the name of the array
        asizes(num)    = real(alen)*8.d0 ! stores the size of the array in bytes
        isalloc(num)   = .true.          ! stores if the array is allocated
        flag = .false.                   ! sets the corresponding flag
      else
c...  In case of an error
        flag     = .true.
        testsize = real(alen)*8.d0/1024.d0/1024.d0/1024.d0
        call pmemo_out(1)
        write(*,1000) name,alen,testsize
        stop
      end if
      
      return
      
1000  format('Error when trying to allocate array ',A,' with size ',i,
     + ' which equals additional ',F10.2,' GB.',/,'If you have enough
     + system memory, keep in mind that a 32-bit version has a limit of
     + 2GB in total.')
      
      end subroutine

      subroutine ialloc(arr,alen,name,flag)
c-----------------------------------------------------------------------
c.... Purpose: allocation of integer array
c-----------------------------------------------------------------------
      USE allocdata
      implicit none
      integer, allocatable, dimension(:) :: arr
      real*8  testsize
      character(*) :: name
      logical flag
      integer alen,num,status
      
c...  setting start values
      status = 0
      num    = 1
      
c...  check for length overflow
      if(alen.lt.0) then
        write(*,*) 'Something wrong in subroutine ialloc'
        write(*,*) 'Tried to allocate array: ', name, 'with length: '
     +            , alen
        write(*,*) 'Maybe integer length of: 2.147.483.647 exceeded'
        stop
      end if

      if(allocated(arr)) then
c...    Get current number of the array
        do while(loc(arr).ne.memloc(num))
          num = num + 1
          if(num.ge.maxnum) then
            write(*,*) 'Something wrong in subroutine ialloc'
            write(*,*) 'Array already allocated but has no number'
            write(*,*) 'Cannot allocate array ', name
            write(*,*) 'Maybe array was not allocated with ialloc'
            stop
          end if  
        end do
c...    Check if array has the correct size
        if(size(arr).ne.alen) then
          deallocate(arr)
          allocate(arr(alen),stat=status)
        end if
c...  Array is not allocated
      else
c...    Calculate new number for the array 
        do while(isalloc(num))
          num = num + 1
          if(num.gt.maxnum) then
            write(*,*) 'Too much arrays allocated array list'
            write(*,*) 'Increase maxnum in module allocdata'
            stop
          end if
        end do
        allocate(arr(alen),stat=status)
      end if
      
c...  Update storage informations for [show,memo]
      if(status.eq.0) then
        memloc(num)    = loc(arr)
        arr = 0
        allocname(num) = name            ! stores the name of the array
        asizes(num)    = real(alen)*4.d0 ! stores the size of the array in bytes
        isalloc(num)   = .true.          ! stores if the array is allocated
        flag = .false.                   ! sets the corresponding flag
      else
c...  In case of an error
        flag     = .true.
        testsize = real(alen)*4.d0/1024.d0/1024.d0/1024.d0
        call pmemo_out(1)
        write(*,1000) name,alen,testsize
        stop
      end if
      
      return
      
1000  format('Error when trying to allocate array ',A,' with size ',i12,
     + ' which equals additional ',F10.2,' GB.',/,'INTEL-Version is a 
     + 64-bit version. The total limit is governed by the core memory.')
      
      end subroutine

      subroutine rdealloc(arr,flag)
c-----------------------------------------------------------------------
c.... Purpose: deallocation of a real array
c-----------------------------------------------------------------------
      Use allocdata
      implicit none
      real*8, allocatable, dimension(:) :: arr
      integer                           :: num
      character*30 name
      logical                           :: flag
      if(allocated(arr)) then 
        num = 1
        do while(loc(arr).ne.memloc(num))
          num = num + 1
          if(num.ge.maxnum) then
            write(*,*) 'Something wrong in subroutine rdealloc'
            write(*,*) 'Array is allocated but has no number'
            write(*,*) 'Maybe array was not allocated with ralloc'
            stop
          end if  
        end do
        deallocate(arr)
        memloc(num)    = 0
        flag           = .true.         ! sets the corresponding flag
        isalloc(num)   = .false.        ! stores if the array ist allocated
        asizes(num)    = 0.d0           ! reset the size
        name           = allocname(num) ! see name of the array
      end if
      return
      end subroutine

      subroutine idealloc(arr,flag)
c-----------------------------------------------------------------------
c.... Purpose: deallocation of a integer array
c-----------------------------------------------------------------------
      Use allocdata
      implicit none
      character*30 name
      integer, allocatable, dimension(:)  :: arr
      integer                             :: num
      logical                             :: flag
      if(allocated(arr)) then 
        num = 1
        do while(loc(arr).ne.memloc(num))
          num = num + 1
          if(num.ge.maxnum) then
            write(*,*) 'Something wrong in subroutine idealloc'
            write(*,*) 'Array is allocated but has no number'
            write(*,*) 'Maybe array was not allocated with ialloc'
            stop
          end if  
        end do
        deallocate(arr)
        memloc(num)    = 0
        flag           = .true.         ! sets the corresponding flag
        isalloc(num)   = .false .       ! stores if the array ist allocated
        asizes(num)    = 0.d0           ! reset the size
        name           = allocname(num) ! see name of the array
      end if
      return      
      end subroutine

      end module
