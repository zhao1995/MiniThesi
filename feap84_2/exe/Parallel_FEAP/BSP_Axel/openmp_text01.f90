program bsp
  use omp_lib
  implicit none
  
  ! fork
  !$omp parallel num_threads(4)

    ! Das nur 1x ausgeben (beim Master Thread)
    if( omp_get_thread_num() == 0 ) then
      write( *, * ) 'Insgesamt gibt es ', omp_get_num_threads(), 'Thread(s)'
    end if      

    ! Das bei jedem Thread ausgeben
    write( *, * ) 'Thread ', omp_get_thread_num(), 'ist aktiv'     
  ! join  
  !$omp end parallel

! Ausgabe:
!    Insgesamt gibt es            3 Thread(s)
!    Thread            0 ist aktiv
!    Thread            1 ist aktiv
!    Thread            2 ist aktiv
  end program bsp
  
