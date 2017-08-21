
MODULE ROUTINES


  INTERFACE Swap
    MODULE PROCEDURE Swap_int
    MODULE PROCEDURE Swap_real
  END INTERFACE Swap

  INTERFACE Sort
    MODULE PROCEDURE Sort_int
    MODULE PROCEDURE Sort_real
  END INTERFACE Sort


  INTERFACE FindMinimum
    MODULE PROCEDURE FindMinimum_int
    MODULE PROCEDURE FindMinimum_real
  END INTERFACE FindMinimum

  PUBLIC Swap, Sort, FindMinimum

CONTAINS



!!!!! the following subroutines are adapted from !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! university lectures available at http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90 !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   INTEGER FUNCTION  FindMinimum_int(x, Start, end)
      IMPLICIT  NONE
      INTEGER, DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, end
      INTEGER                            :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)		! assume the first is the min
      Location = Start			! record its position
      DO i = Start+1, end		! start with next elements
         IF (x(i) < Minimum) THEN	!   if x(i) less than the min?
            Minimum  = x(i)		!      Yes, a new minimum found
            Location = i                !      record its position
         end IF
      end DO
      FindMinimum_int = Location        	! return the position
   end FUNCTION  FindMinimum_int

! --------------------------------------------------------------------
! subroutine  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   subroutine  Swap_int(a, b)
      IMPLICIT  NONE
      INTEGER, INTENT(INOUT) :: a, b
      INTEGER                :: Temp

      Temp = a
      a    = b
      b    = Temp
   end subroutine  Swap_int


! --------------------------------------------------------------------
! subroutine  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   subroutine  Sort_int(x,x_sorted, Size)
      IMPLICIT  NONE
      INTEGER, DIMENSION(1:), intent(in)    :: x
      INTEGER, DIMENSION(1:), intent(out)   :: x_sorted
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location
      x_sorted(:) = x(:) 
      DO i = 1, Size-1			! except for the last
         Location = FindMinimum(x_sorted, i, Size)	! find min from this to last
         CALL  Swap(x_sorted(i), x_sorted(Location))	! swap this and the minimum
      end DO
   end subroutine  Sort_int

! --------------------------------------------------------------------
! the following sorts for real numbers
! -------------------------------------------------------------------
  
   real(kind=8) FUNCTION  FindMinimum_real(x, Start, end)
      IMPLICIT  NONE
      real(kind=8), DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, end
      real(kind=8)                            :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)		! assume the first is the min
      Location = Start			! record its position
      DO i = Start+1, end		! start with next elements
         IF (x(i) < Minimum) THEN	!   if x(i) less than the min?
            Minimum  = x(i)		!      Yes, a new minimum found
            Location = i                !      record its position
         end IF
      end DO
      FindMinimum_real = Location       ! return the position
   end FUNCTION  FindMinimum_real

! --------------------------------------------------------------------
! subroutine  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   subroutine  Swap_real(a, b)
      IMPLICIT  NONE
      real(kind=8), INTENT(INOUT) :: a, b
      real(kind=8)                :: Temp

      Temp = a
      a    = b
      b    = Temp
   end subroutine  Swap_real


! --------------------------------------------------------------------
! subroutine  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   subroutine  Sort_real(x,x_sorted, Size)
      IMPLICIT  NONE
      real(kind=8), DIMENSION(1:), intent(in)    :: x
      real(kind=8), DIMENSION(1:), intent(out)   :: x_sorted
      INTEGER, INTENT(IN)                        :: Size
      INTEGER                                    :: i
      INTEGER                                    :: Location
      x_sorted(:) = x(:) 
      DO i = 1, Size-1			! except for the last
         Location = FindMinimum(x_sorted, i, Size)	! find min from this to last
         CALL  Swap(x_sorted(i), x_sorted(Location))	! swap this and the minimum
      end DO
   end subroutine  Sort_real


END MODULE ROUTINES
