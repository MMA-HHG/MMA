! This module is used for buffering data during the computation
!
! This module was taken directly from http://fortranwiki.org/fortran/show/gen_list

MODULE linked_list
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: list_t
  PUBLIC :: list_data
  PUBLIC :: list_init
  PUBLIC :: list_free
  PUBLIC :: list_insert
  PUBLIC :: list_append
  PUBLIC :: list_put
  PUBLIC :: list_get
  PUBLIC :: list_next
  PUBLIC :: list_end

  ! A PUBLIC variable to use as a MOLD for transfer()
  INTEGER, DIMENSION(:), ALLOCATABLE :: list_data

  ! Linked list node data TYPE
  TYPE :: list_t
     !private
     INTEGER, DIMENSION(:), POINTER :: data => NULL()
     TYPE(list_t), POINTER :: next => NULL()
  END TYPE list_t
  
  
CONTAINS

  ! Initialize a head node SELF and OPTIONALly store the provided DATA.
  SUBROUTINE list_init(self, data)
    TYPE(list_t), POINTER :: self
    INTEGER, DIMENSION(:), INTENT(in), OPTIONAL :: data
    
    ALLOCATE(self)
    NULLIFY(self%next)

    IF (present(data)) THEN
       ALLOCATE(self%data(size(data)))
       self%data = data
    ELSE
       NULLIFY(self%data)
    END IF
  END SUBROUTINE list_init

  ! Free the entire list and all data, beginning at SELF
  SUBROUTINE list_free(self)
    TYPE(list_t), POINTER :: self
    TYPE(list_t), POINTER :: current
    TYPE(list_t), POINTER :: next

    current => self
    DO WHILE (ASSOCIATED(current))
       next => current%next
       IF (ASSOCIATED(current%data)) THEN
          DEALLOCATE(current%data)
          NULLIFY(current%data)
       END IF
       DEALLOCATE(current)
       NULLIFY(current)
       current => next
    END DO
  END SUBROUTINE list_free

  ! Return the next node after SELF
  FUNCTION list_next(self) result(next)
    TYPE(list_t), POINTER :: self
    TYPE(list_t), POINTER :: next
    next => self%next
  END FUNCTION list_next

  ! Loop through the linked list and print the data of the last link
  SUBROUTINE list_end(self)
    TYPE(list_t), POINTER :: self
    TYPE(list_t), POINTER :: current
    TYPE(list_t), POINTER :: next

    current => self
    DO WHILE (ASSOCIATED(current))
       next => current%next
       IF (ASSOCIATED(next)) THEN
         current => next
       ELSE
         print *,current%data
         EXIT
       ENDIF
    END DO
  END SUBROUTINE list_end

  ! Insert a list node after SELF containing DATA (optional)
  SUBROUTINE list_insert(self, data)
    TYPE(list_t), POINTER :: self
    INTEGER, DIMENSION(:), INTENT(in), OPTIONAL :: data
    TYPE(list_t), POINTER :: next

    ALLOCATE(next)
    
    IF (present(data)) THEN
       ALLOCATE(next%data(size(data)))
       next%data = data
    ELSE
       NULLIFY(next%data)
    END IF

    next%next => self%next
    self%next => next
  END SUBROUTINE list_insert

  ! Loop through the list and add a new node at the end of the list
  SUBROUTINE list_append(self, data)
    TYPE(list_t), POINTER :: self
    INTEGER, DIMENSION(:), INTENT(in), OPTIONAL :: data
    TYPE(list_t), POINTER :: current
    TYPE(list_t), POINTER :: next

    current => self
    DO WHILE (ASSOCIATED(current))
       next => current%next
       IF (ASSOCIATED(next)) THEN
         current => next
       ELSE
         ALLOCATE(next)
         IF (present(data)) THEN
           ALLOCATE(next%data(size(data)))
           next%data = data
         ELSE
           NULLIFY(next%data)
         END IF
         next%next => current%next
         current%next => next
         EXIT
       ENDIF
    END DO
  END SUBROUTINE list_append

  ! Store the encoded DATA in list node self
  SUBROUTINE list_put(self, data)
    TYPE(list_t), POINTER :: self
    INTEGER, DIMENSION(:), INTENT(in) :: data

    IF (ASSOCIATED(self%data)) THEN
       DEALLOCATE(self%data)
       NULLIFY(self%data)
    END IF
    self%data = data
  END SUBROUTINE list_put

  ! Return the DATA stored in the node SELF
  FUNCTION list_get(self) result(data)
    TYPE(list_t), POINTER :: self
    INTEGER, DIMENSION(:), POINTER :: data
    data => self%data
  END FUNCTION list_get

END MODULE linked_list

! A derived type for storing data.
module ll_data
  implicit none

  private
  public :: data_t
  public :: data_ptr

  ! Data is stored in data_t
  type :: data_t
     real :: x
  end type data_t

  ! A trick to allow us to store pointers in the list
  type :: data_ptr
     type(data_t), pointer :: p
  end type data_ptr
end module ll_data