module itertools
  implicit none
  private
  public vararray, cproduct

  !!<summary>Provides a structure for variable length arrays that need to have their
  !!cartesian product taken.</summary>
  !!<usage>
  !!type(vararray) single
  !!allocate(single)
  !!single%init(array)
  !!</usage>
  type vararray
     !!<member name="items">The array of items to take cartesian product over.</member>
     !!<member name="length">The number of items in @CREF[this.items].</member>
     integer, pointer :: items(:)
     integer :: length
  contains
    procedure, public :: init => vararray_init
  end type vararray
contains
  !!<summary>Initializes the array items and length property.</summary>
  subroutine vararray_init(self, array)
    class(vararray) :: self
    integer, target, intent(in) :: array(:)
    self%items => array
    self%length = size(self%items, 1)
  end subroutine vararray_init

  !!<summary>Builds the cartesian product of the specified arrays.</summary>
  !!<parameter name="elements" regular="true">A 1-D ragged-array of vararray instances to take the 
  !!cartesian product over.</parameter>
  !!<parameter name="presult" regular="true">A 2-D array with each row being a cartesian product entry
  !!with one item contributed from each of the vararrays in elements.
  !!</parameter>
  subroutine cproduct(elements, presult)
    class(vararray), allocatable, intent(in) :: elements(:)
    integer, allocatable, intent(inout) :: presult(:,:)

    integer :: i, width = 1
    do i=1, size(elements,1)
       width = width*elements(i)%length 
    end do
    call rproduct(elements, presult, 1, 1, width)
  end subroutine cproduct

  !!<summary>Recursively takes the product of the existing subsets in presult with the
  !!new subset specified by depth.</summary>
  !!<parameter name="subsets">A subset whose value need to be cartesian-multiplied with the existing
  !!subsets in presult.</parameter>
  !!<parameter name="presult">The total cproduct result to alter.</parameter>
  !!<parameter name="depth">Recursive depth on the function stack.</parameter>
  !!<parameter name="start">The row number in presult to start altering values for.</parameter>
  !!<parameter name="right">The product of the array sizes from this array to the right.</parameter>
  recursive subroutine rproduct(subsets, presult, depth, start, right)
    class(vararray), intent(in) :: subsets(:)
    integer, allocatable, intent(inout) :: presult(:,:)
    integer, intent(in) :: depth, start, right
    
    !!<local name="cursor">Points to the row in the presult that is currently being updated
    !!with the values of the subset being multiplied.</local>
    !!<local name="kids">The product of array lengths for subsets to the right of this one
    !!*not* including this one.</local>
    integer :: i, j, cursor
    integer :: kids
    integer :: items(subsets(depth)%length)
    items = subsets(depth)%items

    !First we need to get the product of the list sizes to the right of this one.
    cursor = start
    kids = right/subsets(depth)%length
    do i=1, subsets(depth)%length
       do j=1, kids
          presult(cursor+j-1,depth) = items(i)
       end do
       !Now we just fill in the spots to the left of this column in the presult table.
       if (depth < size(subsets, 1)) then
          call rproduct(subsets, presult, depth+1, cursor, kids)
       end if
       cursor = cursor + kids
    end do
  end subroutine rproduct
end module itertools
