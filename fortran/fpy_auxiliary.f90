!!<fortpy version="1.6.2" />
!!<summary>Auto-generated auxiliary module exposing interfaces to save
!!instances of user-derived types.</summary>
module fpy_auxiliary
  use fortpy
  use iso_c_binding, only: c_loc, c_ptr
  use classes
  use itertools
  implicit none
  private
  
  public auxsave

  !!<summary>Provides a single call interface for user-derived data types
  !!for single values and arrays.</summary>
  interface auxsave
     module procedure auxsave_multinomial1d, auxsave_product1d, auxsave_sequence0d, &
                       auxsave_sequence1d, auxsave_vararray1d, auxsave_vararray2d1d, &
                       auxsave_multinomial0d, auxsave_product0d, auxsave_vararray0d, &
                       auxsave_vararray2d0d
  end interface auxsave

contains
  subroutine auxsave_Multinomial1d(variable, filename)
    character(len=*), intent(in) :: filename
    type(Multinomial), intent(in) :: variable(:)
    integer :: fpy0
    character(100) :: pslist

    do fpy0=1, size(variable, 1)
      call fpy_period_join_indices(pslist, (/ fpy0 /), 1)
      call auxsave_Multinomial0d(variable(fpy0), filename//'-'//pslist)
    end do
  end subroutine auxsave_Multinomial1d

  subroutine auxsave_Product1d(variable, filename)
    character(len=*), intent(in) :: filename
    type(Product), intent(in) :: variable(:)
    integer :: fpy0
    character(100) :: pslist

    do fpy0=1, size(variable, 1)
      call fpy_period_join_indices(pslist, (/ fpy0 /), 1)
      call auxsave_Product0d(variable(fpy0), filename//'-'//pslist)
    end do
  end subroutine auxsave_Product1d

  recursive subroutine auxsave_Sequence0d_(variable, prefix, stack)
    type(Sequence), intent(in) :: variable
    integer(fli), allocatable, optional, intent(inout) :: stack(:)
    character(len=*), intent(in) :: prefix

    integer :: i, loc
    integer(fli) :: tadd
    character(len=:), allocatable :: nprefix
    character(100) :: pslist
    integer(fli), allocatable :: stack_(:)

    !V: varseq auto-class support variables
    integer :: variable_kids_ac1
    character(100) :: variable_kids_acps2

    if (present(stack)) then
       stack_ = stack
    else
       allocate(stack_(0))
    end if

    call pysave(variable%varcount, 'filename/varseq/_-varcount')
    call pysave(variable%used, 'filename/varseq/_-used')
    tadd = fpy_address(c_loc(variable%parent))
    loc = address_loc(stack_, tadd)
    nprefix = 'filename/varseq/_-parent'
    if (loc .eq. -1) then
      call append_address(stack_, tadd)
      call auxsave_Sequence0d_(variable%parent, nprefix, stack_)
    else
      call pysave(loc, nprefix)
    end if
    call pysave(variable%kidcount, 'filename/varseq/_-kidcount')
    do variable_kids_ac1=1, size(variable%kids, 1)
      call fpy_period_join_indices(variable_kids_acps2, (/ variable_kids_ac1 /), 1)
      tadd = fpy_address(c_loc(variable%kids(variable_kids_ac1)))
      loc = address_loc(stack_, tadd)
      nprefix = 'filename/varseq/_-kids-'//trim(adjustl(variable_kids_acps2))
      if (loc .eq. -1) then
        call append_address(stack_, tadd)
        call auxsave_Sequence0d_(variable%kids(variable_kids_ac1), nprefix, stack_)
      else
        call pysave(loc, nprefix)
      end if
    end do
    call pysave(variable%root, 'filename/varseq/_-root')
  end subroutine auxsave_Sequence0d_

  subroutine auxsave_Sequence0d(variable, filename)
    character(len=*), intent(in) :: filename
    type(Sequence), intent(in) :: variable
    call auxsave_Sequence0d_(variable, filename//'-')
  end subroutine auxsave_Sequence0d

  subroutine auxsave_Sequence1d(variable, filename)
    character(len=*), intent(in) :: filename
    type(Sequence), intent(in) :: variable(:)
    integer :: fpy0
    character(100) :: pslist

    do fpy0=1, size(variable, 1)
      call fpy_period_join_indices(pslist, (/ fpy0 /), 1)
      call auxsave_Sequence0d(variable(fpy0), filename//'-'//pslist)
    end do
  end subroutine auxsave_Sequence1d

  subroutine auxsave_vararray1d(variable, filename)
    character(len=*), intent(in) :: filename
    type(vararray), intent(in) :: variable(:)
    integer :: fpy0
    character(100) :: pslist

    do fpy0=1, size(variable, 1)
      call fpy_period_join_indices(pslist, (/ fpy0 /), 1)
      call auxsave_vararray0d(variable(fpy0), filename//'-'//pslist)
    end do
  end subroutine auxsave_vararray1d

  subroutine auxsave_vararray2d1d(variable, filename)
    character(len=*), intent(in) :: filename
    type(vararray2d), intent(in) :: variable(:)
    integer :: fpy0
    character(100) :: pslist

    do fpy0=1, size(variable, 1)
      call fpy_period_join_indices(pslist, (/ fpy0 /), 1)
      call auxsave_vararray2d0d(variable(fpy0), filename//'-'//pslist)
    end do
  end subroutine auxsave_vararray2d1d

  subroutine auxsave_Multinomial0d(variable, filename)
    character(len=*), intent(in) :: filename
    type(Multinomial), intent(in) :: variable
    !V: multinoms auto-class support variables

  call pysave(variable%powersum, 'filename/multinoms/_-powersum')
    call pysave(variable%exponent, 'filename/multinoms/_-exponent')
    call pysave(variable%power, 'filename/multinoms/_-power')
    call pysave(variable%possible_powers, 'filename/multinoms/_-possible_powers')
  end subroutine auxsave_Multinomial0d

  subroutine auxsave_Product0d(variable, filename)
    character(len=*), intent(in) :: filename
    type(Product), intent(in) :: variable
    !V: polyndict auto-class support variables
    integer :: variable_multinoms_ac1
    character(100) :: variable_multinoms_acps2

  call pysave(variable%coefficient, 'filename/polyndict/_-coefficient')
    call pysave(variable%multinom_rd, 'filename/polyndict/_-multinom_rd')
    call pysave(variable%targets, 'filename/polyndict/_-targets')
    do variable_multinoms_ac1=1, size(variable%multinoms, 1)
      call fpy_period_join_indices(variable_multinoms_acps2, (/ variable_multinoms_ac1 /), 1)
      call auxsave(variable%multinoms(variable_multinoms_ac1), 'filename/polyndict/_-multinoms-'//trim(adjustl(variable_multinoms_acps2)))
    end do
  end subroutine auxsave_Product0d

  subroutine auxsave_vararray0d(variable, filename)
    character(len=*), intent(in) :: filename
    type(vararray), intent(in) :: variable
    !V: passin auto-class support variables

  call pysave(variable%items, 'filename/passin/_-items')
    call pysave(variable%length, 'filename/passin/_-length')
  end subroutine auxsave_vararray0d

  subroutine auxsave_vararray2d0d(variable, filename)
    character(len=*), intent(in) :: filename
    type(vararray2d), intent(in) :: variable
    !V: mnseq auto-class support variables

  call pysave(variable%items, 'filename/mnseq/_-items')
    call pysave(variable%length, 'filename/mnseq/_-length')
  end subroutine auxsave_vararray2d0d


  !!<summary>Adds the specified value to the end of the integer-valued list.</summary>
  !!<parameter name="list">The integer-valued list to append to.</parameter>
  !!<parameter name="appendage">The extra value to append to the list.</parameter>
  subroutine append_address(list, appendage)
    integer(fli), allocatable, intent(inout) :: list(:)
    integer(fli), intent(in) :: appendage
    integer(fli), allocatable :: templist(:)

    allocate(templist(size(list,1)+1))
    templist(1:size(list,1)) = list
    templist(size(list,1)+1) = appendage

    call move_alloc(templist, list)
  end subroutine append_address
  
  !!<summary>Returns the integer index of the specified address in the stack.</summary>
  !!<parameter name="stack">Array of 64-bit addresses to pointers.</parameter>
  !!<parameter name="value">Address to search for in the stack.</parameter>
  integer function address_loc(stack, avalue)
    integer(fli), intent(in) :: stack(:), avalue
    integer :: i
    address_loc = -1
    do i=1, size(stack, 1)
       if (stack(i) .eq. avalue) then
          address_loc = i
          exit
       end if       
    end do
  end function address_loc
  
  !!<summary>Returns the 64-bit integer address of the pointer at @CREF[param.cloc].</summary>
  !!<parameter name="cloc">The object returned by calling the @CREF[c_loc] interface on
  !!the variable.</parameter>
  integer(fli) function fpy_address(cloc)
    type(c_ptr), intent(in) :: cloc
    character(50) :: saddress
    write(saddress, *) cloc
    read(saddress, *) fpy_address    
  end function fpy_address    
end module fpy_auxiliary
