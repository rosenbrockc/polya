MODULE group_theory
implicit none
private
public grouper

CONTAINS

  !!<summary>Reallocates a 2D integer array to have the new size 'n'x'm'.</summary>
  !!<parameter name="p">The existing 2D array to be reallocated.</parameter>
  !!<parameter name="n, m">The number of rows and columns in the new array.</parameter>
  function ralloc(p, n, m)
    integer, pointer, dimension(:,:) :: p, ralloc
    integer, intent(in) :: n, m
    
    !!<local name="nold, mold">The number of rows and columns in the existing array.</local>
    !!<local name="ierr">Used for checking memory allocation error.</local>
    integer :: nold, mold, ierr

    allocate(ralloc(1:n,1:m), STAT=ierr)
    if (ierr/=0) stop "Allocate error in ralloc_integer_table in utilities module"

    if (.not.associated(p)) return
    nold = min(size(p,1),n); mold = min(size(p,2),m)
    ralloc = 0
    ralloc(1:nold,1:mold) = p(1:nold,1:mold)
    deallocate(p)
  end function ralloc

!!<summary>Generates a group from a list of generators</summary>
!!<usage>Pass in 2D array where the rows are the permutations. The permutations will be taken in all 
!!possible pairs to generate new group elements. The array of elements will be added to until the
!!list of permutations has closure and constitutes a group.</usage>
!!<comments>Permutations are indexed 1..n, not 0..n-1</comments>
!!<parameter name="g" regular="true">On entry contains the generators, on exit contains the whole group </parameter>
subroutine grouper(g)
integer, pointer:: g(:,:) 

!!<local name="nG">Number of elements in the group (on exit).</local>
!!<local name="growing, new">Flags: group is still being added to, new group element was found</local>
!!<local name="mG">Number of elements already checked (don't need to be re-checked)</local>
!!<local name="ig, jg, kg, cg">Generic loop variables</local>
!!<local name="n">Number of permutation elements.</local>
!!<local name="newg, vs">Permuted group element, vector subscript.</local>
integer :: nG
logical :: growing, new
integer :: mG
integer :: ig, jg, kg, cg
integer :: n
integer, allocatable :: newg(:), vs(:)

n = size(g,2) ! Number of permutation elements
if (n<1) stop "ERROR: Empty list passed to 'grouper' subroutine in 'group_theory' module"
allocate(newg(n),vs(n))
if(any(g>n) .or. any(g<1)) then
   write(*,'("n: ",i4)') n
   write(*,'("g: ",5000(i3,1x))') 
   stop "ERROR: generators passed to 'grouper' (group_theory in celib) contain indices outside of 1..n"
endif

growing = .true.
mG = 0
nG = size(g,1)
do while (growing)
   growing = .false.
   g => ralloc(g,nG+nG**2,n) ! Make room for possible new elements when the existing elements are combined
   ! Loop over each pair of group elements...
   cg = 0 ! Keep track of the number of new elements of the group that have been found in this iteration
   do ig = 1, nG
      do jg = 1, nG
         if (ig <= mG .and. jg <= mG) cycle !...but skip those already tested in a previous iteration
         ! Construct a permuted element, that is, multiply two group elements
         vs = g(jg,:) ! Make a vector subscript for the permuted element
         newg = g(ig,vs) ! Permute the element
         ! Is this new g unique in the list of existing group elements?
         new = .true.
         newcheck: do kg = 1, nG+cg
            if (all(newg == g(kg,:))) then
               new = .false.
               exit newcheck
            endif
         enddo newcheck
         if (new) then
            growing = .true.
            cg = cg + 1
            g(nG+cg,:) = newg
         end if
      end do
   end do ! Loops over pairs
   mG = nG
   nG = nG + cg
end do ! Loop for group still growing
g => ralloc(g,nG,n)
end subroutine grouper
END MODULE group_theory
