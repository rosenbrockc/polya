!!<summary>AUTHORS: Conrad W. Rosenbrock, Wiley S. Morgan (October 2014)
!! Classes to support the calculation of coefficients for specific terms in a product
!! of multinomials. Construct a product class by specifying the exponent and target term
!! and then add multinomials using the Product instance's append(). The coefficient is
!! then available from the coeff().

!! EXAMPLE: find the coefficient of x^4.y^4.z^4 in 3*(x+y+z)^3*(x^2+y^2+z^2)^3*(x^3+y^3+z^3).
!! ANSWER: 162 from Mathematica

!! >> p = Product(3, [4,4,4])
!! >> p.append(Multinomial([1,1,1],3))
!! >> p.append(Multinomial([2,2,2],3))
!! >> p.append(Multinomial([3,3,3],1))
!! >> print(p.coeff())
!! 162</summary>

MODULE classes

use itertools

IMPLICIT NONE

  !!<summary>Represents an exponent-limited sequence with a single root. Here, sequence represents a
  !!sequence of integer values k_1, k_2...k_j that are the exponents of a single term in a multinomial.
  !!The root of the sequence is one of the k_i; its children become sets of sequences including
  !!variables to the right of i</summary>
  type, public :: Sequence
     private
     !!<member name="root">The exponent of the variable to the left in the multinomial term.</member>
     !!<member name="parent">A Sequence instance for the variable to the *left* of this one (i.e. has index i-1).</member>
     !!<member name="varcount">The number of variables in the multinomial (also length of targets).</member>
     !!<member name="kids"></member>
     integer :: root, used, varcount
     class(Sequence), pointer :: parent => null()
     class(Sequence), pointer :: kids(:)
     contains
       procedure, public :: initialize => sequence_initialize
       procedure, public :: expand
       procedure, public :: kidcount
       procedure, public :: finalize => sequence_finalize
  end type Sequence

  !!<summary>Represents a multinomial expansion.</summary>
  type Multinomial 
     !!<member name="pawer">the powers on each of the *unexpanded* variables in the multinomial;
     !!of the form (x^2+y^2+z^2) => 2. No support is provided for coefficients on the
     !!variables inside the brackets (only Polya implementation).</member>
     !!<member name="exponent">the exponent of the entire multinomial.</member>
     !!<member name="powersum">The maximum allowed exponent on any term in the expansion.</member>
     !!<member name="possible_powers">The list of possible exponents a single variable can have in
     !!any of the terms in the expansion.</member>
     integer, allocatable :: power
     integer :: exponent, powersum
     integer, allocatable :: possible_powers(:)
     contains
       procedure, public :: initialize => multinomial_initialize
       procedure, public :: nchoosekm
  end type Multinomial

  !!<summary>Represents a product of multinomials for which only a single term is interesting.</summary>
  type Product
     !!<member name="coefficient">The scalar integer multiplying this product of multinomials.</member>
     !!<member name="targets">a list of exponents for the only interesting term in the product. The
     !!list is in the order that the variables appear in each multinomial.</member>
     !!<member name="multinoms">An array containing the multinomials.</member>   
     integer :: coefficient
     integer, allocatable :: targets(:)
     class(Multinomial), allocatable :: multinoms(:)
     integer, allocatable :: multinom_rd(:,:)
     contains
       procedure, public :: initialize => product_initialize
       procedure, public :: coeff 
       procedure, public :: sum_sequence
  end type Product

contains

  !!<summary>Initializes a sequence collector for a variable. 'Term' refers to a product
  !!of variables like x^i.y^j.z^k, then x has index 0, y has 1, etc.</summary>
  !!<parameter name="root">The exponent of the variable to the left in the multinomial term.</parameter>
  !!<parameter name="possibles">A list of possible values for each variable in the multinomial.</parameter>
  !!<parameter name="depth">The index of the variable being sequenced in the term.</parameter>
  !!<parameter name="powersum">The maximum value that the sum of exponents in the sequence 
  !!is allowed to have.</parameter>
  !!<parameter name="targets">An array of the target exponents on each variable.</parameter>
  !!<parameter name="parent">A Sequence instance for the variable to the *left* of this one 
  !!(i.e. has index i-1).</parameter>
  recursive subroutine sequence_initialize(this, root, possibles, depth, powersum, targets, parent)
    class(Sequence) :: this
    integer, intent(in) :: powersum, root, depth
    integer, allocatable, intent(in) :: targets(:)
    integer, intent(in) :: possibles(:)
    class(Sequence), target, optional :: parent
    
    !!<local name="kidroots">An array of the exponents af the children, ie the possibles-root 
    !!that still meet the conditions</local>
    !!<local name="kidcount">The number of kids that each depth has</local>
    integer, allocatable :: kidroots(:)
    integer :: j, p, kidcount

    this%root = root
    if (present(parent)) then
       this%parent => parent
       this%used = this%root + this%parent%used
    else
       this%used = this%root
    end if

    this%varcount = size(possibles, 1)

    allocate(kidroots(this%varcount))
    !We only keep recursively defining sequences until we run out of variables in
    !the term. Possibles is a list of possible exponents for each variable in the
    !term and has the same number of items as variables in the term.
    if (depth .le. this%varcount) then
       !Filter the possible values for the variable being considered based on the
       !exponent of the multinomial. When multinomials are expanded, the sum of
       !the exponents in any term must be less than the exponent on the multinomial
       !times the maximum power of any of its (unexpanded) terms.
       !We find all the possible values for this variable by ensuring that:
       !1) it's exponent is compatible with the exponents of all variables to the left of it.
       !2) the exponent we are suggesting is in the list of possible values for the variable.
       !3) the exponent remains positive.
       kidcount = 0
       do j = 1, this%varcount
          p = possibles(j)
          if ((p-root >= 0) .and. (p-root <= targets(j)) .and. &
              (abs(p-root) <= powersum-this%used) .and. (mod(abs(p-this%used), possibles(1)) == 0)) then
             kidcount = kidcount+1
             kidroots(kidcount)= p-root
          end if
       end do

       allocate(this%kids(kidcount))
       !Construct possible exponents for the variable to the right of this one based on the
       !possible choices we just calculated. This gets done recursively until we have a tree
       !of possible exponent values for the variables from left to right.
       do j = 1, kidcount
          call this%kids(j)%initialize(kidroots(j), possibles, depth+1, powersum, targets, this)
       end do
    else 
       this%kids => null()
    end if
  end subroutine sequence_initialize
  
  !!<summary>Cleans up the local variables so that @CREF[this.initialize()] can be called on the
  !!*same* instance that is already allocated.</summary>
  subroutine sequence_finalize(this)
    class(Sequence) :: this
    deallocate(this%kids)
  end subroutine sequence_finalize

  !!<summary>Counts the total number of kids and grandkids etc. to the latest generation.</summary>
  recursive integer function kidcount(this)
    class(Sequence) :: this
    integer :: i

    kidcount = 0
    do i=1, size(this%kids,1)
       kidcount = kidcount + this%kids(i)%kidcount()
    end do
  end function kidcount

  !!<summary>Recursively generates a list of all relevant sequences for this multinomial term.</summary>
  recursive subroutine expand(this, sequences, start, varindex)
    class(Sequence) :: this 
    integer, intent(inout), allocatable :: sequences(:,:)
    integer, intent(in) :: start, varindex

    !!<local name=""></local>
    !!<local name=""></local>
    integer :: j, k, cursor

    !!<comments>Iterate through the child sequences and add their variable root values if
    !!the total sequence sums to the target.</comments>
    cursor = start
    if (.not. allocated(sequences)) then
       allocate(sequences(this%kidcount(), this%varcount))
    end if
    
    do j = 1, size(this%kids,1)
       !!<comments>Here is where the recursion happens; we add the sequence of this variable's
       !!children to the right of this root.</comments>
       call this%kids(j)%expand(sequences, cursor, varindex+1)
       cursor = cursor + this%kids(j)%kidcount()
       if (varindex .eq. this%varcount) cursor = cursor+1
       do k = start, cursor-1
          sequences(k,varindex) = this%root
       end do
    end do
  end subroutine expand

  !!<summary>Initializes the empty product of multinomials.</summary>
  !!<parameter name="multinomials">A mx2 array defining the power and exponent of each multinomial
  !!in the product. Each row is a multinomial definition. DEALLOCATES this variable.</parameter>
  subroutine product_initialize(this, coefficient, targets, multinomials)
    class(Product) :: this
    integer, intent(in) :: coefficient
    integer, allocatable, intent(in) :: targets(:)
    integer, allocatable, intent(inout) :: multinomials(:,:)

    integer :: i

    this%coefficient = coefficient
    this%targets = targets
    allocate(this%multinoms(size(multinomials,1)))
    do i=1, size(multinomials,1)
       call this%multinoms(i)%initialize(multinomials(i,1), multinomials(i,2))
    end do
    !Change the multinomial array to be embedded inside of the product class.
    call move_alloc(multinomials, this%multinom_rd)

  end subroutine product_initialize

  !!<summary>Returns the coefficient of the term with the target exponents if all the multinomials
  !!in the product were expanded and had their terms collected.</summary>
  integer function coeff(this)
    class(Product), intent(in) :: this

    integer, allocatable :: possibles(:,:)!possibles(size(this%multinoms,1), size(this%targets,1))
    integer :: i, j, z, coeffs, k, t 
    integer, allocatable :: seq0(:,:), expandedseq(:,:)
    integer :: seq(size(this%multinoms, 1))
    class(vararray2d), allocatable, target :: mnseq(:)
    class(Sequence), pointer :: varseq
    class(vararray), allocatable :: passin(:)
    class(vararray2d), pointer :: tempseq

    print *, "in coeff"

    allocate(possibles(size(this%multinoms,1), size(this%targets,1)))
    allocate(mnseq(size(this%multinoms, 1)))

    !If this is an isolated multinomial, we only need to check the coefficient of the target term.
    print *, "here"
    if (size(this%multinoms,1) == 1) then
       if (all(this%multinoms(1)%power-this%targets .gt. 0)) then
          coeff = 0          
       else 
          coeff = this%multinoms(1)%nchoosekm(this%targets)*this%coefficient 
       end if           
    end if

    !Get a list of the possible exponents for each variable in each of the multinomials.
    !We start with the first variable and choose only those combinations of exponents
    !across *all* the multinomials that give the correct target exponent for that variable.
    print *, "here2"
    do i = 1, size(this%multinoms,1)
       possibles(i,:) = this%multinoms(i)%possible_powers
       print *, this%multinoms(i)%possible_powers
    end do

    print *, "here3"
    allocate(passin(size(possibles,1)))
    print *, "here3.1"
    do i=1, size(this%multinoms,1)
       print *, i
       call passin(i)%init(possibles(i,:))
    end do
    print *, "here3.5"
    print *, "here3.55"
!    print *, passin
    print *, size(possibles,2)
    allocate(seq0(,size(possibles,2)))
    call cproduct(passin, seq0)
    print *, "here3.6"
    !Next, we construct Sequence instances for each of the first variable compatible
    !possibilities and follow them through to the other variables.
    print *, "here4"
    coeffs = 0
    allocate(varseq)
    allocate(tempseq)
    do i = 1, size(seq0,1)
       !Each sequence calculated from the first variable has an entry for each multinomial
       !in this product. The Sequence instances construct smart sequences for the remaining
       !variables in each multinomial separately.
       print *,"here4.1"
       seq = seq0(i,:)
       print *,"here4.2"
       do j = 1, size(seq,1)
          allocate(expandedseq())
          call varseq%initialize(seq(i), possibles(i,:), 1, this%multinoms(i)%powersum, this%targets)
          call varseq%expand(expandedseq, 1, 1)
          print *, "here4.4"
          call mnseq(j)%init(expandedseq)
          !Clean up the varsequence and expanded sequence variable.
          call varseq%finalize()
          deallocate(expandedseq)
       end do
       print *, "here4.7"
       coeffs = coeffs + this%sum_sequence(mnseq)
       do j = 1, size(seq,1)    
          tempseq => mnseq(j)
          call tempseq%finalize()
       end do
    end do
    print *, "here5"
    coeff = coeffs*this%coefficient
  end function coeff

  !!<summary>Sums all the possible combinations of relevant sequences based of the variable sequence
  !!lists specified.</summary>
  !!<parameter name="mnseq">a list of possible variable sequences in each multinomial (one for each multinomial)
  !!that might contribute to the correct target variable.</parameter>
  integer function sum_sequence(this, mnseq)
    class(Product), intent(in) :: this
    class(vararray2d), intent(in) :: mnseq(size(this%multinoms,1))
    class(vararray), allocatable :: passin(:)
    
    integer :: i, j, z
    !!<local name="seqprod">Temporary variable for computing the product over adjacent term coefficients.</local>
    !!<local name="allseqs">Cartesian product of the possible sequences from each multinomial in the product.</local>
    !!<local name="seqlj">Re-constructed S_lj for a single result in 'allseqs' using the sequences from 'mnseq'.</local>
    !!<local name="expsum">Temporary sum of the exponents across all adjacent sequences in the product (per variable).</local>
    !!<local name="irange">Temporary range from 1..len(mnesq) for constructing the cartesian product.</local>
    integer :: seqprod
    integer, allocatable :: allseqs(:,:)
    integer :: seqlj(size(this%multinoms, 1), size(this%targets, 1)), expsum(size(this%targets, 1))
    integer, allocatable :: irange(:)

    sum_sequence = 0
    allocate(allseqs(size(mnseq),size(mnseq)))
    allocate(passin(size(mnseq,1)))
    !We could speed this up by keeping the ranges between loop iterations; if a two of the mnsequences have
    !the same number of elements, we could re-use the range. For now its not worth the extra complexity.
    do i=1, size(this%multinoms,1)
       allocate(irange(mnseq(i)%length))
       irange = (/( j, j=1, mnseq(i)%length, 1 )/)
       call passin(i)%init(irange)
       deallocate(irange)
    end do
    call cproduct(passin, allseqs)

    !Each row in allseqs now has a list of integer *indices* from mnseq that need to be used for the
    !coefficient calculations.
    do i = 1, size(allseqs, 1)
       !First we loop over the multinomials and extract the sequence corresponding to the index determined
       !by allseqs.
       do j = 1, size(this%multinoms, 1)
          seqlj(j,:) = mnseq(j)%items(allseqs(i,j),:)
       end do
       !Now that we have the sequences, we sum them up one variable at a time to make sure that they produce the
       !correct target term.
       do j=1, size(this%targets, 1)
          expsum(j) = sum(seqlj(j,:))
       end do
       if (all(expsum .eq. this%targets)) then
          seqprod = 1
          do z = 1, size(seqlj, 1)
             seqprod = seqprod*this%multinoms(z)%nchoosekm(seqlj(z,:))
          end do
          sum_sequence = sum_sequence + seqprod
       end if
    end do

  END function sum_sequence

  !!<summary>Sets up the multinomial.</summary>
  subroutine multinomial_initialize(this, power, exponent)
    class(Multinomial) :: this
    integer, intent(in) :: exponent, power
    integer :: j

    this%power = power
    this%exponent = exponent
    this%powersum = power*exponent
    this%possible_powers = (/(j,j=0,power*exponent+1,power)/)

  end subroutine multinomial_initialize

  !!<summary>Returns the number of different ways to partition an n-element
  !!set into disjoint subsets of sizes k1, ..., km.</summary>
  !!<parameter name="sequencei">an un-normed tuple of form (k1, k2, k3)."</parameter>
  integer function nchoosekm(this, sequencei)
    class(Multinomial), intent(in) :: this
    integer, intent(in) :: sequencei(:)
    integer :: j, nsum
    integer :: normseq(size(sequencei,1))

    if (.not. all(mod(sequencei, this%power) == 0)) then
       nchoosekm = 0
    else
       normseq = sequencei/this%power
       do j = 1, size(sequencei,1)
          nsum = sum(normseq(0:j+1))
          nchoosekm = nchoosekm*nchoosek(nsum, normseq(j))
       end do
    end if
  end function nchoosekm

  !!<summary>This implementation was taken from "Binomial Coefﬁcient Computation: Recursion 
  !!or Iteration?" by Yannis Manolopoulos, ACM SIGCSE Bulletin InRoads, Vol.34, No.4, 
  !!December 2002. http://delab.csd.auth.gr/papers/SBI02m.pdf It is supposed to be robust 
  !!against large, intermediate values and to have optimal complexity.</summary>
  integer function nchoosek(n,k)
    integer, intent(in) :: k, n
    integer :: i 

    if ((k < 0) .or. (k > n)) then
       nchoosek = 0
    end if
    if ((k == 0) .or. (n == 0)) then
       nchoosek = 1
    end if

    nchoosek = 1
    if (k < n-k) then
       do i = n, n-k, -1
          nchoosek = nchoosek*i/(n-i+1)
       end do
    else 
       do i = n, k, -1
          nchoosek = nchoosek*i/(n-i+1)
       end do
    end if
  end function nchoosek

  !!<summary>Uses a group and concentrations to find the number of unique arrangements as 
  !!described by polya..</summary>
  !!<parameter name="concentrations">Specifies a list of integers specifying how many of 
  !!each coloring should be present in each of the enumerated lists.</parameter>
  !!<parameter name="group">Group operations for permuting the colorings.</parameter>
  integer function polya(concentrations, group)
    integer, allocatable, intent(in) :: concentrations(:)
    integer, intent(in) :: group(:,:)

    !!<local name="polyndict">A dictionary of the unique products of multinomials in the polya calculation.</local>
    !!<local name="mult">A mx2 array of the powers and degeneracies of each multinomial in a product.</local>
    !!<local name="operation">Temporary variable for holding a single group operation.</local>
    !!<local name="visited">Temporary array for keeping track of which operation elements have already been
    !!visited when trying to determine the r-cycle structure of the group operation.</local>
    !!<local name="polynomials">A list of integers where the index is the power on the variables in the multinomial
    !!and the value is the degeneracy of that multinomial in the group operation.</local>
    class(Product), allocatable :: polyndict(:)
    integer, allocatable :: mult(:,:)
    integer :: operation(size(group,2)), visited(size(group,2)), polynomials(size(group,2))
    !!<local name="cursor">Points to the current element in the group operation that is being visited.</local>
    !!<local name="vindex">Temporary, used to find the location of the next operation to process when finding
    !!the r-cycles.</local>
    !!<local name="powers">Temporary, the value r for the r-cycle.</local>
    !!<local name="m">The number of unique multinomials in the product.</local>
    !!<local name="mindex">The index of the degeneracy for the lowest valued r-cycle left in the polynomials
    !!array. Corresponds to the value r.</local>
    !!<local name="dupindex">The index of the Product instance in polyndict that has the same multinomial set
    !!as the one for the current group operation.</local>
    !!<parameter name="pused">The number of unique Product instances created so far.</parameter>
    integer :: cursor(1), vindex(1), powers
    integer :: i, j, m, mindex, dupindex, pused = 0, b =0, d = 1
    
    allocate(polyndict(size(group,1)))
    do i=1, size(group,1)
       operation = group(i,:)
       visited = 0
       polynomials = 0
       do while(any(visited == 0))
          cursor = minloc(visited)
          vindex = minloc(visited)
          powers = 1
          visited(cursor(1)) = 1
          cursor = operation(cursor(1))
          do while(cursor(1) /= vindex(1))
             visited(cursor(1)) = 1
             powers = powers + 1
             cursor = operation(cursor(1))
          end do
          
          polynomials(powers) = polynomials(powers) + 1          
       end do
       !We now know the powers and degeneracies of all the r-cycles in the group operation.
       !Create the 2D multinomial array that defines r,d for each M_j^r.
       m = count(polynomials .gt. 0)
       allocate(mult(m,2))
       do j=1, m
          mindex = minval(polynomials, polynomials .gt. 0)
          mult(j,1) = mindex
          mult(j,2) = polynomials(mindex)
          polynomials(mindex) = 0
       end do
       !Find out if any of the existing products has this multinomial set already.
       dupindex = 0
       do j=1, pused
          if (all(polyndict(j)%multinom_rd .eq. mult)) then
             polyndict(j)%coefficient = polyndict(j)%coefficient + 1
             dupindex = j
             exit
          end if   
       end do

       if (dupindex == 0) then
          !We found no duplicate product in the set already, just initialize the next empty slot
          pused = pused + 1
          call polyndict(pused)%initialize(1, concentrations, mult)
       end if
       if (allocated(mult)) then 
          deallocate(mult)
       end if
    end do

    polya = 0
    do i = 1, pused
      polya = polyndict(i)%coeff() + polya
    end do    
    print *, "end"
  end function polya

end MODULE classes
