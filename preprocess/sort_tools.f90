module sort_tools
!
!  Very simple sorting tools. We are using merge sort, which has guaranteed
!  O(N log(N)) scaling and optimal (sequential) memory access pattern.
!  The code is modification of a subroutine from "sparse.f90"
!
  use accuracy
  implicit none
  private
  public sort, sort_unique, order_keys
  public rcsid_sort_tools

  interface sort
     module procedure sort_integer
     module procedure sort_real
     module procedure sort_and_order_integer
     module procedure sort_and_order_real
!*qd module procedure sort_quad
  end interface sort

  interface sort_unique
    module procedure sort_unique_integer
  end interface sort_unique

  interface order_keys
     module procedure order_integer
     module procedure order_real
!*qd module procedure order_quad
  end interface order_keys

  character(len=clen), save :: rcsid_sort_tools = "$Id: sort_tools.f90,v 1.4 2021/09/29 13:43:22 ps Exp $"

  integer(ik), parameter :: verbose = 0

  contains
  !
  !  Externally callable subroutines
  !
  subroutine sort_integer(key)
    integer(ik), intent(inout) :: key(:)  ! List of keys to sort
    !
    integer(ik) :: n_in
    integer(ik) :: order(size(key))       ! The contents of the index array are immaterial
    !
    n_in  = size(key)
    order = 0
    call sort_integer2(n_in,key,order)
  end subroutine sort_integer

  subroutine sort_real(key)
    real(rk), intent(inout) :: key(:)     ! List of keys to sort
    !
    integer(ik) :: n_in
    integer(ik) :: order(size(key))       ! The contents of the index array are immaterial
    !
    n_in  = size(key)
    order = 0
    call sort_real2(n_in,key,order)
  end subroutine sort_real

  subroutine sort_quad(key)
    real(xrk), intent(inout) :: key(:)     ! List of keys to sort
    !
    integer(ik) :: n_in
    integer(ik) :: order(size(key))       ! The contents of the index array are immaterial
    !
    n_in  = size(key)
    order = 0
    call sort_quad2(n_in,key,order)
  end subroutine sort_quad

  function sort_unique_integer(key) result (nkey)
    integer(ik), intent(inout) :: key(:)  ! In:  List of keys to sort
                                          ! Out: List of keys in numerically ascending oreder,
                                          !      with duplicates removed.
    integer(ik)                :: nkey    ! Number of unique keys present in the key() array
    !
    nkey = sort_unique_integer2(size(key),key)
  end function sort_unique_integer

  subroutine order_integer(key,order)
    integer(ik), intent(in)  :: key  (:)  ! List of keys to sort
    integer(ik), intent(out) :: order(:)  ! Sorting order of the keys
    !
    integer(ik) :: n_in, i
    integer(ik) :: key_copy(size(key))
    !
    n_in = size(key)
    if (n_in/=size(order)) stop 'sort_tools%order_integer - called with inconsistent array sizes'
    fill_order: do i=1,n_in
      order(i) = i
    end do fill_order
    key_copy = key
    call sort_integer2(n_in,key_copy,order)
  end subroutine order_integer

  subroutine order_real(key,order)
    real(rk), intent(in)     :: key  (:)  ! List of keys to sort
    integer(ik), intent(out) :: order(:)  ! Sorting order of the keys
    !
    integer(ik) :: n_in, i
    real(rk)    :: key_copy(size(key))
    !
    n_in = size(key)
    if (n_in/=size(order)) stop 'sort_tools%order_real - called with inconsistent array sizes'
    fill_order: do i=1,n_in
      order(i) = i
    end do fill_order
    key_copy = key
    call sort_real2(n_in,key_copy,order)
  end subroutine order_real

  subroutine order_quad(key,order)
    real(xrk), intent(in)    :: key  (:)  ! List of keys to sort
    integer(ik), intent(out) :: order(:)  ! Sorting order of the keys
    !
    integer(ik) :: n_in, i
    real(xrk)   :: key_copy(size(key))
    !
    n_in = size(key)
    if (n_in/=size(order)) stop 'sort_tools%order_quad - called with inconsistent array sizes'
    fill_order: do i=1,n_in
      order(i) = i
    end do fill_order
    key_copy = key
    call sort_quad2(n_in,key_copy,order)
  end subroutine order_quad

  subroutine sort_and_order_integer(key,order)
    integer(ik), intent(inout) :: key  (:)  ! List of keys to sort
    integer(ik), intent(out)   :: order(:)  ! Sorting order of the keys
    !
    integer(ik) :: n_in, i
    !
    n_in = size(key)
    if (n_in/=size(order)) stop 'sort_tools%sort_and_order_integer - called with inconsistent array sizes'
    fill_order: do i=1,n_in
      order(i) = i
    end do fill_order
    call sort_integer2(n_in,key,order)
  end subroutine sort_and_order_integer

  subroutine sort_and_order_real(key,order)
    real(rk), intent(inout)  :: key  (:)  ! List of keys to sort
    integer(ik), intent(out) :: order(:)  ! Sorting order of the keys
    !
    integer(ik) :: n_in, i
    !
    n_in = size(key)
    if (n_in/=size(order)) stop 'sort_tools%sort_and_order_real - called with inconsistent array sizes'
    fill_order: do i=1,n_in
      order(i) = i
    end do fill_order
    call sort_real2(n_in,key,order)
  end subroutine sort_and_order_real
  !
  !  Internal subroutines beyond this point
  !
  recursive subroutine sort_real2(n_in,key,val)
    integer(ik), intent(in)  :: n_in    ! Size of the 
    real(rk), intent(inout)  :: key(:)  ! List of keys to sort
    integer(ik), intent(out) :: val(:)  ! Payload to sort
    !
    real(rk)    :: k_left (2+n_in/2)  ! Temporary buffers for keys
    real(rk)    :: k_right(2+n_in/2)  
    integer(ik) :: v_left (2+n_in/2)  ! Temporary buffers for payload
    integer(ik) :: v_right(2+n_in/2)
    integer(ik) :: n_left, n_right    ! Number of left/right keys
    integer(ik) :: p_left, p_right    ! Positions of the left/right sections
    integer(ik) :: p_out              ! Position of the output key
    !
    if (n_in<=1) then
      return ! Already sorted
    end if
    !
    !  Partition the input array
    !
    p_left  = 1
    n_left  = n_in/2
    p_right = p_left + n_left
    n_right = n_in - n_left
    !
    if (verbose>=2) then
      write (out,"('sort_real2: Asked to sort ',i10,' elements')") n_in
      write (out,"((t4,8(g16.8,1x)))") key(1:n_in)
      write (out,"('sort_real2:  Left: ',i10,' elements from ',i10)") n_left, p_left
      write (out,"('sort_real2: Right: ',i10,' elements from ',i10)") n_right, p_right
    end if
    !
    !  Copy out the keys and values
    !
    k_left (1:n_left ) = key( p_left:p_left+n_left-1  )
    v_left (1:n_left ) = val( p_left:p_left+n_left-1  )
    k_right(1:n_right) = key(p_right:p_right+n_right-1)
    v_right(1:n_right) = val(p_right:p_right+n_right-1)
    !
    !  Merge sort left and right arrays
    !
    call sort_real2(n_left, k_left, v_left )
    call sort_real2(n_right,k_right,v_right)
    !
    !  Add impossibly large, key values as sentinels at the right.
    !  Having the sentinels lets us simplify merging logic.
    !
    k_left (n_left+1)  = huge(key)
    k_right(n_right+1) = huge(key)
    !
    !  Merge pre-sorted arrays. 
    !  At this point, left and right sections are already sorted. 
    !
    p_out   = 0
    p_left  = 1
    p_right = 1
    merge_left_right: do while (p_left<=n_left .or. p_right<=n_right)
      if (k_left(p_left)<=k_right(p_right)) then
        !
        !  Left key is smaller or equal; copy it
        !
        p_out      = p_out + 1
        key(p_out) = k_left(p_left)
        val(p_out) = v_left(p_left)
        p_left     = p_left + 1
      else
        !
        !  Right key is smaller; copy it
        !
        p_out      = p_out + 1
        key(p_out) = k_right(p_right)
        val(p_out) = v_right(p_right)
        p_right    = p_right + 1
      end if
    end do merge_left_right
    !
    if (p_out/=n_in) stop 'sort_tools%sort_real2 - Logical error in merge'
  end subroutine sort_real2

  recursive subroutine sort_quad2(n_in,key,val)
    integer(ik), intent(in)  :: n_in    ! Size of the 
    real(xrk), intent(inout) :: key(:)  ! List of keys to sort
    integer(ik), intent(out) :: val(:)  ! Payload to sort
    !
    real(xrk)   :: k_left (2+n_in/2)  ! Temporary buffers for keys
    real(xrk)   :: k_right(2+n_in/2)  
    integer(ik) :: v_left (2+n_in/2)  ! Temporary buffers for payload
    integer(ik) :: v_right(2+n_in/2)
    integer(ik) :: n_left, n_right    ! Number of left/right keys
    integer(ik) :: p_left, p_right    ! Positions of the left/right sections
    integer(ik) :: p_out              ! Position of the output key
    !
    if (n_in<=1) then
      return ! Already sorted
    end if
    !
    !  Partition the input array
    !
    p_left  = 1
    n_left  = n_in/2
    p_right = p_left + n_left
    n_right = n_in - n_left
    !
    if (verbose>=2) then
      write (out,"('sort_quad2: Asked to sort ',i10,' elements')") n_in
      write (out,"((t4,8(g16.8,1x)))") key(1:n_in)
      write (out,"('sort_quad2:  Left: ',i10,' elements from ',i10)") n_left, p_left
      write (out,"('sort_quad2: Right: ',i10,' elements from ',i10)") n_right, p_right
    end if
    !
    !  Copy out the keys and values
    !
    k_left (1:n_left ) = key( p_left:p_left+n_left-1  )
    v_left (1:n_left ) = val( p_left:p_left+n_left-1  )
    k_right(1:n_right) = key(p_right:p_right+n_right-1)
    v_right(1:n_right) = val(p_right:p_right+n_right-1)
    !
    !  Merge sort left and right arrays
    !
    call sort_quad2(n_left, k_left, v_left )
    call sort_quad2(n_right,k_right,v_right)
    !
    !  Add impossibly large, key values as sentinels at the right.
    !  Having the sentinels lets us simplify merging logic.
    !
    k_left (n_left+1)  = huge(key)
    k_right(n_right+1) = huge(key)
    !
    !  Merge pre-sorted arrays. 
    !  At this point, left and right sections are already sorted. 
    !
    p_out   = 0
    p_left  = 1
    p_right = 1
    merge_left_right: do while (p_left<=n_left .or. p_right<=n_right)
      if (k_left(p_left)<=k_right(p_right)) then
        !
        !  Left key is smaller or equal; copy it
        !
        p_out      = p_out + 1
        key(p_out) = k_left(p_left)
        val(p_out) = v_left(p_left)
        p_left     = p_left + 1
      else
        !
        !  Right key is smaller; copy it
        !
        p_out      = p_out + 1
        key(p_out) = k_right(p_right)
        val(p_out) = v_right(p_right)
        p_right    = p_right + 1
      end if
    end do merge_left_right
    !
    if (p_out/=n_in) stop 'sort_tools%sort_quad2 - Logical error in merge'
  end subroutine sort_quad2

  recursive subroutine sort_integer2(n_in,key,val)
    integer(ik), intent(in)    :: n_in    ! Size of the 
    integer(ik), intent(inout) :: key(:)  ! List of keys to sort
    integer(ik), intent(out)   :: val(:)  ! Payload to sort
    !
    integer(ik) :: k_left (2+n_in/2)  ! Temporary buffers for keys
    integer(ik) :: k_right(2+n_in/2)  
    integer(ik) :: v_left (2+n_in/2)  ! Temporary buffers for payload
    integer(ik) :: v_right(2+n_in/2)
    integer(ik) :: n_left, n_right    ! Number of left/right keys
    integer(ik) :: p_left, p_right    ! Positions of the left/right sections
    integer(ik) :: p_out              ! Position of the output key
    !
    if (n_in<=1) then
      return ! Already sorted
    end if
    !
    !  Partition the input array
    !
    p_left  = 1
    n_left  = n_in/2
    p_right = p_left + n_left
    n_right = n_in - n_left
    !
    if (verbose>=2) then
      write (out,"('sort_integer2: Asked to sort ',i10,' elements')") n_in
      write (out,"((t4,8(i12,1x)))") key(1:n_in)
      write (out,"('sort_integer2:  Left: ',i10,' elements from ',i10)") n_left, p_left
      write (out,"('sort_integer2: Right: ',i10,' elements from ',i10)") n_right, p_right
    end if
    !
    !  Copy out the keys and values
    !
    k_left (1:n_left ) = key( p_left:p_left+n_left-1  )
    v_left (1:n_left ) = val( p_left:p_left+n_left-1  )
    k_right(1:n_right) = key(p_right:p_right+n_right-1)
    v_right(1:n_right) = val(p_right:p_right+n_right-1)
    !
    !  Merge sort left and right arrays
    !
    call sort_integer2(n_left, k_left, v_left )
    call sort_integer2(n_right,k_right,v_right)
    !
    !  Add impossibly large, key values as sentinels at the right.
    !  Having the sentinels lets us simplify merging logic.
    !
    k_left (n_left+1)  = huge(key)
    k_right(n_right+1) = huge(key)
    !
    !  Merge pre-sorted arrays. 
    !  At this point, left and right sections are already sorted. 
    !
    p_out   = 0
    p_left  = 1
    p_right = 1
    merge_left_right: do while (p_left<=n_left .or. p_right<=n_right)
      if (k_left(p_left)<=k_right(p_right)) then
        !
        !  Left key is smaller or equal; copy it
        !
        p_out      = p_out + 1
        key(p_out) = k_left(p_left)
        val(p_out) = v_left(p_left)
        p_left     = p_left + 1
      else
        !
        !  Right key is smaller; copy it
        !
        p_out      = p_out + 1
        key(p_out) = k_right(p_right)
        val(p_out) = v_right(p_right)
        p_right    = p_right + 1
      end if
    end do merge_left_right
    !
    if (p_out/=n_in) stop 'sort_tools%sort_integer2 - Logical error in merge'
  end subroutine sort_integer2
  !
  !  sort_unique_integer2 is a modification of mergeSort2 from sparse.f90, with
  !  payload deleted.
  !
  recursive function sort_unique_integer2(n_in,key) result(n_out)
    integer(ik), intent(in)    :: n_in   ! Number of elements to sort
    integer(ik), intent(inout) :: key(:) ! Keys
    integer(ik)                :: n_out  ! Number of surviving keys
    !
    integer(ik) :: k_left (2+n_in/2)  ! Temporary buffers for keys
    integer(ik) :: k_right(2+n_in/2)  
    integer(ik) :: c_left, c_right    ! Initial number of left/right keys
    integer(ik) :: n_left, n_right    ! Number of surviving keys in left/right arrays
    integer(ik) :: p_left, p_right    ! Positions of the left/right sections
    !
    if (n_in<=1) then
      n_out = n_in
      return ! Already sorted
    end if
    !
    !  Partition the input array
    !
    p_left  = 1
    c_left  = n_in/2
    p_right = p_left + c_left
    c_right = n_in - c_left
    !
    if (verbose>3) then
      write (out,"('sort_unique_integer2: Asked to sort ',i10,' elements')") n_in
      write (out,"((t12,8(i10,1x)))") key(1:n_in)
      write (out,"('sort_unique_integer2:  Left: ',i10,' elements from ',i10)") c_left, p_left
      write (out,"('sort_unique_integer2: Right: ',i10,' elements from ',i10)") c_right, p_right
      call flush(out)
    end if
    !
    !  Copy out the keys and values
    !
    k_left (1:c_left ) = key( p_left:p_left+c_left-1  )
    k_right(1:c_right) = key(p_right:p_right+c_right-1)
    !
    !  Merge sort left and right arrays
    !
    n_left  = sort_unique_integer2(c_left, k_left )
    n_right = sort_unique_integer2(c_right,k_right)
    !
    if (verbose>3) then
      write (out,"('sort_unique_integer2: Survivors left: ',i10,' right: ',i10)") n_left, n_right
      call flush(out)
    end if
    !
    !  Add impossibly large, key values as sentinels at the right
    !
    k_left (n_left+1)  = huge(1_ik)
    k_right(n_right+1) = huge(1_ik)
    !
    !  Merge pre-sorted arrays. At this point, left and right sections are already
    !  sorted, and contain no duplicate keys.
    !
    n_out   = 0
    p_left  = 1
    p_right = 1
    merge_left_right: do 
      if (k_left(p_left)==k_right(p_right)) then
        !
        !  Have we reached the end of both arrays? Thanks to sentinels, only
        !  one index need to be tested.
        !
        if (p_left>n_left) exit merge_left_right
        !
        !  Key are identical - merge and consume both values
        !
        n_out      = n_out + 1
        key(n_out) = k_left(p_left)
        p_left     = p_left + 1
        p_right    = p_right + 1
        cycle merge_left_right
      end if
      if (k_left(p_left)<k_right(p_right)) then
        !
        !  Left key is smaller - copy the left value
        !
        n_out      = n_out + 1
        key(n_out) = k_left(p_left)
        p_left     = p_left + 1
      else
        !
        !  Right key is smaller - copy the right value
        !
        n_out      = n_out + 1
        key(n_out) = k_right(p_right)
        p_right    = p_right + 1
      end if
    end do merge_left_right
    !
    if (verbose>3) then
      write (out,"('sort_unique_integer2: Final merge gave ',i10)") n_out
    end if
  end function sort_unique_integer2

end module sort_tools
