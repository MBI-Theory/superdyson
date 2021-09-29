module block_determinant
!
!  2014 Aug 29 - Initial version
!
  use accuracy
  use timer
  use block_matrices
  use lapack
  use printing

  implicit none
  private
  public block_det
  public rcsid_block_determinant
  !
  character(len=clen), save :: rcsid_block_determinant = "$Id: block_determinant.f90,v 1.5 2021/09/29 13:43:22 ps Exp $"
  !
  integer(ik), parameter :: verbose = 0

  contains

  function block_det(h,bmh) result(det)
    real(rk), intent(in)           :: h(:,:)   ! Matrix for which we need the determinant
    type(bm_structure), intent(in) :: bmh      ! Block structure of matrix h
    real(rk)                       :: det      ! Determinant
    !
    integer(ik) :: sz                    ! Overall size
    integer(ik) :: ib                    ! Block counter
    integer(ik) :: ir                    ! First-row counter
    integer(ik) :: bs                    ! Block size
    real(rk)    :: b_det                 ! Sub-block determinant
    integer(ik) :: row_list(bmh%max_row) ! Permulation list for the rows
    integer(ik) :: col_list(bmh%max_col) ! Permulation list for the columns
    integer(ik) :: r_pos, r_end          ! First and last positions of row indices in %rows
    integer(ik) :: c_pos, c_end          ! ditto, columns
    !
    sz = size(h,dim=1)
    if (sz/=size(h,dim=2)) stop 'block_determinant%block_det - dimensions mismatch'
    !
    if (bmh%n_rows>sz .or. bmh%n_cols>sz) then 
      write (out,"(/'block_det: Determinant of size ',i0)") sz
      call print_matrix(h,10,"e10.4")
      call bm_print(bmh,'block_det received bad descriptor')
      stop 'block_determinant%block_det - bad block descriptor (1)'
    end if
    ! 
    !  Try to make a quick exit without computing anything
    !
    if (bmh%n_rows<sz .or. bmh%n_cols<sz) then
      det = 0 ! Matrix is defective; don't bother making the calculation
      return
    end if
    screen_blocks: do ib=1,bmh%n_blocks
      if (bmh%block_data(1,ib)==bmh%block_data(2,ib)) cycle screen_blocks
      ! Matrix contains a defective sub-block; don't bother any further
      det = 0
      return
    end do screen_blocks
    !
    det = 1
    ir  = 1
    process_blocks: do ib=1,bmh%n_blocks
      bs    = bmh%block_data(1,ib)  ! Guaranteed to be the same as column size because of the loop above
      r_pos = bmh%block_data(3,ib) 
      c_pos = bmh%block_data(4,ib) 
      r_end = r_pos + bs - 1
      c_end = c_pos + bs - 1
      b_det = subblock_det(bs,bmh%rows(r_pos:r_end),bmh%cols(c_pos:c_end))
      det   = det * b_det
      !
      !  Unless we do something at this point, the determinant we calculate is for a block-diagonal
      !  matrix, with block appearing in the order they are listed in (bmh). Since this is not
      !  necessarity true, we need a correction for row/column interchange
      !
      call store_permutation(row_list,ir,bmh%rows(r_pos:r_end))
      call store_permutation(col_list,ir,bmh%cols(c_pos:c_end))
      !
      ir    = ir + bs
    end do process_blocks
    det = det * permutation_parity(row_list) * permutation_parity(col_list)
    ! ???? DEBUG
    ! b_det = linpack_determinant(h)
    ! if (abs(det-b_det)>spacing(1._rk)) then
    !   write (out,"('det= ',g26.16,' block det= ',g26.16,' err= ',g26.16)") b_det, det, b_det-det
    ! end if
    ! ???? DEBUG
    return
    !
    contains

    function permutation_parity(list) result(factor)
      integer(ik), intent(in) :: list(:) ! Elements of the original list, in canonical order
      real(rk)                :: factor
      !
      logical     :: taken(size(list)) ! Element taken care of
      integer(ik) :: head              ! Start of a permutation loop
      integer(ik) :: pos               ! Current position of a permulation loop
      integer(ik) :: loop_len          ! Length of a n interchange loop
      !
      taken  = .false.
      factor = 1
      scan_start: do head=1,size(list)
        if (taken(head)) cycle scan_start
        taken(head) = .true.
        pos         = list(head)
        loop_len    = 1
        follow_loop: do while(pos/=head)
          taken(pos) = .true.
          pos        = list(pos)
          loop_len   = loop_len + 1
        end do follow_loop
        if (mod(loop_len,2)==0) factor = -factor
      end do scan_start
    end function permutation_parity

    subroutine store_permutation(list,pos,elem)
      integer(ik), intent(inout) :: list(:)  ! Permutation list
      integer(ik), intent(in)    :: pos      ! Current position in the list
      integer(ik), intent(in)    :: elem(:)  ! New elements to add
      !
      integer(ik) :: ne
      !
      ne = size(elem)
      if (pos+ne-1>size(list)) stop 'block_determinant%block_det%store_permulation - list() overflows'
      list(pos:pos+ne-1) = elem
    end subroutine store_permutation

    function subblock_det(bs,rows,cols) result(b_det)
      integer(ik), intent(in) :: bs      ! Block size
      integer(ik), intent(in) :: rows(:) ! List of rows
      integer(ik), intent(in) :: cols(:) ! List of columns
      real(rk)                :: b_det
      !
      real(rk) :: tmp_h(bs,bs)
      !
      select case (bs)
        case (1)
          b_det = h(rows(1),cols(1))
        case (2)
          b_det = h(rows(1),cols(1))*h(rows(2),cols(2)) &
                - h(rows(1),cols(2))*h(rows(2),cols(1))
        case (3)
          b_det = h(rows(1),cols(1))*h(rows(2),cols(2))*h(rows(3),cols(3)) &
                + h(rows(1),cols(2))*h(rows(2),cols(3))*h(rows(3),cols(1)) &
                + h(rows(1),cols(3))*h(rows(2),cols(1))*h(rows(3),cols(2)) &
                - h(rows(1),cols(3))*h(rows(2),cols(2))*h(rows(3),cols(1)) &
                - h(rows(1),cols(2))*h(rows(2),cols(1))*h(rows(3),cols(3)) &
                - h(rows(1),cols(1))*h(rows(2),cols(3))*h(rows(3),cols(2))
        case default
          tmp_h = h(rows,cols)
          b_det = linpack_determinant_trash_input(tmp_h)
      end select
    end function subblock_det
  end function block_det

end module block_determinant
