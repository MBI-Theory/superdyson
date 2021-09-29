module block_matrices
!
!  Routines for handling block matrices.
!
!  Matrices are still stored as a standard Fortran rectangular matrix,
!  with no attempt to use block structure to save storage space.
!
!  Our goal is to keep track of non-zero blocks of a matrix, without
!  having to do a full matrix sweep every time an operation is performed
!  on a matrix. Calling routines may use this information if they wish;
!  however, using standard dense matrix operations is always an option.
!
!  We do not make any assumptions on the structure of the blocks; in particular,
!  they do not have to be contiguous in the space of matrix indices.
!
!  Every time a block-matrix descriptor is created (e.g. usign bm_scan_matrix etc.),
!  it must be destroyed with bm_destroy once it is no longer needed. If this
!  is not done, dragons may occur.
!
!  Note that we have a fixed, compile-time limit on the number of rows/columns
!  in a matrix. This is an unfortunate consequence of the rather unreasonable
!  overhead of memory allocation in gfortran: with dynamical structures, we
!  simply spend all the time in malloc() and free().
!
!  Note that a number of routines in superdyson, tr_rdm, and sd_core 
!  make extensive use of the internal structures of bm_structure,
!  including constructing new bm_structure's for some cases.
!
!  Be extra careful if changing semantics of any part of this module!
!
!  Routines:
!
!    bm_destroy        - deletes block-matrix descriptor
!    bm_print          - (debug) reports on the block structure 
!    bm_verify         - (debug) checks descriptor for consistency with the matrix
!    bm_scan_matrix    - analyzes a matrix, and returns a new block descriptor.
!    bm_find_zero_row  - find the first zero row in a given matrix
!    bm_drop_row       - update matrix structure upon row deletion
!    bm_drop_column    - update matrix structure upon column deletion
!    bm_substitute_row - update matrix structure upon row substitution
!    bm_add_block      - service routine; adds another block to the descriptor
!    bm_intersect      - service routine; true if two sorted lists intersect
!    bm_intersection   - service routine; returns intersection of two sorted lists
!
!  2014 Aug 28 - Initial version, using bits extracted from block_diag.f90
!                At some point, we should convert block_diag.f90 to use
!                this module instead of trying to find block structure
!                itself.
!  2014 Sep 01 - Rewritten to eschew dynamical allocation
!
  use accuracy
  use timer
  use printing
  use sort_tools
  implicit none
  private
  public bm_structure
  public bm_scan_matrix, bm_destroy, bm_print, bm_verify
  public bm_find_zero_row, bm_drop_row, bm_drop_column, bm_substitute_row
  public bm_add_block, bm_intersect, bm_intersection
  public rcsid_block_matrices
  !
  character(len=clen), save :: rcsid_block_matrices = "$Id: block_matrices.f90,v 1.9 2021/09/29 13:43:22 ps Exp $"
  !
  integer(ik), parameter :: max_nrows = 1024_ik ! Increase if necessary
  !
  type bm_structure
    integer(ik)    :: n_rows    = 0             ! Number of non-zero rows. Note that the first unused position
                                                ! in rows() is rows(n_rows+1)
    integer(ik)    :: n_cols    = 0             ! Number of non-zero columns
    integer(ik)    :: n_blocks  = -1            ! Number of non-empty blocks;
                                                ! -1 indicates an un-unitialized structure
    integer(ik)    :: max_row   = 0             ! Largest row index in the matrix (including zero blocks!)
    integer(ik)    :: max_col   = 0             ! Largest column index in the matrix
    integer(ik)    :: block_data(4,max_nrows)   ! 1 = number of rows in each block
                                                ! 2 = number of columns in each block
                                                ! 3 = starting position in rows() for the row list
                                                ! 4 = starting position in cols() for the column list
    integer(ik)    :: rows(max_nrows)           ! List of rows within each block
    integer(ik)    :: cols(max_nrows)           ! List of columns within each block
  end type bm_structure
  !
  contains
  !
  subroutine bm_destroy(bms)
    type(bm_structure), intent(inout) :: bms ! Description to destroy
    !
    bms%n_blocks = -1
  end subroutine bm_destroy

  subroutine bm_print(bms,text)
    type(bm_structure), intent(in) :: bms  ! Block descriptor
    character(len=*), intent(in)   :: text ! Descriptive text to associate with the report
    !
    integer(ik) :: ib, i1, ie
    integer(ik) :: rp, cp, rc, cc
    !
    if (bms%n_blocks/=-1) then
      write (out,"(a,': max. number of rows: ',i0)") text, bms%max_row
      write (out,"(a,': max. number of colums: ',i0)") text, bms%max_col
      write (out,"(a,': found ',i0,' blocks')") text, bms%n_blocks
      write (out,"(a,': ',i0,' non-zero rows')") text, bms%n_rows
      write (out,"(a,': ',i0,' non-zero columns')") text, bms%n_cols
      print_blocks: do ib=1,bms%n_blocks
        rc = bms%block_data(1,ib)  ! Rows count
        cc = bms%block_data(2,ib)  ! Columns count
        rp = bms%block_data(3,ib)  ! Rows position
        cp = bms%block_data(4,ib)  ! Columns position
        print_rows: do i1=1,rc,20
          ie = min(i1+19,rc)
          write (out,"(a,': block ',i4,' rows    ',20(1x,i6))") text, ib, bms%rows(rp+i1-1:rp+ie-1)
        end do print_rows
        !
        print_cols: do i1=1,cc,20
          ie = min(i1+19,cc)
          write (out,"(a,': block ',i4,' columns ',20(1x,i6))") text, ib, bms%cols(cp+i1-1:cp+ie-1)
        end do print_cols
      end do print_blocks
    else
      write (out,"(a,': descriptor is not initialized.')") text
    end if
    write (out,"()")
  end subroutine bm_print

  subroutine bm_verify(mat,bms,text,eps)
    real(rk), intent(in)           :: mat(:,:) ! Matrix 
    type(bm_structure), intent(in) :: bms      ! Block structure to be checked for consistency with matrix
    character(len=*), intent(in)   :: text     ! Descriptive text
    real(rk), intent(in), optional :: eps      ! Threshold for treating elements as zero;
                                               ! Missing or negative value will be replaces by the
                                               ! default threshold of 1000*spacing(maxval(abs(a))
    !
    integer(ik)  :: ib                                    ! Block index
    integer(ik)  :: c_pos, r_pos, c_end, r_end            ! Block indices in bms%rows/%cols
    real(rk)     :: thr                                   ! zero threshold
    logical      :: mask(size(mat,dim=1),size(mat,dim=2)) ! Areas of mat which are supposed to be zero
    integer(ik)  :: bigboy(2)
    !
    thr = -1
    if (present(eps)) then
      if (eps>=0) thr = eps
    end if
    if (thr<0) thr = 1000._rk * spacing(maxval(abs(mat)))
    !
    !  We do not really care (too much) if blocks we think are non-zero turn out to be
    !  zero; however, it is unacceptable to have non-zeros in parts of the matrix we
    !  think are nil.
    !
    !  The other cardinal sin is to have overlapping blocks
    !
    mask = .true.
    check_block: do
      assemble_mask: do ib=1,bms%n_blocks
        r_pos = bms%block_data(3,ib)
        c_pos = bms%block_data(4,ib)
        r_end = bms%block_data(1,ib) + r_pos - 1
        c_end = bms%block_data(2,ib) + c_pos - 1
        if (.not.all(mask(bms%rows(r_pos:r_end),bms%cols(c_pos:c_end)))) then
          write (out,"(a,' block ',i0,' overlaps one of the previosly defined blocks.')") text, ib
          exit check_block
        end if
        mask(bms%rows(r_pos:r_end),bms%cols(c_pos:c_end)) = .false.
      end do assemble_mask
      !
      !  All elements of mat for which mask is still true must be zero.
      !
      if (any((abs(mat)>thr).and.mask)) then
        write (out,"(a,': non-zero element(s) in the zero part of matrix')") text
        bigboy = maxloc(abs(mat),mask=mask)
        write (out,"(a,': largest element at indices',2(1x,i0))") text, bigboy
        exit check_block
      end if
      !
      !  If we come here, all is peachy
      !
      return 
    end do check_block
    !
    !  We come here if something is wrong; print matrix structure, matrix itself, and get out
    !
    call bm_print(bms,text)
    call print_matrix(mat,15,"e15.8")
    call flush(out)
    stop 'block_matrix%bm_verify - verify failure'
  end subroutine bm_verify

  subroutine bm_scan_matrix(mat,bms,eps)
    real(rk), intent(in)              :: mat(:,:) ! Matrix to be analyzed
    type(bm_structure), intent(inout) :: bms      ! Block structure of the matrix
    real(rk), intent(in), optional    :: eps      ! Threshold for treating elements as zero;
                                                  ! Missing or negative value will be replaces by the
                                                  ! default threshold of 1000*spacing(maxval(abs(a))
    !
    integer(ik)              :: nrow, ncol     ! Number of rows and columns in matrix a
    integer(ik)              :: max_blk        ! Max number of blocks possible
    integer(ik)              :: blk_nrow       ! Number of rows in a block
    integer(ik)              :: blk_ncol       ! ditto, columns
    integer(ik)              :: ir             ! Row indix
    integer(ik)              :: ib             ! Block index
    integer(ik)              :: alloc
    real(rk)                 :: thr            ! zero threshold
    logical, allocatable     :: row_avail(:)   ! Rows still not processed
    logical, allocatable     :: col_avail(:)   ! ditto, columns
    integer(ik), allocatable :: blk_rows(:)    ! Rows within the current block
    integer(ik), allocatable :: blk_cols(:)    ! ditto, columns
    logical                  :: block_grown    ! Indicates that the block has grown, so that an extra
                                               ! connections pass may be required
    integer(ik)              :: n_block        ! Number of blocks found
    integer(ik), allocatable :: n_row(:)       ! Count of rows in each block
    integer(ik), allocatable :: n_col(:)       ! Count of columns in each block
    integer(ik), allocatable :: rows(:,:)      ! Rows in each block
    integer(ik), allocatable :: cols(:,:)      ! Columns in each block
    !
    call TimerStart('bm_scan_matrix')
    if (bms%n_blocks/=-1) stop 'block_matrices%bm_scan_matrix - output descriptor already in use?!'
    !
    nrow  = size(mat,dim=1)
    ncol  = size(mat,dim=2)
    !
    thr = -1
    if (present(eps)) then
      if (eps>=0) thr = eps
    end if
    if (thr<0) thr = 1000._rk * spacing(maxval(abs(mat)))
    !
    max_blk = min(nrow,ncol)
    allocate (row_avail(nrow),col_avail(ncol),blk_rows(nrow),blk_cols(ncol), &
              n_row(max_blk),n_col(max_blk),rows(nrow,max_blk),cols(ncol,max_blk),stat=alloc)
    if (alloc/=0) stop 'block_matrices%bm_scan_matrix - allocation failed (1)'
    !
    !  In order to simplify block extraction, mark all zero rows and columns as unavailable. 
    !
    row_avail = any(abs(mat(:,:))>thr,dim=2)
    col_avail = any(abs(mat(:,:))>thr,dim=1)
    !
    n_block = 0
    block_loop: do while(any(row_avail) .and. any(col_avail))
      !
      !  The logic of block scanning below is a little complicated: 
      !  since we are working with a general matrix, adding a row and 
      !  connected columns may introduce connection to a row which was 
      !  previously unconnected to the block. To catch this situation, 
      !  we have to restart the scan each time we've added some rows
      !
      blk_nrow    = 0 
      blk_ncol    = 0
      block_grown = .true.
      restart_scan_rows: do while (block_grown)
        block_grown = .false.
        scan_unused_rows: do ir=1,nrow
          if (.not.row_avail(ir)) cycle scan_unused_rows
          if (.not.row_connected(ir)) cycle scan_unused_rows
          call add_row(ir)
          block_grown = .true.
        end do scan_unused_rows
      end do restart_scan_rows
      if (blk_nrow==0 .or. blk_ncol==0) stop 'block_matrices%bm_scan_matrix - logic error: no block found'
      !
      !  At this point, we have the list of rows and columns in the block.
      !  Sort and remember them.
      !
      call sort(blk_rows(1:blk_nrow))
      call sort(blk_cols(1:blk_ncol))
      n_block = n_block + 1
      ! Be paranoid
      if (any((/size(n_row),size(n_col),size(rows,dim=2),size(cols,dim=2)/)<n_block) .or. &
          size(rows,dim=1)<blk_nrow .or. size(cols,dim=1)<blk_ncol ) then
        stop 'block_matrices%bm_scan_matrix - dimensions blown'
      end if
      !
      n_row(n_block) = blk_nrow
      n_col(n_block) = blk_ncol
      rows(1:blk_nrow,n_block) = blk_rows(1:blk_nrow)
      cols(1:blk_ncol,n_block) = blk_cols(1:blk_ncol)
      !
    end do block_loop
    if (any(row_avail) .or. any(col_avail)) stop 'block_matrices%bm_scan_matrix - logic error: unclaimed rows/columns'
    !
    !  Copy results into the output structure. Also count non-zero rows and columns
    !
    ! bms%n_blocks = n_block
    bms%n_blocks = 0
    bms%n_rows   = 0
    bms%n_cols   = 0
    bms%max_row  = nrow
    bms%max_col  = ncol
    copy_blocks: do ib=1,n_block
      blk_nrow = n_row(ib)
      blk_ncol = n_col(ib)
      call bm_add_block(bms,rows(1:blk_nrow,ib),cols(1:blk_ncol,ib))
    end do copy_blocks
    !
    !  Clean up.
    !
    deallocate (row_avail,col_avail,blk_rows,blk_cols,n_row,n_col,rows,cols)
    !
    call TimerStop('bm_scan_matrix')
    !
    contains
    !
    !  Return true if row (ir) is connected to the rows already in the block
    !
    logical function row_connected(ir)
      integer(ik), intent(in) :: ir
      !
      if (blk_nrow==0) then  ! Any row is connected to an empty block!
        row_connected = .true.
      else
        row_connected = any(abs(mat(ir,blk_cols(:blk_ncol)))>thr)
      end if
    end function row_connected
    !
    !  Add row to the block; we may also need to add a few columns
    !
    subroutine add_row(ir)
      integer(ik), intent(in) :: ir
      integer(ik)             :: ic
      !
      if (.not.row_avail(ir)) stop 'block_matrices%bm_scan_matrix%add_row - logical impossibility trap (1)'
      !
      blk_nrow           = blk_nrow + 1
      blk_rows(blk_nrow) = ir
      row_avail(ir)      = .false.
      !
      add_connected_columns: do ic=1,ncol
        if (abs(mat(ir,ic))<=thr) cycle add_connected_columns
        if (any(blk_cols(:blk_ncol)==ic)) cycle add_connected_columns ! Already in the block
        if (.not.col_avail(ic)) stop 'block_matrices%bm_scan_matrix%add_row - logical impossibility trap (2)'
        !
        blk_ncol           = blk_ncol + 1
        blk_cols(blk_ncol) = ic
        col_avail(ic)      = .false.
      end do add_connected_columns
    end subroutine add_row
    !
  end subroutine bm_scan_matrix

  function bm_find_zero_row(bms) result(izero)
    type(bm_structure), intent(in) :: bms
    integer(ik)                    :: izero
    !
    logical     :: nzero(bms%max_row)
    !
    nzero = .false.
    nzero(bms%rows(:bms%n_rows)) = .true.
    !
    find_zero: do izero=1,size(nzero)
      if (.not.nzero(izero)) return
    end do find_zero
    izero = 0
  end function bm_find_zero_row

  subroutine bm_drop_row(bms,irow,bmo)
    type(bm_structure), intent(in)    :: bms      ! Matrix structure to modify
    integer(ik), intent(in)           :: irow     ! Row to delete
    type(bm_structure), intent(inout) :: bmo      ! Structure of the substituted matrix
    !
    integer(ik) :: ib
    integer(ik) :: r_pos, r_cnt, r_end      ! Position, size, and last index of the block's row list in bms%rows 
    integer(ik) :: c_pos, c_cnt, c_end      ! Position, size, and last index of the block's column list in bms%cols 
    !
    !  A bit of sanity checking first
    !
    if (irow<=0 .or. irow>bms%max_row) stop 'block_matrices%bm_drop_row - bad irow'
    if (bmo%n_blocks/=-1) stop 'block_matrices%bm_drop_row - output descriptor already in use?!'
    !
    !  There are three scenarios for the modified matrix:
    !  a) We may hit an empty row, in which care nothing happens
    !  b) We may shrink an existing block
    !  c) We may completely eliminate a block
    !  The total number of blocks will not increase, and may shrink.
    !
    bmo%n_rows   = 0
    bmo%n_cols   = 0
    bmo%n_blocks = 0
    bmo%max_row  = bms%max_row - 1  ! We drop one, remember?
    bmo%max_col  = bms%max_col
    !
    !  Go through all blocks, shifting rows as needed
    !
    scan_input_blocks: do ib=1,bms%n_blocks
      r_cnt = bms%block_data(1,ib)
      c_cnt = bms%block_data(2,ib)
      r_pos = bms%block_data(3,ib)
      c_pos = bms%block_data(4,ib)
      r_end = r_pos + r_cnt - 1
      c_end = c_pos + c_cnt - 1
      call shift_rows_and_add_block(bmo,irow,bms%rows(r_pos:r_end),bms%cols(c_pos:c_end))
    end do scan_input_blocks
    !
    contains

    subroutine shift_rows_and_add_block(bmo,idrop,row,col)
      type(bm_structure), intent(inout) :: bmo            ! Block descriptor
      integer(ik), intent(in)           :: idrop          ! The "bad" row
      integer(ik), intent(in)           :: row(:), col(:) ! row and column table of the block
      !
      integer(ik) :: drow(size(row))   ! List of rows sans idrop
      integer(ik) :: ir, id            ! Input/output row indices
      integer(ik) :: rowx              ! Adjusted row
      !
      id = 0
      scan_rows: do ir=1,size(row)
        rowx = row(ir)
        if (rowx==idrop) cycle scan_rows
        if (rowx> idrop) rowx = rowx - 1
        id = id + 1
        drow(id) = rowx
      end do scan_rows
      call bm_add_block(bmo,drow(:id),col)
    end subroutine shift_rows_and_add_block
  end subroutine bm_drop_row
  !
  !  The bm_drop_column() routine is a trivial variation of bm_drop_row() above
  !
  subroutine bm_drop_column(bms,icol,bmo)
    type(bm_structure), intent(in)    :: bms      ! Matrix structure to modify
    integer(ik), intent(in)           :: icol     ! Column to delete
    type(bm_structure), intent(inout) :: bmo      ! Structure of the substituted matrix
    !
    integer(ik) :: ib
    integer(ik) :: r_pos, r_cnt, r_end      ! Position, size, and last index of the block's row list in bms%rows 
    integer(ik) :: c_pos, c_cnt, c_end      ! Position, size, and last index of the block's column list in bms%cols 
    !
    !  A bit of sanity checking first
    !
    if (icol<=0 .or. icol>bms%max_col) stop 'block_matrices%bm_drop_column - bad icol'
    if (bmo%n_blocks/=-1) stop 'block_matrices%bm_drop_column - output descriptor already in use?!'
    !
    bmo%n_rows   = 0
    bmo%n_cols   = 0
    bmo%n_blocks = 0
    bmo%max_row  = bms%max_row
    bmo%max_col  = bms%max_col - 1  ! We drop one, remember?
    !
    !  Go through all blocks, shifting columns as needed
    !
    scan_input_blocks: do ib=1,bms%n_blocks
      r_cnt = bms%block_data(1,ib)
      c_cnt = bms%block_data(2,ib)
      r_pos = bms%block_data(3,ib)
      c_pos = bms%block_data(4,ib)
      r_end = r_pos + r_cnt - 1
      c_end = c_pos + c_cnt - 1
      call shift_columns_and_add_block(bmo,icol,bms%rows(r_pos:r_end),bms%cols(c_pos:c_end))
    end do scan_input_blocks
    !
    contains

    subroutine shift_columns_and_add_block(bmo,idrop,row,col)
      type(bm_structure), intent(inout) :: bmo            ! Block descriptor
      integer(ik), intent(in)           :: idrop          ! The "bad" row
      integer(ik), intent(in)           :: row(:), col(:) ! row and column table of the block
      !
      integer(ik) :: dcol(size(col))   ! List of columns sans idrop
      integer(ik) :: ic, id            ! Input/output column indices
      integer(ik) :: colx              ! Adjusted column
      !
      id = 0
      scan_columns: do ic=1,size(col)
        colx = col(ic)
        if (colx==idrop) cycle scan_columns
        if (colx> idrop) colx = colx - 1
        id = id + 1
        dcol(id) = colx
      end do scan_columns
      call bm_add_block(bmo,row,dcol(:id))
    end subroutine shift_columns_and_add_block
  end subroutine bm_drop_column

  subroutine bm_substitute_row(bms,irow,links,bmo)
    type(bm_structure), intent(in)    :: bms      ! Matrix structure to modify
    integer(ik), intent(in)           :: irow     ! Row to substitute
    integer(ik), intent(in)           :: links(:) ! Columns this row connects to, must be sorted in ascensing order
    type(bm_structure), intent(inout) :: bmo      ! Structure of the substituted matrix
    !
    integer(ik) :: ib
    integer(ik) :: merge_row(2*bms%max_row) ! List of rows of merged blocks, possibly with duplicates
    integer(ik) :: merge_col(2*bms%max_col) ! ditto columns
    integer(ik) :: n_mr, n_mc               ! Counters for merge_row/merge_col
    integer(ik) :: r_pos, r_cnt, r_end      ! Position, size, and last index of the block's row list in bms%rows 
    integer(ik) :: c_pos, c_cnt, c_end      ! Position, size, and last index of the block's column list in bms%cols 
    !
    !  A bit of sanity checking first
    !
    if (irow<=0 .or. irow>bms%max_row) stop 'block_matrices%bm_substitute_row - bad irow'
    if (bmo%n_blocks/=-1) stop 'block_matrices%bm_substiture_row - output descriptor already in use?!'
    if (any(links<=0) .or. any(links>bms%max_col)) stop 'block_matrices%bm_substitute_row - bad number of links'
    !
    !  There are two scenarios for the modified matrix:
    !  a) We may create a completely new block
    !  b) We may create new connections between existing blocks
    !  The total number of blocks may at most increase by 1.
    !
    bmo%n_rows   = 0
    bmo%n_cols   = 0
    bmo%n_blocks = 0
    bmo%max_row  = bms%max_row
    bmo%max_col  = bms%max_col
    !
    !  Go through all blocks; those which do not overlap with the replaced row
    !  are copied; those which do get merged.
    !
    n_mr = 0 ; n_mc = 0
    call merge_block((/irow/),links)
    scan_input_blocks: do ib=1,bms%n_blocks
      r_cnt = bms%block_data(1,ib)
      c_cnt = bms%block_data(2,ib)
      r_pos = bms%block_data(3,ib)
      c_pos = bms%block_data(4,ib)
      r_end = r_pos + r_cnt - 1
      c_end = c_pos + c_cnt - 1
      if (bm_intersect((/irow/),bms%rows(r_pos:r_end))) then ! Rows overlap
        if(bm_intersect(links,bms%cols(c_pos:c_end))) then 
          !
          !  Both rows and columns overlap; merge the block into the mashup block
          !
          call merge_block(bms%rows(r_pos:r_end),bms%cols(c_pos:c_end))
        else  
          !
          !  Rows overlap, but columns don't; drop the affected row, and copy the remainder (if any)
          !
          call drop_row_and_add_block(bmo,irow,bms%rows(r_pos:r_end),bms%cols(c_pos:c_end))
        end if
      else ! Rows do not overlap
        if(bm_intersect(links,bms%cols(c_pos:c_end))) then 
          !
          !  Columns overlap, but rows don't; merge the block into the mashup block
          !
          call merge_block(bms%rows(r_pos:r_end),bms%cols(c_pos:c_end))
        else  
          !
          !  Nothing overlaps; Copy the unchanged block to the output
          !
          call bm_add_block(bmo,bms%rows(r_pos:r_end),bms%cols(c_pos:c_end))
        end if
      end if
    end do scan_input_blocks
    call sanitize_merged_block
    call bm_add_block(bmo,merge_row(:n_mr),merge_col(:n_mc))
    !
    contains
    !
    subroutine sanitize_merged_block
      n_mr = sort_unique(merge_row(:n_mr))
      n_mc = sort_unique(merge_col(:n_mc))
    end subroutine sanitize_merged_block 
    !
    subroutine merge_block(row,col)
      integer(ik), intent(in) :: row(:), col(:) ! Block indices to merge
      if (n_mr+size(row)>size(merge_row)) stop 'block_matrices%bm_substitute_row%merge_block - merge_row() overflows'
      if (n_mc+size(col)>size(merge_col)) stop 'block_matrices%bm_substitute_row%merge_block - merge_col() overflows'
      merge_row(n_mr+1:n_mr+size(row)) = row
      merge_col(n_mc+1:n_mc+size(col)) = col
      n_mr = n_mr + size(row)
      n_mc = n_mc + size(col)
    end subroutine merge_block
    !
    subroutine drop_row_and_add_block(bmo,idrop,row,col)
      type(bm_structure), intent(inout) :: bmo            ! Block descriptor
      integer(ik), intent(in)           :: idrop          ! The "bad" row
      integer(ik), intent(in)           :: row(:), col(:) ! row and column table of the block
      !
      integer(ik) :: drow(size(row))   ! List of rows sans idrop
      integer(ik) :: ir, id            ! Input/output row indices
      !
      if (size(row)<=1) return ! Block disappears; nothing else to do
      !
      id = 0
      scan_rows: do ir=1,size(row)
        if (idrop==row(ir)) cycle scan_rows
        id = id + 1
        drow(id) = row(ir)
      end do scan_rows
      if (id/=size(row)-1) stop 'block_matrices%bm_substitute_row%bm_drop_row_and_add_block - drop row is not in the list?!'
      !
      call bm_add_block(bmo,drow(:id),col)
    end subroutine drop_row_and_add_block
    !
  end subroutine bm_substitute_row
  !
  !  WARNING: There is absolutely no checking in bm_add_block. It is possible and easy
  !  WARNING: to create corrupt block descriptors using it!
  !
  subroutine bm_add_block(bmo,row,col)
    type(bm_structure), intent(inout) :: bmo            ! Block descriptor to update
    integer(ik), intent(in)           :: row(:), col(:) ! Block indices to copy to the output structure
    !
    integer(ik) :: nb
    !
    if (size(row)==0 .or. size(col)==0) return ! Adding an empty block is a valid no-op
    if (bmo%n_blocks==-1) stop 'block_matrices%bm_add_block - output descriptor is not initialized?!'
    nb = bmo%n_blocks + 1
    bmo%n_blocks = nb
    if (bmo%n_blocks>max_nrows) stop 'block_matrices%bm_add_block - block_data() - increase max_inrows'
    ! Descriptor
    bmo%block_data(1,nb) = size(row)
    bmo%block_data(2,nb) = size(col)
    bmo%block_data(3,nb) = bmo%n_rows + 1
    bmo%block_data(4,nb) = bmo%n_cols + 1
    ! Overflow check
    if (bmo%n_rows+size(row)>max_nrows) stop 'block_matrices%bm_add_block - rows() - increase max_inrows'
    if (bmo%n_cols+size(col)>max_nrows) stop 'block_matrices%bm_add_block - cols() - increase max_inrows'
    ! Data
    bmo%rows(bmo%n_rows+1:bmo%n_rows+size(row)) = row
    bmo%cols(bmo%n_cols+1:bmo%n_cols+size(col)) = col
    ! Advance counters
    bmo%n_rows = bmo%n_rows + size(row)
    bmo%n_cols = bmo%n_cols + size(col)
  end subroutine bm_add_block
  !
  logical function bm_intersect(c1,c2)
    integer(ik), intent(in) :: c1(:), c2(:)   ! Sorted column (or row) lists to compare
    !
    integer(ik) :: i1, i2, n1, n2
    !
    n1 = size(c1) ; n2 = size(c2)
    i1 = 1 ; i2 = 1
    find_common: do while (i1<=n1 .and. i2<=n2)
      if (c1(i1)==c2(i2)) then
        bm_intersect = .true.
        return
      else if (c1(i1)>c2(i2)) then
        i2 = i2 + 1
      else !   c1(i1)<c2(i2))
        i1 = i1 + 1
      end if
    end do find_common
    bm_intersect = .false.
  end function bm_intersect
  !
  subroutine bm_intersection(list,ord,ref,n_inter,ord_inter)
    integer(ik), intent(in)  :: list(:)      ! List of spatial orbitals, sorted (times -1 for betas)
    integer(ik), intent(in)  :: ord (:)      ! Ordinal of each spatial orbital in list()
    integer(ik), intent(in)  :: ref (:)      ! Reference list of spatial orbitals, sorted
    integer(ik), intent(out) :: n_inter      ! Number of entries list() and ref() have in common
    integer(ik), intent(out) :: ord_inter(:) ! Ordinals of the common elements (no sorting)
    !
    integer(ik) :: nl, nr ! Sizes of list() and ref()
    integer(ik) :: il, ir ! Indices within list() and ref()
    !
    nl = size(list)
    nr = size(ref)
    if (size(ord)/=nl) stop 'block_matrix%bm_intersection - mismatched arrays'
    if (size(ord_inter)<min(nl,nr)) stop 'block_matrix%bm_intersection - output too small'
    !
    n_inter = 0
    il = 1 ; ir = 1
    union_loop: do while (il<=nl .and. ir<=nr)
      if      (list(il)==ref(ir)) then
        n_inter = n_inter + 1
        ord_inter(n_inter) = ord(il)
        il = il + 1
        ir = ir + 1
      else if (list(il)> ref(ir)) then
        ir = ir + 1
      else !  (list(il)< ref(ir))
        il = il + 1
      end if
    end do union_loop
    ! write (out,*) 'list = ', list
    ! write (out,*) 'ord  = ', ord
    ! write (out,*) 'ref  = ', ref
    ! write (out,*) 'plus = ', ord_inter(:n_inter)
    ! write (out,"()")
  end subroutine bm_intersection

end module block_matrices
