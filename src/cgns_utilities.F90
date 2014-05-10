! This file contains all the Fortran-back-end functions for working
! with CGNS files. The intent is all high-level functions are actually
! written in python. This interface just facilities actually reading
! and writing the file. 

subroutine openFile(fileName, mode, cg)
  ! This routine opens a file and returns the handle such that it can
  ! be used in other routines. 

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  character*(*), intent(in) :: fileName
  integer, intent(in) :: mode
  integer, intent(out) :: cg

  ! Working 
  integer :: ier, base

  ! Map our "mode" to the CG_MODE
  if (mode == CG_MODE_READ) then
     call cg_open_f(fileName, mode, cg, ier)
  else if (mode == CG_MODE_WRITE) then
     call cg_open_f(fileName, CG_MODE_WRITE, cg, ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     ! Create the base 
     call cg_base_write_f(cg, "BASE#1", 3, 3, base, ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
  else
     print *,'Error: Mode must be 0 for READ or 1 for WRITE!'
     stop
  end if
  if (ier .eq. CG_ERROR) call cg_error_exit_f

end subroutine openFile

subroutine closeFile(cg)
  ! Close the file-handle cg
  implicit none
  include 'cgnslib_f.h'
  
  ! Input/Output
  integer, intent(in) :: cg

  ! Working
  integer :: ier

  ! Close Output File:
  call cg_close_f(cg, ier)
  if (ier .eq. ERROR) call cg_error_exit_f

end subroutine closeFile

subroutine getNBlocks(cg, N)
  ! Interogate an open cgns file at handle cg and return the
  ! number of block N.

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg
  integer, intent(out) :: N

  ! Working
  integer :: ier, base

  base = 1
  call cg_nzones_f(cg, base, N, ier)
  if (ier .eq. ERROR) call cg_error_exit_f

end subroutine getNBlocks

subroutine getGridDimension(cg, cellDim)
  ! Determine the dimension of the grid. Should be one of 2 or 3

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg
  integer, intent(out) :: cellDim

  ! Working
  integer :: ier, base, physDim
  character(len=512) :: baseName
  base = 1
  call cg_base_read_f(cg, base, baseName, cellDim, physDim, ier)
  if (ier .eq. ERROR) call cg_error_exit_f

end subroutine getGridDimension

subroutine getBlockInfo(cg, iBlock, zoneName, dims, nBoco, nB2B)
  ! Determine the critical meta information of the block, iBlock. We
  ! return dims -- the nodal size of the block, nBocos the number of
  ! boundary conditions and nB2B the number of block-to-block
  ! connection information. 

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg, iBlock
  integer, intent(out) :: dims(3)
  character(len=512), intent(out) :: zoneName
  integer, intent(out) :: nBoco, nB2B
  ! Working
  integer :: ier, base, tmp(9)
  
  base = 1
  call cg_goto_f(cg, base, ier, 'end')
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  ! Get block size
  call cg_zone_read_f(cg, base, iBlock, zonename, tmp, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  dims = tmp(1:3)

  ! Get number of bocos
  call cg_nbocos_f(cg, base, iBlock, nBoco, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  ! Get the number of 1to1 connections
  call cg_n1to1_f(cg, base, iBlock, nB2B, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

end subroutine getBlockInfo

subroutine getBCInfo(cg, iBlock, iBC, bocoName, bocoType, ptRange, family)
  ! Get the BCInfor for 'iBC' condition on block 'iBlock' We determine
  ! the bocoName, the botoType, and pointRange which is sufficient to
  ! reproduce the boundary condition. 

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg, iBlock, iBC
  integer, intent(out) :: bocoType, ptRange(3, 2)
  character(len=256), intent(out) :: bocoName, family

  ! Working
  integer :: ier, base
  integer NormalIndex(3), NormalListSize, NormalDataType, nDataSet
  integer ptset_type, npnts, pnts(3,2),pnts_donor(3,2),ncon
  
  base = 1
  call cg_goto_f(cg, base, ier, 'end')
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_boco_info_f(cg, base, iBlock, iBC, bocoName, bocoType, &
       ptset_type, npnts, NormalIndex, NormalListSize, NormalDataType, nDataSet, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  call cg_boco_read_f(cg, base, iBlock, iBC, ptRange, Integer, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_goto_f(cg, base, ier, "Zone_t", iBlock, "ZoneBC_t",1, "BC_t", iBC, "end")
  if (ier == 0) then ! Node exits
     call cg_famname_read_f(family, ier)
  end if

end subroutine getBCInfo

subroutine getB2BInfo(cg, iBlock, iB2B, connectName, donorName, ptRange, donorRange, transform)
  ! Get the block to block connection information for 'iB2B' condition
  ! on block 'iBlock'. We return sufficient information to recreate
  ! this boundary condition. 

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg, iBlock, iB2B
  integer, intent(out) :: transform(3), ptRange(3, 2), donorRange(3, 2)
  character(len=256), intent(out) :: connectName, donorName

  ! Working
  integer :: ier, base
   
  base = 1
  call cg_goto_f(cg, base, ier, 'end')
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_1to1_read_f(cg, base, iBlock, iB2B, connectName, donorName, ptRange, &
       donorRange, transform, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

end subroutine getB2BInfo

subroutine getCoordinates(cg, iBlock, il, jl, kl, X)
  ! Return the subset of block coordinates of block 'iBlock', of dimensions
  ! (il, jl, kl, 3) in X
  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg, iBlock, il, jl, kl
  real(kind=8), intent(out) :: X(il, jl, kl, 3)
  
  ! Working
  integer :: ier, base, tmp(9), blockStart(3), blockEnd(3)
  
  base = 1
  blockStart = (/1,1,1/)
  blockEnd = (/il, jl, kl/)

  call cg_coord_read_f(cg, base, iBlock, 'CoordinateX', RealDouble,&
       blockStart, blockEnd, X(:, :, :, 1), ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_coord_read_f(cg, base, iBlock, 'CoordinateY', RealDouble,&
       blockStart, blockEnd, X(:, :, :, 2), ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_coord_read_f(cg, base, iBlock, 'CoordinateZ', RealDouble,&
       blockStart, blockEnd, X(:, :, :, 3), ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

end subroutine getCoordinates

subroutine writeZone(cg, zoneName, dims, zoneID)

  ! This function writes a zone to the open file name of size
  ! "dims". Note that the cg file must have been open in write mode. 
  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer,intent(in) :: cg, dims(3)
  character*(*), intent(in) :: zoneName
  integer, intent(out) :: zoneID
  ! Working
  integer :: ier, base, sizes(9)
 
  base = 1
  sizes(:) = 0
  sizes(1) = dims(1)
  sizes(2) = dims(2)
  sizes(3) = dims(3)
  sizes(4) = dims(1)-1
  sizes(5) = dims(2)-1
  sizes(6) = dims(3)-1

  call cg_zone_write_f(cg, base, zoneName, sizes, Structured, zoneID, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

end subroutine writeZone

subroutine writeCoordinates(cg, iBlock, il, jl, kl, X)
  ! Return the subset of block coordinates of block 'iBlock', of dimensions
  ! (il, jl, kl, 3) in X
  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg, iBlock, il, jl, kl
  real(kind=8), intent(in) :: X(il, jl, kl, 3)
  
  ! Working
  integer :: ier, base, cordID
  
  base = 1
  
  call cg_coord_write_f(cg, base, iBlock, realDouble, 'CoordinateX', X(:, :, :, 1), &
       cordID, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_coord_write_f(cg, base, iBlock, realDouble, 'CoordinateY', X(:, :, :, 2), &
       cordID, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_coord_write_f(cg, base, iBlock, realDouble, 'CoordinateZ', X(:, :, :, 3), &
       cordID, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

end subroutine writeCoordinates

subroutine writeBC(cg, iBlock, bcName, bcFam, ptRange, bcType)

  ! This function writes a BC to 'iBlockt' with bcName 'bcName' and
  ! family of 'bcFam'. The pt_start and pt_end range defines the range
  ! of the boundary condition and bcType is the boundary condition
  ! type.

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer,intent(in) :: cg, iBlock, ptRange(3,2), bcType
  character*(*), intent(in) :: bcName
  character*(*), intent(in) :: bcFam

  ! Working
  integer :: ier, base, BCOut
  
  base = 1

  call cg_boco_write_f(cg, base, iBlock, trim(bcName), bcType, PointRange, &
       2, ptRange, BCout, ier)

  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_goto_f(cg, base, ier, 'Zone_t', iBlock, "ZoneBC_t", 1, &
       "BC_t", BCOut, "end")
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_famname_write_f(bcFam, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

end subroutine writeBC

subroutine writeB2B(cg, iBlock, connectName, donorName, pts, ptsDonor, transform)

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg, iBlock, pts(3, 2), ptsDonor(3, 2)
  integer, intent(in) :: transform(3)
  character*(*), intent(in) :: connectName, donorName
  
  ! Working
  integer :: ier, base, nCon
  
  base = 1

  call cg_1to1_write_f(cg, base, iBlock, connectName, donorName, pts, &
       ptsDonor, transform, nCon, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  
end subroutine writeB2B

subroutine refine(Xin, Xout, il, jl, kl)
  ! This routine is used to refine a block 'Xin' of size (il, jl, kl,
  ! 3) to form block Xout of size ((il-1)*2+1, (jl-1)*2+1, (kl-1)*2+1,
  ! 3). It uses a local cubic interpolation that attemps to respect
  ! spacing ratios of the original grid. Note that this routine does
  ! not use anything from the cgns library 

  implicit none

  ! Input/Output
  real(kind=8), intent(in), dimension(il, jl, kl, 3) :: Xin
  integer, intent(in) :: il, jl ,kl
  real(kind=8), intent(out), dimension((il-1)*2+1, (jl-1)*2+1, (kl-1)*2+1, 3) :: Xout

  ! Working 
  integer :: i, j, k, ii, jj, kk, idim, ill, jll, kll
  real(kind=8) :: fact26, pt(3) 

  ill = (il-1)*2 + 1
  jll = (jl-1)*2 + 1
  kll = (kl-1)*2 + 1
  fact26 = dble(1.0)/dble(26.0)
  ! First copy the values we have into the new array
  do idim=1,3
     do k=1,kl
        do j=1,jl
           do i=1,il
              Xout((i-1)*2+1, (j-1)*2+1, (k-1)*2+1, idim) = Xin(i, j, k, idim)
           end do
        end do
     end do
  end do

  ! Next interpolate all the edges
  do k=1,kl
     do j=1,jl
        call interpEdge(Xin(:, j, k, :), Xout(:, (j-1)*2+1, (k-1)*2+1, :), il)
     end do
  end do

  do k=1,kl
     do i=1,il
        call interpEdge(Xin(i, :, k, :), Xout((i-1)*2+1, :, (k-1)*2+1, :), jl)
     end do
  end do

  do j=1,jl
     do i=1,il
        call interpEdge(Xin(i, j, :, :), Xout((i-1)*2+1, (j-1)*2+1, :, :), kl)
     end do
  end do

  ! Next interpolate all the faces
  do k=1,kl
     call interpFace(Xin(:, :, k, :), Xout(:, :, (k-1)*2+1, :), il, jl)
  end do

  do j=1,jl

     call interpFace(Xin(:, j, :, :), Xout(:, (j-1)*2+1, :, :), il, kl)
  end do

  do i=1,il
     call interpFace(Xin(i, :, :, :), Xout((i-1)*2+1, :, :, :), jl, kl)
  end do

  ! Finally interpolate the last center points --- by averaging all 26 neighbours
  do k=2, kll, 2
     do j=2, jll, 2
        do i=2, ill, 2
           
           pt = 0.0_8
           ! Get 18 of them from these loops
           do kk=-1,1,2
              do jj=-1,1
                 do ii=-1,1
                    pt = pt + Xout(i+ii, j+jj, k+kk, :)
                 end do
              end do
           end do

           ! And another 6 on the kk=0 sub-plane
           do ii=-1,1
              do jj=-1,1,2
                 pt = pt + Xout(i+ii, j+jj, k, :)
              end do
           end do

           ! And the last 2:
           pt = pt + Xout(i-1, j, k, :)
           pt = pt + Xout(i+1, j, k, :)

           ! Finally set value:
           Xout(i, j, k, :) = pt * fact26
        end do
     end do
  end do

end subroutine refine

subroutine interpEdge(Xcoarse, Xfine, il)

  ! This routine generically interpolates an edge of length il from
  ! Xcoarse to Xfine. 

  implicit none

  ! Input/Output
  real(kind=8), intent(in), dimension(il, 3) :: Xcoarse
  integer, intent(in) :: il
  real(kind=8), intent(inout), dimension((il-1)*2+1, 3) ::Xfine

  ! Working 
  integer :: i, ii

  ! We we really want to do this higher order, but do it linear now:
  do i=1,il-1
     ii = (i-1)*2 + 2
     Xfine(ii, :) = 0.5_8 * (Xcoarse(i, :) + Xcoarse(i+1, :))
  end do

end subroutine interpEdge

subroutine interpFace(Xcoarse, Xfine, il, jl)

  ! This routine generically interpolates an face of size il, jl from
  ! Xcoarse to Xfine. We assume we already have computed where the
  ! edges of the face are: XCoarse is actually unused. 

  implicit none

  ! Input/Output
  real(kind=8), intent(in), dimension(il, jl, 3) :: Xcoarse
  integer, intent(in) :: il, jl
  real(kind=8), intent(inout), dimension((il-1)*2+1, (jl-1)*2+1, 3) ::Xfine

  ! Working 
  integer :: i, j, ill, jll

  ! Average the 8 nodes we already have. 
  ill = (il-1)*2 + 1
  jll = (jl-1)*2 + 1

  do j=2, jll, 2
     do i=2, ill, 2
        Xfine(i, j, :) = &
             0.125_8 * (&
             Xfine(i-1, j-1, :) + &
             Xfine(i-1, j  , :) + &
             Xfine(i-1, j+1, :) + &
             Xfine(i  , j-1, :) + &
             Xfine(i  , j+1, :) + &
             Xfine(i+1, j-1, :) + &
             Xfine(i+1, j  , :) + &
             Xfine(i+1, j+1, :))
     end do
  end do
end subroutine interpFace
