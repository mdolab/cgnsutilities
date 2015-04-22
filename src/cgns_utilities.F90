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

subroutine writeBC(cg, iBlock, bcName, bcFam, ptRange, bcType, bcOut)

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
  integer, intent(out) :: bcOut
  ! Working
  integer :: ier, base
  
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

subroutine writeBCData(cg, iBlock, bcType, bcIn, dataName, dataValue, writeHeader)

  ! This function writes actual BCData. The writeBCDataHeader must
  ! have already been called. 
  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer,intent(in) :: cg, iBlock, bcIn, bcType 
  character*(*),  intent(in) :: dataName
  real(kind=8), intent(in) :: dataValue
  logical :: writeHeader
  ! Working
  integer :: ier, base, i, iDataSet
  
  base = 1
  iDataSet = 1

  if (writeHeader) then 
     call cg_dataset_write_f(cg, base, iBlock, bcIn, "data", BCType, idataset, ier )
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     
     call cg_bcdata_write_f(cg, base, iBlock, bcIn, iDataSet, Dirichlet, ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
  end if
  
  call cg_goto_f(cg, base, ier, 'Zone_t', iBlock, "ZoneBC_t", 1, &
       "BC_t", BCIn, "BCDataSet_t", idataSet, "BCData_t", Dirichlet, &
       "end")
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_array_write_f(trim(dataName), RealDouble, 1, (/1/), dataValue, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

end subroutine writeBCData


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

module zone_vars

  ! This module stores data associated with a single zone and with n_time instances
  implicit none
  save 

  type zoneData

     integer :: nx,ny,nz
     integer :: il,jl,kl
     integer :: n_time, n_vars

     ! Grid has shape: nx by ny by nz ny 3 by n_time
  end type zoneData
end module zone_vars

subroutine time_combine(fileNames, nFiles, outputFile)

  use zone_vars
  implicit none
  include 'cgnslib_f.h'

  ! Input parameters
  integer :: nFiles, lengthMax
  character*256, dimension(nFiles) :: fileNames
  character*256 :: outputFile

  ! Working variables
  integer, dimension(3*3) :: isize
  integer, dimension(6) :: rinde

  ! File handles for cgns files: cg_current is the current step,
  ! cg_unsteady is the new unsteady file
  integer :: cg_current, cg_unsteady

  ! CGNS counters for base
  integer current_base,unsteady_base

  ! CGNS Names for base, zone and iterative data
  character*32 basename,  baseitername
  character*32 fieldname
  character*36 zonename
  ! CGNS Data type
  integer datatype

  ! CGNS Name for unsteady grid
  character*32 gridName,solName
  character*32, dimension(:),allocatable :: motionNames,solNames
  character*32, dimension(:),allocatable :: gridNames
  character*32, dimension(:),allocatable :: zoneNames
  character*32, dimension(:),allocatable :: fieldNames
  character*32, dimension(:),allocatable :: var_names

  ! CGNS cell dimension and physical dimension
  integer  CellDim, PhysDim

  ! CGNS number of bases and zones in current file
  integer nbases, nzones

  ! CGNS error flag
  integer ier

  ! Integer Counters etc
  integer :: N,i_start,i_end,n_steps,n_time,i,j,nn
  integer :: ii,jj,kk,ind(3)
  integer :: counter,nsols,sol,nfields,field,izone
  logical, dimension(5) :: EV_Exist
  integer :: temp_shape(6),rmax(3),farStart

  ! Time step
  real(kind=8) :: dt

  ! Character arrays for names and I/O
  character*100 :: char_temp
  character*128 base_name,unsteady_filename
  character*256 current_filename
  character*32, dimension(3) :: coordNames
  ! Array of all zone data
  type(zoneData), dimension(:), allocatable :: zone

  ! Misc Data Arrays
  real(kind=4), allocatable,dimension(:,:,:  ) :: data3d
  real(kind=8), allocatable,dimension(:,:,:,:) :: gridTemp
  real(kind=8), allocatable,dimension(:,:,:  ) :: fieldData

  integer         , allocatable,dimension(:)       :: int1d
  integer :: dummy_int

  ! Tecplot Data Stuff
  ! Working 
  ! These Never Change
  integer VIsDouble                /1/
  integer debug                    /0/
  integer FileType                 /0/ 
  INTEGER ZoneType                 /0/
  INTEGER TECINI112,tecdat112,teczne112,tecend112

  INTEGER ICellMax                 / 0/
  INTEGER JCellMax                 / 0/
  INTEGER KCellMax                 / 0/
  INTEGER DIsDouble                / 1/
  INTEGER StrandID                 / 0/      !/* StaticZone */
  INTEGER ParentZn                 / 0/
  INTEGER IsBlock                  / 1/      !/* Block */
  INTEGER NFConns                  / 0/
  INTEGER FNMode                   / 0/
  INTEGER TotalNumFaceNodes        / 0/
  INTEGER TotalNumBndryFaces       / 0/
  INTEGER TotalNumBndryConnections / 0/
  INTEGER ShrConn                  / 0/
  INTEGER,dimension(:),allocatable :: ValueLocation,passiveVarList,ShareVarFromZone

  ! Start of Code
  coordNames(1) = "CoordinateX"
  coordNames(2) = "CoordinateY"
  coordNames(3) = "CoordinateZ"

  call cg_open_f(trim(fileNames(1)), MODE_READ, cg_current, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_nbases_f(cg_current, nbases, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  if (nbases .lt. 1) then
     print *, 'Error: No bases found in first CGNS file'
     stop
  end if

  ! We will only deal with 1 base in each zone so ALWAYS only use 1 base
  current_base = 1
  unsteady_base = 1

  ! Get the cell and phys dim 
  call cg_base_read_f(cg_current, current_base, basename, &
       CellDim, PhysDim, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_nzones_f(cg_current, current_base, nzones, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  allocate(zoneNames(nzones))
  ! ------------------------------------------------------------------
  ! Allocate the master data array:

  allocate(zone(nzones))
  ! ------------------------------------------------------------------
  ! Loop over the zones to get the required size info:
  do i=1,nzones
     call cg_zone_read_f(cg_current, current_base, i,&
          zoneNames(i), isize,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     call cg_nsols_f(cg_current, current_base, i, nsols , ier )
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     if (nsols >= 1) then
        call cg_nfields_f(cg_current, current_base, i, 1, nfields, ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f
     else
        nfields = 0
     end if
     
     if (cellDim == 2) then
        zone(i)%nx = isize(1)
        zone(i)%ny = isize(2)
        zone(i)%nz = 1
        zone(i)%il = isize(3)
        zone(i)%jl = isize(4)
        zone(i)%kl = 1
     else
        zone(i)%nx = isize(1)
        zone(i)%ny = isize(2)
        zone(i)%nz = isize(3)
        zone(i)%il = isize(4)
        zone(i)%jl = isize(5)
        zone(i)%kl = isize(6)
     end if
  end do

  ! Allocate some names
  allocate(fieldNames(nfields))
  allocate(var_names(nFields+4))
  allocate(valueLocation(nFields+physDim),&
       passiveVarList(nFields+physDim),&
       ShareVarFromZone(nFields+physDim))

  ! Loop over to get field names
  do j=1,nFields
     call cg_field_info_f(cg_current,current_base,1,1,j,&
          datatype , fieldname , ier )

     fieldNames(j) = fieldName
  end do

  call cg_close_f(cg_current,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  ! ------------------------------------------------------------------
  ! Load in all the data:
  print *,'Reading Data from File: '

  ! No Passive Var List or Share from Zone
  passiveVarList(:) = 0
  ShareVarFromZone(:) = 0

  ! All Data is at Nodes
  valueLocation(:) = 1

  ! Accumulate Names
  do j=1,physDim
     var_names(j) = coordNames(j)
  end do

  do j=1,nFields
    var_names(j+physDim) = fieldNames(j)
  end do
  var_names(nFields+3+1) = char(0)

  ! Open Ouput Tecplot File
  ier = TECINI112("Unsteady Data"//char(0),var_names,trim(outputFile)//char(0),"."//char(0),&
       FileType,Debug,VIsDouble)

  nn = 0
  dt = 1.0
  do n=1,nFiles
     print *,'Processing file:', n
     ! Open Current File
     call cg_open_f(trim(fileNames(n)) ,MODE_READ, cg_current, ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     ! Goto Base Node in Current File
     call cg_goto_f(cg_current, current_base, ier, 'end')
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     ! Loop over the zones
     do i=1,nzones

        ! Write Tecplot Zone Data
        ier = TECZNE112(trim(zoneNames(i))//char(0),ZoneType,zone(i)%nx,zone(i)%ny,zone(i)%nz,&
             ICellMax,JCellMax,KCellMax, &
             dt*dble(n-1),i,ParentZn,IsBlock,NFConns,FNMode,&
             TotalNumFaceNodes,TotalNumBndryFaces,TotalNumBndryConnections,&
             PassiveVarList, valueLocation, ShareVarFromZone,ShrConn)
        
        !-------------------------------------------------------------
        !                         Grid Coodinates 
        !-------------------------------------------------------------
           
        ! Load the grid coordinates
        call cg_zone_read_f(cg_current,current_base,i,zoneNames(i),isize,ier)
           
        rmax(1) = zone(i)%nx
        rmax(2) = zone(i)%ny
        rmax(3) = zone(i)%nz

        allocate(gridTemp(rmax(1),rmax(2),rmax(3),physDim))

        do j=1,physDim
           call cg_coord_read_f(cg_current,current_base,i,&
                coordNames(j),RealDouble,(/1,1,1/),&
                rmax,gridTemp(:,:,:,j),ier)
           if (ier .eq. CG_ERROR) call cg_error_exit_f
        end do

        ! Write Grid Coordinates
        do j=1,physDim
           ier = TECDAT112(zone(i)%nx*zone(i)%ny*zone(i)%nz,gridTemp(:,:,:,j),&
                DIsDouble)
        end do

        deallocate(gridTemp)

        !-------------------------------------------------------------
        !                     Field Values
        !-------------------------------------------------------------

        if (nfields > 0) then
           ! Read the Rind Data
           call cg_goto_f(cg_current,current_base,ier,'Zone_t',i, &
                'FlowSolution_t',1,'end')
           rinde(:) = 0
           call cg_rind_read_f(rinde , ier )
           if (ier .eq. CG_ERROR) call cg_error_exit_f


           ! Load the Variables
           do j=1,nFields

              ! Allocate node based field data
              allocate(fieldData(zone(i)%nx,zone(i)%ny,zone(i)%nz))

              call cg_field_info_f(cg_current,current_base,i,1,j,&
                   datatype , fieldname , ier )

              if (ier .eq. CG_ERROR) call cg_error_exit_f
              temp_shape(:) = 0
              temp_shape = (/1 - rinde(1), zone(i)%il + rinde(2), &
                   1 - rinde(3), zone(i)%jl + rinde(4), &
                   1 - rinde(5), zone(i)%kl + rinde(6)/)

              ! Allocate Temporary Array for Cell Data, possibly with rind cells
              allocate(data3d(temp_shape(1):temp_shape(2),&
                   temp_shape(3):temp_shape(4),&
                   temp_shape(5):temp_shape(6)))

              rmax(1) = zone(i)%il+rinde(1)+rinde(2)
              rmax(2) = zone(i)%jl+rinde(3)+rinde(4)
              rmax(3) = zone(i)%kl+rinde(5)+rinde(6)

              call cg_field_read_f(cg_current,current_base,i,1, &
                   fieldName, datatype, (/1,1,1/),&
                   rmax,data3d,ier)
              if (ier .eq. CG_ERROR) call cg_error_exit_f

              ! Call the interpolate function to reconstruct node data
              call interpolate(data3d,temp_shape,&
                   fieldData,zone(i)%nx,zone(i)%ny,zone(i)%nz)
              
              ! Deallocate Temporary Data Array
              deallocate(data3d)
              
              ! Write out Field Values
              ier   = TECDAT112(zone(i)%nx*zone(i)%ny*zone(i)%nz,&
                   fieldData, DIsDouble)
              deallocate(fieldData)
           end do ! Field Loop
        end if ! If we have fields
     end do ! Zone Loop
     ! Close the current CGNS file
     call cg_close_f(cg_current,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
  end do

  ! Close Tecplot File
  ier = TECEND112()
  ! Deallocate data and quit
  deallocate(zone)  
  deallocate(fieldNames,var_Names,valueLocation)
  deallocate(passiveVarList,ShareVarFromZone)
  deallocate(zoneNames)
end subroutine time_combine

subroutine interpolate(cell_data,cell_shape,node_data,inode,jnode,knode)

  ! Take cell centered data which possibly contain rind cells and
  ! perform a linear reconstruction back to the nodes. cell_shape
  ! contains the start and end values for each index in
  ! cell_data. Retruned data comes back in node_data. It may not be
  ! the most efficient possible, but it works.

  implicit none

  integer, intent(in) :: cell_shape(6),inode,jnode,knode

  real(kind=4), intent(in ), dimension(&
       cell_shape(1):cell_shape(2),&
       cell_shape(3):cell_shape(4),&
       cell_shape(5):cell_shape(6)) :: cell_data
  real(kind=8), intent(out), dimension(inode,jnode,knode) :: node_data
  integer :: ic_s,ic_e,jc_s,jc_e,kc_s,kc_e
  integer, dimension(inode,jnode,knode) :: node_cnts
  integer :: i,j,k

  ic_s = cell_shape(1)
  ic_e = cell_shape(2)
  jc_s = cell_shape(3)
  jc_e = cell_shape(4)
  kc_s = cell_shape(5)
  kc_e = cell_shape(6)

  node_data(:,:,:) = 0.0
  node_cnts(:,:,:) = 0

  do k=kc_s,kc_e
     do j=jc_s,jc_e
        do i=ic_s,ic_e

           if (i >= 1 .and. i <= inode .and. j >= 1 .and. j <= jnode .and. k >=1 .and. k <= knode) then
              node_data(i  ,j  ,k  ) = node_data(i  ,j  ,k  ) + cell_data(i,j,k)
              node_cnts(i  ,j  ,k  ) = node_cnts(i  ,j  ,k  ) + 1
           end if

           if (i+1 >= 1 .and. i+1 <= inode .and. j >= 1 .and. j <= jnode .and. k >=1 .and. k <= knode) then
              node_data(i+1,j  ,k  ) = node_data(i+1,j  ,k  ) + cell_data(i,j,k)
              node_cnts(i+1,j  ,k  ) = node_cnts(i+1,j  ,k  ) + 1
           end if

           if (i >= 1 .and. i <= inode .and. j+1 >= 1 .and. j+1 <= jnode .and. k >=1 .and. k <= knode) then
              node_data(i  ,j+1,k  ) = node_data(i  ,j+1,k  ) + cell_data(i,j,k)
              node_cnts(i  ,j+1,k  ) = node_cnts(i  ,j+1,k  ) + 1
           end if

           if (i+1 >= 1 .and. i+1 <= inode .and. j+1 >= 1 .and. j+1 <= jnode .and. k >=1 .and. k <= knode) then
              node_data(i+1,j+1,k  ) = node_data(i+1,j+1,k  ) + cell_data(i,j,k)
              node_cnts(i+1,j+1,k  ) = node_cnts(i+1,j+1,k  ) + 1
           end if

           if (i >= 1 .and. i <= inode .and. j >= 1 .and. j <= jnode .and. k+1 >=1 .and. k+1 <= knode) then
              node_data(i  ,j  ,k+1) = node_data(i  ,j  ,k+1) + cell_data(i,j,k)
              node_cnts(i  ,j  ,k+1) = node_cnts(i  ,j  ,k+1) + 1
           end if

           if (i+1 >= 1 .and. i+1 <= inode .and. j >= 1 .and. j <= jnode .and. k+1>=1 .and. k+1 <= knode) then
              node_data(i+1,j  ,k+1) = node_data(i+1,j  ,k+1) + cell_data(i,j,k)
              node_cnts(i+1,j  ,k+1) = node_cnts(i+1,j  ,k+1) + 1
           end if

           if (i >= 1 .and. i <= inode .and. j+1 >= 1 .and. j+1 <= jnode .and. k+1>=1 .and. k+1<= knode) then
              node_data(i  ,j+1,k+1) = node_data(i  ,j+1,k+1) + cell_data(i,j,k)
              node_cnts(i  ,j+1,k+1) = node_cnts(i  ,j+1,k+1) + 1
           end if

           if (i+1 >= 1 .and. i+1 <= inode .and. j+1 >= 1 .and. j+1 <= jnode .and. k+1 >=1 .and. k+1 <= knode) then
              node_data(i+1,j+1,k+1) = node_data(i+1,j+1,k+1) + cell_data(i,j,k)
              node_cnts(i+1,j+1,k+1) = node_cnts(i+1,j+1,k+1) + 1
           end if
        end do
     end do
  end do

  ! Divide through by the number of times a value was added to a node

  do k=1,knode
     do j=1,jnode
        do i=1,inode
           node_data(i,j,k) = node_data(i,j,k) / dble(node_cnts(i,j,k))
        end do
     end do
  end do

end subroutine interpolate
