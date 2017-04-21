
! This file contains all the Fortran-back-end functions for working
! with CGNS files. The intent is all high-level functions are actually
! written in python. This interface just facilities actually reading
! and writing the file. 

subroutine openFile(fileName, mode, cellDim, cg)
  ! This routine opens a file and returns the handle such that it can
  ! be used in other routines. 
  !
  ! The available modes are:
  ! CG_MODE_READ = 0
  ! CG_MODE_WRITE = 1

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  character*(*), intent(in) :: fileName
  integer, intent(in) :: mode, cellDim
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
     call cg_base_write_f(cg, "BASE#1", cellDim, 3, base, ier)
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

subroutine getBCInfo(cg, iBlock, iBC, cellDim, bocoName, bocoType, ptRange, family, nDataSet)
  ! Get the BCInfor for 'iBC' condition on block 'iBlock' We determine
  ! the bocoName, the botoType, and pointRange which is sufficient to
  ! reproduce the boundary condition. 

  ! Note that ptRange can be a 2D i.e. (2,2) but we force it 3D (3,2)
  ! tmpPtRange contains the actual size 

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg, iBlock, iBC, cellDim
  integer, intent(out) :: bocoType, ptRange(3, 2), nDataSet
  character(len=256), intent(out) :: bocoName, family

  ! Working
  integer :: ier, base, i
  integer :: NormalIndex(3), NormalListSize, NormalDataType
  integer :: ptset_type, npnts, tmpPtRange(cellDim,2),pnts_donor(cellDim,2),ncon
  integer :: nuserdata
  character(len=256) :: name

  base = 1
  call cg_goto_f(cg, base, ier, 'end')
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_boco_info_f(cg, base, iBlock, iBC, bocoName, bocoType, &
       ptset_type, npnts, NormalIndex, NormalListSize, NormalDataType, nDataSet, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  call cg_boco_read_f(cg, base, iBlock, iBC, tmpPtRange, Integer, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  if (cellDim == 2) then
     ptRange(1:2,1:2) = tmpPtRange(1:2,1:2)
  else
     ptRange = tmpPtRange
  end if 
  
  call cg_goto_f(cg, base, ier, "Zone_t", iBlock, "ZoneBC_t",1, "BC_t", iBC, "end")
  if (ier == 0) then ! Node exits

     ! Get family name
     call cg_famname_read_f(family, ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     ! Check for overset BC
     ! First get number of user defined data
     call cg_nuser_data_f(nuserdata, ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     ! Loop over all user defined data to check for overset bc
     do i = 1,nuserdata

        ! Read user defined data
        call cg_user_data_read_f(i, name, ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f
        
     end do
  end if
  
end subroutine getBCInfo

subroutine getBCDataSetInfo(cg, iBlock, iBC, iBCDataSet, bocoDatasetName, &
        bocoType, nDirichletArrays, nNeumannArrays)
  
  ! This subroutine returns information about any BC datasets 
  ! that might be associated to this BC and the number of each type.

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg, iBlock, iBC, iBCDataSet
  integer, intent(out) :: nDirichletArrays, nNeumannArrays, bocoType
  character(len=256), intent(out) :: bocoDatasetName

  ! Working
  integer :: ier, base
  integer :: DirichletFlag, NeumannFlag

  base = 1
  call cg_dataset_read_f(cg, base, iBlock, iBC, iBCDataSet, &
       bocoDatasetName, bocoType, DirichletFlag, NeumannFlag, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  
  if (DirichletFlag == 1) then
    call cg_goto_f(cg, base, ier, 'Zone_t', iBlock, 'ZoneBC_t', 1, 'BC_t', &
          iBC, 'BCDataSet_t', iBCDataSet, 'BCData_t', Dirichlet, 'end')
    if (ier .eq. CG_ERROR) call cg_error_exit_f

    ! Read number of variables specified here.
    call cg_narrays_f(nDirichletArrays, ier)
  end if
  
  if (NeumannFlag == 1) then
    call cg_goto_f(cg, base, ier, 'Zone_t', iBlock, 'ZoneBC_t', 1, 'BC_t', &
          iBC, 'BCDataSet_t', iBCDataSet, 'BCData_t', Neumann, 'end')
    if (ier .eq. CG_ERROR) call cg_error_exit_f

    ! Read number of variables specified here.
    call cg_narrays_f(nNeumannArrays, ier)
  end if
  
end subroutine getBCDataSetInfo


subroutine getBCDataArrayInfo(cg, iBlock, iBC, iBCDataSet, iDataArr, &
        flagDirNeu, dataArrayName, dataType, nDataDimensions, &
        dataDimensionVector)

  ! This subroutine returns information about data arrays in a particular BC dataset.

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg, iBlock, iBC, iBCDataSet, iDataArr, flagDirNeu
  integer, intent(out) :: dataType, nDataDimensions
  integer, intent(out), dimension(3) :: dataDimensionVector
  character(len=256), intent(out) :: dataArrayName

  ! Working
  integer :: ier, base, i
  
  ! Initialize from zeros to ones to make sure nothing bad happens
  ! when we get the total number of elements
  dataDimensionVector = 1

  base = 1  
  call cg_goto_f(cg, base, ier, 'Zone_t', iBlock, 'ZoneBC_t', 1, 'BC_t', &
        iBC, 'BCDataSet_t', iBCDataSet, 'BCData_t', flagDirNeu, 'end')
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  
  call cg_array_info_f(iDataArr, &
        dataArrayName, dataType, nDataDimensions, dataDimensionVector, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  
end subroutine getBCDataArrayInfo

subroutine getBCDataArray(cg, iBlock, iBC, iBCDataSet, iDataArr, flagDirNeu, dataArr, nDataArr)
  ! This subroutine retrieves and returns BC dataset array
  
  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg, iBlock, iBC, iBCDataSet, iDataArr, flagDirNeu, nDataArr  
  real(kind=8), intent(inout), dimension(nDataArr) :: dataArr

  ! Working
  integer :: ier, base
  
  ! Call the goto to make sure we are at the right location in the tree
  base = 1  
  call cg_goto_f(cg, base, ier, 'Zone_t', iBlock, 'ZoneBC_t', 1, 'BC_t', &
        iBC, 'BCDataSet_t', iBCDataSet, 'BCData_t', flagDirNeu, 'end')
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  ! Note we read all data as real double even though it is integer or single.
  call cg_array_read_as_f(iDataArr, RealDouble, dataArr, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

end subroutine getBCDataArray

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

subroutine getConvInfo(cg, niterations, narrays)
  ! Get the convegence information stored under the "GlobalConvergenceHistory" node.
  ! This function only get the size of the data.
  ! Ney Secco, 2016

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg
  integer, intent(out) :: niterations, narrays

  ! Working
  character(len=256) :: NormDefinitions
  integer :: ier, base

  ! Navigate to base 1, which usually contains convergence data
  base = 1
  call cg_goto_f(cg, base, ier, 'end')
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  ! Get the number of iterations in the convergence history
  call cg_convergence_read_f(niterations, NormDefinitions, ier)

  ! If no history was found, set the number of iterations to zero and exit
  if (ier .ne. CG_OK) then

     niterations = 0
     narrays = 0

  else ! Continue reading

     print *,'Found convergence history...'

     ! Navigate to the convergence data node
     call cg_goto_f(cg, base, ier, 'ConvergenceHistory_t', 1, 'end')
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     ! Get the number of entries (arrays) in the convergence history
     call cg_narrays_f(narrays, ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

  end if

  ! The we can call the getConvArray subroutine from python to get data for each array

end subroutine getConvInfo

subroutine getConvArray(cg, niterations, arrayID, arrayName, arrayData)
  ! Get the convegence data from a the arrayID data node.
  ! arrayID folows 1-based ordering (unlike Python).
  ! This subroutine only reads 1D arrays with real numbers.
  ! Ney Secco, 2016

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg, niterations, arrayID
  character(len=256), intent(out) :: arrayName
  real(kind=8), intent(out) :: arrayData(niterations)

  ! Working
  integer :: ier, base, DataType, DataDimension
  integer :: DimensionVector(12) ! 12 is the max size allowed to DimensionVector according to the manual

  ! We will only read info from base 1
  base = 1

  ! Navigate to the convergence data node
  call cg_goto_f(cg, base, ier, 'ConvergenceHistory_t', 1, 'end')
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  ! Initialize DimensionVector
  DImensionVector = 0

  ! Get the data array info
  call cg_array_info_f(arrayID, arrayName, DataType, DataDimension, DimensionVector, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  ! Print name of the array
  print *,'reading array: ',arrayName(:33)

  ! Initialize array data
  arrayData = 0

  ! Check if we have 1D array with no chars
  if (DimensionVector(2) .gt. 0) then
     print *,'   found 2D array in convergence history. Unable to read'

  else if (DataType .eq. Character) then
     print *,'   found char in convergence history. Unable to read'

  else ! We can proceed

     ! Read Data
     call cg_array_read_as_f(arrayID, RealDouble, arrayData, ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

  end if

end subroutine getConvArray

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

  if (dims(3) .ne. 1) then

     sizes(1) = dims(1)
     sizes(2) = dims(2)
     sizes(3) = dims(3)
     sizes(4) = dims(1)-1
     sizes(5) = dims(2)-1
     sizes(6) = dims(3)-1
     
  else

     sizes(1) = dims(1)
     sizes(2) = dims(2)
     sizes(3) = dims(1)-1
     sizes(4) = dims(2)-1

  end if

  call cg_zone_write_f(cg, base, zoneName, sizes, Structured, zoneID, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

end subroutine writeZone

subroutine writeCoordinates(cg, iBlock, il, jl, kl, X)
  ! Writes the subset of block coordinates 'X' of dimensions
  ! (il, jl, kl, 3) into the block 'iBlock' of  the file 'cg',
  ! (which should be open in write mode)
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

  !Add an user-defined data node for the overset BC case
  if (bcType == 1) call cg_user_data_write_f('BCOverset', ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

end subroutine writeBC



subroutine writeBCDataHeader(cg, iBlock, bcType, iBC, datasetName, iDataSet)

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer,intent(in) :: cg, iBlock, bcType, iBC
  character*(*),  intent(in) :: datasetName
  integer,intent(out) :: iDataSet
  
  ! Working
  integer :: ier, base

  base = 1
  
  call cg_dataset_write_f(cg, base, iBlock, iBC, datasetName, bcType, iDataSet, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f


end subroutine writeBCDataHeader


subroutine writeBCData(cg, iBlock, iBC, iDataSet, flagDirNeu, writeBCDataHeader,&
           dataArrayName, dataType, nDataDimensions, dataDimensionVector, &
           dataArr, nDataArr)

  ! This function writes actual BCData. The writeBCDataHeader must
  ! have already been called. 
  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer, intent(in) :: cg, iBlock, iBC, iDataSet, flagDirNeu
  integer, intent(in) :: dataType, nDataDimensions, nDataArr
  logical, intent(in) :: writeBCDataHeader
  character*(*),  intent(in) :: dataArrayName
  integer, intent(in), dimension(3) :: dataDimensionVector
  real(kind=8), intent(in), dimension(nDataArr) :: dataArr
  
  ! Working
  integer :: ier, base

  base = 1

  if (writeBCDataHeader) then 
     call cg_bcdata_write_f(cg, base, iBlock, iBC, iDataSet, flagDirNeu, ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
  end if

  call cg_goto_f(cg, base, ier, 'Zone_t', iBlock, "ZoneBC_t", 1, &
       "BC_t", iBC, "BCDataSet_t", idataSet, "BCData_t", flagDirNeu, "end")
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  ! Note that dataType is ignored here since the data is stored 
  ! as RealDouble in python. Everything is thus written as double.
  call cg_array_write_f(trim(dataArrayName), RealDouble, nDataDimensions, &
        dataDimensionVector, dataArr, ier)
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

module dataTypes

  ! This module stores data associated with a single zone and with n_time instances
  implicit none
  save 

  integer, parameter :: iMin=1
  integer, parameter :: iMax=2
  integer, parameter :: jMin=3
  integer, parameter :: jMax=4
  integer, parameter :: kMin=5
  integer, parameter :: kMax=6
  integer, parameter :: BC=0
  integer, parameter :: B2B=1

  type zoneData

     integer :: nx,ny,nz
     integer :: il,jl,kl
     integer :: n_time, n_vars

     ! Grid has shape: nx by ny by nz ny 3 by n_time
  end type zoneData

  type patch
     ! This derivied data type is used to store either boundary
     ! condition or block to block information for a given mesh. Since
     ! we don't know how many patches will be on a given multiblock
     ! mesh before we start, we use a (singly) linked list to
     ! incrementally add patches as we progress through the grid blocks. 

     ! One of BC or B2B
     integer :: type

     ! Data common to BC and B2B
     integer :: pointRange(3,2)
     integer :: myID
     real(kind=8), dimension(3) :: faceAvg
     real(kind=8), dimension(3) :: faceNormal

     ! Data for B2B
     integer :: pointRangeDonor(3, 2)
     integer :: transform(3)
     integer :: donorID

     ! Pointer for next BC
     type(patch), pointer :: next

  end type patch

  type(patch), pointer :: patches
  integer :: nPatches

end module dataTypes

subroutine time_combine(fileNames, nFiles, outputFile)

  use dataTypes
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
  INTEGER StrandID                 / 1/      !/* StaticZone */
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
             dt*(n-1),i,ParentZn,IsBlock,NFConns,FNMode,&
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

subroutine computeConnectivity(inCoords, nCoords, sizes, nBlock, tol)

  ! This routine will compute all bock to block connections as well as
  ! identify the subfaces with boundary conditions. 
  ! 
  ! The basic outline of the routine is as follows: We are given a
  ! flattened list of coordinates containing all coordinates in all
  ! blocks along with the sizes of each blocks. Using the size
  ! information (and the total number of blocks) we can go through
  ! inCoords and temporarily reconstruct the coordinates of each block
  ! to be a 3d array (actually 4D due to the face there are three
  ! x-y-z coordinates). 
  !
  ! The key observation is that it is easier to deal with connection
  ! *faces* instead of nodes. The reason is that every face (a side of
  ! a hexahedreal element) must be connected to precisely one other
  ! element (a block to block connection) or not connected to any
  ! other element (a boundary condition).  To that end, for each block
  ! we've reconstructed, we compute the center of each face on each of
  ! the 6 sides of a block (iMin, iMax, jMin, jMax, kMin, kMax). Given
  ! the block sizes, we can easily determine the total number of these
  ! faces. These spatial coordinates are gathered into an (flatted)
  ! 2-array (3xN) called 'coords'. While we are gathering the
  ! coordinates we also record the block index, the faceID and the
  ! index of the face into the 'info' array (size 5xN).  Once we have
  ! these two arrays, we create a KD-tree using the spatial
  ! coordinates.  This enables fast spatial searchs. 

  ! The main part of the algorithm begins at the lowest face index of
  ! a given block (iBlock (1->nBlocks)) and face (iFace (i->6)). We
  ! know the spatial coordinate of the face center, and can then
  ! search the tree for the two closest other points. We will
  ! obviously will always find the point we searched for since it is
  ! already in the tree. However, if both distances are zero, that
  ! means we found another face on aother block (or potentially the
  ! same block) that is attached the face I searched for. Using the
  ! index returned from the KD-tree we can index into the gobal info
  ! array to determine the blockID, faceID and i, j, k values of the
  ! connection. After checking the (1,1) index on my block we proceed
  ! to the (2, 1). If the other blockID and other faceID remain the
  ! same and one index on the other block only changed by 1, we are on
  ! the same connection (or boundary condition). We keep incrementing
  ! the I index until we hit either the end of the block or we have a
  ! change in blockID/faceID. That determines the I-extent of the
  ! sub-patch on my block. We then run the same algorithm in the
  ! Y-direction....tracing along the j-direction to find the extent of
  ! the sub-patch. Again, once we hit the end of the block or a change
  ! in blockID/faceID, we stop. That now has determined the range of
  ! my subface. Also, when we are tracing out the I and J-directions
  ! we determine which index has changed in the other block. This is
  ! necessary to form the cgns "transform" array. 

  ! Once we have identified the sub-patch on my block, and determined
  ! two of the three values of hte transform array, the third value
  ! can be determine by considering the faceID's of the two connected
  ! blocks. 

  ! The information for a B2B or boundary is then saved in a singely
  ! linked list of patch information. An auxilary array called
  ! faceConsumed keeps track of how much of the face we have already
  ! used. Once a patch is found with range (iStart:iEnd, jStart:jEnd)
  ! we flag the values in the faceConsumed array. Then we loop through
  ! the faceConsumed array to find the lowest (i,j) index that has not
  ! already been consumed. The entire operation then restarts at the
  ! (i,j) faceID that has not already been "consumed". An aribtrary
  ! number of subpatches can exist. Once all of faceConsumed is
  ! .True. this block's face is complete and we move to the next face.

  ! A little ascii art may help explain how subfaces work:


  ! +------------------------+--------------------+------------------------+
  ! |                        |                    |                        |
  ! |                        |                    |                        |
  ! |                        |                    |                        |
  ! |                        |                    |                        |
  ! |                        |    subPatch 4      |      subPatch 5        |
  ! |                        |                    |                        |
  ! |      subPatch 3        |                    |                        |
  ! |                        |                    |                        |
  ! |                        |                    |                        |
  ! |                        |                    |                        |
  ! |                        +--------------------+------------------------+
  ! |                        |                                             |
  ! |                        |                                             |
  ! |                        |                                             |
  ! |                        |                                             |
  ! +------------------------+                                             |
  ! |           (iEnd, jEnd) |                                             |
  ! |                      ^ |                subPatch 2                   |
  ! |                      | |                                         ^   |
  ! |                      | |                                         |   |
  ! |      subPatch 1      | |                                         |   |
  ! |                      ^ |                                         ^   |
  ! |                      | |                                         |   |
  ! |                      | |                                         |   |
  ! | (1,1)---->  (iEnd, j)| | (iStart, jStart)------------->------>   |   |
  ! +------------------------+---------------------------------------------+

  ! Starting at (1,) we trace along until the bottom until we hit
  ! subpatch 2 (becuase the blockID or faceID changed) and then move
  ! in the J direction until again we hit a change in
  ! blockID/faceID. SubPatch 2 is then found next becuase (iStart,
  ! iEnd) is the lowest index that has not been consumed. Likewise,
  ! the remaining patches are found in the order labelled. 


  use kdtree2_module
  use dataTypes
  implicit none

  ! Input/Output
  real(kind=8), dimension(3*nCoords), intent(in):: inCoords
  integer, intent(in) :: nCoords, nBlock
  integer, dimension(3, nBlock), intent(in) :: sizes
  real(kind=8), intent(in) :: tol
  ! Working
  integer :: i, j, k, ii, jj, iDim, iBlock, faceCount, iFace, oFace, faceIndex
  integer :: masterCount
  type(kdtree2), pointer :: tree
  real(kind=8), dimension(:, :), allocatable :: coords
  real(kind=8), dimension(:, :, :, :), allocatable :: X
  integer, dimension(:, :), allocatable :: info
  logical, dimension(:, :), allocatable :: faceConsumed
  type(patch), pointer :: curPatch
  integer :: il, jl ,kl
  integer :: curOtherBlock, curOtherFace, curOtherI, curOtherJ, curOtherK
  integer :: nextOtherBlock, nextOtherFace, nextOtherI, nextOtherJ, nextOtherK
  integer :: delI, delJ, delK, iStart, iEnd, jStart, jEnd
  integer :: IDirIndex, jDirIndex, transform(3), T(3,3), a, b, c
  integer, dimension(3) :: index1, index2, begin1, begin2
  logical :: complete, redo
  real(kind=8), dimension(3) :: n1, n2, n3, n4, v1, v2, ss
  ! Setup the linked list for the patches
  nullify(patches)
  nPatches = 0

  ! Figure out the total number of block faces
  ii = 0
  do iBlock=1, nBlock
     il = sizes(1, iBlock)
     jl = sizes(2, iBlock)
     kl = sizes(3, iBlock)

     ii = ii + 2*(  (il-1)*(jl-1) + (il-1)*(kl-1) + (jl-1)*(kl-1))
  end do

  ! Allocate the space for the nodes on all faces
  allocate(coords(3, ii), info(5, ii))

  ii = 0
  jj = 0
  do iBlock=1, nBlock
     ! Extract block sizes to make things easier to read
     il = sizes(1, iBlock)
     jl = sizes(2, iBlock)
     kl = sizes(3, iBlock)

     ! Reconstruct the blockd ata
     allocate(X(3, il, jl, kl))

     ! Copy out the data...note the order...the data was in python (C)
     ! ordering so we transform it while we copy it out.
     do i=1, il
        do j=1, jl
           do k=1, kl
              do iDim=1,3
                 jj = jj + 1
                 X(iDim, i, j, k) = inCoords(jj)
              end do
           end do
        end do
     end do

     ! iMin Face
     do k=1, kl-1
        do j=1, jl-1
           ii = ii + 1
           coords(:, ii) = X(:, 1, j, k) + X(:, 1, j, k+1) + X(:, 1, j+1, k) + X(:, 1, j+1, k+1)
           info(:, ii) = (/iBlock, iMin, 1, j, k/)
        end do
     end do

     ! iMax Face
     do k=1, kl-1
        do j=1, jl-1
           ii = ii + 1
           coords(:, ii) = X(:, il, j, k) + X(:, il, j, k+1) + &
                X(:, il, j+1, k) + X(:, il, j+1, k+1)
           info(:, ii) = (/iBlock, iMax,  il, j, k/)
        end do
     end do

     ! jMin Face
     do k=1, kl-1
        do i=1, il-1
           ii = ii + 1
           coords(:, ii) = X(:, i, 1, k) + X(:, i, 1, k+1) + &
                X(:, i+1, 1, k) + X(:, i+1, 1, k+1)
           info(:, ii) = (/iBlock, jMin, i, 1, k/)
        end do
     end do

     ! jMax Face
     do k=1, kl-1
        do i=1, il-1
           ii = ii + 1
           coords(:, ii) = X(:, i, jl, k) + X(:, i, jl, k+1) + &
                X(:, i+1, jl, k) + X(:, i+1, jl, k+1)
           info(:, ii) = (/iBlock, jMax, i, jl, k/)
        end do
     end do

     ! kMin Face
     do j=1, jl-1
        do i=1, il-1
           ii = ii + 1
           coords(:, ii) = X(:, i, j, 1) + X(:, i, j+1, 1) + &
                X(:, i+1, j, 1) + X(:, i+1, j+1, 1)
           info(:, ii) = (/iBlock, kMin, i, j, 1/)
        end do
     end do

     ! kMax Face
     do j=1, jl-1
        do i=1, il-1
           ii = ii + 1
           coords(:, ii) = X(:, i, j, kl) + X(:, i, j+1, kl) + &
                X(:, i+1, j, kl) + X(:, i+1, j+1, kl)
           info(:, ii) = (/iBlock, kMax, i, j, kl/)
        end do
     end do
     deallocate(X)
  end do

  ! Divide everything by 4 since we didn't above
  coords = 0.25*coords

  ! Next thing we do is build the KD-tree
  tree => kdtree2_create(coords)

  ! Next we loop over each of the blocks:
  faceCount = 0
  blockLoop: do iBlock=1, nBlock
     faceLoop: do iFace=1,6

        select case(iFace)
        case (iMin, iMax)
           il = sizes(2, iBlock)-1
           jl = sizes(3, iBlock)-1
        case (jMin, jMax)
           il = sizes(1, iBlock)-1
           jl = sizes(3, iBlock)-1
        case (kMin, kMax)
           il = sizes(1, iBlock)-1
           jl = sizes(2, iBlock)-1
        end select

        ! Allocate space for the information for my face, and the
        ! faces we are potentially connected to. 
        allocate(faceConsumed(il, jl))
        faceConsumed = .False. 
        masterCount = 0

        ! Initialize the ranges
        iStart = 1
        jStart = 1
        iEnd = iStart
        jEnd = jStart
        complete = .False.
        do while (.not. complete) 

           call otherInfo(iStart, jStart, curOtherBlock, curOtherFace, &
                curOtherI, curOtherJ, curOtherk, tol)

           begin2 = (/curOtherI, curOtherJ, curOtherk/)
           iDirIndex = 0
           jDirIndex = 0

           ! Find out how far we can go in the I-direction
           iLoop: do while (iEnd < il)
              iEnd = iEnd + 1

              call otherInfo(iEnd, jStart, nextOtherBlock, nextOtherFace, &
                   nextOtherI, nextOtherJ, nextOtherk, tol)

              delI = abs(nextOtherI - curOtherI)
              delJ = abs(nextOtherJ - curOtherJ)
              delK = abs(nextOtherK - curOtherK)
              
              if (  (nextOtherBlock == curOtherBlock  .and. &
                   nextOtherBlock == -1)  .or.  & ! BC
                   
                   (nextOtherBlock == curOtherBlock .and. & 
                    nextOtherFace  == curOtherFace .and. &
                   (delI == 1 .or. delJ == 1 .or. delK == 1))) then  ! B2B
                 
                 ! Determine iDirIndex if not done so already

                 if (iDirIndex == 0 .and. nextOtherBlock /= -1) then 
                    if (delI == 1) then 
                       iDirIndex = sign(1, nextOtherI - curOtherI)
                    else if (delJ == 1) then 
                       iDirIndex = sign(2, nextOtherJ - curOtherJ)
                    else if (delK == 1) then 
                       iDirIndex = sign(3, nextOtherK - curOtherk)
                    else
                       print *,'Something is screwed up 1...'
                       stop
                    end if
                 end if
                 curOtherBlock = nextOtherBlock
                 curOtherFace = nextOtherFace
                 curOtherI = nextOtherI
                 curOtherJ = nextOtherJ
                 curOtherK = nextOtherK
                
              else
                 ! we've stepped to far so back off one.
                 iEnd = iEnd - 1
                 exit 
              end if
           end do iLoop
           
           call otherInfo(iEnd, jStart, curOtherBlock, curOtherFace, &
                curOtherI, curOtherJ, curOtherk, tol)

           ! Find out how far we can go in the J-direction
           jDirIndex = 0
           
           jLoop: do while (jEnd < jl)
              jEnd = jEnd + 1
              
              call otherInfo(iEnd, jend, nextOtherBlock, nextOtherFace, &
                   nextOtherI, nextOtherJ, nextOtherk, tol)

              delI = abs(nextOtherI - curOtherI)
              delJ = abs(nextOtherJ - curOtherJ)
              delK = abs(nextOtherK - curOtherK)

              if (  (nextOtherBlock == curOtherBlock  .and. &
                   nextOtherBlock == -1)  .or.  & ! BC
                   
                   (nextOtherBlock == curOtherBlock .and. & 
                    nextOtherFace == curOtherFace .and. &
                   (delI == 1 .or. delJ == 1 .or. delK == 1))) then  ! B2B
                 
                 ! Determine jDirIndex if not done so already
                 if (jDirIndex == 0 .and. nextOtherBlock /= -1) then
                    if (delI == 1) then 
                       jDirIndex = sign(1, nextOtherI - curOtherI)
                    else if (delJ == 1) then 
                       jDirIndex = sign(2, nextOtherJ - curOtherJ)
                    else if (delK == 1) then 
                       jDirIndex = sign(3, nextOtherK - curOtherK)
                    else
                       print *,'Something is screwed up 2...'
                       stop
                    end if
                 end if

                 curOtherBlock = nextOtherBlock
                 curOtherFace = nextOtherFace
                 curOtherI = nextOtherI
                 curOtherJ = nextOtherJ
                 curOtherK = nextOtherK
              else
                 ! We've stepped too far so back off 1
                 jEnd = jEnd - 1
                 exit 
              end if
           end do jLoop
           
           ! Now the range of the patch we just found is simply
           ! iStart->iEnd, jStart->jEnd. We need to check the faceID to
           ! determine the other index
           
           ! First allocation
           if (.not. associated(patches)) then 
              allocate(patches)
              patches%next => patches
              curPatch => patches
           else
              ! Add the next patch
              allocate(curPatch%next)
              curPatch%next%next => patches
              curPatch => curPatch%next
           end if
           nPatches = nPatches + 1

           select case (iFace)
              
           case (iMin)
              curPatch%pointRange = reshape(&
                   (/1, iStart, jStart, 1, iEnd+1, jEnd+1/), (/3,2/))
           case (iMax)
              curPatch%pointRange = reshape(&
                   (/sizes(1, iBlock), iStart, jStart, sizes(1, iBLock), iEnd+1, jEnd+1/), (/3,2/))
           case (jMin)
              curPatch%pointRange = reshape(&
                   (/iStart, 1, jStart, iEnd+1, 1, jEnd+1/), (/3,2/))
           case (jMax)
              curPatch%pointRange = reshape(&
                   (/iStart, sizes(2, iBlock), jStart, iEnd+1, sizes(2, iBlock), jEnd+1/), (/3,2/))
           case (kMin)
              curPatch%pointRange = reshape(&
                   (/iStart, jStart, 1, iEnd+1, jEnd+1, 1/), (/3,2/))
           case (kMax)
              curPatch%pointRange = reshape(&
                   (/iStart, jStart, sizes(3, iBlock), iEnd+1, jEnd+1, sizes(3, iBlock)/), (/3,2/))
           end select
           
           ! Default type to BC
           curPatch%type = BC
           curPatch%myID = iBlock

           if (curOtherBlock /= -1) then 

              oFace = curOtherFace
              
              ! This is a B2B patch so we need a bit more info
              curPatch%type = B2B
              curPatch%donorID = curOtherBlock
              
              ! The first two if blocks here are for the special case
              ! when there is only once face in teh i, j, or i and j
              ! directions. Becuase of these we couldn't determine the
              ! {i,j}DirIndex. If there is only 1 missing, there is
              ! only one choise that is left and we arbitrarily make
              ! it positive. When both are missing, it is simplier, we
              ! just assume the positve indices based on the faceID of
              ! the other face. 

              if (iDirIndex == 0 .and. jDirIndex /= 0) then 
                 select case(oFace)
                 case (iMin, iMax)
                    if (abs(jDirIndex) == 2) then 
                       iDirIndex = 3
                    else
                       iDirIndex = 2
                    end if
                 case (jMin, jMax)
                    if (abs(jDirIndex) == 1) then
                       iDirIndex = 3
                    else
                       iDirIndex = 1
                    end if
                    
                 case (kMin, kMax)
                    if (abs(jDirIndex) == 1) then 
                       iDirIndex = 2
                    else
                       iDirIndex = 1
                    end if
                 end select
                 
              else if (iDirIndex /=0 .and. jDirIndex == 0) then 
                 select case(oFace)
                 case (iMin, iMax)
                    if (abs(iDirIndex) == 2) then 
                       jDirIndex = 3
                    else
                       jDirIndex = 2
                    end if
                 case (jMin, jMax)
                    if (abs(iDirIndex) == 1) then
                       jDirIndex = 3
                    else
                       jDirIndex = 1
                    end if

                 case (kMin, kMax)
                    if (abs(iDirIndex) == 1) then 
                       jDirIndex = 2
                    else
                       jDirIndex = 1
                    end if
                 end select
              else if (iDirIndex == 0 .or. jDirIndex == 0) then 

                 ! The ordering of iDirIndex and jDirIndex doesn't matter
                 select case(oFace)
                 case (iMin, iMax)
                    iDirIndex = 2
                    jDirIndex = 3
                 case (jMin, jMax)
                    iDirIndex = 1
                    jDirIndex = 3
                 case (kMin, kMax)
                    iDirIndex = 1
                    jDirIndex = 2
                 end select
              end if
              
              ! We need to compute the stupid transform array. Once we
              ! know the transform array, we can use that to correctly
              ! determine the donorPointRange.

              select case(iFace)
              case (iMin, iMax)
                 transform(2) = iDirIndex
                 transform(3) = jDirIndex
              case(jMin, jMax)
                 transform(1) = iDirIndex
                 transform(3) = jDirIndex
              case (kMin, kMax)
                 transform(1) = iDirIndex
                 transform(2) = jDirIndex
              end select
              
              ! Do the final normal face transformation
              call setNormalTransform(transform, iFace, oFace)
              
              curPatch%transform = transform

              ! According to the cgns documentation once T is known,
              ! the following holds: 

              ! index2 = T.(index1 - beging1) + begin2

              a = transform(1)
              b = transform(2)
              c = transform(3)

              T(1, 1) = sgn(a)*del(a, 1)
              T(1, 2) = sgn(b)*del(b, 1)
              T(1, 3) = sgn(c)*del(c, 1)

              T(2, 1) = sgn(a)*del(a, 2)
              T(2, 2) = sgn(b)*del(b, 2)
              T(2, 3) = sgn(c)*del(c, 2)

              T(3, 1) = sgn(a)*del(a, 3)
              T(3, 2) = sgn(b)*del(b, 3)
              T(3, 3) = sgn(c)*del(c, 3)

              index1 = curPatch%pointRange(:, 2)
              begin1 = curPatch%pointRange(:, 1)
              index2 = matmul(T, (index1 - begin1)) + begin2

              ! We may have to do a correction here: Since begin2 was
              ! based on a face value, if the index2(1..3) is *lower*
              ! than begin2(1..3) we have to correct the upper index
              ! by adding 1 and then redoing the transformation
              ! again. This is a bit of a hack, but I don't think there
              ! is a way around it since you don't know at the
              ! beginnign if begin2 will be a lower bound or an upper
              ! bound.
              
              redo = .False.
              do idim=1,3
                 if (index2(idim) < begin2(idim)) then 
                    begin2(idim) = begin2(idim) + 1
                    redo = .True.
                 end if
              end do
              
              ! Update the upper range (index2) for the donor
              if (redo) then 
                 index2 = matmul(T, (index1 - begin1)) + begin2
              end if

              curPatch%pointRangeDonor(:, 1) = begin2
              curPatch%pointRangeDonor(:, 2) = index2
           end if
       
           ! Before we're completely done with this patch, extract the
           ! coordinates of the faces at the start and end:

           faceIndex = faceCount + (jStart-1)*il + iStart
           n1 = coords(:, faceIndex)

           faceIndex = faceCount + (jStart-1)*il + iEnd
           n2 = coords(:, faceIndex)

           faceIndex = faceCount + (jEnd-1)*il + iStart
           n3 = coords(:, faceIndex)

           faceIndex = faceCount + (jEnd-1)*il + iEnd
           n4 = coords(:, faceIndex)

           ! Average of these 4 nodes
           curPatch%faceAvg = (n1 + n2 + n3 + n4)/4

           ! And now the (normalized) normal
           v1 = n4 - n1
           v2 = n3 - n2

           ss(1) = (v1(2)*v2(3) - v1(3)*v2(2))
           ss(2) = (v1(3)*v2(1) - v1(1)*v2(3))
           ss(3) = (v1(1)*v2(2) - v1(2)*v2(1))

           ! Set the faceNormal
           curPatch%faceNormal = ss / sqrt(ss(1)**2 + ss(2)**2 + ss(3)**2 + 1e-14)

           ! Indiscrimentaly flag everything in the range as consumed
           faceConsumed(iStart:iEnd, jStart:jEnd) = .True.

           ! Determine the next availble starting index. If we used
           ! everything up, we set complete to .True. to exit the
           ! while loop

           do while (masterCount < il*jl)

              masterCount = masterCount + 1

              i = mod(masterCount-1, il) + 1
              j = (masterCount-1)/il + 1

              if (.not. faceConsumed(i, j) ) then 

                 ! Set the ranges for the next face. 
                 iStart = i
                 jStart = j
                 iEnd = iStart
                 jEnd = jStart
                 exit
              end if
           end do

           if (masterCount == il*jl) then 
              ! We're done with the face
              complete = .True.
           end if

        end do ! While loop for sub-faces

        deallocate(faceConsumed)

        ! Increment the face counter by the total number of faces gone
        ! through
        faceCount = faceCount + il*jl
     end do faceLoop
  end do blockLoop

contains 

  subroutine setNormalTransform(transform, iFace, oFace)
    implicit none
    integer, intent(inout) :: transform(3)
    integer, intent(in) :: iFace, oFace
    integer ::  normalIndexSign

    ! To get the sign of the third face you sum the two face
    ! ids. If that is even, then index is negative, other wise
    ! it is positive. Why does this work? {i,j,k}Min faces are
    ! odd (1,3,5) while {i,j,k}Max faces are even. A min-min
    ! connection or a max-max connection has as the normal
    ! index increase away from the face on each side, so the
    ! index must be negative (opposite direction). Summing and
    ! min+min ID or max+max ID must be position. Thefore for a
    ! min attached to a max has an odd sum and the sign of the
    ! transform is positive. 

    normalIndexSign = 1
    if (mod( iFace + oFace, 2) == 0) then 
       normalIndexSign = -1
    end if

    select case(iFace)
    case (iMin, iMax)
       select case(oFace)
       case (iMin, iMax)
          transform(1) = 1*normalIndexSign
       case(jMin, jMax)
          transform(1) = 2*normalIndexSign
       case (kMin, kMax) 
          transform(1) = 3*normalIndexSign
       end select

    case (jMin, jMax)
       select case(oFace)
       case (iMin, iMax)
          transform(2) = 1*normalIndexSign
       case(jMin, jMax)
          transform(2) = 2*normalIndexSign
       case (kMin, kMax) 
          transform(2) = 3*normalIndexSign
       end select

    case (kMin, kMax)
       select case(oFace)
       case (iMin, iMax)
          transform(3) = 1*normalIndexSign
       case(jMin, jMax)
          transform(3) = 2*normalIndexSign
       case (kMin, kMax) 
          transform(3) = 3*normalIndexSign
       end select
    end select
  end subroutine setNormalTransform
  
  function sgn(x)
    implicit none
    integer :: sgn, x
    if (x >= 0) then 
       sgn = 1
    else
       sgn = -1
    end if
  end function sgn

  function del(x, y)
    implicit none
    integer :: del, x,y
    if (abs(x) == abs(y)) then 
       del = 1
    else
       del = 0
    end if
  end function del
  
  subroutine otherInfo(i, j, otherBlock, otherFace, otherI, otherJ, otherK, tol)
    implicit none
    
    integer, intent(in) :: i, j
    integer, intent(out) :: otherBlock, otherFace, otherI, otherJ, otherK
    integer :: faceIndex
    real(kind=8) :: qv(3)
    type(kdtree2_result) :: results(2)
    real(kind=8), intent(in) :: tol

    ! Reconstruct the index
    faceIndex = faceCount + (j-1)*il + i
    
    ! Now search for a neighbour:
    qv = coords(:, faceIndex)
    call kdtree2_n_nearest(tree, qv, 2, results)
    
    ! Check the distances
    if (sqrt(results(1)%dis) < tol .and. sqrt(results(2)%dis) < tol) then
       ! We found a match: Now determine if it was the first
       ! or second one that isn't myself:
       
       if (results(1)%idx == faceIndex) then  ! Must be the second one
          otherBlock = info(1, results(2)%idx)
          otherFace  = info(2, results(2)%idx)
          otherI     = info(3, results(2)%idx)
          otherJ     = info(4, results(2)%idx)
          otherK     = info(5, results(2)%idx)
       else
          otherBlock = info(1, results(1)%idx)
          otherFace  = info(2, results(1)%idx)
          otherI     = info(3, results(1)%idx)
          otherJ     = info(4, results(1)%idx)
          otherK     = info(5, results(1)%idx)
       end if
    else
       ! We have a BC so flag 'otherInfo' so we know not to
       ! deal with it later.
       otherBlock = -1
       otherFace  = -1
       otherI     = -1
       otherJ     = -1
       otherK     = -1
    end if
  end subroutine otherInfo
end subroutine computeConnectivity

subroutine getPatchInfo(n, types, pointRanges, myIDs, pointRangeDonors, &
     transforms, donorIDs, faceAvgs, faceNormals)

  use dataTypes
  implicit none

  ! Input variables
  integer, intent(in) :: n

  ! Output variables
  integer, dimension(n), intent(out) :: myIDs, donorIDs, types
  integer, dimension(3, 2, n), intent(out) :: pointRanges, pointRangeDonors
  integer, dimension(3, n), intent(out) :: transforms
  real(kind=8), dimension(3, n), intent(out) :: faceAvgs, faceNormals
  ! Working variables
  type(patch), pointer :: curPatch
  integer :: i

  i = 0 
  curPatch => patches
  do while (i < n)
     i = i + 1
     myIDs(i) = curPatch%myID
     types(i) = curPatch%type
     pointRanges(:, :, i) = curPatch%pointRange
     pointRangeDonors(:, :, i) = curPatch%pointRangeDonor
     transforms(:, i) = curPatch%transform
     donorIDs(i) = curPatch%donorID
     faceAvgs(:, i) = curPatch%faceAvg
     faceNormals(:, i) = curPatch%faceNormal
     curPatch => curPatch%next
  end do
end subroutine getPatchInfo

subroutine getnPatches(n)
  use dataTypes
  implicit none
  integer , intent(out) :: N
  n = nPatches
end subroutine getnPatches

subroutine deallocpatches()
  use dataTypes
  implicit none
  integer :: i
  type(patch), pointer :: curPatch, tmp
  
  ! Deallocate all the memory associated with the list of patches

  curPatch => patches
  i = 0
  do while (i < nPatches)
     i = i + 1
     tmp => curPatch%next
     deallocate(curPatch)
     curPatch => tmp
  end  do
  nPatches= 0
end subroutine deallocpatches

!=============================================!
!                                             !
!         CARTESIAN MESH SUBROUTINES          !
!                                             !
!=============================================!
! Ney Secco, 2015
! neysecco@umich.edu

subroutine findBounds(x, xBounds, il, jl, kl)

  ! xBounds has 6 elements that represents the maximum and minumim values of
  ! coordinates in this block:
  ! xBounds = ((xmin, ymin, zmin),
  !            (xmax, ymax, zmax))
  
  implicit none
  
  ! Subroutine inputs
  real(kind=8), dimension(il,jl,kl,3), intent(in) :: x
  integer, intent(in) :: il, jl, kl

  ! Subroutine inputs/outputs.
  real(kind=8), dimension(2,3), intent(inout) :: xBounds

  !
  ! BEGIN EXECUTION
  !

  ! Find the maximum bounds in each face of the block. xBounds will
  ! be updated within each call
  
  ! Face imin
  call checkBounds(1,1,1,jl,1,kl,x,xBounds)
  ! Face imax
  call checkBounds(il,il,1,jl,1,kl,x,xBounds)
  ! Face jmin
  call checkBounds(1,il,1,1,1,kl,x,xBounds)
  ! Face jmax
  call checkBounds(1,il,jl,jl,1,kl,x,xBounds)
  ! Face kmin
  call checkBounds(1,il,1,jl,1,1,x,xBounds)
  ! Face kmax
  call checkBounds(1,il,1,jl,kl,kl,x,xBounds)
  
  ! xBounds was updated and will be returned by this subroutine
  
contains
  
  !        ================================================================
  
  subroutine checkBounds(imin,imax,jmin,jmax,kmin,kmax,x,xBounds)
    
    ! This function will compute the maximum values within the bounds
    ! specified by imin, imax, jmin, jmax, kmin, kmax.
    ! This is useful when we want to loop over each surface.
    ! xBounds will be updated
    
    implicit none
    
    ! Subroutine Inputs
    integer, intent(in) :: imin, imax, jmin, jmax, kmin, kmax
    real(kind=8), dimension(:,:,:,:), intent(in) :: x
    
    ! Subroutine Outputs
    real(kind=8), dimension(2,3), intent(inout) :: xBounds
    
    ! Working variables
    integer :: i, j, k
    
    !
    ! BEGIN EXECUTION
    !
    
    do k=kmin, kmax
       do j=jmin, jmax
          do i=imin, imax
             
             ! Check the coordinates bounds
             if (x(i,j,k,1) < xBounds(1,1)) xBounds(1,1) = x(i,j,k,1)
             if (x(i,j,k,2) < xBounds(1,2)) xBounds(1,2) = x(i,j,k,2)
             if (x(i,j,k,3) < xBounds(1,3)) xBounds(1,3) = x(i,j,k,3)
             if (x(i,j,k,1) > xBounds(2,1)) xBounds(2,1) = x(i,j,k,1)
             if (x(i,j,k,2) > xBounds(2,2)) xBounds(2,2) = x(i,j,k,2)
             if (x(i,j,k,3) > xBounds(2,3)) xBounds(2,3) = x(i,j,k,3)
             
          end do
       end do
    end do
    
  end subroutine checkBounds

end subroutine findBounds

!        ================================================================
!        ================================================================
!        ================================================================
!        ================================================================

subroutine computeVolumes(x, xBounds, binVolX, binVolY, binVolZ, &
     binCellsX, binCellsY, binCellsZ, il, jl, kl, nBinX, nBinY, nBinZ)

  ! This subroutine will loop over all the faces of the block, compute their
  ! volumes, and if this volume is greater than the reference volume stored
  ! in maxVol for the corresponding bin, then the reference value will be
  ! updated.
  ! xBounds has 6 elements that represents the maximum and minumim values of
  ! coordinates in this block:
  ! xBounds = ((xmin, ymin, zmin),
  !            (xmax, ymax, zmax))
  ! maxVol = maximum volume in each bin (will be updated then returned)

  implicit none

  ! Subroutine inputs
  real(kind=8), dimension(il,jl,kl,3), intent(in) :: x
  real(kind=8), dimension(2,3), intent(in) :: xBounds
  integer, intent(in) :: il, jl, kl, nBinX, nBinY, nBinZ

  ! Subroutine inputs/outputs.
  real(kind=8), dimension(nBinX), intent(inout) :: binVolX
  real(kind=8), dimension(nBinY), intent(inout) :: binVolY
  real(kind=8), dimension(nBinZ), intent(inout) :: binVolZ
  integer, dimension(nBinX), intent(inout) :: binCellsX
  integer, dimension(nBinY), intent(inout) :: binCellsY
  integer, dimension(nBinZ), intent(inout) :: binCellsZ

  !
  ! BEGIN EXECUTION
  !

  ! Now that we have information about the total bounds, we can
  ! compute the volumes of the cells and assign them to their
  ! respective bins. maxVol will be updated after each call
  
  ! DOING PER FACE
  ! Face imin
  !call checkVol(2,2,2,jl,2,kl,x,xBounds,maxVol)
  ! Face imax
  !call checkVol(il,il,2,jl,2,kl,x,xBounds,maxVol)
  ! Face jmin
  !call checkVol(2,il,2,2,2,kl,x,xBounds,maxVol)
  ! Face jmax
  !call checkVol(2,il,jl,jl,2,kl,x,xBounds,maxVol)
  ! Face kmin
  !call checkVol(2,il,2,jl,2,2,x,xBounds,maxVol)
  ! Face kmax
  !call checkVol(2,il,2,jl,kl,kl,x,xBounds,maxVol)

  ! Doing full block
  call checkVol(2, il, 2, jl, 2, kl, x, xBounds, &
       binVolX, binVolY, binVolZ, binCellsX, binCellsY, binCellsZ, &
       nBinX, nBinY, nBinZ)

contains
  
  !        ================================================================

  subroutine checkVol(imin, imax, jmin, jmax, kmin, kmax, x, xBounds, &
       binVolX, binVolY, binVolZ, binCellsX, binCellsY, binCellsZ, &
       nBinX, nBinY, nBinZ)

    ! This function will compute the maximum volumes within the bounds
    ! specified by imin, imax, jmin, jmax, kmin, kmax and assign to the correct bin.
    ! This is useful when we want to loop over each surface.
    ! binVol and binCells will be updated
    
    implicit none

    ! Subroutine Inputs
    integer, intent(in) :: imin, imax, jmin, jmax, kmin, kmax
    real(kind=8), dimension(:,:,:,:), intent(in) :: x
    real(kind=8), dimension(2,3), intent(in) :: xBounds
    integer, intent(in) :: nBinX, nBinY, nBinZ

    ! Subroutine  Outputs
    real(kind=8), dimension(:), intent(inout) :: binVolX
    real(kind=8), dimension(:), intent(inout) :: binVolY
    real(kind=8), dimension(:), intent(inout) :: binVolZ
    integer, dimension(:), intent(inout) :: binCellsX
    integer, dimension(:), intent(inout) :: binCellsY
    integer, dimension(:), intent(inout) :: binCellsZ

    ! Working variables
    real(kind=8) :: vol, xp, yp, zp, eighth
    real(kind=8) :: vp1, vp2, vp3, vp4, vp5, vp6
    integer :: i, j, k, l, m, n, iBin, jBin, kBin

    ! Compute the volumes. The hexahedron is split into 6 pyramids
    ! whose volumes are computed. The volume is positive for a
    ! right handed block.
    ! Initialize the volumes to zero. The reasons is that the second
    ! level halo's must be initialized to zero and for convenience
    ! all the volumes are set to zero.

    do k=kmin, kmax
       n = k -1
       
       do j=jmin, jmax
          m = j -1
          
          do i=imin, imax
             l = i -1

             ! Compute the coordinates of the center of gravity.
             
             eighth = 1.0/8.0

             xp = eighth*(x(i,j,k,1) + x(i,m,k,1) &
                  +         x(i,m,n,1) + x(i,j,n,1) &
                  +         x(l,j,k,1) + x(l,m,k,1) &
                  +         x(l,m,n,1) + x(l,j,n,1))
             yp = eighth*(x(i,j,k,2) + x(i,m,k,2) &
                  +         x(i,m,n,2) + x(i,j,n,2) &
                  +         x(l,j,k,2) + x(l,m,k,2) &
                  +         x(l,m,n,2) + x(l,j,n,2))
             zp = eighth*(x(i,j,k,3) + x(i,m,k,3) &
                  +         x(i,m,n,3) + x(i,j,n,3) &
                  +         x(l,j,k,3) + x(l,m,k,3) &
                  +         x(l,m,n,3) + x(l,j,n,3))
             
             ! Compute the volumes of the 6 sub pyramids. The
             ! arguments of volpym must be such that for a (regular)
             ! right handed hexahedron all volumes are positive.
             
             vp1 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                  x(i,j,n,1), x(i,j,n,2), x(i,j,n,3), &
                  x(i,m,n,1), x(i,m,n,2), x(i,m,n,3), &
                  x(i,m,k,1), x(i,m,k,2), x(i,m,k,3),xp,yp,zp)
             
             vp2 = volpym(x(l,j,k,1), x(l,j,k,2), x(l,j,k,3), &
                  x(l,m,k,1), x(l,m,k,2), x(l,m,k,3), &
                  x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                  x(l,j,n,1), x(l,j,n,2), x(l,j,n,3),xp,yp,zp)
               
             vp3 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                  x(l,j,k,1), x(l,j,k,2), x(l,j,k,3), &
                  x(l,j,n,1), x(l,j,n,2), x(l,j,n,3), &
                  x(i,j,n,1), x(i,j,n,2), x(i,j,n,3),xp,yp,zp)
               
             vp4 = volpym(x(i,m,k,1), x(i,m,k,2), x(i,m,k,3), &
                  x(i,m,n,1), x(i,m,n,2), x(i,m,n,3), &
                  x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                  x(l,m,k,1), x(l,m,k,2), x(l,m,k,3),xp,yp,zp)
             
             vp5 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                  x(i,m,k,1), x(i,m,k,2), x(i,m,k,3), &
                  x(l,m,k,1), x(l,m,k,2), x(l,m,k,3), &
                  x(l,j,k,1), x(l,j,k,2), x(l,j,k,3),xp,yp,zp)
             
             vp6 = volpym(x(i,j,n,1), x(i,j,n,2), x(i,j,n,3), &
                  x(l,j,n,1), x(l,j,n,2), x(l,j,n,3), &
                  x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                  x(i,m,n,1), x(i,m,n,2), x(i,m,n,3),xp,yp,zp)
             
             ! Set the volume to 1/6 of the sum of the volumes of the
             ! pyramid. Remember that volpym computes 6 times the
             ! volume.
             
             ! Avoid negative volumes
             vp1 = abs(vp1)
             vp2 = abs(vp2)
             vp3 = abs(vp3)
             vp4 = abs(vp4)
             vp5 = abs(vp5)
             vp6 = abs(vp6)

             vol = (vp1 + vp2 + vp3 + vp4 + vp5 + vp6)/6.0

             ! Now we need to find the bin which the current element belongs
             ! x coordinate
             iBin = findBin(xBounds(1,1), xBounds(2,1), nBinX, xp)
             ! y coordinate
             jBin = findBin(xBounds(1,2), xBounds(2,2), nBinY, yp)
             ! z coordinate
             kBin = findBin(xBounds(1,3), xBounds(2,3), nBinZ, zp)
             
             ! Increment the cells counter of the bins
             binCellsX(iBin) = binCellsX(iBin) + 1
             binCellsY(jBin) = binCellsY(jBin) + 1
             binCellsZ(kBin) = binCellsZ(kBin) + 1

             ! Update the average volumes of the bins
             binVolX(iBin) = ((binCellsX(iBin)-1)*binVolX(iBin) + vol)/binCellsX(iBin)
             binVolY(jBin) = ((binCellsY(jBin)-1)*binVolY(jBin) + vol)/binCellsY(jBin)
             binVolZ(kBin) = ((binCellsZ(kBin)-1)*binVolZ(kBin) + vol)/binCellsZ(kBin)
             
          end do
       end do
    end do
    
  end subroutine checkVol
             
  !        ================================================================
  
  !        ================================================================
  
  function volpym(xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd,xp,yp,zp)
    !
    !        ****************************************************************
    !        *                                                              *
    !        * volpym computes 6 times the volume of a pyramid. Node p,     *
    !        * whose coordinates are set in the subroutine metric itself,   *
    !        * is the top node and a-b-c-d is the quadrilateral surface.    *
    !        * It is assumed that the cross product vCa * vDb points in     *
    !        * the direction of the top node. Here vCa is the diagonal      *
    !        * running from node c to node a and vDb the diagonal from      *
    !        * node d to node b.                                            *
    !        *                                                              *
    !        ****************************************************************
    !

    implicit none
    !
    !        Function type.
    !
    real(kind=8) :: volpym
    !
    !        Function arguments.
    !
    real(kind=8), intent(in) :: xa, ya, za, xb, yb, zb
    real(kind=8), intent(in) :: xc, yc, zc, xd, yd, zd
    real(kind=8), intent(in) :: xp, yp, zp
    !
    !        Working
    !
    real(kind=8) :: fourth

    !
    !        ****************************************************************
    !        *                                                              *
    !        * Begin execution                                              *
    !        *                                                              *
    !        ****************************************************************
    !
    
    fourth = 1.0/4.0

    volpym = (xp - fourth*(xa + xb  + xc + xd))              &
         * ((ya - yc)*(zb - zd) - (za - zc)*(yb - yd))   + &
         (yp - fourth*(ya + yb  + yc + yd))              &
         * ((za - zc)*(xb - xd) - (xa - xc)*(zb - zd))   + &
         (zp - fourth*(za + zb  + zc + zd))              &
         * ((xa - xc)*(yb - yd) - (ya - yc)*(xb - xd))

  end function volpym

  !        ================================================================
  
  !        ================================================================

  function findBin(xmin,xmax,numBins,x)
    
    ! This function returns the bin index where the coordinate x
    ! belongs when the interval [xmin,xmax] is split in numBins bins
    
    implicit none
    
    ! Function type
    integer :: findBin
    
    ! Function inputs
    real(kind=8), intent(in) :: xmin, xmax, x
    integer, intent(in) :: numBins
    
    ! Working variables
    real(kind=8) :: dx
    
    !
    ! BEGIN EXECUTION
    !
    
    ! Find the bin size
    dx = (xmax - xmin)/numBins
    
    ! Find the bin index
    findBin = floor((x-xmin)/dx) + 1
    
  end function findBin
  
end subroutine computeVolumes


!=============================================!
!                                             !
!     END OF CARTESIAN MESH SUBROUTINES       !
!                                             !
!=============================================!


subroutine calcGridRatio(N, s0, S, ratio)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: calcGridRatio() calculates the exponential
  !     distribution Turns out we need to solve a transendental
  !     equation here. We will do this with a bisection search
  !
  !     Parameters
  !     ----------
  !     N : integer
  !         The number of nodes in sequence
  !     s0 : real
  !         The initial grid spacing
  !     S : real
  !         The total integrated length
  !     
  !     Returns
  !     -------
  !     ratio : real
  !         The computed grid ratio that satifies, N, s0, and S.


  implicit none

  ! Input Parameters
  integer(kind=4), intent(in) :: N
  real(kind=8), intent(in) :: s0, S

  ! Output Parameters
  real(kind=8), intent(out) :: ratio

  ! Working Parameters
  integer(kind=4) :: i, M
  real(kind=8) ::  r, a,b, c, f, fa, fb

  ! function 'f' is S - s0*(1-r^n)/(1-r) where S is total length, s0 is
  ! initial ratio and r is the grid ratio. 

  M = N-1

  ! Do a bisection search
  ! Max and min bounds...root must be in here...
  a = 1.0_8 + 1e-8
  b = 4.0_8

  fa = S - s0*(1-a**M)/(1-a)
  fb = S - s0*(1-b**M)/(1-b)
  do i=1, 100
     c = 0.5_8*(a + b)
     f = S - s0*(1-c**M)/(1-c)
     if (abs(f) < 1e-6) then ! Converged
        exit
     end if

     if (f * fa > 0) then 
        a = c
     else
        b = c
     end if
  end do

  ! Finally set the ratio variable to r
  ratio = c

end subroutine calcGridRatio
