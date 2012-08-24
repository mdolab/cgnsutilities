! This program is used to convert surface files from SUmb into binary
! tecplot files. This script optionally strips off the symmetry and
! far field zones which are often not needed. 


program cgns_surf_convert

  implicit none
  include 'cgnslib_f.h'

  integer, dimension(3*3) :: isize
  integer, dimension(6) :: rinde

  ! File handles for cgns files
  integer :: cg_current

  ! CGNS counters for base
  integer current_base

  ! CGNS Names for base, zone and iterative data
  character*32 basename
  character*1024 input_name
  character*1024 output_name
  character*1024 max_name

  character*32 fieldname
  character*36 zonename
  ! CGNS Data type
  integer datatype

  ! CGNS Names 
  character*32 gridName,solName
  character*32, dimension(:),allocatable :: fieldNames
  character*32, dimension(:),allocatable :: var_names

  ! CGNS cell dimension and physical dimension
  integer  CellDim, PhysDim

  ! CGNS number of bases and zones in current file
  integer nbases, nzones

  ! CGNS error flag
  integer ier

  ! Integer Counters etc
  integer :: N,i,j
  integer :: nx,ny,nz,il,jl,kl
  integer :: nsols,sol,nfields,field
  integer :: temp_shape(6),rmax(3)
  integer :: farStart, symStart, cgnsStart

  ! Character arrays for names and I/O
  character*32, dimension(3) :: coordNames

  ! Misc Data Arrays
  double precision, allocatable,dimension(:,:,:  ) :: data3d
  double precision, allocatable,dimension(:,:,:,:) :: gridTemp
  double precision, allocatable,dimension(:,:,:  ) :: fieldData
  double precision :: xmin(3), xmax(3)

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
  xmin(:) = 0d0
  xmax(:) = 0d0
  N = IARGC ()

  if (.not. ( N == 1 )) then
     print *,'Error: cgns_surf_convert must be called with ONE argument:'
     print *,'./cgns_surf_convert <infile.cgns>'
     print *,'Output file will be called <infile.plt>'
     stop
  end if
  CALL GETARG(1 , input_name)

  cgnsStart = INDEX (trim(input_name),'.cgns')
  if (cgnsStart == 0) then
     print *,'Error: no .cgns extension found'
     stop
  end if
 10 format(A,A)
  write(output_name,10),input_name(1:cgnsStart),'plt'
  write(max_name,10),input_name(1:cgnsStart-1),'_max.txt'
  print *,'Converting file:',trim(input_name)
  ! ------------------------------------------------------------------
  ! Open a file to get the number of zones:

  call cg_open_f(input_name,MODE_READ,cg_current,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_nbases_f(cg_current, nbases, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  if (nbases .lt. 1) then
     print *, 'Error: No bases found in first CGNS file'
     stop
  end if

  ! We will only deal with 1 base in each zone so ALWAYS only use 1 base
  current_base = 1

  ! Get the cell and phys dim 
  call cg_base_read_f(cg_current, current_base, basename, &
       CellDim, PhysDim, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_nzones_f(cg_current, current_base, nzones, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_nsols_f(cg_current, current_base, 1, nsols , ier )
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  if (nsols >= 1) then
     call cg_nfields_f(cg_current,current_base,1,1,nfields,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
  else
     nfields = 0
  end if

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
  ier = TECINI112("CGNS Surf Data"//char(0),var_names,trim(output_name)//char(0),"."//char(0),&
       FileType,Debug,VIsDouble)

  ! ------------------------------------------------------------------
  ! Loop over the zones to get the required size info:
  do i=1,nzones

     call cg_zone_read_f(cg_current,current_base,i,&
          zoneName,isize,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f


     farStart = INDEX (trim(zoneName),'Far')
     symStart = INDEX (trim(zoneName),'Sym')

     if (farStart == 0 .and. symStart == 0) then 

        ! Load the zone sizes:
        call cg_zone_read_f(cg_current,current_base,i,zoneName,isize,ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f

        if (cellDim == 2) then
           nx = isize(1)
           ny = isize(2)
           nz = 1
           il = isize(3)
           jl = isize(4)
           kl = 1
        else
           nx = isize(1)
           ny = isize(2)
           nz = isize(3)
           il = isize(4)
           jl = isize(5)
           kl = isize(6)
        end if

        ! Write Tecplot Zone Data
        ier = TECZNE112(trim(zoneName),ZoneType,nx,ny,nz,&
             ICellMax,JCellMax,KCellMax, &
             0.0,0,ParentZn,IsBlock,NFConns,FNMode,&
             TotalNumFaceNodes,TotalNumBndryFaces,TotalNumBndryConnections,&
             PassiveVarList, valueLocation, ShareVarFromZone,ShrConn)

        !-------------------------------------------------------------
        !                         Grid Coodinates 
        !-------------------------------------------------------------

        allocate(gridTemp(nx,ny,nz,physDim))

        do j=1,physDim
           call cg_coord_read_f(cg_current,current_base,i,&
                coordNames(j),RealDouble,(/1,1,1/),(/nx,ny,nz/),&
                gridTemp(:,:,:,j),ier)
        end do

        ! Get the min and max vals:
        if (i == 1) then
           do j=1,physDim
              xmin(j) = minval(gridTemp(:,:,:,j))
              xmax(j) = maxval(gridTemp(:,:,:,j))
           end do
        else
           do j=1,physDim
              xmin(j) = minval((/minval(gridTemp(:,:,:,j)),xmin(j)/))
              xmax(j) = maxval((/maxval(gridTemp(:,:,:,j)),xmax(j)/))
           end do
        end if
        

        ! Write Grid Coordinates
        do j=1,physDim
           ier   = TECDAT112(nx*ny*nz,gridTemp(:,:,:,j), DIsDouble)
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
              allocate(fieldData(nx,ny,nz))

              call cg_field_info_f(cg_current,current_base,i,1,j,&
                   datatype , fieldname , ier )
              if (ier .eq. CG_ERROR) call cg_error_exit_f
              temp_shape(:) = 0
              temp_shape = (/1 - rinde(1), il + rinde(2), &
                   1 - rinde(3), jl + rinde(4), &
                   1 - rinde(5), kl + rinde(6)/)

              ! Allocate Temporary Array for Cell Data, possibly with rind cells
              allocate(data3d(temp_shape(1):temp_shape(2),&
                   temp_shape(3):temp_shape(4),&
                   temp_shape(5):temp_shape(6)))

              rmax(1) = il+rinde(1)+rinde(2)
              rmax(2) = jl+rinde(3)+rinde(4)
              rmax(3) = kl+rinde(5)+rinde(6)

              call cg_field_read_f(cg_current,current_base,i,1, &
                   fieldName, datatype, (/1,1,1/),&
                   rmax,data3d,ier)
              if (ier .eq. CG_ERROR) call cg_error_exit_f

              ! Call the interpolate function to reconstruct node data
              call interpolate(data3d,temp_shape,&
                   fieldData,nx,ny,nz)

              ! Deallocate Temporary Data Array
              deallocate(data3d)

              ! Write out Field Values
              ier   = TECDAT112(nx*ny*nz,fieldData, DIsDouble)
              deallocate(fieldData)
           end do ! Field Loop
        end if ! If we have fields
     end if
  end do ! Zone Loop

  call cg_close_f(cg_current,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  ! Close Tecplot File
  ier = TECEND112()

  ! Open the max file and write the 6 values
  OPEN(UNIT=11, FILE=max_name, STATUS="replace")
  write(11, '(15f)') xmin(1)
  write(11, '(15f)') xmin(2)
  write(11, '(15f)') xmin(3)
  write(11, '(15f)') xmax(1)
  write(11, '(15f)') xmax(2)
  write(11, '(15f)') xmax(3)
  close(11)

end program cgns_surf_convert

subroutine interpolate(cell_data,cell_shape,node_data,inode,jnode,knode)

  ! Take cell centered data which possibly contain rind cells and
  ! perform a linear reconstruction back to the nodes. cell_shape
  ! contains the start and end values for each index in
  ! cell_data. Retruned data comes back in node_data. It may not be
  ! the most efficient possible, but it works.

  implicit none

  integer, intent(in) :: cell_shape(6),inode,jnode,knode

  double precision, intent(in ), dimension(cell_shape(1):cell_shape(2),&
       cell_shape(3):cell_shape(4),&
       cell_shape(5):cell_shape(6)) :: cell_data
  double precision, intent(out), dimension(inode,jnode,knode) :: node_data
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
