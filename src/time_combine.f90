module zone_vars

  ! This module stores data associated with a single zone and with n_time instances
  implicit none
  save 

  type zoneData

     integer :: nx,ny,nz
     integer :: il,jl,kl
     integer :: n_time, n_vars

     ! Grid has shape: nx by ny by nz ny 3 by n_time

     double precision, dimension(:,:,:,:,:), pointer :: grid

     ! Vars has shape: il by jl by kl by n_vars by n_time

     double precision, dimension(:,:,:,:,:),pointer :: vars

  end type zoneData
end module zone_vars

program time_combine

  use zone_vars
  implicit none
  include 'cgnslib_f.h'

  integer, dimension(3*3) :: isize
  integer, dimension(6) :: rinde

  ! File handles for cgns files: cg_current is the current step,
  ! cg_unsteady is the new unsteady file
  integer :: cg_current, cg_unsteady

  !   ! CGNS Zone Counter
  !   integer  zone

  ! CGNS counters for base
  integer current_base,unsteady_base

  ! CGNS Names for base, zone and iterative data
  character*32 basename,  baseitername
  character*32 output_name
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
  integer :: temp_shape(6),rmax(3)

  ! Time step
  double precision :: dt

  ! Character arrays for names and I/O
  character*100 :: char_temp
  character*128 base_name,unsteady_filename
  character*132 current_filename
  character*32, dimension(3) :: coordNames
  ! Array of all zone data
  type(zoneData), dimension(:), allocatable :: zone

  ! Misc Data Arrays
  double precision, allocatable,dimension(:      ) :: data1d
  double precision, allocatable,dimension(:,:    ) :: data2d
  double precision, allocatable,dimension(:,:,:  ) :: data3d
  double precision, allocatable,dimension(:,:,:,:) :: data4d
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

  N = IARGC ()

  if (.not. ( N == 4 .or. N == 5)) then
     print *,'Error: cgns_combine must be called with FOUR OR FIVE arguments:'
     print *,'./cgns_combine <base_name> <start> <stop> <dt> [output_name]'
     print *,'<base_name> is the file name of the time step solutions'
     print *,'use %d, %1.1d,%2.2d,%3.3d or %4.4d to specify counters'
     print *,'<start> is starting timestep number'
     print *,'<end> is last timestep number'
     print *,'<dt> Is the physical time step'
     print *,'[output_name] is the output Tecplot File'
     stop
  end if
  CALL GETARG(1 , base_name)
  CALL GETARG(2 , char_temp)
  Read(char_temp, '(i4)' )  i_start
  CALL GETARG(3 , char_temp)
  Read(char_temp, '(i4)' )  i_end
  CALL GETARG(4 , char_temp)
  Read(char_temp, '(f16.12)' )  dt
  if (N==5) then
     CALL GETARG(5 , output_name)
  else
     output_name = 'data.plt'
  end if

  print *,'Check of Input:'
  print *,'Base Name     :',trim(base_name)
  print *,'starting step :',i_start
  print *,'ending step   :',i_end
  print *,'dt            :',dt
  print *,'Output Name   :',output_name
  n_steps = i_end-i_start + 1

  ! ------------------------------------------------------------------
  ! Open a file to get the number of zones:
  call getFilename(base_name,i_start,current_filename)

  call cg_open_f(current_filename,MODE_READ,cg_current,ier)
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

     call cg_zone_read_f(cg_current,current_base,i,&
          zoneNames(i),isize,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     call cg_nsols_f(cg_current, current_base, i, nsols , ier )
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     if (nsols >= 1) then
        call cg_nfields_f(cg_current,current_base,i,1,nfields,ier)
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

     ! Allocate Grid Coords
     allocate(zone(i)%grid(zone(i)%nx,zone(i)%ny,zone(i)%nz,3,n_steps))
     zone(i)%grid = 0.0
     ! Allocate Var Coords

     allocate(zone(i)%vars(zone(i)%nx,zone(i)% ny,zone(i)%nz,nfields,n_steps))
     zone(i)%vars(:,:,:,:,:) = 0.0

  end do

  call cg_close_f(cg_current,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f



  ! ------------------------------------------------------------------
  ! Load in all the data:
  print *,'Reading Data from File: '

  allocate(fieldNames(nfields))

  nn = 0
  do n=i_start,i_end
10   format(A,I4)
11   format(A,I4)
     write(6,11,advance='yes'),'File ',n
     nn = nn + 1
     !  ! Open Current File
     call getFilename(base_name,n,current_filename)
     call cg_open_f(current_filename,MODE_READ,cg_current,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     !      ! Goto Base Node in Current File
     call cg_goto_f(cg_current, current_base, ier, 'end')
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     ! Loop over the zones
     do i=1,nzones
        call cg_zone_read_f(cg_current,current_base,i,zoneNames(i),isize,ier)
        ! Load the grid coordinates

        rmax(1) = zone(i)%nx
        rmax(2) = zone(i)%ny
        rmax(3) = zone(i)%nz

        do j=1,physDim
           call cg_coord_read_f(cg_current,current_base,i,&
                coordNames(j),RealDouble,(/1,1,1/),&
                rmax,zone(i)%grid(:,:,:,j,nn),ier)
           if (ier .eq. CG_ERROR) call cg_error_exit_f
        end do

        if (nfields > 0) then
           ! Read the Rind Data
           call cg_goto_f(cg_current,current_base,ier,'Zone_t',i, &
                'FlowSolution_t',1,'end')
           rinde(:) = 0
           call cg_rind_read_f(rinde , ier )
           if (ier .eq. CG_ERROR) call cg_error_exit_f

           ! Load the Variables
           do j=1,nFields

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
                   zone(i)%vars(:,:,:,j,nn),&
                   zone(i)%nx,zone(i)%ny,zone(i)%nz)

              ! Deallocate Temporary Data Array
              deallocate(data3d)

              ! DO NOT put fieldNames(j) in the cg_field_read_f call!!!!
              fieldNames(j) = fieldName
           end do
        end if

     end do ! Zone Loop
     ! Close the current file

     call cg_close_f(cg_current,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
  end do
  write(*,*),' '

  ! ------------------------------------------------------------------
  ! Create the new Tecplot File:
  print *,'Writing unsteady Tecplot file ...'

  allocate(var_names(nFields+4))
  allocate(valueLocation(nFields+physDim),&
       passiveVarList(nFields+physDim),&
       ShareVarFromZone(nFields+physDim))

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

  ! Open Tecplot File
  ier = TECINI112("Unsteady Data"//char(0),var_names,trim(output_name)//char(0),"."//char(0),&
       FileType,Debug,VIsDouble)

  ! Loop over timesteps and zones and Write to File
  do n=1,n_steps
     do i=1,nzones
        ier = TECZNE112(trim(zoneNames(i)),ZoneType,zone(i)%nx,zone(i)%ny,zone(i)%nz,&
             ICellMax,JCellMax,KCellMax, &
             dt*dble(n-1),i,ParentZn,IsBlock,NFConns,FNMode,&
             TotalNumFaceNodes,TotalNumBndryFaces,TotalNumBndryConnections,&
             PassiveVarList, valueLocation, ShareVarFromZone,ShrConn)

        ! Write Grid Coordinates
        do j=1,physDim
           ier   = TECDAT112(zone(i)%nx*zone(i)%ny*zone(i)%nz,&
                zone(i)%grid(:,:,:,j,n), DIsDouble)
        end do

        ! Write out Field Values

        do j=1,nFields
           ier   = TECDAT112(zone(i)%nx*zone(i)%ny*zone(i)%nz,&
                zone(i)%vars(:,:,:,j,n), DIsDouble)
        end do
     end do
  end do

  ! Close Tecplot File
  ier = TECEND112()
  ! Deallocate data and quit
  deallocate(zone)  
  deallocate(fieldNames,var_Names,valueLocation)
  deallocate(passiveVarList,ShareVarFromZone)
  deallocate(zoneNames)
end program time_combine

subroutine getFileName(base_name,i,current_filename)

  implicit none
  character*128, intent(in) ::  base_name
  character*128 :: head_str
  character*128 :: tail_str
  integer :: tail_len,head_len,str_len,spec_len
  character*132, intent(out) ::  current_filename
  integer, intent(in) :: i 
  integer :: split_ind,split_ind1,split_ind2,split_ind3,split_ind4,split_ind5
  integer :: s_ind
  CHARACTER(LEN=128) :: FMT

  ! We're going to try to "smart" here...we will use a %d specifier to
  ! determine how to make the filename. Accpeted specifiers are:
  ! %d, %1.1d, %2.2d, %3.3d, %4.4d and %5.5d

  ! Find one of the specifiers in the string:

  split_ind5 = INDEX (base_name, '%5.5d', BACK = .TRUE.)
  split_ind4 = INDEX (base_name, '%4.4d', BACK = .TRUE.)
  split_ind3 = INDEX (base_name, '%3.3d', BACK = .TRUE.)
  split_ind2 = INDEX (base_name, '%2.2d', BACK = .TRUE.)
  split_ind1 = INDEX (base_name, '%1.1d', BACK = .TRUE.)
  split_ind  = INDEX (base_name, '%d', BACK = .TRUE.)
  str_len = len(trim(base_name))

  if (split_ind5 > 0) then
     s_ind = split_ind5
     spec_len = 5
     write(FMT,800),'(A,','I5.5,','A)'
  else if(split_ind4 > 0) then
     s_ind = split_ind4
     spec_len = 5
     write(FMT,800),'(A,','I4.4,','A)'
  else if(split_ind3 > 0) then
     s_ind = split_ind3
     spec_len = 5
     write(FMT,800),'(A,','I3.3,','A)'
  else if(split_ind2 > 0) then
     s_ind = split_ind2
     spec_len = 5
     write(FMT,800),'(A,','I2.2,','A)'
  else if(split_ind1 > 0) then
     s_ind = split_ind1
     spec_len = 5
     write(FMT,800),'(A,','I1.1,','A)'
  else if(split_ind  > 0) then
     s_ind = split_ind
     spec_len = 2
     if (i < 10) then
        write(FMT,800),'(A,','I1,','A)'
     else if (i < 100) then
        write(FMT,800),'(A,','I2,','A)'
     else if (i < 1000) then
        write(FMT,800),'(A,','I3,','A)'
     else if (i < 10000) then
        write(FMT,800),'(A,','I4,','A)'
     else
        print *,'Error: Index too large!'
        stop
     end if
  else
     print *,'a %d specified was not found!'
     stop
  end if

  head_len = s_ind-1
  if (head_len>0) then
     head_str(1:head_len) = base_name(1:head_len)
  end if
  
  tail_len = str_len - spec_len - s_ind + 1
  if (tail_len > 0) then
     tail_str(1:tail_len) = base_name(s_ind+spec_len:str_len)
  end if
  write(current_filename,FMT),head_str(1:head_len),i,tail_str(1:tail_len)

800 format(A,A,A)

end subroutine getFileName

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
