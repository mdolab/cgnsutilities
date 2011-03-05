module zone_vars

  ! This module stores data associated with a single zone and with n_time instances
  implicit none
  save 

  type zoneData

     integer :: nx,ny,nz
     integer :: il,jl,kl
     integer :: n_time, n_vars

     ! Grid has shape: nx by ny by nz ny 3 by n_time
     
     double precision, dimension(:,:,:,:,:), pointer :: ts_grid
     double precision, dimension(:,:,:,:,:), pointer :: grid

     ! Vars has shape: il by jl by kl by n_vars by n_time
     
     double precision, dimension(:,:,:,:,:),pointer :: ts_vars
     double precision, dimension(:,:,:,:,:),pointer :: ev_in
     double precision, dimension(:,:,:,:,:),pointer :: ev_out
     double precision, dimension(:,:,:,:,:),pointer :: vars
     
  end type zoneData
end module zone_vars

program ts_combine

  use zone_vars
  implicit none
  include 'cgnslib_f.h'

  integer, dimension(3*3) :: isize

  ! File handles for cgns files: cg_current is the current step,
  ! cg_unsteady is the new unsteady file
  integer :: cg_current, cg_unsteady

!   ! CGNS Zone Counter
!   integer  zone

  ! CGNS counters for base
  integer current_base,unsteady_base

  ! CGNS Names for base, zone and iterative data
  character*32 basename, current_zonename, baseitername
  character*37 ArbitraryGridMotionName,ZoneIterName
  character*32 fieldname

  ! CGNS Data type
  integer datatype

  ! CGNS Name for unsteady grid
  character*32 gridName,solName
  character*32, dimension(:),allocatable :: motionNames,solNames
  character*32, dimension(:),allocatable :: gridNames
  character*32, dimension(:),allocatable :: zoneNames
  character*32, dimension(:),allocatable :: fieldNames


  ! CGNS cell dimension and physical dimension
  integer  CellDim, PhysDim

  ! CGNS number of bases and zones in current file
  integer nbases, nzones

  ! CGNS error flag
  integer ier

  ! Integer Counters etc
  integer :: N,i_start,i_end,n_steps,n_time,i,j
  integer :: ii,jj,kk,ind(3)
  integer :: counter,nsols,sol,nfields,field,izone
  logical, dimension(5) :: EV_Exist
  

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

  ! Start of Code
  coordNames(1) = "CoordinateX"
  coordNames(2) = "CoordinateY"
  coordNames(3) = "CoordinateZ"

  N = IARGC ()

  if (N .ne. 3) then
     print *,'Error: cgns_combine must be called with FOUR arguments:'
     print *,'./ts_combine <base_name> <n_files> <n_steps>'
     print *,'<base_name> is the file name of the TS solutions'
     print *,'without the number at the end'
     print *,'<n_files> is number of time spectral solutions'
     print *,'<n_steps> is the number of steps to interpolate'
     stop
  end if
  CALL GETARG(1 , base_name)
  CALL GETARG(2 , char_temp)
  Read(char_temp, '(i4)' )  n_time
  CALL GETARG(3 , char_temp)
  Read(char_temp, '(i4)' )  n_steps

  print *,'Check of Input:'
  print *,'Base Name     :',trim(base_name)
  print *,'n_files       :',n_time
  print *,'n_steps       :',n_steps

  ! Allocate some names
  allocate(gridNames(n_steps),motionNames(n_steps),solNames(n_steps))
  ! ------------------------------------------------------------------
  ! Open a file to get the number of zones:
  call getFilename(base_name,1,current_filename)

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
   
     call cg_nfields_f(cg_current,current_base,i,1,nfields,ier)
     if (i==1) then
        allocate(fieldNames(nfields))
     end if
     if (ier .eq. CG_ERROR) call cg_error_exit_f
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
     allocate(zone(i)%ts_grid(zone(i)%nx,zone(i)%ny,zone(i)%nz,3,n_time),stat=ier)
     allocate(zone(i)%grid(zone(i)%nx,zone(i)%ny,zone(i)%nz,3,n_steps))

     ! Allocate Var Coords

     ! Input and Ouput Fields
     allocate(zone(i)%ts_vars(zone(i)%il,zone(i)%jl,zone(i)%kl,nfields,n_time))
     !allocate(zone(i)%vars(zone(i)%il,zone(i)%jl,zone(i)%kl,nfields,n_steps))

     ! Conservative Input and Output Fields (ev => Euler Variables)
     allocate(zone(i)%ev_in(zone(i)%il,zone(i)%jl,zone(i)%kl,5,n_time))
     allocate(zone(i)%ev_out(zone(i)%il,zone(i)%jl,zone(i)%kl,5,n_steps))

     zone(i)%grid = 0.0
     zone(i)%ev_in = 0.0
     zone(i)%ev_out = 0.0
  end do

  call cg_close_f(cg_current,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
 
  ! ------------------------------------------------------------------
  ! Load in all the data:
  print *,'Reading Data from File: '
 
  do n=1,n_time
10   format(A,I4)
11   format(A,I4)
     write(*,11,advance='no'),' ',n
     
     ! Open Current File
     call getFilename(base_name,n,current_filename)
     call cg_open_f(current_filename,MODE_READ,cg_current,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     
     ! Goto Base Node in Current File
     call cg_goto_f(cg_current, current_base, ier, 'end')
     
     ! Loop over the zones
     do i=1,nzones
        call cg_zone_read_f(cg_current,current_base,i,zoneNames(i),isize,ier)
        ! Load the grid coordinates
        do j=1,physDim
           call cg_coord_read_f(cg_current,current_base,i,&
                coordNames(j),RealDouble,(/1,1,1/),&
                (/zone(i)%nx,zone(i)%ny,zone(i)%nz/),&
                zone(i)%ts_grid(:,:,:,j,n),ier)
           if (ier .eq. CG_ERROR) call cg_error_exit_f
        end do

        ! Load the Variables
        do j=1,nFields
            call cg_field_info_f(cg_current,current_base,i,1,j,&
                 datatype , fieldname , ier )
            if (ier .eq. CG_ERROR) call cg_error_exit_f

            call cg_field_read_f(cg_current,current_base,i,1, &
                 fieldName, datatype, (/1,1,1/),&
            (/zone(i)%il,zone(i)%jl,zone(i)%kl/),&
                 zone(i)%ts_vars(:,:,:,j,n) , ier )
            fieldNames(j) = fieldName
            if (ier .eq. CG_ERROR) call cg_error_exit_f
         end do
      end do ! Zone Loop
     
     ! Close the current file
     call cg_close_f(cg_current,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
  end do
  write(*,*),' '

  ! Next we have to figure out what the actual conservative variabes
  ! are, and spectrally interpolate those varibales, and NOT Velocity,
  ! Cp, Mach etc. 
  EV_Exist(:) = .False.
  ! ~~~~~~~~~~~~~~~~~~~~~~ DENSITY ~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! First check for rho, we basically need it for everything:
  call getFieldIndex(fieldNames,nFields,"Density",ind(1))

  if (ind(1) == 0) then
     print *,'Error: "Density" field not found. Must have Density for spectral interpolation'
     stop
  end if
  
  EV_Exist(1) = .True.
  ! Copy over Denstiy
  do i=1,nzones
     zone(i)%ev_in(:,:,:,1,:) = zone(i)%ts_vars(:,:,:,ind(1),:)
  end do
  
  ! ~~~~~~~~~~~~~~~~~ Velocity / Momentum ~~~~~~~~~~~~~~~~~~~~

  call getFieldIndex(fieldNames,nFields,"VelocityX",ind(1))
  call getFieldIndex(fieldNames,nFields,"VelocityY",ind(2))
  call getFieldIndex(fieldNames,nFields,"VelocityZ",ind(3))
  do j = 1,3 ! Loop over the three dimensions
     if (ind(j) .ne. 0) then
        EV_Exist(j+1) = .True.
        do i=1,nzones
           do kk=1,zone(i)%kl
              do jj=1,zone(i)%jl
                 do ii=1,zone(i)%il
                    ! Set the momentum = rho * vel. They are indices 2,3 and 4
                    zone(i)%ev_in(ii,jj,kk,j+1,:) = &
                         zone(i)%ts_vars(ii,jj,kk,ind(j),:)*zone(i)%ev_in(ii,jj,kk,1,:)
                 end do
              end do
           end do
        end do
     end if
  end do

  ! ~~~~~~~~~~~~~~~~~ Pressure  ~~~~~~~~~~~~~~~~~~~~

  call getFieldIndex(fieldNames,nFields,"Pressure",ind(1))

  ! We need to back out rhoE from Pressure
  
  if(ind(1) .ne. 0) then
     EV_Exist(5) = .True.
     ! Copy over Denstiy
     do i=1,nzones
        zone(i)%ev_in(:,:,:,5,:) = zone(i)%ts_vars(:,:,:,ind(1),:)
     end do
  end if
  
  ! Loop over the Euler Data and perform the spectral interpolation
  
  print *,'Generating spectral interpolation data for zone...'

  do i=1,nzones
     write(*,11,advance='no'),' ',i
     
     ! Do the grid coordinates
     do j=1,physDim
        do kk=1,zone(i)%nz
           do jj=1,zone(i)%ny
              do ii=1,zone(i)%nx             
                 call spectralInterpolate(n_time,n_steps,&
                      zone(i)%ts_grid(ii,jj,kk,j,:), &
                      zone(i)%grid   (ii,jj,kk,j,:))
              end do
           end do
        end do
     end do

     ! Dothe Field Data
     do j=1,5
        do kk=1,zone(i)%kl
           do jj=1,zone(i)%jl
              do ii=1,zone(i)%il
                 call spectralInterpolate(n_time,n_steps,&
                      zone(i)%ev_in (ii,jj,kk,j,:), &
                      zone(i)%ev_out(ii,jj,kk,j,:))
                 
              end do
           end do
        end do
     end do
  end do
  write(*,*),' '

  ! ------------------------------------------------------------------
  ! Create the new unsteady file:
  print *,'Writing unsteady CGNS file ...'

  unsteady_filename = base_name

  ! Open the new cgns file
  call cg_open_f(unsteady_filename,CG_MODE_WRITE, cg_unsteady, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  ! Write the first base
  if (celldim == 2) then
     basename = 'Unsteady Surface Data'
  else
     basename = 'Unsteady Volume Data'
  end if
  call cg_base_write_f(cg_unsteady,basename,celldim,physdim,unsteady_base,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  baseitername = 'Base Iterative Data'
  call cg_biter_write_f(cg_unsteady,unsteady_base, BaseIterName, n_steps, ier )
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  allocate(data1d(n_steps))
  allocate(int1d(n_steps))
  dt = 1.0
  do i=1,n_steps
     data1d(i) = (dble(i)-1)/n_steps * dt
     int1d(i) = i
  end do
  call cg_goto_f(cg_unsteady,unsteady_base,ier,'BaseIterativeData_t',1,'end')
  call cg_array_write_f('TimeValues',RealDouble,1,n_steps,data1d,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  call cg_array_write_f('IterationValues',Integer,1,n_steps,int1d,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  deallocate(data1d,int1d)
  call cg_simulation_type_write_f(cg_unsteady,unsteady_base,TimeAccurate,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  ! We now have all the data generated, and can just dump it out into the file

  do i=1,nzones
     isize(:) = 0     
     if (cellDim == 2) then
        isize(1) = zone(i)%nx
        isize(2) = zone(i)%ny
        isize(3) = zone(i)%il
        isize(4) = zone(i)%jl
     else
        isize(1) = zone(i)%nx
        isize(2) = zone(i)%ny
        isize(3) = zone(i)%nz
        isize(4) = zone(i)%il
        isize(5) = zone(i)%jl
        isize(6) = zone(i)%kl
     end if

     call cg_zone_write_f(cg_unsteady,unsteady_base,zoneNames(i),&
          isize,Structured,dummy_int,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

   !   ! Write Iterative Data
      call cg_ziter_write_f(cg_unsteady,unsteady_base,i,'Iter Zone Data', ier )
      if (ier .eq. CG_ERROR) call cg_error_exit_f

     ! Get the unique names
     do j=1,n_steps
        write(motionNames(j),12),'Motion',j
        write(gridNames(j),12),'Grid',j
        write(solNames(j),12),'Solution',j
12      format(A,I4.4)
     end do

     ! Write out the coordinates
     do j=1,physDim
        call cg_coord_write_f(cg_unsteady,unsteady_base,i,RealDouble,&
             coordNames(j),zone(i)%ts_grid(:,:,:,j,1),dummy_int,ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f
     end do

     ! Write the links
      call cg_goto_f(cg_unsteady,unsteady_base,ier,'Zone_t',i,'ZoneIterativeData_t',1,'end')
      call cg_array_write_f('GridCoordinatesPointers',Character,2,(/32,n_steps/),GridNames,ier)  
      if (ier .eq. CG_ERROR) call cg_error_exit_f
      call cg_array_write_f('ArbitraryGridMotionPointers',Character,2,(/32,n_steps/),MotionNames,ier)  
      if (ier .eq. CG_ERROR) call cg_error_exit_f
      call cg_array_write_f('FlowSolutionPointers',Character,2,(/32,n_steps/),SolNames,ier)  
      if (ier .eq. CG_ERROR) call cg_error_exit_f
 
     ! Write out the arbitrary motion data, grid and sol handles
     do j=1,n_steps
        call cg_grid_write_f(cg_unsteady,unsteady_base,i,gridNames(j),&
             dummy_int , ier )
        if (ier .eq. CG_ERROR) call cg_error_exit_f

        call cg_arbitrary_motion_write_f(&
             cg_unsteady,unsteady_base,i,motionNames(j),&
             DeformingGrid,dummy_int,ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f

        call cg_sol_write_f(cg_unsteady,unsteady_base,i,solNames(j),&
             CellCenter, dummy_int,ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f
     end do

     ! Write out the coordinates
    
     do n=1,n_steps
        do j=1,physDim
           call cg_goto_f(cg_unsteady,unsteady_base, ier, 'Zone_t', i,&
                'GridCoordinates_t',n+1,'end')
           call cg_array_write_f(coordNames(j), RealDouble, cellDim,&
           (/zone(i)%nx,zone(i)%ny,zone(i)%nz/),zone(i)%grid(:,:,:,j,n),ier)
           if (ier .eq. CG_ERROR) call cg_error_exit_f
        end do
     end do

!     Write out the variables
     do n=1,n_steps
        !do j=1,nfields
        ! Density
       call cg_field_write_f(cg_unsteady,unsteady_base,i,n,&
            datatype, "Denstiy", zone(i)%ev_out(:,:,:,1,n), dummy_int,ier)
        
        ! VelX
        call cg_field_write_f(cg_unsteady,unsteady_base,i,n,datatype,"VelocityX",&
             zone(i)%ev_out(:,:,:,2,n)/zone(i)%ev_out(:,:,:,1,n), dummy_int,ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f

        ! VelY
        call cg_field_write_f(cg_unsteady,unsteady_base,i,n,datatype,"VelocityY",&
             zone(i)%ev_out(:,:,:,3,n)/zone(i)%ev_out(:,:,:,1,n), dummy_int,ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f

        ! VelZ
        call cg_field_write_f(cg_unsteady,unsteady_base,i,n,datatype,"VelocityZ",&
             zone(i)%ev_out(:,:,:,4,n)/zone(i)%ev_out(:,:,:,1,n), dummy_int,ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f
     end do
  end do
  call cg_close_f(cg_unsteady,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  ! Deallocate data and quit
  deallocate(gridNames,solNames,motionNames)
  deallocate(zoneNames)
  deallocate(zone)  
end program ts_combine

 
subroutine spectralInterpolate(n_in,n_out,data_in,data_out)

  ! Take data, data_in of length n_in and spectrally interpolate to
  ! get data_out of length n_out
  implicit none
  integer, intent(in) :: n_in,n_out
  double precision, intent(in) ::  data_in(n_in) 
  double precision, intent(out) ::  data_out(n_out)
  double precision :: t_out(n_out)
  complex*16 :: X(n_in)
  double precision :: ix,rx,pi
  integer :: n_break,k,n

  pi = 3.14159265358979
  ! Generate the DFT of the input data...
  X(:) = cmplx(0.0,0.0)
  do k=0,n_in-1
     do n=0,n_in-1
        X(k+1) = X(k+1) + cmplx(data_in(n+1),0)* &
             exp(-2*cmplx(pi,0)*cmplx(0,1)*cmplx(k,0)*cmplx(n,0)/cmplx(n_in,0))
     end do
  end do

  ! Generate the output time:
  do n=1,n_out
     t_out(n) = 2*pi*(dble(n)-1)/(n_out)
  end do

  data_out(:) = 0.0
  ! Now do the spectral interpolation
  n_break = (n_in-1)/2
  do n=1,n_out
     do k=0,n_break
        rx = real(X(k+1))
        ix = aimag(X(k+1))
        data_out(n) = data_out(n) + &
             rx*cos(k*t_out(n)) - ix*sin(k*t_out(n))
       
     end do
     do k=n_break+1,n_in-1
        rx = real(X(k+1))
        ix = aimag(X(k+1))
        data_out(n) = data_out(n) + &
             rx*cos((n_in-k)*t_out(n)) + ix*sin((n_in-k)*t_out(n))
     end do
  end do
  data_out = data_out * (1/dble(n_in))

end subroutine spectralInterpolate


subroutine getFileName(base_name,i,current_filename)

  implicit none
  character*128, intent(in) ::  base_name
  character*132, intent(out) ::  current_filename
  integer, intent(in) :: i 

  if (i < 10) then
     write(current_filename,1000),trim(base_name),i
  else if (i < 100) then
     write(current_filename,1001),trim(base_name),i
  else if (i < 1000) then
     write(current_filename,1002),trim(base_name),i
  else if (i < 10000) then
     write(current_filename,1003),trim(base_name),i
  else
     print *,'Error: Index too large!'
     stop
  end if

  current_filename = trim(current_filename)

1000 format(A,I1)
1001 format(A,I2)
1002 format(A,I3)
1003 format(A,I4)

end subroutine getFileName

subroutine getFieldIndex(fieldList,nFields,name,ind)

  ! Get the index, "ind" of "name" in "fieldList" of length
  ! "nFields". ind returns 0 if it isnt' foudn

  implicit none
  integer, intent(in) :: nFields
  character*32, dimension(nFields), intent(in) :: fieldList
  character(*), intent(in) :: name
  integer, intent(out) :: ind

  integer i

  ind = 0

  do i=1,nFields
     if (trim(fieldList(i)) == name) then
        ind = i
     end if
  end do
end subroutine getFieldIndex

