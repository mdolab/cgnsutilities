program cgns_scale

  implicit none
  include 'cgnslib_f.h'

!   integer Ndim,boco
!   parameter (Ndim = 3)

!  , zonetype
!   character*32 filename
!   integer  base, zone,j
! size(Ndim*3),zonesize(3)
!   character*32 basename, zonename,boconame
!   integer nbocos,bocotype
!   integer NormalIndex(3), NormalListFlag, ndataset,datatype
!   integer ptset_type, npnts, pnts(6)
!   double precision data_double(6)
!   integer nfamilies,nFamBC,nGeo,famID,ifam,BC,coordID,start(3)


  integer, dimension(3*3) :: isize

  ! File handles for cgns files: cg_current is the current step,
  ! cg_unsteady is the new unsteady file
  integer :: cg_current, cg_unsteady

  ! CGNS Zone Counter
  integer  zone

  ! CGNS counters for bases and zones
  integer current_base,current_zonem,unsteady_base,unsteady_zone
  
  ! CGNS Names for base, zone and iterative data
  character*32 basename, current_zonename, baseitername
  character*37 ArbitraryGridMotionName,ZoneIterName

  ! CGNS Name for unsteady grid
  character*32 gridName
  character*32, dimension(:),allocatable :: gridNames

  ! CGNS cell dimension and physical dimension
  integer  CellDim, PhysDim

  ! CGNS number of bases and zones in current file
  integer nbases, nzones
  
  ! CGNS error flag
  integer ier
  integer :: N,i_start,i_end,n_steps,i,j,counter
  double precision :: dt
  character*100 :: char_temp
  character*128 base_name,unsteady_filename
  character*132 current_filename
  character*32, dimension(3) :: coordNames
  ! Misc Data Arrays
  double precision, allocatable,dimension(:      ) :: data1d
  double precision, allocatable,dimension(:,:    ) :: data2d
  double precision, allocatable,dimension(:,:,:  ) :: data3d
  double precision, allocatable,dimension(:,:,:,:) :: data4d
  integer         , allocatable,dimension(:)       :: int1d
  integer :: dummy_int

  coordNames(1) = "CoordinateX"
  coordNames(2) = "CoordinateY"
  coordNames(3) = "CoordinateZ"

  ! Start of Code
  N = IARGC ()

  if (N .ne. 4) then
     print *,'Error: cgns_combine must be called with THREE arguments:'
     print *,'./cgns_combine <base_name> <start> <stop> <dt>'
     print *,'<base_name> is the file name of the time step solutions'
     print *,'without the number at the end'
     print *,'<start> is starting timestep number'
     print *,'<end> is last timestep number'
     print *,'<dt> Is the physical time step'
     stop
  end if
  CALL GETARG(1 , base_name)
  CALL GETARG(2 , char_temp)
  Read(char_temp, '(i)' )  i_start
  CALL GETARG(3 , char_temp)
  Read(char_temp, '(i)' )  i_end
  CALL GETARG(4 , char_temp)
  Read(char_temp, '(f)' )  dt

  print *,'Check of Input:'
  print *,'Base Name     :',trim(base_name)
  print *,'starting step :',i_start
  print *,'ending step   :',i_end
  print *,'dt            :',dt
  n_steps = i_end-i_start + 1
  allocate(gridNames(n_steps))
  ! ------------------------------------------------------------------
  ! Open the first time-step solution to get the physical and cell dim

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
  call cg_base_read_f(cg_current, current_base, basename, &
       CellDim, PhysDim, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_close_f(cg_current,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  ! ------------------------------------------------------------------
  ! Create the new unsteady file:

  ! First create the new filename:
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
  do i=1,n_steps
     data1d(i) = (dble(i)-1)/n_steps * dt
     int1d(i) = i
  end do
  call cg_goto_f(cg_unsteady,unsteady_base,ier,'BaseIterativeData_t',1,'end')
  !call cg_array_write_f('TimeValues',RealDouble,1,n_steps,data1d,ier)
  !if (ier .eq. CG_ERROR) call cg_error_exit_f
  call cg_array_write_f('IterationValues',Integer,1,n_steps,int1d,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  deallocate(data1d)
  call cg_simulation_type_write_f(cg_unsteady,unsteady_base,TimeAccurate,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f


  ! ------------------------------------------------------------------
  ! Master Loop Over the Unsteady Files
  counter = 0
 do i=i_start,i_end
    counter = counter + 1
!  do i=1,1
10   format(A,I4)
     write(*,10),'Processing File:',i

     ! Open Current File
     call getFilename(base_name,i,current_filename)
     call cg_open_f(current_filename,MODE_READ,cg_current,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     
     ! Get the number of zones in current file
     ! Goto Base Node in Current File

     call cg_goto_f(cg_current, current_base, ier, 'end')
     call cg_nzones_f(cg_current, current_base, nzones, ier)

     do zone=1,nzones

      
        ! Read Zone from current file and write
        call cg_zone_read_f(cg_current,current_base,zone,current_zonename,isize,ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f

11 format(I4.4,'_',A)
        write(ZoneIterName,11), counter,current_zonename

!        if (i == i_start) then
           call cg_zone_write_f(cg_unsteady,unsteady_base,ZoneIterName,isize, &
                Structured,dummy_int,ier)
           if (ier .eq. CG_ERROR) call cg_error_exit_f

        !    zoneIterName = 'Iterative Zone Data'
!            call cg_ziter_write_f(cg_unsteady,unsteady_base,zone, ZoneIterName, ier )
!            if (ier .eq. CG_ERROR) call cg_error_exit_f

           

           if (cellDim == 2) then
              allocate(data2d(isize(1),isize(2)))
              do j=1,physDim
                 call cg_coord_read_f(cg_current,current_base,zone,&
                      coordNames(j),RealDouble,(/1,1/),(/isize(1),isize(2)/),&
                      data2d,ier)
                 if (ier .eq. CG_ERROR) call cg_error_exit_f
                 call cg_coord_write_f(cg_unsteady,unsteady_base,zone+(counter-1)*nzones,RealDouble,&
                      coordNames(j),data2d,dummy_int,ier)
              end do
              deallocate(data2d)
           end if
!        end if

! 11      format(A,I4.4)
!         write(ArbitraryGridMotionName,11),'Motion',i

!         write(gridName,11),'Grid',i
!         gridNames(counter) = gridName
    
!         if (cellDim == 2) then
!            allocate(data2d(isize(1),isize(2)))
!            do j=1,physDim
!               data2d = 0.0
!               call cg_coord_read_f(cg_current,current_base,zone,&
!                    coordNames(j),RealDouble,(/1,1/),(/isize(1),isize(2)/),&
!                    data2d,ier)
!               if (ier .eq. CG_ERROR) call cg_error_exit_f

!               ! Only needs to written once
!               if (j==1) then
!                  call cg_grid_write_f(cg_unsteady,unsteady_base,zone,gridName, dummy_int , ier )
!                  if (ier .eq. CG_ERROR) call cg_error_exit_f
!                  call cg_arbitrary_motion_write_f(&
!                       cg_unsteady,unsteady_base,zone, ArbitraryGridMotionName,&
!                       DeformingGrid,dummy_int,ier)
!                  if (ier .eq. CG_ERROR) call cg_error_exit_f
!               END if
!               call cg_goto_f(cg_unsteady,unsteady_base, ier, 'Zone_t', zone,'GridCoordinates_t',counter + 1,'end')
            
!               call cg_array_write_f(coordNames(j), RealDouble, 2,(/isize(1),isize(2)/),&
!                    data2d,ier)
!               if (ier .eq. CG_ERROR) call cg_error_exit_f

              
                          
!            end do
!            deallocate(data2d)
!         end if

     end do
  
  ! Close the current file
  call cg_close_f(cg_current,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
end do
  


do zone=1,nzones
   call cg_goto_f(cg_unsteady,unsteady_base,ier,'Zone_t',zone,'ZoneIterativeData_t',1,'end')
   call cg_array_write_f('GridCoordinatesPointers',Character,2,(/32,n_steps/),GridNames,ier)  
   if (ier .eq. CG_ERROR) call cg_error_exit_f
end do

  ! Close the unsteady file
  call cg_close_f(cg_unsteady,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  deallocate(gridNames)
end program cgns_scale

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

1000  format(A,I1)
1001  format(A,I2)
1002  format(A,I3)
1003  format(A,I4)

end subroutine getFileName
