program cgns_mirror

  !use cgnsGrid

  implicit none
  include 'cgnslib_f.h'

  !type(cgnsBlock), dimension(:),allocatable :: blocks
  integer  CellDim, PhysDim
  integer nbases, nzones, in_dims(9),out_dims(9)

  integer ier, zonetype
  integer n,nn,mm,i,j,k,idim

  integer cg_in, cg_out, base, zone, zoneCounter

  integer :: ptSetType, normalDataType
  integer :: sizes(9),coordID,blockStart(3),blockEnd(3)

  integer donor_range(6),transform(3)
  integer nbocos,n1to1,bocotype,BCout
  integer NormalIndex(3), NormalListFlag, ndataset,datatype
  integer ptset_type, npnts, pnts(3,2),pnts_donor(3,2),ncon

  character*32 in_filename,out_filename,mirror
  character*32 basename, zonename,boconame
  character*32 connectname,donorname, famname
  double precision, allocatable,dimension(:,:,:) :: data3d
  double precision, allocatable,dimension(:,:,:) :: data3d_mirr
  double precision data_double(6)

  N = IARGC ()
  if (N .ne. 3) then
     print *,'Error: cgns_mirror must be called with THREE arguments:'
     print *,'./cgns_mirror cgns_in_file.cgns <x,y,z> cgns_out_file.cgns'
     stop
  end if

  call getarg(1, in_filename)
  call getarg(3, mirror)
  call getarg(2, out_filename)

  call cg_open_f(in_filename,CG_MODE_READ, cg_in, ier)

  if (ier .eq. CG_ERROR) call cg_error_exit_f
  print *,'Reading input file: ',in_filename

  call cg_nbases_f(cg_in, nbases, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  if (nbases .lt. 1) then
     print *, 'Error: No bases found in CGNS file'
  else if(nbases .gt. 1) then
     print *, 'Warning: More than one base found. Using only the first base'
  end if

  ! Explicitly use ONLY the first base
  base = 1

  call cg_base_read_f(cg_in, base, basename, CellDim, PhysDim, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  if (cellDim .ne. 3 .or. PhysDim .ne. 3) then
     print *,'The Cells must be hexahedreal in 3 dimensions'
     stop
  end if

  ! Goto Base Node
  call cg_goto_f(cg_in, base, ier, 'end')
  call cg_nzones_f(cg_in, base, nzones, ier)

  ! Open the output file:
  call cg_open_f(out_filename,CG_MODE_WRITE, cg_out, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_base_write_f(cg_out,"Base#1", Celldim,Physdim, base, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  
  do nn=1,nzones
     call cg_zone_read_f(cg_in, base, nn, zonename, in_dims, ier)

     ! ------------------------ Regular Zone ----------------
999  FORMAT('zone',I4.4)
     write(zonename,999) 2*nn-1
     call cg_zone_write_f(cg_out,base,zonename,in_dims,Structured,zoneCounter,ier)
     write(zonename,999) 2*nn
     call cg_zone_write_f(cg_out,base,zonename,in_dims,Structured,zoneCounter,ier)
     allocate(data3d(in_dims(1),in_dims(2),in_dims(3)))
     allocate(data3d_mirr(in_dims(1),in_dims(2),in_dims(3)))

     blockStart = (/1,1,1/)
     blockEnd   = in_dims(1:3)

     ! X coordinate
     call cg_coord_read_f(cg_in,base,nn,'CoordinateX',RealDouble,&
          blockStart,blockEnd,data3d,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     
     if (trim(mirror) == 'x') then
        data3d_mirr = -data3d
     else
        data3d_mirr = data3d
     end if
     
     call cg_coord_write_f(cg_out,base,zoneCounter-1,realDouble,&
          'CoordinateX',data3d, coordID,ier)

     call cg_coord_write_f(cg_out,base,zoneCounter  ,realDouble,&
          'CoordinateX',data3d_mirr, coordID,ier)

     ! Y coordinate
     call cg_coord_read_f(cg_in,base,nn,'CoordinateY',RealDouble,&
          blockStart,blockEnd,data3d,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     
     if (trim(mirror) == 'y') then
        data3d_mirr = -data3d
     else
        data3d_mirr = data3d
     end if
     
     call cg_coord_write_f(cg_out,base,zoneCounter-1,realDouble,&
          'CoordinateY',data3d, coordID,ier)

     call cg_coord_write_f(cg_out,base,zoneCounter  ,realDouble,&
          'CoordinateY',data3d_mirr, coordID,ier)
     

     ! Z coordinate
     call cg_coord_read_f(cg_in,base,nn,'CoordinateZ',RealDouble,&
          blockStart,blockEnd,data3d,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     
     if (trim(mirror) == 'z') then
        data3d_mirr = -data3d
     else
        data3d_mirr = data3d
     end if
     
     call cg_coord_write_f(cg_out,base,zoneCounter-1,realDouble,&
          'CoordinateZ',data3d, coordID,ier)

     call cg_coord_write_f(cg_out,base,zoneCounter  ,realDouble,&
          'CoordinateZ',data3d_mirr, coordID,ier)
     
     deallocate(data3d,data3d_mirr)


  end do ! Zone Loop

  ! Close both CGNS files
  call cg_close_f(cg_in, ier)
  call cg_close_f(cg_out,ier)

end program cgns_mirror
