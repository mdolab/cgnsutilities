program cgns_flip

  implicit none
  include 'cgnslib_f.h'

  integer Ndim,boco
  parameter (Ndim = 3)
  integer  CellDim, PhysDim
  integer ier, n, zonetype
  character*32 filename
  integer cg, base, zone,j
  integer nbases, nzones, size(Ndim*3),zonesize(3)
  character*32 basename, zonename,boconame
  integer nbocos,bocotype
  integer NormalIndex(3), NormalListFlag, ndataset,datatype
  integer ptset_type, npnts, pnts(6)
  double precision data_double(6)
  integer nfamilies,nFamBC,nGeo,famID,ifam,BC,coordID,start(3)
  character *100 FamilyName,famName,char_flip
  double precision, allocatable,dimension(:,:,:,:) :: coorX,coorY,coorZ
  
  N = IARGC ()
  if (N .ne. 2) then
     print *,'Error: cgns_flip must be called with TWO arguments:'
     print *,'./cgns_flip <cgns_file.cgns> <x,y, or z>'
     stop
  end if
  CALL GETARG(1 , filename)
  CALL GETARG(2 , char_flip)
  call cg_open_f(filename,CG_MODE_MODIFY, cg, ier)

  if (ier .eq. CG_ERROR) call cg_error_exit_f
  print *,'Processing file: ',filename

  call cg_nbases_f(cg, nbases, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  if (nbases .lt. 1) then
     print *, 'Error: No bases found in CGNS file'
     stop
  end if

  base = 1

  call cg_base_read_f(cg, base, basename, CellDim, PhysDim, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  if (cellDim .ne. 3 .or. PhysDim .ne. 3) then
     print *,'The Cells must be hexahedreal in 3 dimensions'
     stop
  end if

  ! Goto Base Node
  call cg_goto_f(cg, base, ier, 'end')
  call cg_nzones_f(cg, base, nzones, ier)

  do zone=1, nzones
     call cg_zone_read_f(cg, base, zone, zonename, zonesize, ier)
     call cg_zone_type_f(cg, base, zone, zonetype, ier)
     
     allocate(coorX(zonesize(1),zonesize(2),zonesize(3),3))
     allocate(coorY(zonesize(1),zonesize(2),zonesize(3),3))
     allocate(coorZ(zonesize(1),zonesize(2),zonesize(3),3))

     start(:) = 1
     call cg_coord_read_f(cg,base,zone,'CoordinateX',RealDouble,start,zonesize,coorX,ier)
     call cg_coord_read_f(cg,base,zone,'CoordinateY',RealDouble,start,zonesize,coorY,ier)
     call cg_coord_read_f(cg,base,zone,'CoordinateZ',RealDouble,start,zonesize,coorZ,ier)
     
     if (char_flip == 'x') then
        coorX = -coorX
     else if (char_flip == 'y') then
        coorY = -coorY
     else if (char_flip == 'z') then
        coorZ = -coorZ
     else
        print *,'Error: Flip was not understood. Use x,y or z.'
        stop
     end if

     call cg_coord_write_f(cg,base,zone,RealDouble,"CoordinateX",coorX,coordID,ier)
     call cg_coord_write_f(cg,base,zone,RealDouble,"CoordinateY",coorY,coordID,ier)
     call cg_coord_write_f(cg,base,zone,RealDouble,"CoordinateZ",coorZ,coordID,ier)

     deallocate(coorX,coorY,coorZ)

  end do
  call cg_close_f(cg, ier)
end program cgns_flip
