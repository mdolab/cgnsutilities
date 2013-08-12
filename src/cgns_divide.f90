program cgns_divide

  !use cgnsGrid

  implicit none
  include 'cgnslib_f.h'

  integer  CellDim, PhysDim
  integer nbases, nzones, in_dims(9),out_dims(9)

  integer ier, zonetype
  integer n,nn,mm,i,j,k,idim,ii,jj,kk
  integer istart,iend,jstart,jend,kstart,kend
  integer cg_in, cg_out, base, zone, zoneCounter
  integer tempInt

  integer :: ptSetType, normalDataType
  integer :: sizes(9),coordID,blockStart(3),blockEnd(3)

  integer donor_range(6),transform(3)
  integer nbocos,n1to1,bocotype,BCout
  integer NormalIndex(3), NormalListFlag, ndataset,datatype
  integer ptset_type, npnts, pnts(3,2),pnts_donor(3,2),ncon

  character*32 in_filename,out_filename
  character*32 basename, zonename,boconame
  character*32 connectname,donorname, famname
  double precision, allocatable,dimension(:,:,:) :: data3d
  double precision, allocatable,dimension(:,:,:) :: data3dtmp
  double precision data_double(6)

  N = IARGC ()
  if (N .ne. 2) then
     print *,'Error: cgns_divide must be called with TWO arguments:'
     print *,'./cgns_divide cgns_in_file.cgns cgns_out_file.cgns'
     stop
  end if

  call getarg(1, in_filename)
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

  zoneCounter = 0
  do nn=1,nzones
     call cg_zone_read_f(cg_in, base, nn, zonename, in_dims, ier)

     ! All blocks are the same size??
     do i=1,3
        out_dims(i) = (in_dims(i)-1)/2 + 1
     end do
     do i=4,6
        out_dims(i) = in_dims(i)/2
     end do

     allocate(data3d(in_dims(1),in_dims(2),in_dims(3)))

     blockStart = (/1,1,1/)
     blockEnd   = in_dims(1:3)

     ! Split the block into 8:
     do ii=1,2
        do jj=1,2
           do kk=1,2
              
              ! We need a unique zone name:
              zoneCounter = zoneCounter + 1
999           FORMAT('domain.',I5.5)
              write(zonename,999) zoneCounter

              call cg_zone_write_f(cg_out,base,zonename,out_dims,Structured,tempInt,ier)
              if (ier .eq. CG_ERROR) call cg_error_exit_f

              ! Determine the start and stop for each i,j,k index based on ii,jj,kk
              istart = (ii-1)*(out_dims(1)-1) + 1
              iend   = istart + out_dims(1)-1

              jstart = (jj-1)*(out_dims(2)-1) + 1
              jend   = jstart + out_dims(2)-1

              kstart = (kk-1)*(out_dims(3)-1) + 1
              kend   = kstart + out_dims(3)-1

              allocate(data3dTmp(istart:iend, jstart:jend, kstart:kend))
          
              ! X coordinate
              call cg_coord_read_f(cg_in,base,nn,'CoordinateX',RealDouble,&
                   blockStart,blockEnd,data3d,ier)
              if (ier .eq. CG_ERROR) call cg_error_exit_f
     
              ! Copy
              do k=kstart,kend
                 do j=jstart,jend
                    do i=istart,iend
                       data3dtmp(i,j,k) = data3d(i,j,k)
                    end do
                 end do
              end do

              call cg_coord_write_f(cg_out,base,zoneCounter,realDouble,&
                   'CoordinateX',data3dtmp, coordID,ier)
              
              ! Y coordinate
              call cg_coord_read_f(cg_in,base,nn,'CoordinateY',RealDouble,&
                   blockStart,blockEnd,data3d,ier)
              if (ier .eq. CG_ERROR) call cg_error_exit_f

              ! Copy
              do k=kstart,kend
                 do j=jstart,jend
                    do i=istart,iend
                       data3dtmp(i,j,k) = data3d(i,j,k)
                    end do
                 end do
              end do

              call cg_coord_write_f(cg_out,base,zoneCounter,realDouble,&
                   'CoordinateY',data3dtmp, coordID,ier)
              
              ! Z coordinate
              call cg_coord_read_f(cg_in,base,nn,'CoordinateZ',RealDouble,&
                   blockStart,blockEnd,data3d,ier)
              if (ier .eq. CG_ERROR) call cg_error_exit_f

              ! Copy
              do k=kstart,kend
                 do j=jstart,jend
                    do i=istart,iend
                       data3dtmp(i,j,k) = data3d(i,j,k)
                    end do
                 end do
              end do
     
              call cg_coord_write_f(cg_out,base,zoneCounter,realDouble,&
                   'CoordinateZ',data3dtmp, coordID,ier)

              deallocate(data3dTmp)

           end do
        end do
     end do
     deallocate(data3d)
  end do ! Zone Loop

  ! Close both CGNS files
  call cg_close_f(cg_in, ier)
  call cg_close_f(cg_out,ier)

end program cgns_divide
