! This file contains functions (to be wrapped with f2py) which are
! used with pyCreateCGNS to produce a cgns file that is actually
! useful

subroutine getNBlocks(f_in,N)
  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  character*(*), intent(in) :: f_in
  integer, intent(out) :: N

  ! Working
  integer :: ier,cg,base,nbases

  call cg_open_f(f_in, MODE_READ, cg, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  base = 1
  call cg_goto_f(cg, base, ier, 'end')
  if (ier .eq. CG_ERROR) call cg_error_exit_f
 
  call cg_nzones_f(cg, base, N, ier)
  if (ier .eq. ERROR) call cg_error_exit_f

  call cg_close_f(cg,ier)
  if (ier .eq. ERROR) call cg_error_exit_f

end subroutine getNBlocks

subroutine preprocess(f_in,f_out,N,block_dims,coords)

  ! This function takes f_in and copies all the grid data to f_out,
  ! and returns to python the block dimensions and the coordindates of
  ! the corners, midpoints of edges and midpoints of faces

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  character*(*), intent(in) :: f_in,f_out
  integer,intent(in) :: N
  integer, dimension(N,3), intent(out) :: block_dims
  double precision, dimension(N,26,3), intent(out) :: coords
  
  ! Working
  integer :: ier,cg_in,cg_out,base,nbases,nn,dims(9)
  character*32  zonename, basename
  double precision, allocatable, dimension(:,:,:,:) :: tempx
  integer  CellDim, PhysDim,blockStart(3),blockEnd(3)
  integer :: zoneCounter,coordID,il,jl,kl
  integer :: midu(2),midv(2),midw(2)

  ! Open and read zone sizes:
  call cg_open_f(f_in, MODE_READ, cg_in, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

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

  call cg_goto_f(cg_in, base, ier, 'end')
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  do nn=1,N
     call cg_zone_read_f(cg_in, base, nn, zonename, dims, ier)
     block_dims(nn,:) = dims(1:3)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
  end do

  ! Open Output File:
  call cg_open_f(f_out,CG_MODE_WRITE, cg_out, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_base_write_f(cg_out,"Base#1", Celldim,Physdim, base, ier)

  do nn=1,N
     call cg_zone_read_f(cg_in, base, nn, zonename, dims, ier)

999  FORMAT('domain.',I5.5)
     write(zonename,999) nn

     write(zonename,999) nn
     call cg_zone_write_f(cg_out,base,zonename,dims,Structured,zoneCounter,ier)
     
     allocate(tempx(block_dims(nn,1),block_dims(nn,2),block_dims(nn,3),3))
     blockStart = (/1,1,1/)
     blockEnd   = block_dims(nn,:)
     ! Read Grid Section
     call cg_coord_read_f(cg_in,base,nn,'CoordinateX',RealDouble,&
          blockStart,blockEnd,tempx(:,:,:,1),ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     call cg_coord_read_f(cg_in,base,nn,'CoordinateY',RealDouble,&
          blockStart,blockEnd,tempx(:,:,:,2),ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     call cg_coord_read_f(cg_in,base,nn,'CoordinateZ',RealDouble,&
          blockStart,blockEnd,tempx(:,:,:,3),ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     
     ! Write Grid Section
     call cg_coord_write_f(cg_out,base,zoneCounter,realDouble,&
          'CoordinateX',tempx(:,:,:,1), coordID,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     call cg_coord_write_f(cg_out,base,zoneCounter,realDouble,&
          'CoordinateY',tempx(:,:,:,2), coordID,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     call cg_coord_write_f(cg_out,base,zoneCounter,realDouble,&
          'CoordinateZ',tempx(:,:,:,3), coordID,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
     
     ! Get the required coordinates while we have tempx in memory:

     il = block_dims(nn,1)
     jl = block_dims(nn,2)
     kl = block_dims(nn,3)

     coords(nn,1,:) = tempx(1 ,1 ,1,:)
     coords(nn,2,:) = tempx(il,1 ,1,:)
     coords(nn,3,:) = tempx(1 ,jl,1,:)
     coords(nn,4,:) = tempx(il,jl,1,:)
     coords(nn,5,:) = tempx(1 ,1 ,kl,:)
     coords(nn,6,:) = tempx(il,1 ,kl,:)
     coords(nn,7,:) = tempx(1 ,jl,kl,:)
     coords(nn,8,:) = tempx(il,jl,kl,:)

     ! Mid points of edges

     if (mod(il,2) == 1) then
        midu(1) = (il-1)/2+1
        midu(2) = (il-1)/2+1
     else
        midu(1) = il/2+1
        midu(2) = il/2
     end if

     if (mod(jl,2) == 1) then
        midv(1) = (jl-1)/2+1
        midv(2) = (jl-1)/2+1
     else
        midv(1) = jl/2+1
        midv(2) = jl/2
     end if

     if (mod(kl,2) == 1) then
        midw(1) = (kl-1)/2+1
        midw(2) = (kl-1)/2+1
     else
        midw(1) = kl/2+1
        midw(2) = kl/2
     end if

     coords(nn,9 ,:) = 0.5*(tempx(midu(1) ,1 ,1,:) + tempx(midu(2) ,1 ,1,:))
     coords(nn,10,:) = 0.5*(tempx(midu(1) ,jl,1,:) + tempx(midu(2) ,jl,1,:))
     coords(nn,11,:) = 0.5*(tempx(1,midv(1),1,:)   + tempx(1,midv(2),1,:))
     coords(nn,12,:) = 0.5*(tempx(il,midv(1),1,:)  + tempx(il,midV(2),1,:))

     coords(nn,13,:) = 0.5*(tempx(midu(1) ,1 ,kl,:) + tempx(midu(2) ,1 ,kl,:))
     coords(nn,14,:) = 0.5*(tempx(midu(1) ,jl,kl,:) + tempx(midu(2) ,jl,kl,:))
     coords(nn,15,:) = 0.5*(tempx(1,midv(1),kl,:)   + tempx(1,midv(2),kl,:))
     coords(nn,16,:) = 0.5*(tempx(il,midv(1),kl,:)  + tempx(il,midv(2),kl,:))

     coords(nn,17,:) = 0.5*(tempx(1,1,midw(1),:) + tempx(1,1,midw(2),:))  
     coords(nn,18,:) = 0.5*(tempx(il,1,midw(1),:)+ tempx(il,1,midw(2),:))
     coords(nn,19,:) = 0.5*(tempx(1,jl,midw(1),:)+ tempx(1,jl,midw(2),:))
     coords(nn,20,:) = 0.5*(tempx(il,jl,midw(1),:)+tempx(il,jl,midw(2),:))

     ! Mid Points of faces

     coords(nn,21,:) = .25*(tempx(midu(1),midv(1),1,:) + &
          tempx(midu(2),midv(1),1,:) + &
          tempx(midu(1),midv(2),1,:) + &
          tempx(midu(2),midv(2),1,:))
     coords(nn,22,:) = .25*(tempx(midu(1),midv(1),kl,:) + &
          tempx(midu(2),midv(1),kl,:) + &
          tempx(midu(1),midv(2),kl,:) + &
          tempx(midu(2),midv(2),kl,:))
     coords(nn,23,:) = .25*(tempx(1,midv(1),midw(1),:) + &
          tempx(1,midv(2),midw(1),:) + &
          tempx(1,midv(1),midw(2),:) + &
          tempx(1,midv(2),midw(2),:))
     coords(nn,24,:) = .25*(tempx(il,midv(1),midw(1),:) + &
          tempx(il,midv(2),midw(1),:) + &
          tempx(il,midv(1),midw(2),:) + &
          tempx(il,midv(2),midw(2),:))
     coords(nn,25,:) = .25*(tempx(midu(1),1,midw(1),:) + &
          tempx(midu(2),1,midw(1),:) + &
          tempx(midu(1),1,midw(2),:) + &
          tempx(midu(2),1,midw(2),:))
     coords(nn,26,:) = .25*(tempx(midu(1),jl,midw(1),:) + &
          tempx(midu(2),jl,midw(1),:) + &
          tempx(midu(1),jl,midw(2),:) + &
          tempx(midu(2),jl,midw(2),:))

     deallocate(tempx)

  end do

  call cg_close_f(cg_in,ier)
  if (ier .eq. ERROR) call cg_error_exit_f

  call cg_close_f(cg_out,ier)
  if (ier .eq. ERROR) call cg_error_exit_f
 
end subroutine preprocess

subroutine openFile(f_out,cg_out)
  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  character*(*), intent(in) :: f_out
  integer, intent(out) :: cg_out

  ! Working 
  integer :: ier

  ! Open Output File:
  call cg_open_f(f_out,CG_MODE_MODIFY, cg_out, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

end subroutine openFile

subroutine writeBC(cg_out,nn,bcName,bcFam,pt_start,pt_end,bcType)

  ! This function takes f_in and copies all the grid data to f_out,
  ! and returns to python the block dimensions and the coordindates of
  ! the corners, midpoints of edges and midpoints of faces

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer,intent(in) :: cg_out
  integer,intent(in) :: nn
  integer,intent(in) :: pt_start(3)
  integer,intent(in) :: pt_end(3)
  integer,intent(in) :: bcType
  character*(*), intent(in) :: bcName
  character*(*), intent(in) :: bcFam

  ! Working
  integer :: ier,base,bcout,pnts(3,2),cgns_bc
  
  base =1
  pnts(:,1) = pt_start
  pnts(:,2) = pt_end

  if (bcType == 0) then
     cgns_bc = BCSymmetryPlane
  else if(bcType == 1) then
     cgns_bc = BCWallViscous
  else
     cgns_bc = BCFarfield
  end if
    
  call cg_boco_write_f(cg_out,base,nn+1,bcName, cgns_bc, PointRange, &
       2, pnts, BCout ,  ier )
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_goto_f(cg_out,base,ier,'Zone_t', nn+1,"ZoneBC_t", 1,&
       "BC_t", BCOut, "end")
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_famname_write_f(bcFam, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

end subroutine writeBC

subroutine writeB2B(cg_out,nn,conName,donorName,pt_start,pt_end,pt_start_donor,&
     pt_end_donor,transform)

  ! This function takes f_in and copies all the grid data to f_out,
  ! and returns to python the block dimensions and the coordindates of
  ! the corners, midpoints of edges and midpoints of faces

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output
  integer,intent(in) :: cg_out
  integer,intent(in) :: nn
  character*(*), intent(in) :: conName
  character*(*), intent(in) :: donorName
  integer,intent(in) :: pt_start(3)
  integer,intent(in) :: pt_end(3)
  integer,intent(in) :: pt_start_donor(3)
  integer,intent(in) :: pt_end_donor(3)
  integer,intent(in) :: transform(3)
  
  ! Working
  integer :: ier,base,pnts(3,2),pnts_donor(3,2),nCon
  
  base = 1

  pnts(:,1) = pt_start
  pnts(:,2) = pt_end

  pnts_donor(:,1) = pt_start_donor
  pnts_donor(:,2) = pt_end_donor

  call cg_1to1_write_f(cg_out,base,nn+1,conName,donorName,pnts,pnts_donor, &
       transform,nCon,ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  
end subroutine writeB2B

subroutine closeFile(cg_out)
  implicit none
  include 'cgnslib_f.h'
  
  ! Input/Output
  integer, intent(in) :: cg_out

  ! Working
  integer :: ier

  ! Close Output File:
  
  call cg_close_f(cg_out,ier)
  if (ier .eq. ERROR) call cg_error_exit_f

end subroutine closeFile
