# =============================================================================
# Standard Python modules
# =============================================================================

import os, sys, string, copy, pdb, time

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================

from mdo_import_helper import *
exec(import_modules('geo_utils'))
import cgns_create


# This file creates a proper CGNS file form an unconnected CGNS file
# with only coordinates

N = len(sys.argv)

in_file = sys.argv[1]
out_file = sys.argv[2]
sym = sys.argv[3]
R = float(sys.argv[4])

try:
    family_file = sys.argv[5]
except:
    family_file = None
# end if


offset = [0,0,0]
# Preprocess:
timeA = time.time()
print 'Copying Grid and Preprocessing...'
nVol = cgns_create.getnblocks(in_file)
block_dims,coords = cgns_create.preprocess(in_file,out_file,nVol)

# Create the Topology:
print 'Creating Topology...'
topo = geo_utils.BlockTopology(coords)

class BCInfo(object):
    
    def __init__(self,ID,bcType,bcFam,pt_start,pt_end):
        self.bcName = 'BC_on_SF%d'%(ID)
        self.bcFam = bcFam
        self.pt_start = pt_start
        self.pt_end = pt_end
        self.bcType = bcType
        
        return
    
class B2BInfo(object):

    def __init__(self,ID,Vol1,Vol2,pt_start,pt_end,
                 pt_start_donor,pt_end_donor,transform):
        # Note the +1 is for the 1 based ordering in the cgns file
        self.conName = 'SF%d (domains %d and %d)'%(ID,Vol1+1,Vol2+1)
        self.donorName = 'domain.%5.5d'%(Vol2+1)
        self.pt_start = pt_start
        self.pt_end = pt_end
        self.pt_start_donor = pt_start_donor
        self.pt_end_donor = pt_end_donor
        self.transform = transform
        
        return

def generalized_coord_dir(iFace):
    if iFace in [0,1]:
        return [0,1,2]
    elif iFace in [2,3]:
        return [1,2,0]
    elif iFace in [4,5]:
        return [0,2,1]
    # end if

def isodd(num):
        return num & 1 and True or False

def point_range(iFace,blk_dim):
    # Return the correct point range for face iFace on a block with
    # dimensions given in blk_dim:
    il = blk_dim[0]
    jl = blk_dim[1]
    kl = blk_dim[2]
    if iFace == 0:
        pt_start = [1,1,1]
        pt_end = [il,jl,1]
    elif iFace == 1:
        pt_start = [1,1,kl]
        pt_end = [il,jl,kl]
    elif iFace == 2:
        pt_start = [1,1,1]
        pt_end = [1,jl,kl]
    elif iFace == 3:
        pt_start = [il,1,1]
        pt_end = [il,jl,kl]
    elif iFace == 4:
        pt_start = [1,1,1]
        pt_end = [il,1,kl]
    elif iFace == 5:
        pt_start = [1,jl,1]
        pt_end = [il,jl,kl]
    # end if

    return pt_start,pt_end

def normal_direction(iFace1,iFace2):
    # Normal direction is positive if iFace1 and iFace two are of
    # opposite odditiy, even if they are the same oddity
    isOdd1 = isodd(iFace1)
    isOdd2 = isodd(iFace2)

    if isOdd1==True and isOdd2==True:
        return -1
    if isOdd1==False and isOdd2==False:
        return -1

    # otherwise:
    return 1

# Generate symmetry data
sym_normal = [0,0,0]
if sym == 'x':
    sym_normal_index = 0
elif sym == 'y':
    sym_normal_index = 1
elif sym == 'z':
    sym_normal_index = 2
# end if

sym_normal[sym_normal_index] = 1

# Generate the face_ptr: This is an array of length unique faces, that
# contains the (iVol,iFace) that this unique faces points to. This is
# used to determine if a face is a boundary condition, (only 1 face)
# or a block-to-block connection (2 faces)
face_ptr = []
for iFace in xrange(topo.nFace):
    face_ptr.append([])

for iVol in xrange(nVol):
    for iFace in xrange(6):
        iUniqueFace = topo.face_link[iVol,iFace]
        face_ptr[iUniqueFace].append([iVol,iFace])
    # end for
# end for

# Now we have to go though and fill up the faceInfo
face_info = numpy.zeros((topo.nVol,6),numpy.object)

for iUFace in xrange(topo.nFace):
    nFace_connected = len(face_ptr[iUFace])
    
    if nFace_connected == 1: # Boudnary Condition

        # The is only one face here
        iVol = face_ptr[iUFace][0][0]
        iFace = face_ptr[iUFace][0][1]

        # By default we don't know what the family is
        bcType = None
        bcFam = None

        # First check for a symmetry condition:

        # -> Get the coordinates for this face:
        nodes = geo_utils.nodesFromFace(iFace)

        # -> Extract the coordiantes for this face from coords:
        face_coords = numpy.take(coords[iVol],nodes,axis=0)

        # -> Take the cross-product to get the normal:
        v1 = face_coords[3]-face_coords[0]
        v2 = face_coords[2]-face_coords[1]
        normal = numpy.cross(v1,v2)

        # -> Normalize with a zero-lenght check:
        normal /= max(numpy.linalg.norm(normal),1e-14)

        coor_check = abs(numpy.average(face_coords[:,sym_normal_index]))<1e-3
        dp_check = numpy.dot(normal,sym_normal) > 0.98 or numpy.dot(normal,sym_normal) < -0.98
      
        if dp_check and coor_check:
            bcType = 0
            bcFam = 'sym'
        # end if

        # Next check for a wall-type boundary condition:

        inside_sphere = True
        for i in xrange(4):
            if not numpy.linalg.norm(face_coords[i]-offset) < R:
                inside_sphere = False
            # end if
        # end for

        if inside_sphere and bcType is None:
            bcType = 1 # BCWall Viscous
            bcFam = 'wall'
        # end if

        # Finally, the default boudnary condition if either of
        # these two conditions are met, is the farfield condition

        if bcType is None:
            bcType = 2 # Farfield Condition
            bcFam = 'far'
        # end if

        # Get the starting and ending point ranges
        pt_start,pt_end = point_range(iFace,block_dims[iVol])

        # Add to faceInfo:
        face_info[iVol,iFace] = BCInfo(iUFace,bcType,bcFam,pt_start,pt_end)

    elif nFace_connected == 2: # Block to Block

        # The is two faces:
        iVol1  = face_ptr[iUFace][0][0]
        iFace1 = face_ptr[iUFace][0][1]

        iVol2  = face_ptr[iUFace][1][0]
        iFace2 = face_ptr[iUFace][1][1]

        # Get the pt ranges:
        pt_start,pt_end = point_range(iFace1,block_dims[iVol1])
        pt_start_donor,pt_end_donor = point_range(iFace2,block_dims[iVol2])

        # Now for the fun part...the transformation matrix:
        transform1 = [0,0,0]

        g_dir_1 = generalized_coord_dir(iFace1)
        g_dir_2 = generalized_coord_dir(iFace2)

        fl = topo.face_dir[iVol2,iFace2]

        # Do the generalized coordinate direction 1

        if fl in [0,2]:
            transform1[g_dir_1[0]] =  (g_dir_2[0] + 1)
        if fl in [1,3]:
            transform1[g_dir_1[0]] = -(g_dir_2[0] + 1)
        if fl in [4,6]:
            transform1[g_dir_1[0]] =  (g_dir_2[1] + 1)
        if fl in [5,7]:
            transform1[g_dir_1[0]] = -(g_dir_2[1] + 1)

        # Do the generalized coordinate direction 2
        if fl in [0,1]:
            transform1[g_dir_1[1]] =  (g_dir_2[1] + 1)
        if fl in [2,3]:
            transform1[g_dir_1[1]] = -(g_dir_2[1] + 1)
        if fl in [4,5]:
            transform1[g_dir_1[1]] =  (g_dir_2[0] + 1)
        if fl in [6,7]:
            transform1[g_dir_1[1]] = -(g_dir_2[0] + 1)
        # end if

        # Do the generalized coordiante direction 3 (or normal)

        transform1[g_dir_1[2]] = normal_direction(iFace1,iFace2)*(g_dir_2[2]+1)
        
        # Now that we have transform1 we can fairly easily get transform2:
        transform2 = [0,0,0]
        for ii in xrange(3):
            ind = abs(transform1[ii])-1
            transform2[ind] = numpy.sign(transform1[ii])*(ii+1)

        # One last thing to do: if the transform is -ve swap the
        # ranges on the donor
        for ii in xrange(3):
            if transform2[ii] < 0:
                temp = pt_start_donor[ii]
                pt_start_donor[ii] = pt_end_donor[ii]
                pt_end_donor[ii] = temp
            # end if
        # end for

        # Finally add the two face connectivities:

        face_info[iVol1,iFace1] = B2BInfo(iUFace,iVol1,iVol2,pt_start,pt_end,
                                          pt_start_donor,pt_end_donor,
                                          transform1)

        face_info[iVol2,iFace2] = B2BInfo(iUFace,iVol2,iVol1,pt_start_donor,
                                          pt_end_donor,pt_start,pt_end,
                                          transform2)
    else:
        print 'Error: More than two faces ocnnected to a unique face. This results in a non-physical topology.'
        print nFace_connected 
        sys.exit(-1)
    # end if
# end for

# Finally write out the BC's to the cgns file
cg_out = cgns_create.openfile(out_file)
print 'Writing Boundary and Connectivity Info...'
for iVol in xrange(topo.nVol):
    for iFace in xrange(6):
        if isinstance(face_info[iVol,iFace],BCInfo):
            cgns_create.writebc(cg_out,iVol,
                                face_info[iVol,iFace].bcName,
                                face_info[iVol,iFace].bcFam,
                                face_info[iVol,iFace].pt_start,
                                face_info[iVol,iFace].pt_end,
                                face_info[iVol,iFace].bcType)
        # end if
        elif isinstance(face_info[iVol,iFace],B2BInfo):
            cgns_create.writeb2b(cg_out,iVol,
                                face_info[iVol,iFace].conName,
                                face_info[iVol,iFace].donorName,
                                face_info[iVol,iFace].pt_start,
                                face_info[iVol,iFace].pt_end,
                                face_info[iVol,iFace].pt_start_donor,
                                face_info[iVol,iFace].pt_end_donor,
                                face_info[iVol,iFace].transform)
        # end if
    # end for
# end for

# Close CGNS File
cgns_create.closefile(cg_out)

print 'Done! Time was %5.3f seconds'%(time.time()-timeA)
