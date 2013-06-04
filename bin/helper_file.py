# =============================================================================
# Utility Functions for Use in pyCGNSCreate
# =============================================================================

import numpy as np

def e_dist(x1, x2):
    '''Get the eculidean distance between two points'''
    if len(x1) == 3:
        return np.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2)
    elif len(x1) == 2:
        return np.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2)
    elif len(x1) == 1:
        return np.abs(x1[0]-x2[0])

def unique_index(s, s_hash=None):
    '''
    This function is based on unique

    The idea is to take a list s, and reduce it as per unique.

    However, it additionally calculates a linking index arrary that is
    the same size as the original s, and points to where it ends up in
    the the reduced list

    if s_hash is not specified for sorting, s is used

    '''
    if s_hash != None:
        ind = np.argsort(np.argsort(s_hash))
    else:
        ind = np.argsort(np.argsort(s))
    # end if
    n = len(s)
    t = list(s)
    t.sort()
    
    diff = np.zeros(n, 'bool')

    last = t[0]
    lasti = i = 1
    while i < n:
        if t[i] != last:
            t[lasti] = last = t[i]
            lasti += 1
        else:
            diff[i] = True
        # end if
        i += 1
    # end while
    b = np.where(diff == True)[0]
    for i in xrange(n):
        ind[i] -= b.searchsorted(ind[i], side='right')
    # end for

    return t[:lasti], ind

def pointReduce(points, node_tol=1e-4):
    '''Given a list of N points in ndim space, with possible
    duplicates, return a list of the unique points AND a pointer list
    for the original points to the reduced set'''

    # First 
    points = np.array(points)
    N = len(points)
    dists = []
    for ipt in xrange(N): 
        dists.append(np.sqrt(np.dot(points[ipt], points[ipt])))
    # end for
    temp = np.array(dists)
    temp.sort()
    ind = np.argsort(dists)
    i = 0
    cont = True
    new_points = []
    link = np.zeros(N, 'intc')
    link_counter = 0
   
    while cont:
        cont2 = True
        temp_ind = []
        j = i
        while cont2:
            if abs(dists[ind[i]]-dists[ind[j]])<node_tol:
                temp_ind.append(ind[j])
                j = j + 1
                if j == N: # Overrun check
                    cont2 = False
                # end if
            else:
                cont2 = False
            #end if
        # end while
        sub_points = [] # Copy of the list of sub points with the dists
        for ii in xrange(len(temp_ind)):
            sub_points.append(points[temp_ind[ii]])

        # Brute Force Search them 
        sub_unique_pts, sub_link = pointReduceBruteForce(sub_points, node_tol)
        new_points.extend(sub_unique_pts)

        for ii in xrange(len(temp_ind)):
            link[temp_ind[ii]] = sub_link[ii] + link_counter
        # end if
        link_counter += max(sub_link) + 1

        i = j - 1 + 1
        if i == N:
            cont = False
        # end if
    # end while
    return np.array(new_points), np.array(link)

def pointReduceBruteForce(points,  node_tol=1e-4):
    '''Given a list of N points in ndim space, with possible
    duplicates, return a list of the unique points AND a pointer list
    for the original points to the reduced set

    BRUTE FORCE VERSION

    '''
    N = len(points)
    unique_points = [points[0]]
    link = [0]
    for i in xrange(1, N):
        found_it = False
        for j in xrange(len(unique_points)):
            if e_dist(points[i], unique_points[j]) < node_tol:
                link.append(j)
                found_it = True
                break
            # end if
        # end for
        if not found_it:
            unique_points.append(points[i])
            link.append(j+1)
        # end if
    # end for
    return np.array(unique_points), np.array(link)

def edgeOrientation(e1, e2):
    '''Compare two edge orientations. Basically if the two nodes are
    in the same order return 1 if they are in opposite order, return
    1'''
    
    if [e1[0], e1[1]] == [e2[0], e2[1]]:
        return 1
    elif [e1[1], e1[0]] == [e2[0], e2[1]]:
        return -1
    else:
        mpiPrint('Error with edgeOrientation: Not possible.')
        mpiPrint('Orientation 1 [%d %d]'%(e1[0], e1[1]))
        mpiPrint('Orientation 2 [%d %d]'%(e2[0], e2[1]))
        sys.exit(0)
    # end if

def faceOrientation(f1, f2):
    '''Compare two face orientations f1 and f2 and return the
    transform to get f1 back to f2'''
    
    if [f1[0], f1[1], f1[2], f1[3]] == [f2[0], f2[1], f2[2], f2[3]]:
        return 0
    elif [f1[0], f1[1], f1[2], f1[3]] == [f2[1], f2[0], f2[3], f2[2]]:
        return 1
    elif [f1[0], f1[1], f1[2], f1[3]] == [f2[2], f2[3], f2[0], f2[1]]:
        return 2
    elif [f1[0], f1[1], f1[2], f1[3]] == [f2[3], f2[2], f2[1], f2[0]]:
        return 3
    elif [f1[0], f1[1], f1[2], f1[3]] == [f2[0], f2[2], f2[1], f2[3]]:
        return 4
    elif [f1[0], f1[1], f1[2], f1[3]] == [f2[2], f2[0], f2[3], f2[1]]:
        return 5
    elif [f1[0], f1[1], f1[2], f1[3]] == [f2[1], f2[3], f2[0], f2[2]]:
        return 6
    elif [f1[0], f1[1], f1[2], f1[3]] == [f2[3], f2[1], f2[2], f2[0]]:
        return 7
    else:
        mpiPrint('Error with faceOrientation: Not possible.')
        mpiPrint('Orientation 1 [%d %d %d %d]'%(f1[0], f1[1], f1[2], f1[3]))
        mpiPrint('Orientation 2 [%d %d %d %d]'%(f2[0], f2[1], f2[2], f2[3]))
        sys.exit(0)

# --------------------------------------------------------------
#                     Node/Edge Functions
# --------------------------------------------------------------

def nodesFromEdge(edge):
    '''Return the nodes on either edge of a standard edge'''
    if edge == 0:
        return 0, 1
    elif edge == 1:
        return 2, 3
    elif edge == 2:
        return 0, 2
    elif edge == 3:
        return 1, 3
    elif edge == 4:
        return 4, 5
    elif edge == 5:
        return 6, 7
    elif edge == 6:
        return 4, 6
    elif edge == 7:
        return 5, 7
    elif edge == 8:
        return 0, 4
    elif edge == 9:
        return 1, 5
    elif edge == 10:
        return 2, 6
    elif edge == 11:
        return 3, 7

def nodesFromFace(face):
    if face == 0:
        return [0, 1, 2, 3]
    elif face == 1:
        return [4, 5, 6, 7]
    elif face == 2:
        return [0, 2, 4, 6]
    elif face == 3:
        return [1, 3, 5, 7]
    elif face == 4:
        return [0, 1, 4, 5]
    elif face == 5:
        return [2, 3, 6, 7]

class topology(object):
    '''
    The base topology class from which the BlockTopology,
    SurfaceTology and CuveTopology classes inherit from
    
    The topology object contains all the info required for the block
    topology (most complex) however, simpiler topologies are handled
    accordingly.

    Class Attributes:
        nVol : The number of volumes in the topology (may be 0)
        nFace: The number of unique faces on the topology (may be 0)
        nEdge: The number of uniuqe edges on the topology 
        nNode: The number of unique nodes on the topology

        nEnt: The number of "entities" in the topology class. This may
        be curves, faces or volumes

        mNodeEnt: The number of NODES per entity. For curves it's 2, for
        surfaces 4 and for volumes 8.

        mEdgeEnt: The number of EDGES per entity. For curves it's 1,
        for surfaces, 4 and for volumes, 12

        mFaceEnt: The number of faces per entity. For curves its's 0,
        for surfaces, 1 and for volumes,6

        mVolEnt: The number of volumes per entity. For curves it's 0,
        for surfaces, 0 and for volumnes, 1

        node_link: The array of size nEnt x mNodesEnt which points
                   to the node for each entity
        edge_link: The array of size nEnt x mEdgeEnt which points
                   to the edge for each edge of entity
        face_link: The array of size nEnt x mFaceEnt which points to 
                   the face of each face on an entity

        edge_dir:  The array of size nEnt x mEdgeEnt which detrmines
                   if the intrinsic direction of this edge is
                   opposite of the direction as recorded in the
                   edge list. edge_dir[entity#][#] = 1 means same direction;
                   -1 is opposite direction.
                  
        face_dir:  The array of size nFace x 6 which determines the 
                   intrinsic direction of this face. It is one of 0->7
                   
        l_index:   The local->global list of arrays for each volue
        g_index:   The global->local list points for the entire topology
        edges:     The list of edge objects defining the topology
        simple    : A flag to determine of this is a "simple" topology 
                   which means there are NO degernate Edges, 
                   NO multiple edges sharing the same nodes and NO 
                   edges which loop back and have the same nodes
                   MUST BE SIMPLE
    '''

    def __init__(self):
        # Not sure what should go here...
        return
    def _calcDGs(self, edges, edge_link, edge_link_sorted, edge_link_ind):

        dg_counter = -1
        for i in xrange(self.nEdge):
            if edges[i][2] == -1: # Not set yet
                dg_counter += 1
                edges[i][2] = dg_counter
                self._addDGEdge(i, edges, edge_link, 
                                edge_link_sorted, edge_link_ind)
            # end if
        # end for
        self.nDG = dg_counter + 1
   
    def _addDGEdge(self, i, edges, edge_link, edge_link_sorted, edge_link_ind):
        left  = edge_link_sorted.searchsorted(i, side='left')
        right = edge_link_sorted.searchsorted(i, side='right')
        res   = edge_link_ind[slice(left, right)]

        for j in xrange(len(res)):
            ient = res[j]/self.mEdgeEnt #Integer Division
            iedge = np.mod(res[j], self.mEdgeEnt)

            pEdges = self._getParallelEdges(iedge)
            oppositeEdges = []
            for iii in xrange(len(pEdges)):
                oppositeEdges.append(
                    edge_link[self.mEdgeEnt*ient + pEdges[iii]])
            
            for ii in xrange(len(pEdges)):
                if edges[oppositeEdges[ii]][2] == -1:
                    edges[oppositeEdges[ii]][2] = edges[i][2]
                    if not edges[oppositeEdges[ii]][0] == \
                            edges[oppositeEdges[ii]][1]:
                        self._addDGEdge(oppositeEdges[ii], edges, 
                                        edge_link, edge_link_sorted, 
                                        edge_link_ind)
                # end if
            # end if
        # end for

    def _getParallelEdges(self, iedge):
        '''Return parallel edges for surfaces and volumes'''

        if self.topo_type == 'surface':
            if iedge == 0: return [1]
            if iedge == 1: return [0]
            if iedge == 2: return [3]
            if iedge == 3: return [2]

        if self.topo_type == 'volume':
            if iedge == 0: 
                return [1, 4, 5]
            if iedge == 1: 
                return [0, 4, 5]
            if iedge == 2: 
                return [3, 6, 7]
            if iedge == 3: 
                return [2, 6, 7]
            if iedge == 4: 
                return [0, 1, 5]
            if iedge == 5: 
                return [0, 1, 4]
            if iedge == 6: 
                return [2, 3, 7]
            if iedge == 7: 
                return [2, 3, 6]
            if iedge == 8: 
                return [9, 10, 11]
            if iedge == 9: 
                return [8, 10, 11]
            if iedge == 10: 
                return [8, 9, 11]
            if iedge == 11: 
                return [8, 9, 10]
        if self.topo_type == 'curve':
            return None

    def _getDGList(self):
        '''After calcGlobalNumbering is called with the size
        parameters, we can now produce a list of length ndg with the
        each entry coorsponing to the number N associated with that DG'''

        # This can be run in linear time...just loop over each edge
        # and add to dg list
        N_list = np.zeros(self.nDG, 'intc')
        for iedge in xrange(self.nEdge):
            N_list[self.edges[iedge].dg] = self.edges[iedge].N
        # end for
            
        return N_list

class BlockTopology(topology):
    '''
    See Topology base class for more information
    '''

    def __init__(self, coords=None, node_tol=1e-4, edge_tol=1e-4, file=None):
        '''Initialize the class with data required to compute the topology'''
        
        topology.__init__(self)
        self.mNodeEnt = 8
        self.mEdgeEnt = 12
        self.mFaceEnt = 6
        self.mVolEnt  = 1
        self.topo_type = 'volume'
        self.g_index = None
        self.l_index = None
        self.nGlobal = None
        if file != None:
            self.readConnectivity(file)
            return
        # end if

        coords = np.atleast_2d(coords)
        nVol = len(coords)

        if coords.shape[1] == 8: # Just the corners are given --- Just
                                 # put in np.zeros for the edge and face
                                 # mid points
            temp = np.zeros((nvol, (8 + 12 + 6), 3))
            temp[:, 0:8, :] = coords
            coords = temp.copy()
        # end if

        # ----------------------------------------------------------
        #                     Unique Nodes
        # ----------------------------------------------------------

        # Do the pointReduce Agorithm on the corners
        un, node_link = pointReduce(coords[:, 0:8, :].reshape((nVol*8, 3)))
        node_link = node_link.reshape((nVol, 8))
         
        # ----------------------------------------------------------
        #                     Unique Edges
        # ----------------------------------------------------------
 
        # Now determine the unique edges:
        edge_objs = []
        orig_edges = []
        for ivol in xrange(nVol):
            for iedge in xrange(12):
                # Node number on volume
                n1, n2 = nodesFromEdge(iedge)

                # Actual Global Node Number
                n1 = node_link[ivol][n1]
                n2 = node_link[ivol][n2]

                # Midpoint
                midpoint = coords[ivol][iedge + 8]

                # Sorted Nodes:
                ns = sorted([n1, n2])
        
                # Append the new edge_cmp Object
                edge_objs.append(edge_cmp_object(
                        ns[0], ns[1], n1, n2, midpoint, edge_tol))

                # Keep track of original edge orientation---needed for
                # face direction
                orig_edges.append([n1, n2])
            # end for
        # end for
                    
        # Generate unique set of edges
        unique_edge_objs,  edge_link = unique_index(edge_objs)

        edge_dir = []
        for i in xrange(len(edge_objs)): # This is nVol * 12
            edge_dir.append(edgeOrientation(
                    orig_edges[i], unique_edge_objs[edge_link[i]].nodes))
        # end for

        # ----------------------------------------------------------
        #                     Unique Faces
        # ----------------------------------------------------------

        face_objs = []
        orig_faces = []
        for ivol in xrange(nVol):
            for iface in xrange(6):
                # Node number on volume
                n1, n2, n3, n4 = nodesFromFace(iface)

                # Actual Global Node Number
                n1 = node_link[ivol][n1]
                n2 = node_link[ivol][n2] 
                n3 = node_link[ivol][n3]
                n4 = node_link[ivol][n4] 
                
                # Midpoint --> May be [0, 0, 0] -> This is OK
                midpoint = coords[ivol][iface + 8 + 12]
                
                # Sort the nodes before they go into the faceObject
                ns = sorted([n1, n2, n3, n4])
                face_objs.append(face_cmp_object(ns[0], ns[1], ns[2], ns[3], 
                                                 n1, n2, n3, n4, 
                                                 midpoint, 1e-4))
                # Keep track of original face orientation---needed for
                # face direction
                orig_faces.append([n1, n2, n3, n4])
            # end for
        # end for
                    
        # Generate unique set of faces
        unique_face_objs, face_link = unique_index(face_objs)

        face_dir = []
        for i in xrange(len(face_objs)): # This is nVol * 12
            face_dir.append(faceOrientation(
                    unique_face_objs[face_link[i]].nodes, orig_faces[i]))
            uEdge = face_link[i]
        # end for

        # --------- Set the Requried Data for this class ------------
        self.nNode = len(un)
        self.nEdge = len(unique_edge_objs)
        self.nFace = len(unique_face_objs)
        self.nVol  = len(coords)
        self.nEnt  = self.nVol

        self.node_link = node_link
        self.edge_link = np.array(edge_link).reshape((nVol, 12))
        self.face_link = np.array(face_link).reshape((nVol, 6))

        self.edge_dir  = np.array(edge_dir).reshape((nVol, 12))
        self.face_dir  = np.array(face_dir).reshape((nVol, 6))

        # Next Calculate the Design Group Information
        edge_link_sorted = np.sort(edge_link.flatten())
        edge_link_ind    = np.argsort(edge_link.flatten())

        ue = []
        for i in xrange(len(unique_edge_objs)):
            ue.append([unique_edge_objs[i].nodes[0], 
                       unique_edge_objs[i].nodes[1], -1, 0, 0])
        # end for

        self._calcDGs(ue, edge_link, edge_link_sorted, edge_link_ind)
        
        # Set the edge ojects
        self.edges = []
        for i in xrange(self.nEdge): # Create the edge objects
            self.edges.append(edge(
                    ue[i][0], ue[i][1], 0, 0, 0, ue[i][2], ue[i][3]))
        # end for

        return

class edge(object):
    '''A class for edge objects'''

    def __init__(self, n1, n2, cont, degen, intersect, dg, N):
        self.n1        = n1        # Integer for node 1
        self.n2        = n2        # Integer for node 2
        self.cont      = cont      # Integer: 0 for c0 continuity, 1
                                   # for c1 continuity
        self.degen     = degen     # Integer: 1 for degenerate, 0 otherwise
        self.intersect = intersect # Integer: 1 for an intersected
                                   # edge, 0 otherwise
        self.dg        = dg        # Design Group index
        self.N         = N         # Number of control points for this edge
        
    def write_info(self, i, handle):
        handle.write('  %5d        | %5d | %5d | %5d | %5d | %5d |\
  %5d |  %5d |\n'\
                     %(i, self.n1, self.n2, self.cont, self.degen, 
                       self.intersect, self.dg, self.N))

class edge_cmp_object(object):
    '''A temporary class for sorting edge objects'''

    def __init__(self, n1, n2, n1o, n2o, mid_pt, tol):
        self.n1 = n1
        self.n2 = n2
        self.nodes = [n1o, n2o]
        self.mid_pt = mid_pt
        self.tol = tol

    def __repr__(self):
        return 'Node1: %d Node2: %d Mid_pt: %f %f %f'% (
            self.n1, self.n2, self.mid_pt[0], self.mid_pt[1], self.mid_pt[2])

    def __cmp__(self, other):
        # This function should return :
        # -1 if self < other
        #  0 if self == other
        #  1 if self > other

        # Basically we want to make three comparisons: n1, n2 and the
        # mid_pt Its (substantially) faster if we break before all 3
        # comparisons are done
        
        n1_cmp = cmp(self.n1, other.n1)
        if n1_cmp: # n1_cmp is non-zero so return with the result
            return n1_cmp

        n2_cmp = cmp(self.n2, other.n2)

        if n2_cmp: # n2_cmp is non-zero so return 
            return n2_cmp

        x_cmp = cmp(self.mid_pt[0], other .mid_pt[0])
        y_cmp = cmp(self.mid_pt[1], other .mid_pt[1])
        z_cmp = cmp(self.mid_pt[2], other .mid_pt[2])
        
        if e_dist(self.mid_pt, other.mid_pt) < self.tol:
            mid_cmp = 0
        else:
            mid_cmp = x_cmp or y_cmp or z_cmp
        # end if

        return mid_cmp

class face_cmp_object(object):
    '''A temporary class for sorting edge objects'''

    def __init__(self, n1, n2, n3, n4, n1o, n2o, n3o, n4o, mid_pt, tol):
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
        self.n4 = n4
        self.nodes = [n1o, n2o, n3o, n4o]
        self.mid_pt = mid_pt
        self.tol = tol

    def __repr__(self):
        return 'Node1: %d Node2: %d Mid_pt: %f %f %f'% (
            self.n1, self.n2, self.mid_pt[0], self.mid_pt[1], self.mid_pt[2])

    def __cmp__(self, other):
        # This function should return :
        # -1 if self < other
        #  0 if self == other
        #  1 if self > other

        # Basically we want to make three comparisons: n1, n2, n3, n4 and the
        # mid_pt Its (substantially) faster if we break before all 
        # comparisons are done
        
        n1_cmp = cmp(self.n1, other.n1)
        if n1_cmp: # n1_cmp is non-zero so return with the result
            return n1_cmp

        n2_cmp = cmp(self.n2, other.n2)

        if n2_cmp: # n2_cmp is non-zero so return 
            return n2_cmp

        n3_cmp = cmp(self.n3, other.n3)
        if n3_cmp: # n3_cmp is non-zero so return
            return n3_cmp

        n4_cmp = cmp(self.n4, other.n4)
        if n4_cmp: # n4_cmp is non-zero so return
            return n4_cmp

        # Finally do mid-pt calc
        x_cmp = cmp(self.mid_pt[0], other .mid_pt[0])
        y_cmp = cmp(self.mid_pt[1], other .mid_pt[1])
        z_cmp = cmp(self.mid_pt[2], other .mid_pt[2])
        
        if e_dist(self.mid_pt, other.mid_pt) < self.tol:
            mid_cmp = 0
        else:
            mid_cmp = x_cmp or y_cmp or z_cmp
        # end if

        return mid_cmp
