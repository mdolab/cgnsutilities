"""
This is the new gateway program to all of the cgns_utils.

Run cgns_utils -help to get a list of all available options. The basic
idea is as follows:

                                              | write new file
read cngs file -> Do some operations on it -> |     .or.
                                              | write modified file
Developed by Dr. Gaetan K. W. Kenway
"""
import os
import copy
import tempfile
import numpy
from . import libcgns_utils

# These are taken from the CGNS include file (cgnslib_f.h in your cgns library folder)
BC = {
    "bcfarfield": 7,
    "bcsymmetryplane": 16,
    "bcwall": 20,
    "bcwallinviscid": 21,
    "bcwallviscous": 22,
    "bcwallviscousheatflux": 23,
    "bcwallviscousisothermal": 24,
    "bcoutflow": 13,
    "bcoutflowsubsonic": 14,
    "bcoutflowsupersonic": 15,
    "bcinflow": 9,
    "bcinflowsubsonic": 10,
    "bcinflowssupersonic": 11,
    "bcoverset": 1,
}  # The Overset BC will be considered as a CG_USERDEFINED option ()

BCDATATYPE = {"Dirichlet": 2, "Neumann": 3}

CGNSDATATYPES = {"Integer": 2, "RealSingle": 3, "RealDouble": 4, "Character": 5, "LongInteger": 6}

CG_MODE_READ = 0
CG_MODE_WRITE = 1


class Grid(object):
    """Represent a complete 3D multiblock grid"""

    def __init__(self):
        self.blocks = []
        self.convArray = {}
        self.topo = None
        self.name = "domain"
        self.cellDim = 3

    def getTotalCellsNodes(self):
        """
        Returns the total number of Cells and Nodes in the grid.

        Example
        -------
        To determine the number of cells/nodes in a grid use:
            >>> from cgnsutilities.cgnsutilities import Grid, readGrid
            >>> grid = readGrid("gridfilename.cgns")
            >>> totalCells, totalNodes = grid.getTotalCellsNodes()
        """
        totalCells = 0
        totalNodes = 0
        for blk in self.blocks:
            totalCells += blk.getNumCells()
            totalNodes += blk.getNumNodes()

        return totalCells, totalNodes

    def getWallCellsNodes(self):
        """
        Returns the number of Cells and Nodes on wall boundaries

        Example
        -------
        To determine the number of cells/nodes in a grid use:
            >>> from cgnsutilities.cgnsutilities import Grid, readGrid
            >>> grid = readGrid("gridfilename.cgns")
            >>> nWallCells, nWallNodes = grid.getWallCellsNodes()
        """

        boundaryNodes = 0
        boundaryCells = 0
        wallBCs = [
            BC["bcwallviscous"],
            BC["bcwall"],
            BC["bcwallinviscid"],
            BC["bcwallviscousheatflux"],
            BC["bcwallviscousisothermal"],
        ]
        for blk in self.blocks:
            for boco in blk.bocos:
                if boco.type in wallBCs:
                    ptRange = boco.ptRange
                    if ptRange[0, 0] == ptRange[0, 1]:
                        boundaryCells += (ptRange[1, 1] - ptRange[1, 0]) * (ptRange[2, 1] - ptRange[2, 0])

                        boundaryNodes += (ptRange[1, 1] - ptRange[1, 0] + 1) * (ptRange[2, 1] - ptRange[2, 0] + 1)

                    elif ptRange[1, 0] == ptRange[1, 1]:
                        boundaryCells += (ptRange[0, 1] - ptRange[0, 0]) * (ptRange[2, 1] - ptRange[2, 0])

                        boundaryNodes += (ptRange[0, 1] - ptRange[0, 0] + 1) * (ptRange[2, 1] - ptRange[2, 0] + 1)

                    elif ptRange[2, 0] == ptRange[2, 1]:
                        boundaryCells += (ptRange[0, 1] - ptRange[0, 0]) * (ptRange[1, 1] - ptRange[1, 0])

                        boundaryNodes += (ptRange[0, 1] - ptRange[0, 0] + 1) * (ptRange[1, 1] - ptRange[1, 0] + 1)
        return boundaryCells, boundaryNodes

    def printInfo(self):
        """
        Prints the number of zones, total number of cells and nodes,
        and the number of cells and nodes on wall boundaries.
        """
        totalCells, totalNodes = self.getTotalCellsNodes()

        print("Total Zones:", len(self.blocks))
        print("Total Cells:", totalCells)
        print("Total Nodes:", totalNodes)

        boundaryCells, boundaryNodes = self.getWallCellsNodes()

        print("Wall Boundary Cells:", boundaryCells)
        print("Wall Boundary Nodes:", boundaryNodes)

    def getBlockInfo(self):
        """Get the number of nodes, number of cells, BCs, and
        the dimensions for each block. This info can be helpful
        for assessing overset meshes."""

        counter = 1
        allBlocksInfo = {}

        for blk in self.blocks:
            blockInfo = {}
            blockInfo["nCells"] = blk.getNumCells()
            blockInfo["nNodes"] = blk.getNumNodes()
            blockInfo["dims"] = list(blk.dims)
            blockInfo["BCs"] = [boco.type for boco in blk.bocos]
            allBlocksInfo[f"{counter}"] = blockInfo
            counter += 1

        totalCells, totalNodes = self.getTotalCellsNodes()

        allBlocksInfo["totalZones"] = len(self.blocks)
        allBlocksInfo["totalCells"] = totalCells
        allBlocksInfo["totalNodes"] = totalNodes

        return allBlocksInfo

    def printBlockInfo(self):
        """Print the number of nodes, number of cells, and
        the dimensions for each block. This info can be helpful
        for assessing overset meshes."""

        allBlocksInfo = self.getBlockInfo()

        for i in range(len(self.blocks)):
            blockNumber = str(i + 1)
            print("Block Number:", blockNumber)
            blockInfo = allBlocksInfo[blockNumber]
            print("Number of Cells:", blockInfo["nCells"])
            print("Number of Nodes:", blockInfo["nNodes"])
            print("Block dimensions:", blockInfo["dims"])

        print("Total Zones:", allBlocksInfo["totalZones"])
        print("Total Cells:", allBlocksInfo["totalCells"])
        print("Total Nodes:", allBlocksInfo["totalNodes"])

    def addBlock(self, blk):

        """Add a block to the grid"""
        self.blocks.append(blk)

    def removeBlocks(self, blockIDs):

        """
        This function will remove certain blocks from the grid.
        The user should ensure that the final mesh is still valid
        in terms of boundary conditions and connectivities.

        ATTENTION: blockIDs should be 1-indexed
        """

        # Remove the blocks in reverse order
        for ID in sorted(blockIDs, reverse=True):
            del self.blocks[ID - 1]

    def writeToCGNS(self, fileName):
        """Write what is in this grid tree to the fileName provided"""
        self.renameBCs()
        outFile = libcgns_utils.utils.openfile(fileName, CG_MODE_WRITE, self.cellDim)
        for blk in self.blocks:
            blk.writeToCGNS(outFile)
        libcgns_utils.utils.closefile(outFile)

    def writeToCGNSSelected(self, fileName, toWrite):
        """Write what is in this grid tree to the fileName provided"""
        outFile = libcgns_utils.utils.openfile(fileName, CG_MODE_WRITE, self.cellDim)
        for iblk in toWrite:
            self.blocks[iblk - 1].writeToCGNS(outFile)
        libcgns_utils.utils.closefile(outFile)

    def writePlot3d(self, fileName):
        """Write what is in this grid tree to the plot3d filename
        provided. This is mostly done in python so will be slow-ish."""
        f = open(fileName, "w")
        f.write("%d\n" % len(self.blocks))
        for blk in self.blocks:
            blk.writeDimsPlot3d(f)
        for blk in self.blocks:
            blk.writeCoordsPlot3d(f)
        f.close()

    def scale(self, scaleFact):
        """Scale blocks in this grid"""
        for blk in self.blocks:
            blk.scale(scaleFact)

    def flip(self, axis):
        """Flip the grid about a axis, 'x', 'y' or 'z'"""
        for blk in self.blocks:
            blk.flip(axis)

    def coarsen(self):
        """Coarsen the block by taking every-other grid line"""
        for blk in self.blocks:
            blk.coarsen()

    def refine(self, axes):
        """Refine the block by interpolating every-other grid line"""
        for blk in self.blocks:
            blk.refine(axes)

    def renameBlocks(self, actualName=False):
        """Rename all blocks in a consistent fashion"""
        i = 1
        for blk in self.blocks:

            # If we the actualName flag is true, then we use the name stored
            # in the block. Otherwise, we use 'domain' as the base of the name.
            # This is to keep the behavior consistent with previous
            # cgns_utils operations while allowing for different naming
            # for use in pyWarpMulti.
            if actualName:
                blk.name = self.name + ".%5.5d" % i
            else:
                blk.name = "domain.%5.5d" % i
            i += 1

    def renameBCs(self):
        """Rename all block boundary conditions in a consistent fashion"""
        i = 1
        for blk in self.blocks:
            for boco in blk.bocos:
                boco.name = "BC%d" % i
                i += 1

    def extractSurface(self, fileName):
        """Extract wall surfaces and write to plot3d file"""
        patches = []
        for blk in self.blocks:
            patches.extend(blk.extractWallSurfaces())
        if len(patches) > 0:
            f = open(fileName, "w")
            f.write("%d\n" % len(patches))
            for i in range(len(patches)):

                f.write("%d %d 1\n" % (patches[i].shape[0], patches[i].shape[1]))
            for i in range(len(patches)):
                patches[i][:, :, 0].flatten("F").tofile(f, sep="\n", format="%20.15g")
                f.write("\n")
                patches[i][:, :, 1].flatten("F").tofile(f, sep="\n", format="%20.15g")
                f.write("\n")
                patches[i][:, :, 2].flatten("F").tofile(f, sep="\n", format="%20.15g")
                f.write("\n")
            f.close()
        else:
            print("Warning: No wall surfaces found!")

    def extractSpecifiedSurface(self, fileName, blkid, imin, imax, jmin, jmax, kmin, kmax):
        """Extract Specified surfaces and write to plot3d file"""
        patches = []
        blk = self.blocks[int(blkid)]
        patches.extend(blk.extractSpecifiedSurfaces(int(imin), int(imax), int(jmin), int(jmax), int(kmin), int(kmax)))
        if len(patches) > 0:
            f = open(fileName, "w")
            f.write("%d\n" % len(patches))
            for i in range(len(patches)):

                f.write("%d %d 1\n" % (patches[i].shape[0], patches[i].shape[1]))
            for i in range(len(patches)):
                patches[i][:, :, 0].flatten("F").tofile(f, sep="\n", format="%20.15g")
                f.write("\n")
                patches[i][:, :, 1].flatten("F").tofile(f, sep="\n", format="%20.15g")
                f.write("\n")
                patches[i][:, :, 2].flatten("F").tofile(f, sep="\n", format="%20.15g")
                f.write("\n")
            f.close()
        else:
            print("Warning: No surfaces found!")

    def overwriteFamilies(self, familyFile):
        """Overwrite families of BC with information given in the
        family file"""

        with open(familyFile, "r") as f:
            for line in f:
                aux = line.split()

                # Check if the read line has the correct number of arguments, if not we stop
                if len(aux) != 3:
                    raise ValueError(
                        f"FamilyFile incorrectly formatted. Line '{line.strip()}' has incorrect number of parameters, got {len(aux)}, expected 3"
                    )

                blockID = int(aux[0]) - 1
                face = aux[1].lower()
                family = aux[2]

                self.blocks[blockID].overwriteFamily(face, family)

    def writeSubfaceFamily(self, familyFile):
        """Add a number of subface Bocos to replace one full-face boco"""
        # Note that this function could easily be expanded to change
        # other information, like bcDataSets() on subfaces as well
        f = open(familyFile, "r")
        blockID = int(f.readline()) - 1
        face = f.readline().lower()[:-1]
        count = 0

        # Locate the Boco we're replacing
        boco = None
        for i in range(len(self.blocks[blockID].bocos)):
            r = self.blocks[blockID].bocos[i].ptRange  # get the point range for existing boco
            if (
                (r[0][0] == r[0][1] == 1 and face == "ilow")
                or (r[0][0] == r[0][1] == self.blocks[blockID].dims[0] and face == "ihigh")
                or (r[1][0] == r[1][1] == 1 and face == "jlow")
                or (r[1][0] == r[1][1] == self.blocks[blockID].dims[1] and face == "jhigh")
                or (r[2][0] == r[2][1] == 1 and face == "klow")
                or (r[2][0] == r[2][1] == self.blocks[blockID].dims[2] and face == "khigh")
            ):
                boco = i
                break

        oldBoco = self.blocks[blockID].bocos[boco]

        # Write the new bocos on the subfaces
        for line in f:
            # Parse out the familyFile info
            aux = line.split()
            if len(aux) == 2:
                ptRanges = numpy.array(aux[0].split(","), dtype=float).reshape(3, 2)
                famName = aux[1]
            else:
                raise ValueError("familyFile is incorrectly formatted.")

            self.blocks[blockID].addBoco(
                Boco(oldBoco.name + "_" + str(count), oldBoco.type, ptRanges, famName, bcDataSets=oldBoco.dataSets)
            )
            count = count + 1

        self.blocks[blockID].bocos.remove(oldBoco)
        f.close()

    def copyFamilyInfo(self, otherGrid):
        """Copy family information out of another grid"""
        for i in range(len(self.blocks)):
            for j in range(len(self.blocks[i].bocos)):
                self.blocks[i].bocos[j].family = otherGrid.blocks[i].bocos[j].family

    def removeBCs(self):
        """Remove any BC's there may be"""
        for i in range(len(self.blocks)):
            self.blocks[i].bocos = []

    def overwriteBCs(self, bcFile):
        """Overwrite BCs with information given in the file"""

        def isFloat(string):
            try:
                float(string)
                return True
            except ValueError:
                return False

        with open(bcFile, "r") as f:
            for line in f:
                if line.strip():
                    aux = line.split()
                    blockID = int(aux[0]) - 1
                    face = aux[1].lower()
                    bocoType = aux[2].lower()
                    family = aux[3]

                    dataSet = []
                    # Check if we have possible datasets specified
                    if len(aux) > 4:
                        bocoSetName = aux[4]
                        bocoDataSetType = aux[5]
                        DirNeu = aux[6]
                        bocoDataSet = BocoDataSet(bocoSetName, BC[bocoDataSetType.lower()])

                        i = 7
                        while i < len(aux):
                            arrayName = aux[i]
                            i += 1
                            dType = CGNSDATATYPES["RealDouble"]
                            nDims = 1

                            dataArr = []

                            for j in range(i, len(aux)):
                                if isFloat(aux[j]):
                                    dataArr.append(aux[j])
                                    i += 1
                                else:
                                    break

                            dataArr = numpy.array(dataArr, dtype=numpy.float64)

                            nDims = 1
                            dataDims = numpy.ones(3, dtype=numpy.int32, order="F")
                            dataDims[0] = dataArr.size

                            bcDataArr = BocoDataSetArray(arrayName, dType, nDims, dataDims, dataArr)
                            if DirNeu == "Dirichlet":
                                bocoDataSet.addDirichletDataSet(bcDataArr)
                            elif DirNeu == "Neumann":
                                bocoDataSet.addNeumannDataSet(bcDataArr)
                            else:
                                raise ValueError(f"DirNeu must be Dirichlet or Neumann. Input was {DirNeu}.")

                        dataSet.append(bocoDataSet)

                    self.blocks[blockID].overwriteBCs(face, bocoType, family, dataSet)

    def writeBCs(self, bcFile):
        """write BCs to file"""
        if bcFile is None:
            bcFile = self.name + ".bc"

        with open(bcFile, "w") as fid:

            for idx_block, block in enumerate(self.blocks):
                blockID = idx_block + 1
                block.writeBCs(blockID, fid)

    def autoOversetBC(self, sym, connectSelf, tol):
        """This is essentially a simplified version of autoBC that flags all
        kMin faces as walls and all kMax faces as BCOverset"""

        # Remove any BCinfo/B2B info we may have.
        for blk in self.blocks:
            blk.bocos = []
            blk.B2Bs = []

        checkSym = True
        if sym == "x":
            symAxis = 0
        elif sym == "y":
            symAxis = 1
        elif sym == "z":
            symAxis = 2
        else:
            symAxis = 0  # doesn't matter
            checkSym = False

        symNormal = [0.0, 0.0, 0.0]
        symNormal[symAxis] = 1.0

        # Do the b2b by running connect:
        if connectSelf:
            types, pointRanges, myIDs, faceAvg, faceNormal = self.connectSelfOnly(tol)
        else:
            types, pointRanges, myIDs, faceAvg, faceNormal = self.connect(tol)

        # Loop over all subfaces and deal with the BCs
        for i in range(len(types)):
            blockID = myIDs[i] - 1

            if types[i] == 0:  # Boco
                coor_check = abs(faceAvg[symAxis, i]) < 1e-3
                dp_check = abs(numpy.dot(faceNormal[:, i], symNormal)) > 0.98
                if dp_check and coor_check and checkSym:
                    bocoType = BC["bcsymmetryplane"]
                    famName = "sym"
                else:
                    # Next check for a wall-type boundary condition if
                    # we have a kMin face
                    if pointRanges[2, 0, i] == pointRanges[2, 1, i] == 1:
                        bocoType = BC["bcwallviscous"]
                        famName = "wall"
                    else:
                        # Must be a overset outer bound
                        bocoType = BC["bcoverset"]
                        famName = "overset"

                # Now simply add the boco
                self.blocks[blockID].addBoco(Boco("dummy", bocoType, pointRanges[:, :, i], famName))

        # Lastly rename the BCs to be consistent
        self.renameBCs()

    def autoNearfieldBC(self, sym):
        """This is essentially a simplified version of autoBC that flags all
        boundaries as BCOverset except for possible symmetry planes."""

        # Remove any BCinfo/B2B info we may have.
        for blk in self.blocks:
            blk.bocos = []
            blk.B2Bs = []

        if sym == "x":
            symAxis = 0
        elif sym == "y":
            symAxis = 1
        else:
            symAxis = 2

        symNormal = [0.0, 0.0, 0.0]
        symNormal[symAxis] = 1.0

        # Do the b2b by running connect:
        types, pointRanges, myIDs, faceAvg, faceNormal = self.connect()

        # Loop over all subfaces and deal with the BCs
        for i in range(len(types)):
            blockID = myIDs[i] - 1

            if types[i] == 0:  # Boco
                coor_check = abs(faceAvg[symAxis, i]) < 1e-3
                dp_check = abs(numpy.dot(faceNormal[:, i], symNormal)) > 0.98
                if dp_check and coor_check:
                    bocoType = BC["bcsymmetryplane"]
                    famName = "sym"
                else:
                    # Flag as overset
                    bocoType = BC["bcoverset"]
                    famName = "overset"

                # Now simply add the boco
                self.blocks[blockID].addBoco(Boco("dummy", bocoType, pointRanges[:, :, i], famName))

        # Lastly rename the BCs to be consistent
        self.renameBCs()

    def autoFarfieldBC(self, sym):
        """This is essentially a simplified version of autoBC that flags all
        boundaries as BCFarfield except for possible symmetry planes."""

        # Remove any BCinfo/B2B info we may have.
        for blk in self.blocks:
            blk.bocos = []
            blk.B2Bs = []

        if sym == "x":
            symAxis = 0
        elif sym == "y":
            symAxis = 1
        else:
            symAxis = 2

        symNormal = [0.0, 0.0, 0.0]
        symNormal[symAxis] = 1.0

        # Do the b2b by running connect:
        types, pointRanges, myIDs, faceAvg, faceNormal = self.connect()

        # Loop over all subfaces and deal with the BCs
        for i in range(len(types)):
            blockID = myIDs[i] - 1

            if types[i] == 0:  # Boco
                coor_check = abs(faceAvg[symAxis, i]) < 1e-3
                dp_check = abs(numpy.dot(faceNormal[:, i], symNormal)) > 0.98
                if dp_check and coor_check:
                    bocoType = BC["bcsymmetryplane"]
                    famName = "sym"
                else:
                    # Flag as farfield
                    bocoType = BC["bcfarfield"]
                    famName = "far"

                # Now simply add the boco
                self.blocks[blockID].addBoco(Boco("dummy", bocoType, pointRanges[:, :, i], famName))

        # Lastly rename the BCs to be consistent
        self.renameBCs()

    def double2D(self):
        """Doubles a mesh in the "2d" direction. Ie the direction with one
        cell"""
        for blk in self.blocks:
            blk.double2D()

    def simpleCart(self, dh, hExtra, nExtra, sym, mgcycle, outFile):
        """Generates a Cartesian mesh around the provided grid"""

        # Get the bounds of each grid.
        xMin = 1e20 * numpy.ones(3)
        xMax = -1.0 * numpy.ones(3)

        for blk in self.blocks:
            tmp1 = numpy.min(blk.coords, axis=(0, 1, 2))
            tmp2 = numpy.max(blk.coords, axis=(0, 1, 2))
            for iDim in range(3):
                xMin[iDim] = min(xMin[iDim], tmp1[iDim])
                xMax[iDim] = max(xMax[iDim], tmp2[iDim])

        # Call the generic routine
        return simpleCart(xMin, xMax, dh, hExtra, nExtra, sym, mgcycle, outFile)

    def simpleOCart(self, dh, hExtra, nExtra, sym, mgcycle, outFile, userOptions=None):
        """Generates a Cartesian mesh around the provided grid, surrounded by an O-mesh.
        This function requires pyHyp to be installed. If this function is run with MPI,
        pyHyp will be run in parallel.
        """

        # First run simpleCart with no extension:
        X, dx = self.simpleCart(dh, 0.0, 0, sym, mgcycle, outFile=None)

        # Pull out the patches. Note that we have to pay attention to
        # the symmetry and the ordering of the patches to make sure
        # that all the normals are pointing out.
        patches = []

        # First take patches that are opposite from the origin planes
        if "xmax" not in sym:
            patches.append(X[-1, :, :, :])
        if "ymax" not in sym:
            patches.append(X[:, -1, :, :][::-1, :, :])
        if "zmax" not in sym:
            patches.append(X[:, :, -1, :])

        if "x" not in sym and "xmin" not in sym:
            patches.append(X[0, :, :, :][::-1, :, :])
        if "y" not in sym and "ymin" not in sym:
            patches.append(X[:, 0, :, :])
        if "z" not in sym and "zmin" not in sym:
            patches.append(X[:, :, 0, :][::-1, :, :])

        # Set up the generic input for pyHyp
        hypOptions = {
            "patches": patches,
            "unattachedEdgesAreSymmetry": True,
            "outerFaceBC": "farfield",
            "autoConnect": True,
            "BC": {},
            "N": nExtra,
            "s0": numpy.average(dx),
            "marchDist": hExtra,
            "cmax": 3.0,
        }

        # Use user-defined options if provided
        if userOptions is not None:
            hypOptions.update(userOptions)

        # Run pyHyp
        from pyhyp import pyHyp

        hyp = pyHyp(options=hypOptions)
        hyp.run()

        from mpi4py import MPI

        fName = None
        if MPI.COMM_WORLD.rank == 0:
            dirpath = tempfile.mkdtemp()
            fName = os.path.join(dirpath, "tmp.cgns")

        hyp.writeCGNS(MPI.COMM_WORLD.bcast(fName))

        # Reset symmetry to single axis
        if "x" in sym or "xmin" in sym or "xmax" in sym:
            sym = "x"
        elif "y" in sym or "ymin" in sym or "ymax" in sym:
            sym = "y"
        elif "z" in sym or "zmin" in sym or "zmax" in sym:
            sym = "z"

        if MPI.COMM_WORLD.rank == 0:
            # Read the pyhyp mesh back in and add our additional "X" from above.
            grid = readGrid(fName)
            dims = X.shape[0:3]
            grid.addBlock(Block("interiorBlock", dims, X))
            grid.renameBlocks()
            grid.connect()
            grid.BCs = []
            grid.autoFarfieldBC(sym)
            grid.writeToCGNS(outFile)

            # Delete the temp file
            os.remove(fName)

    def cartesian(self, cartFile, outFile):
        """Generates a Cartesian mesh around the provided grid"""

        # PARAMETERS
        inLayer = 2  # How many layers of the overset interpolation
        # faces will be used for volume computation

        print("Running Cartesian grid generator")

        # Preallocate arrays
        extensions = numpy.zeros((2, 3), order="F")
        nNodes = numpy.zeros(3, order="F")
        weightGR = numpy.zeros(3, order="F")
        numBins = numpy.zeros(3, order="F")

        # Read four lines of the cartesian specs file
        with open(cartFile, "r") as f:
            lines = list(f)
        extensions[0, :] = lines[0].split()
        extensions[1, :] = lines[1].split()
        nNodes[:] = lines[2].split()
        weightGR[:] = lines[3].split()

        # Specify number of bins
        numBins[:] = 1  # The tangent law only works for single bin

        # Initialize bounding box coordinates using the first point of the first zone
        xBounds = numpy.zeros((2, 3), order="F")
        xBounds[0, 0] = self.blocks[0].coords[0, 0, 0, 0]  # Using the first point for initialization
        xBounds[1, 0] = self.blocks[0].coords[0, 0, 0, 0]  # because I can't use 0
        xBounds[0, 1] = self.blocks[0].coords[0, 0, 0, 1]
        xBounds[1, 1] = self.blocks[0].coords[0, 0, 0, 1]
        xBounds[0, 2] = self.blocks[0].coords[0, 0, 0, 2]
        xBounds[1, 2] = self.blocks[0].coords[0, 0, 0, 2]
        binVolX = numpy.zeros(numBins[0], order="F")  # Assign zeroes to all bins
        binVolY = numpy.zeros(numBins[1], order="F")
        binVolZ = numpy.zeros(numBins[2], order="F")
        binCellsX = numpy.zeros(numBins[0], order="F", dtype=int)  # Initialize cells counter for each bin
        binCellsY = numpy.zeros(numBins[1], order="F", dtype=int)
        binCellsZ = numpy.zeros(numBins[2], order="F", dtype=int)

        # Loop over all blocks to find the bounding box coordinates
        for index in range(len(self.blocks)):
            # Loop over all BCs of this block
            for boco in self.blocks[index].bocos:
                # Check if we have an overset boundary condition
                if boco.type == BC["bcoverset"]:
                    # Find overset BC face and select some inner layers to compute volume
                    r = boco.ptRange
                    if r[0][0] == r[0][1] == 1:  # ilow detected
                        imin = 0
                        imax = min(0 + inLayer, self.blocks[index].dims[0])
                        jmin = r[1][0] - 1
                        jmax = r[1][1]
                        kmin = r[2][0] - 1
                        kmax = r[2][1]
                    elif r[0][0] == r[0][1] == self.blocks[index].dims[0]:  # ihigh detected
                        imin = max(self.blocks[index].dims[0] - inLayer, 0)
                        imax = self.blocks[index].dims[0]
                        jmin = r[1][0] - 1
                        jmax = r[1][1]
                        kmin = r[2][0] - 1
                        kmax = r[2][1]
                    elif r[1][0] == r[1][1] == 1:  # jlow detected
                        imin = r[0][0] - 1
                        imax = r[0][1]
                        jmin = 0
                        jmax = min(0 + inLayer, self.blocks[index].dims[1])
                        kmin = r[2][0] - 1
                        kmax = r[2][1]
                    elif r[1][0] == r[1][1] == self.blocks[index].dims[1]:  # jhigh detected
                        imin = r[0][0] - 1
                        imax = r[0][1]
                        jmin = max(self.blocks[index].dims[1] - inLayer, 0)
                        jmax = self.blocks[index].dims[1]
                        kmin = r[2][0] - 1
                        kmax = r[2][1]
                    elif r[2][0] == r[2][1] == 1:  # klow detected
                        imin = r[0][0] - 1
                        imax = r[0][1]
                        jmin = r[1][0] - 1
                        jmax = r[1][1]
                        kmin = 0
                        kmax = min(0 + inLayer, self.blocks[index].dims[2])
                    elif r[2][0] == r[2][1] == self.blocks[index].dims[2]:  # khigh detected
                        imin = r[0][0] - 1
                        imax = r[0][1]
                        jmin = r[1][0] - 1
                        jmax = r[1][1]
                        kmin = max(self.blocks[index].dims[2] - inLayer, 0)
                        kmax = self.blocks[index].dims[2]
                    # Use the range to compute average volume
                    libcgns_utils.utils.findbounds(
                        self.blocks[index].coords[imin:imax, jmin:jmax, kmin:kmax, :], xBounds
                    )

        # Loop over all blocks to find the bin volumes
        for index in range(len(self.blocks)):
            # Loop over all BCs of this block
            for boco in self.blocks[index].bocos:
                # Check if we have an overset boundary condition
                if boco.type == BC["bcoverset"]:
                    # Find overset BC face and select some inner layers to compute volume
                    r = boco.ptRange
                    if r[0][0] == r[0][1] == 1:  # ilow detected
                        imin = 0
                        imax = min(0 + inLayer, self.blocks[index].dims[0])
                        jmin = r[1][0] - 1
                        jmax = r[1][1]
                        kmin = r[2][0] - 1
                        kmax = r[2][1]
                    elif r[0][0] == r[0][1] == self.blocks[index].dims[0]:  # ihigh detected
                        imin = max(self.blocks[index].dims[0] - inLayer, 0)
                        imax = self.blocks[index].dims[0]
                        jmin = r[1][0] - 1
                        jmax = r[1][1]
                        kmin = r[2][0] - 1
                        kmax = r[2][1]
                    elif r[1][0] == r[1][1] == 1:  # jlow detected
                        imin = r[0][0] - 1
                        imax = r[0][1]
                        jmin = 0
                        jmax = min(0 + inLayer, self.blocks[index].dims[1])
                        kmin = r[2][0] - 1
                        kmax = r[2][1]
                    elif r[1][0] == r[1][1] == self.blocks[index].dims[1]:  # jhigh detected
                        imin = r[0][0] - 1
                        imax = r[0][1]
                        jmin = max(self.blocks[index].dims[1] - inLayer, 0)
                        jmax = self.blocks[index].dims[1]
                        kmin = r[2][0] - 1
                        kmax = r[2][1]
                    elif r[2][0] == r[2][1] == 1:  # klow detected
                        imin = r[0][0] - 1
                        imax = r[0][1]
                        jmin = r[1][0] - 1
                        jmax = r[1][1]
                        kmin = 0
                        kmax = min(0 + inLayer, self.blocks[index].dims[2])
                    elif r[2][0] == r[2][1] == self.blocks[index].dims[2]:  # khigh detected
                        imin = r[0][0] - 1
                        imax = r[0][1]
                        jmin = r[1][0] - 1
                        jmax = r[1][1]
                        kmin = max(self.blocks[index].dims[2] - inLayer, 0)
                        kmax = self.blocks[index].dims[2]
                    # Use the range to compute average volume
                    libcgns_utils.utils.computevolumes(
                        self.blocks[index].coords[imin:imax, jmin:jmax, kmin:kmax, :],
                        xBounds,
                        binVolX,
                        binVolY,
                        binVolZ,
                        binCellsX,
                        binCellsY,
                        binCellsZ,
                    )

        # DEFINE UNIDIMENSIONAL GRID GENERATION ROUTINES

        # Define tangent bunching law
        def tanDist(Sp1, Sp2, N):

            # This is the tangential spacing developed by Ney Secco
            # This bunching law is coarse at the ends an fine at the middle
            # of the interval, just like shown below:
            # |    |   |  | || |  |   |    |

            # Sp1: initial spacing (within the [0,1] interval)
            # Sp2: final spacing (within the [0,1] interval)
            # N: number of nodes

            # IMPORTS
            from numpy import tan, arange, pi
            from scipy.optimize import minimize

            # Convert number of nodes to number of cells, because I derived the equations using
            # N the as number of cells =P.
            N = N - 1

            # Define objective function
            def func(P):
                # Split variables
                a = P[0]
                e = P[1]
                c = P[2]
                # Find b
                b = e - c
                # Equations
                Eq1 = a * (tan(b + c) - tan(c)) - 1
                Eq2 = a * (tan(b / N + c) - tan(c)) - Sp1
                Eq3 = a * (tan(b + c) - tan(b * (1 - 1 / N) + c)) - Sp2
                # Cost function
                J = Eq1 ** 2 + Eq2 ** 2 + Eq3 ** 2
                # Return
                return J

            # Define bounds for the problem
            a_bounds = [(0, None)]
            e_bounds = [(0, pi / 2)]
            c_bounds = [(-pi / 2, 0)]
            bounds = a_bounds + e_bounds + c_bounds

            # Define initial guess
            a_start = 1.0
            e_start = pi / 4
            c_start = -pi / 4
            x_start = [a_start, e_start, c_start]

            # Optimize
            res = minimize(
                func, x_start, method="SLSQP", bounds=bounds, options={"maxiter": 1000, "disp": False, "ftol": 1e-12}
            )

            # Split variables
            a = res.x[0]
            e = res.x[1]
            c = res.x[2]

            # Find other parameters
            b = e - c
            d = -a * tan(c)

            # Generate spacing
            index = arange(N + 1)
            S = a * tan(b * index / N + c) + d

            # Force the extremes to 0 and 1 so that we always meet the bounds
            # (this is to avoid numerical issues with symmetry planes)
            S[0] = 0.0
            S[-1] = 1.0

            # Return spacing
            return S

        # Define function that optimizes bunching law to match grid resolution

        def generateGrid(xmin, xmax, extension1, extension2, nNodes, binVol, weightGR):

            # xmin: float -> position where the bounding box begins
            # xmax: float -> position where the bounding box ends
            # extension1: float > 0 -> ratio between the negative farfield distance and the bounding box length:
            #                          extension1 = (xmin-negative_farfield_position)/(xmax-xmin)
            # extension2: float > 0 -> ratio between the positive farfield distance and the bounding box length:
            #                          extension2 = (positive_farfield_position-xmax)/(xmax-xmin)
            # nNodes: integer > 0 -> Number of nodes along the edge
            # binVol: float > 0 -> Average volume of the bounding box cells (foreground mesh)
            # weightGR: 0 < float < 1 -> Weight used to balance growth ratio and cell volume during the optimization.
            #                            If weightGR = 0, the optimizer will not care about the growth ratios at the
            #                            farfield and will just try to match the bounding box resolution.
            #                            If weightGR = 1, the optimizer will not care about the bounding box resolution
            #                            and will just try to get an uniform growth ratio. This results in an uniform mesh.

            # IMPORTS
            from numpy import array, mean
            from scipy.optimize import minimize

            # Compute farfield coordinates
            x0 = xmin - (xmax - xmin) * extension1
            xf = xmax + (xmax - xmin) * extension2

            # Get number of bins and bin size
            nBins = len(binVol)
            dxBin = (xmax - xmin) / nBins

            # Get bin edges
            binEdge = binVol ** (1.0 / 3.0)

            # Define objective function
            def func(P):
                # Split variables
                Sp1 = P[0]
                Sp2 = P[1]

                # Generate grid coordinates with tangent bunching law
                S = tanDist(Sp1, Sp2, nNodes)

                # Rescale the interval
                S = x0 + S * (xf - x0)

                # Compute edge size of each cell
                E = S[1:] - S[:-1]

                # Initialize edge error
                edgeError = 0

                # Find cells that are inside each bin and check the edge difference
                for binIndex in range(nBins):
                    # Find bin interval
                    x0bin = xmin + dxBin * binIndex
                    xfbin = xmin + dxBin * (binIndex + 1)
                    # Find cells that touch this interval and get their edges
                    bol = -(((S[:-1] < x0bin) * (S[1:] < x0bin)) + ((S[:-1] > xfbin) * (S[1:] > xfbin)))
                    bolEdges = E[bol]
                    # print bol
                    # Compute edge mismatch and increment variable
                    edgeError = edgeError + mean((bolEdges - binEdge[binIndex]) ** 2) / 2

                # Compute term regarding growing ratios at the ends
                if nNodes > 3:
                    growthRatio = ((S[1] - S[0]) / (S[2] - S[1]) - 1.0) ** 2 + (
                        (S[-1] - S[-2]) / (S[-2] - S[-3]) - 1
                    ) ** 2
                else:  # There's no way to define growth ratio when we have less than 3 cells
                    growthRatio = 0

                # Return objective function
                return (1 - weightGR) * edgeError / mean(binEdge) + weightGR * growthRatio
                # Note that the edgeError is normalized so that the weighed average makes sense

            # Define initial guess based on uniform spacing
            Sp1_start = 1 / (nNodes - 1)
            Sp2_start = 1 / (nNodes - 1)
            x_start = array([Sp1_start, Sp2_start])

            # Optimize
            res = minimize(
                func, x_start, method="Nelder-Mead", options={"maxiter": 2000, "disp": True, "xtol": 1e-8, "ftol": 1e-8}
            )

            # Split variables
            Sp1 = res.x[0]
            Sp2 = res.x[1]

            # Generate grid
            S = tanDist(Sp1, Sp2, nNodes)
            S = x0 + S * (xf - x0)

            # Return grid
            return S

        # Generate grid for each dimension
        Sx = generateGrid(
            xBounds[0, 0], xBounds[1, 0], extensions[0, 0], extensions[1, 0], nNodes[0], binVolX[0:1], weightGR[0]
        )
        Sy = generateGrid(
            xBounds[0, 1], xBounds[1, 1], extensions[0, 1], extensions[1, 1], nNodes[1], binVolY[0:1], weightGR[1]
        )
        Sz = generateGrid(
            xBounds[0, 2], xBounds[1, 2], extensions[0, 2], extensions[1, 2], nNodes[2], binVolZ[0:1], weightGR[2]
        )

        # Compute growth ratios
        if nNodes[0] > 3:
            gx = max((Sx[1] - Sx[0]) / (Sx[2] - Sx[1]), (Sx[-1] - Sx[-2]) / (Sx[-2] - Sx[-3]))
        else:
            gx = None
        if nNodes[1] > 3:
            gy = max((Sy[1] - Sy[0]) / (Sy[2] - Sy[1]), (Sy[-1] - Sy[-2]) / (Sy[-2] - Sy[-3]))
        else:
            gy = None
        if nNodes[2] > 3:
            gz = max((Sz[1] - Sz[0]) / (Sz[2] - Sz[1]), (Sz[-1] - Sz[-2]) / (Sz[-2] - Sz[-3]))
        else:
            gz = None

        # Print growth ratios
        print("")
        print("Maximum growth ratios along each direction:")
        print("X: ", gx)
        print("Y: ", gy)
        print("Z: ", gz)
        if max(gx, gy, gz) > 1.2:
            print("You may bring weightGR closer to 1 to decrease ratios")
        print("")

        # Allocate coordinates block
        X = numpy.zeros((nNodes[0], nNodes[1], nNodes[2], 3))

        # Write grid coordinates
        Xx, Xy, Xz = numpy.meshgrid(Sx, Sy, Sz, indexing="ij")
        X[:, :, :, 0] = Xx
        X[:, :, :, 1] = Xy
        X[:, :, :, 2] = Xz

        # Open a new CGNS file
        cg = libcgns_utils.utils.openfile(outFile, CG_MODE_WRITE, 3)

        # Write a Zone to it
        zoneID = libcgns_utils.utils.writezone(cg, "cartesian", nNodes)

        # Write mesh coordinates
        libcgns_utils.utils.writecoordinates(cg, zoneID, X)

        # CLose file
        libcgns_utils.utils.closefile(cg)

        # Print
        print("Mesh successfully generated and stored in: " + outFile)

    def split(self, extraSplits):

        """Recursively propagate splits due to boundary conditions or
        B2B information"""

        # First generate a mapping between block name and its index:
        mapping = {}
        for iBlock in range(len(self.blocks)):
            mapping[self.blocks[iBlock].name] = iBlock

        for iBlock in range(len(self.blocks)):
            splits = []
            for boco in self.blocks[iBlock].bocos:
                splits.extend(getSplits(boco.ptRange))
            for b2b in self.blocks[iBlock].B2Bs:
                splits.extend(getSplits(b2b.ptRange))

            # Now just add the (unique) splits for this block: DON't
            # USE numpy.unique it doesn't actually work for tuples.
            newSplits = []
            for split in splits:
                if split not in newSplits:
                    newSplits.append(split)
            splits = newSplits

            for split in splits:
                self._addSplit(iBlock, split[0] + 1, split[1], mapping)

        # And Add the extra splits:
        for split in extraSplits:
            self._addSplit(split[0] - 1, split[1], split[2], mapping)

    def _addSplit(self, iBlock, iDim, index, mapping):
        """Recursive routine to add a split to block 'iBlock', on
        dimension 'iDim' at index 'index'. NOTE IDIM is 1 based!"""
        if index in self.blocks[iBlock].splits[iDim - 1]:
            return  # This is the main recursive return
        else:
            # Add the split and call any others we need
            self.blocks[iBlock].splits[iDim - 1].append(index)
            self.blocks[iBlock].splits[iDim - 1].sort()

            for b2b in self.blocks[iBlock].B2Bs:
                low = min(b2b.ptRange[iDim - 1, :])
                high = max(b2b.ptRange[iDim - 1, :])

                # Index must be fully contained:
                if index > low and index < high:
                    newBlock = mapping[b2b.donorName]
                    iDim_new = b2b.transform[iDim - 1]
                    offset = index - low
                    abs_idim = abs(iDim_new)
                    donor_high = max(b2b.donorRange[abs_idim - 1, :])
                    donor_low = min(b2b.donorRange[abs_idim - 1, :])

                    if iDim_new >= 0:
                        index_new = donor_low + offset
                    else:
                        index_new = donor_high - offset

                    # Finally recursively call itself for the new
                    # block, dimension and index
                    self._addSplit(newBlock, abs_idim, index_new, mapping)

    def connect(self, tol=1e-12):
        """Generate block-to-block connectivity information for a grid. It
        does not need to be face matched, only point matched"""
        isize = 0
        for i in range(len(self.blocks)):
            blk = self.blocks[i]
            isize += blk.dims[0] * blk.dims[1] * blk.dims[2]

        # Allocate space for all coordinates
        coords = numpy.zeros(isize * 3)
        sizes = []
        istart = 0
        for i in range(len(self.blocks)):
            blk = self.blocks[i]
            iend = istart + blk.dims[0] * blk.dims[1] * blk.dims[2]
            coords[istart * 3 : 3 * iend] = blk.coords.flatten()
            sizes.append(blk.dims)
            istart = iend

        # Get our list of sizes
        sizes = numpy.vstack(sizes)

        # Run the fortran code to generate all the connectivities
        libcgns_utils.utils.computeconnectivity(coords, sizes.T, tol)
        nPatches = libcgns_utils.utils.getnpatches()
        (
            types,
            pointRanges,
            myIDs,
            pointRangeDonors,
            transforms,
            donorIDs,
            faceAvgs,
            faceNormals,
        ) = libcgns_utils.utils.getpatchinfo(nPatches)
        libcgns_utils.utils.deallocpatches()

        # Remove all existing B2B info
        for blk in self.blocks:
            blk.B2Bs = []

        for i in range(nPatches):
            blockID = myIDs[i] - 1

            if types[i] == 1:  # B2B
                connectName = "SF%d" % i
                donorName = self.blocks[donorIDs[i] - 1].name

                self.blocks[blockID].B2Bs.append(
                    B2B(connectName, donorName, pointRanges[:, :, i], pointRangeDonors[:, :, i], transforms[:, i])
                )

        # Return most of the information we computed since other
        # routines (autobc for example) may need this.

        return types, pointRanges, myIDs, faceAvgs, faceNormals

    def connectSelfOnly(self, tol=1e-12):
        """Generate block-to-block connectivity information for a grid, but
        only for b2b connections within a given block. Ie only periodic conditions
        """

        types = []
        pointRanges = []
        myIDs = []
        faceAvgs = []
        faceNormals = []
        for i in range(len(self.blocks)):
            blk = self.blocks[i]
            coords = blk.coords.flatten()
            sizes = numpy.array([blk.dims])

            # Run the fortran code to generate all the connectivities
            libcgns_utils.utils.computeconnectivity(coords, sizes.T, tol)
            nPatches = libcgns_utils.utils.getnpatches()
            (
                t,
                pointRange,
                myID,
                pointRangeDonor,
                transform,
                donorID,
                faceAvg,
                faceNormal,
            ) = libcgns_utils.utils.getpatchinfo(nPatches)
            libcgns_utils.utils.deallocpatches()

            # Remove all existing B2B info
            blk.B2Bs = []

            for j in range(nPatches):
                if t[j] == 1:  # B2B
                    connectName = "SF%d_%d" % (i, j)
                    donorName = blk.name  # Has to be the same block
                    blk.B2Bs.append(
                        B2B(connectName, donorName, pointRange[:, :, j], pointRangeDonor[:, :, j], transform[:, j])
                    )

                # Also append this information to return the same way
                # that connect does:
                types.append(t[j])
                myIDs.append(myID[j])
                pointRanges.append(pointRange[:, :, j])
                faceAvgs.append(faceAvg[:, j])
                faceNormals.append(faceNormal[:, j])

        # Return the information we computed since other
        # routines (autobc for example) may need this.
        pointRanges = numpy.moveaxis(numpy.array(pointRanges), 0, -1)
        faceNormals = numpy.moveaxis(numpy.array(faceNormals), 0, -1)
        faceAvgs = numpy.moveaxis(numpy.array(faceAvgs), 0, -1)
        return (numpy.array(types), pointRanges, numpy.array(myIDs), faceAvgs, faceNormals)

    def autoBC(self, radius, sym, offset):
        """This function will try to generate boundary condition
        information for all patches that are not part of a
        block-to-block connection. If a surface is inside the sphere,
        it gets counted as a wall, if it is outside it is a farfield
        condition. If the surface is flat and a coordinate is zero, it
        gets treated as a symmetry plane."""

        # Remove any BCinfo/B2B info we may have.
        for blk in self.blocks:
            blk.bocos = []
            blk.B2Bs = []

        if sym == "x":
            symAxis = 0
        elif sym == "y":
            symAxis = 1
        else:
            symAxis = 2

        symNormal = [0.0, 0.0, 0.0]
        symNormal[symAxis] = 1.0

        # Do the b2b by running connect:
        types, pointRanges, myIDs, faceAvg, faceNormal = self.connect()

        # Loop over all subfaces and deal with the BCs
        for i in range(len(types)):
            blockID = myIDs[i] - 1

            if types[i] == 0:  # Boco
                coor_check = abs(faceAvg[symAxis, i]) < 1e-3
                dp_check = abs(numpy.dot(faceNormal[:, i], symNormal)) > 0.98
                if dp_check and coor_check:
                    bocoType = BC["bcsymmetryplane"]
                    famName = "sym"
                else:
                    # Next check for a wall-type boundary condition if
                    # the face avg is inside the sphere:

                    if numpy.linalg.norm(faceAvg[:, i] - offset) < radius:
                        bocoType = BC["bcwallviscous"]
                        famName = "wall"
                    else:
                        # Must be a farfield
                        bocoType = BC["bcfarfield"]
                        famName = "far"

                # Now simply add the boco
                self.blocks[blockID].addBoco(Boco("dummy", bocoType, pointRanges[:, :, i], famName))

        # Lastly rename the BCs to be consistent
        self.renameBCs()

    def fillOpenBCs(self, bocoType, famName):
        """This function will add the desired BC for all faces that are not
        block to block and also do not have any previously assigned BC."""

        # Remove any B2B info we may have.
        for blk in self.blocks:
            blk.B2Bs = []

        # Do the b2b by running connect:
        types, pointRanges, myIDs, faceAvg, faceNormal = self.connect()

        # Loop over all subfaces and deal with the BCs
        for i in range(len(types)):
            # Get reference to block
            blockID = myIDs[i] - 1

            if types[i] == 0:  # Boco

                # Check if face already has a BC
                has_bc = False
                for boco in self.blocks[blockID].bocos:
                    # Get norm of difference of point range
                    diff_norm = numpy.linalg.norm(pointRanges[:, :, i] - boco.ptRange)
                    # Check if BC already exists
                    if diff_norm < 1e-10:
                        has_bc = True

                # Add new BC if necessary
                if not has_bc:
                    self.blocks[blockID].addBoco(Boco("dummy", bocoType, pointRanges[:, :, i], famName))

        # Lastly rename the BCs to be consistent
        self.renameBCs()

    def rebunch(self, spacing, extraCells, nStar):
        """Perform rebunching on offwall-directions. The user should
        be *VERY* careful with this function. It will *only* work for
        grids that that have 'O-type' topologies around the
        surface. This is typical of viscous grids. The main
        application is to rebunch nodes in the boundary layer to adapt
        an existing grid for a different reynolds number"""

        for blk in self.blocks:
            blk.rebunch(spacing, extraCells, nStar)
            blk.B2Bs = []
            blk.BCs = []
        self.connect()

    def randomize(self, seed, keepRHS):
        """Perform random reording of grid orientation and block numbers. This
        method destroys *ALL* boundary condition information. Grid
        connectivity is recomputed on the reorginized grid. Actual
        nodal locations remain precisely unchanged; only the grid
        ordering changes. Block handendness is not necessairly
        preserved.
        """
        numpy.random.seed(seed)
        for blk in self.blocks:
            blk.bocos = []
            blk.B2Bs = []
            blk.randomize(keepRHS)

        # Finally reconnect
        self.connect()

    def reorder(self, intDigits):
        """
        When CGNSlib generates a CGNS file (when converting from a plot3d file, for instance),
        it does not add extra digits to the integers when naming zones. This becomes a problem
        when you have more than 10 zones because the ordering will be:
        Zone1, Zone11, Zone12, ..., Zone19, Zone2, Zone21, ...
        This method will add extra digits to the zone names to give the correct ordering.
        """

        # IMPORTS
        import re

        # Initialize list of names
        nameList = []

        # Loop over the blocks to add more significant digits to the last integer
        for blk in self.blocks:
            # Find last integer in the current name
            last_int = re.findall(rb"\d+", blk.name)[-1]

            # Apply modifications only if we have found an integer
            if last_int:
                # Crop the integer from the name
                blk.name = blk.name[: -len(last_int)]
                # Add zeros to get the necessary number of digits
                last_int = last_int.zfill(intDigits)
                # Append integer with more significant digits back to the name
                blk.name = blk.name + last_int

            # Append the name to the names list
            nameList.append(blk.name)

        # Reorder blocks based on their new names
        self.blocks = [blk for (n, blk) in sorted(zip(nameList, self.blocks))]

    def symmZero(self, sym):
        """Zero nodes along axis 'sym'"""
        if sym == "x":
            idir = 0
        elif sym == "y":
            idir = 1
        elif sym == "z":
            idir = 2
        for blk in self.blocks:
            blk.symmZero(idir)

    def symmZeroNoBC(self, sym, tol):
        """Zero nodes below tol distance from symmetry plane"""
        if sym == "x":
            idir = 0
        elif sym == "y":
            idir = 1
        elif sym == "z":
            idir = 2
        for blk in self.blocks:
            blk.symmZeroNoBC(idir, tol)

    def translate(self, dx, dy, dz):
        for blk in self.blocks:
            blk.coords[:, :, :] += [dx, dy, dz]

    def rotate(self, vx, vy, vz, theta):

        """
        This rotates the grid around an axis that passes through the origin.
        vx, vy, vz are the components of the rotation vector
        theta is the rotation angle, in degrees.

        Ney Secco 2016-11
        """

        # Normalize the components of the rotation vector
        normV = numpy.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
        uu = vx / normV
        vv = vy / normV
        ww = vz / normV

        # Compute sines and cosines of the rotation angle
        ss = numpy.sin(theta * numpy.pi / 180.0)
        cc = numpy.cos(theta * numpy.pi / 180.0)

        # Build rotation matrix
        rotMat = numpy.zeros((3, 3))
        rotMat[0, 0] = uu * uu + (1.0 - uu * uu) * cc
        rotMat[0, 1] = uu * vv * (1.0 - cc) - ww * ss
        rotMat[0, 2] = uu * ww * (1.0 - cc) + vv * ss
        rotMat[1, 0] = uu * vv * (1.0 - cc) + ww * ss
        rotMat[1, 1] = vv * vv + (1.0 - vv * vv) * cc
        rotMat[1, 2] = vv * ww * (1.0 - cc) - uu * ss
        rotMat[2, 0] = uu * ww * (1.0 - cc) - vv * ss
        rotMat[2, 1] = vv * ww * (1.0 - cc) + uu * ss
        rotMat[2, 2] = ww * ww + (1.0 - ww * ww) * cc

        for blk in self.blocks:
            blk.coords[:, :, :] = numpy.dot(blk.coords[:, :, :], rotMat)

    def extrude(self, direction):
        """
        Takes a planar grid in 2D and extrudes into the third
        dimension making a 3D that is single cell wide. This routine
        maintains the BCs and creates 2 new symm BCs for each side.

        direction: "str" {x,y,z}
        """

        # Extrude all blocks
        for blk in self.blocks:
            blk.extrude(direction)

        # Rebuild B2B connectivity
        self.connect()

    def revolve(self, normalDirection, axis, startAngle, endAngle, nThetas):
        """
        Takes a planar grid in 2D and revolves about specified axis to
        make a 3D axisymmetric mesh. This routine maintains the BCs and
        creates 2 new symm BCs for each side.

        normalDirection: "str" {x,y,z}
        axis: "str" {x,y,z}
        angle: "float" degrees
        nThetas: "int" number of points in the theta direction
        """

        newGrid = Grid()  # need a dummy 1-cell wide grid to get the connectives from

        # revolve the blocks
        for blk in self.blocks:
            new_blk = copy.deepcopy(blk)
            newGrid.addBlock(new_blk)

            new_blk.revolve(normalDirection, axis, startAngle, endAngle, 2)
            blk.revolve(normalDirection, axis, startAngle, endAngle, nThetas)

        # Rebuild B2B connectivity
        newGrid.connect()

        for blk, new_blk in zip(self.blocks, newGrid.blocks):
            # empty the connectivities from the current grid
            blk.B2Bs = []

            # grab the connectivities from the 1-cell wide,
            # modify them, then add them to the original grid

            for b2b in new_blk.B2Bs:
                pt_rng = b2b.ptRange
                pt_rng[pt_rng == 2] = nThetas
                # print(b2b.ptRange)

                dnr_rng = b2b.donorRange
                dnr_rng[dnr_rng == 2] = nThetas

                blk.addB2B(b2b)

    def addConvArray(self, arrayName, arrayData):
        # This method just appends a new array data to the convergence history dictionary
        self.convArray[arrayName] = arrayData


class Block(object):
    """Class for storing information related to a single block
    structured zone"""

    def __init__(self, zoneName, dims, coords):
        self.name = zoneName.strip()
        self.dims = dims
        self.coords = coords
        self.bocos = []
        self.B2Bs = []
        self.splits = [[1, dims[0]], [1, dims[1]], [1, dims[2]]]
        self.bocoCounter = 0

    def addBoco(self, boco):
        """A add a boundary condition to this block"""
        self.bocos.append(boco)

    def addB2B(self, b2b):
        """A  block-2-block connection to this block"""
        self.B2Bs.append(b2b)

    def writeToCGNS(self, cg):
        """Write all information in this block to the cg file handle"""
        zoneID = libcgns_utils.utils.writezone(cg, self.name, self.dims)
        libcgns_utils.utils.writecoordinates(cg, zoneID, self.coords)
        for boco in self.bocos:
            iBC = libcgns_utils.utils.writebc(cg, zoneID, boco.name, boco.family, boco.ptRange, boco.type)
            for dataSet in boco.dataSets:
                # Write the header for the BCDataSet
                iDataSet = libcgns_utils.utils.writebcdataheader(cg, zoneID, dataSet.type, iBC, dataSet.name)

                # Loop over all Dirichlet and Neumann sets
                writeBCDataHeader = True
                for dirArr in dataSet.dirichletArrays:
                    libcgns_utils.utils.writebcdata(
                        cg,
                        zoneID,
                        iBC,
                        iDataSet,
                        BCDATATYPE["Dirichlet"],
                        writeBCDataHeader,
                        dirArr.name,
                        dirArr.dataType,
                        dirArr.nDimensions,
                        dirArr.dataDimensions,
                        dirArr.dataArr,
                        dirArr.dataArr.shape,
                    )
                    writeBCDataHeader = False

                writeBCDataHeader = True
                for neuArr in dataSet.neumannArrays:
                    libcgns_utils.utils.writebcdata(
                        cg,
                        zoneID,
                        iBC,
                        iDataSet,
                        BCDATATYPE["Neumann"],
                        writeBCDataHeader,
                        neuArr.name,
                        neuArr.dataType,
                        neuArr.nDimensions,
                        neuArr.dataDimensions,
                        neuArr.dataArr,
                        neuArr.dataArr.shape,
                    )
                    writeBCDataHeader = False

        for b2b in self.B2Bs:
            libcgns_utils.utils.writeb2b(
                cg, zoneID, b2b.name, b2b.donorName, b2b.ptRange, b2b.donorRange, b2b.transform
            )

    def writeDimsPlot3d(self, f):
        """Write dimensions to a plot3d file"""
        f.write("%d %d %d\n" % (self.dims[0], self.dims[1], self.dims[2]))

    def writeCoordsPlot3d(self, f):
        """Write coordinates to plot3d file"""
        for iDim in range(3):
            self.coords[:, :, :, iDim].flatten("F").tofile(f, sep="\n", format="%20.15g")
            f.write("\n")

    def scale(self, scaleFact):
        """Scale the coordinates"""
        self.coords *= scaleFact

    def flip(self, axis):
        """Flip coordinates by plane defined by 'axis'"""
        if axis.lower() == "x":
            index = 0
        elif axis.lower() == "y":
            index = 1
        elif axis.lower() == "z":
            index = 2
        self.coords[:, :, :, index] = -self.coords[:, :, :, index]

        # HOWEVER just doing this results in a left-handed block (if
        # the original block was right handed). So we have to also
        # reverse ONE of the indices
        for k in range(self.dims[2]):
            for j in range(self.dims[1]):
                for idim in range(3):
                    self.coords[:, j, k, idim] = self.coords[::-1, j, k, idim]
        # AND we now have to flip the BC's on i-faces since they will
        # now be on the other side:
        for boco in self.bocos:
            if boco.ptRange[0, 0] == boco.ptRange[0, 1] and boco.ptRange[0, 0] == 1:
                boco.ptRange[0, 0] = self.dims[0]
                boco.ptRange[0, 1] = self.dims[0]
            elif boco.ptRange[0, 0] == boco.ptRange[0, 1] and boco.ptRange[0, 0] == self.dims[0]:
                boco.ptRange[0, 0] = 1
                boco.ptRange[0, 1] = 1

    def coarsen(self):
        """Coarsen the block uniformly. We will update the boundary
        conditions and B2B if necessary"""
        # We will coarsen one direction at a time. We do this to check if the block
        # is already 1-cell wide, which can't be coarsened any further
        
        # the new dimensions are half rounded up of the old dimensions
        new_dims = copy.deepcopy(self.dims)
        for i in range(3):
            if self.dims[i] > 2:
                
                if self.dims[i] % 2 == 0:
                    print(f"INFO: unevenly coarsing block {self.name.decode('utf-8', 'ignore')} along dimension {i} (size {self.dims[i]}) ")
                
                # The RHS takes odd numbers to *1/2 rounded up and even numbers to *1/2 
                # for example 
                # old dim: 0 1 2 3 4 5 6 7 8 9
                # new dim: 1 1 1 2 2 3 3 4 4 5
                new_dims[i] = int(numpy.ceil((self.dims[i]) / 2))
                
                
        new_coords = numpy.zeros((new_dims[0], new_dims[1], new_dims[2], 3))
        
        # Loop over all directions
        s = slice(None)
        fine_slicer = [s] * 3
        for idx_dim in range(3):
            # can this direction be coarsened?
            if self.dims[idx_dim] > 2:
                # set the slice in that direction to take every other point
                fine_slicer[idx_dim] = slice(None, None, 2)

        new_coords = self.coords[tuple(fine_slicer)]

        # set the last point to be the same so we don't create any gaps if
        # the number of points in that direction isn't odd
        for idx_dim in range(3):
            end_fine_slicer = copy.deepcopy(fine_slicer)
            end_coarse_slicer = [s] * 3

            end_fine_slicer[idx_dim] = -1
            end_coarse_slicer[idx_dim] = -1

            new_coords[tuple(end_coarse_slicer)] = self.coords[tuple(end_fine_slicer)]

        # coarsen the bc and connectivity for the blk as well
        for boco in self.bocos:
            boco.coarsen()
        for b2b in self.B2Bs:
            b2b.coarsen()

        self.coords = new_coords

        self.dims = new_dims

    def refine(self, axes):
        """Refine the block uniformly. We will also update the
        boundary conditions and B2Bs if necessary"""
        axes = "".join(axes)
        self.coords = libcgns_utils.utils.refine(self.coords, "i" in axes, "j" in axes, "k" in axes)
        self.dims[0] = self.coords.shape[0]
        self.dims[1] = self.coords.shape[1]
        self.dims[2] = self.coords.shape[2]
        for boco in self.bocos:
            boco.refine(axes)

        for b2b in self.B2Bs:
            b2b.refine(axes)

    def section(self, iStart, iEnd, jStart, jEnd, kStart, kEnd):
        self.bocos = []
        self.B2Bs = []
        self.coords = self.coords[iStart - 1 : iEnd, jStart - 1 : jEnd, kStart - 1 : kEnd, :]
        self.dims[0] = self.coords.shape[0]
        self.dims[1] = self.coords.shape[1]
        self.dims[2] = self.coords.shape[2]

    def double2D(self):
        """Double in just the 2D direction"""
        # First find the 2D direction
        for dim_index in range(3):
            if self.dims[dim_index] == 2:
                # Increase the size of the 2D dimension
                new_dimensions = self.dims[:]
                new_dimensions[dim_index] = new_dimensions[dim_index] + 1
                newCoords = numpy.zeros((new_dimensions[0], new_dimensions[1], new_dimensions[2], 3))

                if dim_index == 0:
                    for i in range(self.dims[1]):
                        for j in range(self.dims[2]):
                            newCoords[0, i, j, :] = self.coords[0, i, j, :]
                            newCoords[2, i, j, :] = self.coords[1, i, j, :]
                            newCoords[1, i, j, :] = 0.5 * (newCoords[0, i, j, :] + newCoords[2, i, j, :])
                elif dim_index == 1:
                    for i in range(self.dims[0]):
                        for j in range(self.dims[2]):
                            newCoords[i, 0, j, :] = self.coords[i, 0, j, :]
                            newCoords[i, 2, j, :] = self.coords[i, 1, j, :]
                            newCoords[i, 1, j, :] = 0.5 * (newCoords[i, 0, j, :] + newCoords[i, 2, j, :])
                elif dim_index == 2:
                    for i in range(self.dims[0]):
                        for j in range(self.dims[1]):
                            newCoords[i, j, 0, :] = self.coords[i, j, 0, :]
                            newCoords[i, j, 2, :] = self.coords[i, j, 1, :]
                            newCoords[i, j, 1, :] = 0.5 * (newCoords[i, j, 0, :] + newCoords[i, j, 2, :])

                for boco in self.bocos:
                    for j in range(2):
                        if boco.ptRange[dim_index, j] == 2:
                            boco.ptRange[dim_index, j] = 3
                for b2b in self.B2Bs:
                    for j in range(2):
                        if b2b.ptRange[dim_index, j] == 2:
                            b2b.ptRange[dim_index, j] = 3
                        if b2b.donorRange[dim_index, j] == 2:
                            b2b.donorRange[dim_index, j] = 3

                # Replace previous coordinates
                self.coords = newCoords
                self.dims = new_dimensions[:]

    def _extrudeGetDataOrderAndDIms(self, directionNormal, nSteps):
        """This is a support function that member functions extrude and revolve call"""

        # Note that the self.dims always has data in the first and second
        # slot like it is a xy plane dataset. The third slot always has ones
        # set in readgrid() function. This will we updated.

        if directionNormal == "x":
            # Data given is in yz-plane
            order = [2, 0, 1]
            newDims = [nSteps, self.dims[0], self.dims[1], 3]
        elif directionNormal == "y":
            # Data given is in xz-plane
            order = [0, 2, 1]
            newDims = [self.dims[0], nSteps, self.dims[1], 3]
        elif directionNormal == "z":
            # Data given is in xy-plane
            order = [0, 1, 2]
            newDims = [self.dims[0], self.dims[1], nSteps, 3]
        else:
            raise ValueError(f"Direction normal must be x, y, or z. Input was {directionNormal}.")

        return order, numpy.array(newDims)

    def _extrudeBocoAndAddSymmBoco(self, order, nSteps=2):
        """This is a support function that member functions extrude and revolve call"""

        # Update current BCs
        for boco in self.bocos:
            # Find where we have zeros. That will indicate dimension that has not been updated
            # We only need to check the last row of ptRange because the data actual data is always
            # in the first two rows
            if boco.ptRange[2, 0] == 0:
                boco.ptRange[2, 0] = 1
                boco.ptRange[2, 1] = nSteps

            # Sort based on which dimension we want to extrude in
            boco.ptRange = boco.ptRange[order]

        # for b2b in self.B2Bs:
        #     print(b2b.ptRange)
        #     print(b2b.donorRange)
        #     if b2b.ptRange[2, 0] == 0:
        #           b2b.ptRange[2, 0] = 1
        #           b2b.ptRange[2, 1] = nSteps
        #           b2b.donorRange[2, 0] = 1
        #           b2b.donorRange[2, 1] = nSteps

        #     # Sort based on which dimension we want to extrude in
        #     b2b.ptRange = b2b.ptRange[order]
        #     b2b.donorRange = b2b.donorRange[order]

        # Create 2 new SYMM BCs for this block (each side). This is the plane which the grid was created in
        bocoType = BC["bcsymmetryplane"]
        family = "sym"

        bocoName = "SYMM-{0}".format(0)
        ptRange = numpy.ones((3, 2))
        ptRange[0, 1] = self.dims[0]
        ptRange[1, 1] = self.dims[1]
        ptRange[2, :] = 1

        # Sort based on which dimension we want to extrude in
        ptRange = ptRange[order]

        # Create and add the BC
        self.addBoco(Boco(bocoName, bocoType, ptRange, family))

        bocoName = "SYMM-{0}".format(1)
        ptRange = numpy.ones((3, 2))
        ptRange[0, 1] = self.dims[0]
        ptRange[1, 1] = self.dims[1]
        ptRange[2, :] = nSteps

        # Sort based on which dimension we want to extrude in
        ptRange = ptRange[order]

        # Create and add the BC
        self.addBoco(Boco(bocoName, bocoType, ptRange, family))

    def extrude(self, direction):
        """Extrudes from 2D panar grid to 3D"""

        # Get the data order and new dims
        order, newDims = self._extrudeGetDataOrderAndDIms(direction)

        # Allocate memory for new coordinates
        newCoords = numpy.zeros(newDims)

        # Now copy current coords into new coord array.

        # As for the dims above the coordinates have only data in the first two slots i,j.
        # The actual coordinates stored however are given in the plane specified by the user they are
        # therefore not updated/changed.
        for i in range(self.dims[0]):
            for j in range(self.dims[1]):
                if direction == "x":
                    newCoords[0, i, j, :] = self.coords[i, j, 0, :]
                    newCoords[1, i, j, :] = self.coords[i, j, 0, :]
                    # Update the x-dimension coord with unit length
                    newCoords[1, i, j, 0] = 1.0
                elif direction == "y":
                    newCoords[i, 0, j, :] = self.coords[i, j, 0, :]
                    newCoords[i, 1, j, :] = self.coords[i, j, 0, :]
                    # Update the y-dimension coord with unit length
                    newCoords[i, 1, j, 1] = 1.0
                elif direction == "z":
                    newCoords[i, j, 0, :] = self.coords[i, j, 0, :]
                    newCoords[i, j, 1, :] = self.coords[i, j, 0, :]
                    # Update the z-dimension coord with unit length
                    newCoords[i, j, 1, 2] = 1.0

        # Update the coordinates
        self.coords = newCoords

        # Update current BCs
        self._extrudeBocoAndAddSymmBoco(order)

        # Update the dims. This is done last since the original dims are used above to simplify and reduce code
        self.dims = newDims[:-1]

    def revolve(self, normalDirection, rotationAxis, startAngle, endAngle, nThetas):
        """Revolves a 2D planar grid to create a 3D axisymmetric grid"""

        wedgeAngleRad = numpy.deg2rad(endAngle - startAngle)

        startAngleRad = numpy.deg2rad(startAngle)
        angleRadStep = wedgeAngleRad / (nThetas - 1)

        # Get the data order and new dims
        order, newDims = self._extrudeGetDataOrderAndDIms(normalDirection, nThetas)

        # Allocate memory for new coordinates
        newCoords = numpy.zeros(newDims)

        # Now copy current coords into new coord array.
        for i in range(self.dims[0]):
            for j in range(self.dims[1]):
                for k in range(nThetas):

                    tc = self.coords[i, j, 0, :].copy()

                    angleRad = startAngleRad + angleRadStep * k

                    if normalDirection == "x":
                        if rotationAxis == "y":
                            r = numpy.linalg.norm(tc[[0, 2]])
                            tc[0] = numpy.sin(angleRad) * r
                            tc[2] = numpy.cos(angleRad) * r
                        elif rotationAxis == "z":
                            r = numpy.linalg.norm(tc[0:2])
                            tc[0] = numpy.sin(angleRad) * r
                            tc[1] = numpy.cos(angleRad) * r

                        newCoords[k, i, j, :] = tc

                    elif normalDirection == "y":
                        if rotationAxis == "x":
                            r = numpy.linalg.norm(tc[1:])
                            tc[1] = numpy.sin(angleRad) * r
                            tc[2] = numpy.cos(angleRad) * r
                        elif rotationAxis == "z":
                            r = numpy.linalg.norm(tc[0:2])
                            tc[0] = numpy.cos(angleRad) * r
                            tc[1] = numpy.sin(angleRad) * r

                        newCoords[i, k, j, :] = tc

                    elif normalDirection == "z":
                        if rotationAxis == "x":
                            r = numpy.linalg.norm(tc[1:])
                            tc[2] = numpy.sin(angleRad) * r
                            tc[1] = numpy.cos(angleRad) * r
                        elif rotationAxis == "y":
                            r = numpy.linalg.norm(tc[0, 2])
                            tc[0] = numpy.sin(angleRad) * r
                            tc[2] = numpy.cos(angleRad) * r

                        newCoords[i, j, k, :] = tc

        # Update the coordinates
        # newCoords_swap = newCoords[order]

        self.coords = newCoords

        # Update current BCs
        self._extrudeBocoAndAddSymmBoco(order, nThetas)

        # Update the dims. This is done last since the original dims
        # are used above to simplify and reduce code

        self.dims = newDims[:-1]

    def getSplitBlocks(self):
        """Return a list of blocks that have been split according to
        the self.splits array. This is used for the 'split' operation
        as well as for the 'divide' operation. Boundary information is
        kept but connectivity information is removed"""
        blkList = []
        s = self.splits  # For cleaner code below

        for i in range(len(s[0]) - 1):
            for j in range(len(s[1]) - 1):
                for k in range(len(s[2]) - 1):
                    newCoords = self.coords[
                        s[0][i] - 1 : s[0][i + 1], s[1][j] - 1 : s[1][j + 1], s[2][k] - 1 : s[2][k + 1], :
                    ].copy()

                    dims = [newCoords.shape[0], newCoords.shape[1], newCoords.shape[2]]
                    blk = Block("dummy", dims, newCoords)

                    # Loop over the 6 faces and determine what BC they
                    # came from:

                    for boco in self.bocos:
                        # iLow
                        chkRange = [[s[0][i], s[0][i]], [s[1][j], s[1][j + 1]], [s[2][k], s[2][k + 1]]]

                        if inRange(boco.ptRange, chkRange):
                            blk.addBoco(Boco(boco.name, boco.type, [[1, 1], [1, dims[1]], [1, dims[2]]], boco.family))

                        # iHigh
                        chkRange = [[s[0][i + 1], s[0][i + 1]], [s[1][j], s[1][j + 1]], [s[2][k], s[2][k + 1]]]

                        if inRange(boco.ptRange, chkRange):
                            blk.addBoco(
                                Boco(
                                    boco.name, boco.type, [[dims[0], dims[0]], [1, dims[1]], [1, dims[2]]], boco.family
                                )
                            )

                        # jLow
                        chkRange = [[s[0][i], s[0][i + 1]], [s[1][j], s[1][j]], [s[2][k], s[2][k + 1]]]

                        if inRange(boco.ptRange, chkRange):
                            blk.addBoco(Boco(boco.name, boco.type, [[1, dims[0]], [1, 1], [1, dims[2]]], boco.family))

                        # jHigh
                        chkRange = [[s[0][i], s[0][i + 1]], [s[1][j + 1], s[1][j + 1]], [s[2][k], s[2][k + 1]]]

                        if inRange(boco.ptRange, chkRange):
                            blk.addBoco(
                                Boco(
                                    boco.name, boco.type, [[1, dims[0]], [dims[1], dims[1]], [1, dims[2]]], boco.family
                                )
                            )

                        # kLow
                        chkRange = [[s[0][i], s[0][i + 1]], [s[1][j], s[1][j + 1]], [s[2][k], s[2][k]]]

                        if inRange(boco.ptRange, chkRange):
                            blk.addBoco(Boco(boco.name, boco.type, [[1, dims[0]], [1, dims[1]], [1, 1]], boco.family))

                        # kHigh
                        chkRange = [[s[0][i], s[0][i + 1]], [s[1][j], s[1][j + 1]], [s[2][k + 1], s[2][k + 1]]]

                        if inRange(boco.ptRange, chkRange):
                            blk.addBoco(
                                Boco(
                                    boco.name, boco.type, [[1, dims[0]], [1, dims[1]], [dims[2], dims[2]]], boco.family
                                )
                            )

                    blkList.append(blk)
        return blkList

    def divide(self):
        """Return a list of 8 blocks split derivied from the current
        block. Boundary condition information is kept, but
        connectivity information is removed"""

        # Just add the splits and run getSplitBlocks
        for iDim in range(3):
            self.splits[iDim].append((self.dims[iDim] - 1) // 2 + 1)
            self.splits[iDim].sort()

        return self.getSplitBlocks()

    def removeSymBCs(self):
        """Remove any sym BC's there may be"""
        self.bocos = [boco for boco in self.bocos if not boco.type == BC["bcsymmetryplane"]]

    def extractWallSurfaces(self):
        """Return patches for any surfaces that have BCViscous on them"""
        patches = []
        for boco in self.bocos:
            if isWall(boco.type):
                ptRange = boco.ptRange - 1  # Convert to python ordering
                patches.append(
                    self.coords[
                        ptRange[0, 0] : ptRange[0, 1] + 1,
                        ptRange[1, 0] : ptRange[1, 1] + 1,
                        ptRange[2, 0] : ptRange[2, 1] + 1,
                        :,
                    ].squeeze()
                )
                # Make sure the patch is correctly orientated since we
                # might have left-handed faces. Essentially we have to
                # flip an index on any "high" boundary condition:

                if (
                    (ptRange[0, 0] == ptRange[0, 1] and ptRange[0, 0] + 1 == self.dims[0])
                    or (ptRange[1, 0] == ptRange[1, 1] and ptRange[1, 0] + 1 == 1)
                    or (ptRange[2, 0] == ptRange[2, 1] and ptRange[2, 0] + 1 == self.dims[2])
                ):
                    patches[-1] = patches[-1][::-1, :, :]

        return patches

    def extractSpecifiedSurfaces(self, imin, imax, jmin, jmax, kmin, kmax):
        """Return patches for slices at the six specified indices"""

        # check the indices against the block dimensions and cap as neccessary
        if imin < 0:
            imin = 0
        if jmin < 0:
            jmin = 0
        if kmin < 0:
            kmin = 0
        if imax > self.dims[0] - 1:
            imax = self.dims[0] - 1
        if jmax > self.dims[1] - 1:
            jmax = self.dims[1] - 1
        if kmax > self.dims[2] - 1:
            kmax = self.dims[2] - 1

        patches = []
        # Setup the slice dimensions
        ptRanges = [
            numpy.array([[imin, imin], [jmin, jmax], [kmin, kmax]]),
            numpy.array([[imax, imax], [jmin, jmax], [kmin, kmax]]),
            numpy.array([[imin, imax], [jmin, jmin], [kmin, kmax]]),
            numpy.array([[imin, imax], [jmax, jmax], [kmin, kmax]]),
            numpy.array([[imin, imax], [jmin, jmax], [kmin, kmin]]),
            numpy.array([[imin, imax], [jmin, jmax], [kmax, kmax]]),
        ]

        for i in range(len(ptRanges)):
            ptRange = ptRanges[i]

            patches.append(
                self.coords[
                    ptRange[0, 0] : ptRange[0, 1] + 1,
                    ptRange[1, 0] : ptRange[1, 1] + 1,
                    ptRange[2, 0] : ptRange[2, 1] + 1,
                    :,
                ].squeeze()
            )

            # Make sure the patch is correctly orientated since we
            # might have left-handed faces. Essentially we have to
            # flip an index on any "high" boundary condition:

            if (
                (ptRange[0, 0] == ptRange[0, 1] and ptRange[0, 0] == imax)
                or (ptRange[1, 0] == ptRange[1, 1] and ptRange[1, 0] == jmin)
                or (ptRange[2, 0] == ptRange[2, 1] and ptRange[2, 0] == kmax)
            ):
                patches[-1] = patches[-1][::-1, :, :]

            # Flip all the normals
            patches[-1] = patches[-1][::-1, :, :]

        return patches

    def overwriteFamily(self, face, family):
        """Possibly overwrite the family in the bocos if possible"""
        for boco in self.bocos:
            if self.isFaceInPtRange(face, boco.ptRange):
                boco.family = family

    def overwriteBCs(self, face, bocoType, family, dataSet):
        """Find any BCs on this face and toast them. Note that we *ONLY ALLOW
        ONE BC per face*
        """

        # Check for existing boco and pop if necessary
        face = face.lower()
        pop_list = []
        for index, boco in enumerate(self.bocos):
            # Check if this boco point range matches the face
            if self.isFaceInPtRange(face, boco.ptRange):
                pop_list = pop_list + [index]

        # Pop all bcs in the face
        pop_list.reverse()  # We have to remove the higher indices first
        for index in pop_list:
            self.bocos.pop(index)

        d = self.dims
        if face == "ilow":
            ptRange = [[1, 1, 1], [1, d[1], d[2]]]
        elif face == "ihigh":
            ptRange = [[d[0], 1, 1], [d[0], d[1], d[2]]]
        elif face == "jlow":
            ptRange = [[1, 1, 1], [d[0], 1, d[2]]]
        elif face == "jhigh":
            ptRange = [[1, d[1], 1], [d[0], d[1], d[2]]]
        elif face == "klow":
            ptRange = [[1, 1, 1], [d[0], d[1], 1]]
        elif face == "khigh":
            ptRange = [[1, 1, d[2]], [d[0], d[1], d[2]]]
        else:
            raise ValueError(f"Face must be one of iLow, iHigh, jLow, jHigh, kLow, or kHigh. Input was {face}")

        ptRange = numpy.array(ptRange).T
        self.addBoco(Boco("boco_%d" % self.bocoCounter, BC[bocoType.lower()], ptRange, family, dataSet))
        self.bocoCounter += 1

    def writeBCs(self, blk_num, file_handle):
        """write the bc data to a file"""

        d = self.dims

        for boco in self.bocos:
            # check the point range to see what face it is on
            ptRange = boco.ptRange
            if (ptRange[0] == [1, 1]).all():
                face = "ilow"
            elif (ptRange[0] == [d[0], d[0]]).all():
                face = "ihigh"

            elif (ptRange[1] == [1, 1]).all():
                face = "jlow"
            elif (ptRange[1] == [d[1], d[1]]).all():
                face = "jhigh"

            elif (ptRange[2] == [1, 1]).all():
                face = "klow"
            elif (ptRange[2] == [d[2], d[2]]).all():
                face = "khigh"
            else:
                raise ValueError("Face could not be determined to be one of (iLow, iHigh, jLow, jHigh, kLow, or kHigh)")

            data_arr_str = ""
            for data_arr in boco.dataSets:

                # use the BC dictionary in reverse to find the bc type string
                for bctype in BC:
                    if data_arr.type == BC[bctype]:
                        bctype_str = bctype
                        break

                data_arr_str += " " + data_arr.name.decode("utf-8", "ignore")
                data_arr_str += " " + bctype_str

                if data_arr.dirichletArrays:
                    data_arr_str += " Dirichlet"

                    for d_arr in data_arr.dirichletArrays:
                        data_arr_str += " " + d_arr.name.decode("utf-8", "ignore")
                        data_arr_str += " " + " ".join([f"{data}" for data in d_arr.dataArr])

                if data_arr.neumannArrays:
                    data_arr_str += " Neumann"

                    for n_arr in data_arr.neumannArrays:
                        data_arr_str += " " + n_arr.name
                        data_arr_str += " " + " ".join([f"{data}" for data in n_arr.dataArr])

            # use the BC dictionary in reverse to find the bc type string
            for bctype in BC:
                if boco.type == BC[bctype]:
                    bctype_str = bctype
                    break

            fam_name = boco.family.strip().decode("utf-8", "ignore")
            file_handle.write(f"{blk_num} {face} {bctype_str} {fam_name} {data_arr_str}\n")

    def isFaceInPtRange(self, face, ptRange):
        """
        Identifies if the provided face matches a given point range.

        Parameters
        ----------
        face : str
            Should be one of iLow, iHigh, etc.
        ptRange : array (3,2)
            Point range

        Returns
        -------
        isFound : bool
            Returns True if face is found in the given point range.
        """

        face = face.lower()
        r = ptRange
        isFound = (
            (r[0][0] == r[0][1] == 1 and face == "ilow")
            or (r[0][0] == r[0][1] == self.dims[0] and face == "ihigh")
            or (r[1][0] == r[1][1] == 1 and face == "jlow")
            or (r[1][0] == r[1][1] == self.dims[1] and face == "jhigh")
            or (r[2][0] == r[2][1] == 1 and face == "klow")
            or (r[2][0] == r[2][1] == self.dims[2] and face == "khigh")
        )
        return isFound

    def rebunch(self, spacing, extraCells, nStar):
        """Perform rebunching for this block"""
        from pyspline import Curve

        # ********* WARNING THIS HARD CODED TO K-MIN PLANES *********
        self.dims[2] += extraCells
        if nStar == -1:
            nStar = self.dims[2]

        newNodes = numpy.zeros((self.dims[0], self.dims[1], self.dims[2], 3))
        for i in range(self.dims[0]):
            for j in range(self.dims[1]):

                xx = self.coords[i, j, :, :]
                c = Curve(X=xx, localInterp=True)
                # First get the distance off-wall:
                d = numpy.linalg.norm(self.coords[i, j, 0, :] - self.coords[i, j, 1, :])

                # This is the segment of S we are dealing with:
                sSegment = c.s[0:nStar]

                # Compute the new S0
                s0 = (spacing / d) * c.s[1]
                # Get the newS.
                newS = getS(len(sSegment) + extraCells + 1, s0, sSegment[-1])
                # The final 's' for evaluation
                newS = numpy.hstack([newS, c.s[nStar + 1 :]])
                newNodes[i, j, :, :] = c(newS)

        self.coords = newNodes

    def randomize(self, keepRHS):
        """Randomly reorder the indices in the block. No attempt is made to
        change BCs or B2Bs since these should be deleted already
        """
        flipCount = 0
        if numpy.random.random() > 0.5:
            # We will flip the i-index
            flipCount += 1
            for k in range(self.dims[2]):
                for j in range(self.dims[1]):
                    for idim in range(3):
                        self.coords[:, j, k, idim] = self.coords[::-1, j, k, idim]

        if numpy.random.random() > 0.5:
            flipCount += 1
            # We will flip the j-index
            for k in range(self.dims[2]):
                for i in range(self.dims[0]):
                    for idim in range(3):
                        self.coords[i, :, k, idim] = self.coords[i, ::-1, k, idim]

        if numpy.random.random() > 0.5:
            flipCount += 1
            # We will flip the k-index
            for j in range(self.dims[1]):
                for i in range(self.dims[0]):
                    for idim in range(3):
                        self.coords[i, j, :, idim] = self.coords[i, j, ::-1, idim]

        # So that filps the order of the axis. We can also perform
        # axis swapping.
        if numpy.random.random() > 0.5:

            # Swap X and Y axis
            newCoords = numpy.zeros((self.dims[1], self.dims[0], self.dims[2], 3))
            for k in range(self.dims[2]):
                for idim in range(3):
                    newCoords[:, :, k, idim] = numpy.rot90(self.coords[:, :, k, idim].copy())

            self.dims = list(newCoords.shape[0:3])
            self.coords = newCoords.copy()

        if numpy.random.random() > 0.5:
            # Swap Z and X axis
            newCoords = numpy.zeros((self.dims[2], self.dims[1], self.dims[0], 3))
            for j in range(self.dims[1]):
                for idim in range(3):
                    newCoords[:, j, :, idim] = numpy.rot90(self.coords[:, j, :, idim])

            self.dims = list(newCoords.shape[0:3])
            self.coords = newCoords.copy()

        if numpy.random.random() > 0.5:
            # Swap Y and Z axis
            newCoords = numpy.zeros((self.dims[0], self.dims[2], self.dims[1], 3))
            for i in range(self.dims[0]):
                for idim in range(3):
                    newCoords[i, :, :, idim] = numpy.rot90(self.coords[i, :, :, idim])

            self.dims = list(newCoords.shape[0:3])
            self.coords = newCoords.copy()

        # if the flip count is odd, do one final flip of the j-index
        # to keep the same handed-ness
        if numpy.mod(flipCount, 2) == 1 and keepRHS:
            for k in range(self.dims[2]):
                for j in range(self.dims[1]):
                    for idim in range(3):
                        self.coords[:, j, k, idim] = self.coords[::-1, j, k, idim]

    def symmZero(self, idir):
        for bc in self.bocos:
            if bc.type == BC["bcsymmetryplane"]:
                # 'r' is the range. We need to subtract off -1 from
                # the low end since it was in fortran 1-based ordering
                r = bc.ptRange.copy()
                self.coords[r[0, 0] - 1 : r[0, 1], r[1, 0] - 1 : r[1, 1], r[2, 0] - 1 : r[2, 1], idir] = 0.0

    def symmZeroNoBC(self, idir, tol):

        # Find which nodes are closer than the tolerance from the symmetry plane
        nodeIDs = numpy.where(self.coords[:, :, :, idir] < tol)

        # Zero those nodes
        self.coords[:, :, :, idir][nodeIDs] = 0.0

    def getFaceCoords(self, blockID):
        """Return the list of coordinates on the face as well as its index info"""

        il = self.dims[0]
        jl = self.dims[1]
        kl = self.dims[2]
        nFace = 2 * ((il - 1) * (jl - 1) + (il - 1) * (kl - 1) + (jl - 1) * (kl - 1))

        return libcgns_utils.utils.computefacecoords(self.coords, nFace, blockID)

    def getNumCells(self):
        """Computes and returns the number of cells for this block"""
        return (self.dims[0] - 1) * (self.dims[1] - 1) * (self.dims[2] - 1)

    def getNumNodes(self):
        """Computes and returns the number of nodes for this block"""
        return self.dims[0] * self.dims[1] * self.dims[2]


class Boco(object):

    """Class for storing information related to a boundary condition"""

    def __init__(self, bocoName, bocoType, ptRange, family, bcDataSets=None):
        self.name = bocoName.strip()
        self.type = bocoType
        self.ptRange = ptRange

        if bcDataSets is None:
            self.dataSets = []
        else:
            self.dataSets = bcDataSets

        if family is None or family.strip() == "":
            self.family = "default"
        else:
            self.family = family

    def addBocoDataSet(self, bocoDataSet):
        """Add a boundary condition dataset to this bc"""
        self.dataSets.append(bocoDataSet)

    def coarsen(self):
        """Coarsen the range of the BC"""

        for idim in range(3):

            self.ptRange[idim, 0] = int(numpy.floor((self.ptRange[idim, 1]) / 2)) + 1
            if self.ptRange[idim, 1] > 2:

                # coarsen the data set if it is an array
                if self.dataSets:
                    for data_set in self.dataSets:
                        for dir_arr in data_set.dirichletArrays:
                            if numpy.prod(dir_arr.dataDimensions) == 1:
                                # one value is being used for all points
                                # thus, there is no need to coarsen the data
                                continue

                            slicer = [slice(None)] * 3
                            slicer[idim] = slice(None, None, 2)

                            data_mat = dir_arr.dataArr.reshape(self.ptRange[:, 1])
                            new_data_mat = data_mat[tuple(slicer)]

                            # make sure the last data point is always copied
                            slicer[idim] = -1

                            new_data_mat[tuple(slicer)] = data_mat[tuple(slicer)]

                            dir_arr.dataDimensions[0] = numpy.prod(new_data_mat.shape)

                            dir_arr.dataArr = new_data_mat.flatten()

                self.ptRange[idim, 1] = int(numpy.ceil((self.ptRange[idim, 1]) / 2))

    def refine(self, axes):
        """refine the range of the BC"""
        for i, axis in enumerate(["i", "j", "k"]):
            for j in range(2):
                self.ptRange[i, j] = (self.ptRange[i, j] - 1) * 2 ** (axis in axes) + 1


class BocoDataSet(object):
    """Container class that contains list of data arrays that are associated to a boundary condition"""

    def __init__(self, bocoSetName, bocoDataSetType):
        self.name = bocoSetName.strip()
        self.type = bocoDataSetType  # BC type
        self.dirichletArrays = []
        self.neumannArrays = []

    def addDirichletDataSet(self, dirDataSet):
        self.dirichletArrays.append(dirDataSet)

    def addNeumannDataSet(self, neuDataSet):
        self.neumannArrays.append(neuDataSet)


class BocoDataSetArray(object):
    """Class that contains the actual dataset associated to a boundary condition"""

    def __init__(self, arrayName, dType, nDims, dataDims, dataArr):
        self.name = arrayName.strip()
        self.dataType = dType  # This is the CGNS datatype that was read from the CGNS file.
        self.nDimensions = nDims  # Number of dimensions that the data has
        self.dataDimensions = dataDims  # Number of data points of every dimensions. Array of ints
        self.dataArr = dataArr  # Note that this is a flat 1D numpy array and is float64


class B2B(object):
    """
    Class for storing information related to a Block-to-block or
    (1to1 in cgns speak) connection. More details at http://cgns.github.io/CGNS_docs_current/sids/cnct.html#GridConnectivity1to1.

    Parameters
    ----------
    connectName : str
        Name of the surface patch.

    donorName : str
        Name of the adjacent block (that sits on the other side of the block-to-
        block interface).

    ptRange : array (3,2)
        ptRange contains the subrange of indices that makes up the interface
        patch in the current zone.

    donorRange : array (3,2)
        donorRange contains the interface patch subrange of indices for the
        adjacent block (whose identifier is given by donorName). By
        convention the indices contained in ptRange and donorRange
        refer to vertices.

    transform : array (3)
        Information to produce transformation matrix between ijk axes from one
        block to the other. Each entry, transform[i], gives the axis in the
        donor block that corresponds to the ith axis in the owner block. If the
        blocks are perfectly aligned, transform = [1, 2, 3].
    """

    def __init__(self, connectName, donorName, ptRange, donorRange, transform):
        self.name = connectName.strip()
        self.donorName = donorName.strip()
        self.ptRange = ptRange
        self.donorRange = donorRange
        self.transform = transform

    def coarsen(self):
        """Coarsen the range of the B2B along the specified direction"""
        for idim in range(3):

            donorDir = abs(self.transform[idim]) - 1

            for j in range(2):

                if self.ptRange[idim, j] > 2:
                    self.ptRange[idim, j] = (self.ptRange[idim, j] - 1) // 2 + 1

                if self.donorRange[donorDir, j] > 2:
                    self.donorRange[donorDir, j] = (self.donorRange[donorDir, j] + 1) // 2

    def refine(self, axes):
        """refine the range of the B2B"""
        for i, axis in enumerate(["i", "j", "k"]):
            for j in range(2):
                self.ptRange[i, j] = (self.ptRange[i, j] - 1) * 2 ** (axis in axes) + 1
                self.donorRange[i, j] = (self.donorRange[i, j] - 1) * 2 ** (axis in axes) + 1


# ----------------------------------------
# These are miscellaneous helper functions
# ----------------------------------------
def isWall(bc):
    """Determine if a bc is a wall-type boundary condition"""
    if (
        bc == BC["bcwall"]
        or bc == BC["bcwallinviscid"]
        or bc == BC["bcwallviscous"]
        or bc == BC["bcwallviscousheatflux"]
        or bc == BC["bcwallviscousisothermal"]
    ):
        return True
    else:
        return False


def getS(N, s0, S):
    """Determine the new set of parameters that geometrically fit N
    nodes with the last distance S"""

    # function 'f' is 1 - s0*(1-r^n)/(1-r), s0 is initial ratio and r
    # is the grid ratio.

    # Bisection search:
    a = 1.0 + 1e-8
    b = 4.0

    def f(r):
        s = numpy.zeros(N)
        s[1] = s0
        for i in range(2, N):
            s[i] = s[i - 1] + r * (s[i - 1] - s[i - 2])

        return s[-1]

    fa = S - f(a)

    for _i in range(100):
        c = 0.5 * (a + b)
        ff = S - f(c)
        if abs(ff) < 1e-6:
            break

        if ff * fa > 0:
            a = c
        else:
            b = c
    s = numpy.zeros(N)
    s[1] = s0

    for i in range(2, N):
        s[i] = s[i - 1] + c * (s[i - 1] - s[i - 2])

    return s


def getSplits(ptRange):
    """Return info required to split this face to make it face
    matched"""
    if ptRange[0][0] == ptRange[0][1]:
        splits = [(1, ptRange[1][0]), (1, ptRange[1][1]), (2, ptRange[2][0]), (2, ptRange[2][1])]
    elif ptRange[1][0] == ptRange[1][1]:
        splits = [(0, ptRange[0][0]), (0, ptRange[0][1]), (2, ptRange[2][0]), (2, ptRange[2][1])]
    elif ptRange[2][0] == ptRange[2][1]:
        splits = [(0, ptRange[0][0]), (0, ptRange[0][1]), (1, ptRange[1][0]), (1, ptRange[1][1])]
    return splits


def generalizedCoordDir(iFace):
    """Not really sure how this works..."""
    if iFace in [0, 1]:
        return [0, 1, 2]
    elif iFace in [2, 3]:
        return [1, 2, 0]
    elif iFace in [4, 5]:
        return [0, 2, 1]


def isodd(num):
    """check if a number is odd"""
    return num & 1 and True or False


def getPointRange(iFace, dims):
    """Return the correct point range for face iFace on a block with
    dimensions given in dims"""
    il = dims[0]
    jl = dims[1]
    kl = dims[2]
    if iFace == 0:
        return [[1, il], [1, jl], [1, 1]]
    elif iFace == 1:
        return [[1, il], [1, jl], [kl, kl]]
    elif iFace == 2:
        return [[1, 1], [1, jl], [1, kl]]
    elif iFace == 3:
        return [[il, il], [1, jl], [1, kl]]
    elif iFace == 4:
        return [[1, il], [1, 1], [1, kl]]
    elif iFace == 5:
        return [[1, il], [jl, jl], [1, kl]]


def inRange(ptRange, chkRange):
    """Determine if 'chkRange' fully overlaps with 'ptRange'"""
    val = True
    for iDim in range(3):
        if not (chkRange[iDim][0] >= ptRange[iDim][0] and chkRange[iDim][1] <= ptRange[iDim][1]):
            val = False

    return val


def simpleCart(xMin, xMax, dh, hExtra, nExtra, sym, mgcycle, outFile):
    """
    Generates a Cartesian mesh

    Parameters
    ----------
    xMin : array (3)
        Minimum along each coordinate axis.

    xMax : array (3)
        Maximum along each coordinate axis.

    dh : float OR array (3)
        Approximate edge length of each cell.

    hExtra : float

    nExtra : float

    sym : str OR list
        Axis of symmetry plane, one or more of ('x', 'y', 'z', 'xmin', 'xmax',
        'ymin', 'ymax', 'zmin', 'zmax').

    mgcycle : int OR array (3)
        Number of times mesh should be able to be coarsened for multigrid cycles.

    outFile : str
        Output file name (optional).
    """
    assert len(xMin) == 3
    assert len(xMax) == 3

    if isinstance(dh, float) or isinstance(dh, int):
        dh = [dh] * 3
    else:
        assert len(dh) == 3

    if isinstance(sym, str):
        sym = [sym]

    if isinstance(mgcycle, int):
        mgcycle = [mgcycle] * 3
    else:
        assert len(mgcycle) == 3

    # Now determine how many nodes we need on the inside
    N = numpy.zeros(3, "intc")
    dx = numpy.zeros(3)
    r = numpy.zeros(3)
    Xcart = []
    for iDim in range(3):
        assert isinstance(mgcycle[iDim], int)
        MGFact = 2 ** (mgcycle[iDim] - 1)
        n = int((xMax[iDim] - xMin[iDim]) / dh[iDim])
        n = n // MGFact + 1
        N[iDim] = n * MGFact  # Number of CELLS

        # Compute the *actual* dx's
        dx[iDim] = (xMax[iDim] - xMin[iDim]) / N[iDim]

        # Next we need to find the grid stretch ratios for each
        # direction to satify our requested extra distance.
        r[iDim] = libcgns_utils.utils.calcgridratio(nExtra, dx[iDim], hExtra)

        # Determine if this direction should have a sym plane:
        pos = True
        neg = True

        if ("x" in sym or "xmin" in sym) and iDim == 0:
            neg = False
        if ("y" in sym or "ymin" in sym) and iDim == 1:
            neg = False
        if ("z" in sym or "zmin" in sym) and iDim == 2:
            neg = False
        if ("xmax" in sym) and iDim == 0:
            pos = False
        if ("ymax" in sym) and iDim == 1:
            pos = False
        if ("zmax" in sym) and iDim == 2:
            pos = False

        # Now fill up the cartesian direction
        n = N[iDim]
        iStart = 0
        if neg:
            n += nExtra
            iStart = nExtra

        if pos:
            n += nExtra

        # cordinates for this dimension
        x = numpy.zeros(n + 1)

        # First coordinate is at iStart:
        x[iStart] = xMin[iDim]

        # Add remainder of the uniform part:
        for i in range(N[iDim]):
            x[i + 1 + iStart] = x[i + iStart] + dx[iDim]

        # Add neg part if necessary:
        if neg:
            for i in range(nExtra):
                x[iStart - 1 - i] = x[iStart - i] - r[iDim] * (x[iStart - i + 1] - x[iStart - i])
        if pos:
            iStart = iStart + N[iDim]
            for i in range(nExtra):
                x[iStart + i + 1] = x[iStart + i] + r[iDim] * (x[iStart + i] - x[iStart + i - 1])

        Xcart.append(x)

    # Allocate coordinates block
    shp = [Xcart[0].shape[0], Xcart[1].shape[0], Xcart[2].shape[0]]
    X = numpy.zeros((shp[0], shp[1], shp[2], 3))

    print("Grid Dimensions:", shp)
    print("Grid Ratios:", r)
    # Write grid coordinates
    Xx, Xy, Xz = numpy.meshgrid(Xcart[0], Xcart[1], Xcart[2], indexing="ij")
    X[:, :, :, 0] = Xx
    X[:, :, :, 1] = Xy
    X[:, :, :, 2] = Xz

    if outFile is not None:
        # Open a new CGNS file and write if necessary:
        cg = libcgns_utils.utils.openfile(outFile, CG_MODE_WRITE, 3)

        # Write a Zone to it
        zoneID = libcgns_utils.utils.writezone(cg, "cartesian", shp)

        # Write mesh coordinates
        libcgns_utils.utils.writecoordinates(cg, zoneID, X)

        # CLose file
        libcgns_utils.utils.closefile(cg)

    return X, dx


# def normal_direction(iFace1, iFace2):
#     """Normal direction is positive if iFace1 and iFace two are of
#     opposite oddity, even if they are the same oddity"""
#     isOdd1 = isodd(iFace1)
#     isOdd2 = isodd(iFace2)

#     if isOdd1 is True and isOdd2 is True:
#         return -1
#     if isOdd1 is False and isOdd2 is False:
#         return -1

#     # otherwise:
#     return 1

# -----------------------------------------------------------------
# These functions perform operations that return new 'Grid' objects
# -----------------------------------------------------------------
def readGrid(fileName):
    """Internal routine to return a 'grid' object that contains all
    the information that is in the file 'fileName'"""

    inFile = libcgns_utils.utils.openfile(fileName, CG_MODE_READ, 3)
    cellDim = libcgns_utils.utils.getgriddimension(inFile)
    nBlock = libcgns_utils.utils.getnblocks(inFile)
    nIterations, nArrays = libcgns_utils.utils.getconvinfo(inFile)

    newGrid = Grid()

    # Assign the fileName as the grid name. We need to remove that path
    # and the file extension.
    newGrid.name = os.path.splitext(os.path.basename(fileName))[0]

    for iBlock in range(1, nBlock + 1):
        zoneName, dims, nBoco, nB2B = libcgns_utils.utils.getblockinfo(inFile, iBlock)

        if cellDim == 2:
            dims[2] = 1
        coords = libcgns_utils.utils.getcoordinates(inFile, iBlock, dims[0], dims[1], dims[2])
        blk = Block(zoneName, dims, coords)

        for iBoco in range(1, nBoco + 1):
            # Get the BCs
            bocoName, bocoType, ptRange, family, nDataSets = libcgns_utils.utils.getbcinfo(
                inFile, iBlock, iBoco, cellDim
            )
            bc = Boco(bocoName, bocoType, ptRange, family)

            # Get the BCDataSets
            if nDataSets != 0:
                # Loop over all the datasets for this BC
                for iBocoDataSet in range(1, nDataSets + 1):

                    (
                        bocoDatasetName,
                        bocoDataSetType,
                        nDirichletArrays,
                        nNeumannArrays,
                    ) = libcgns_utils.utils.getbcdatasetinfo(inFile, iBlock, iBoco, iBocoDataSet)
                    bcDSet = BocoDataSet(bocoDatasetName, bocoType)

                    def getBocoDataSetArray(flagDirNeu, iDir):
                        # Get data information
                        (
                            dataArrayName,
                            dataType,
                            nDimensions,
                            dataDimensionVector,
                        ) = libcgns_utils.utils.getbcdataarrayinfo(
                            inFile, iBlock, iBoco, iBocoDataSet, iDir, flagDirNeu
                        )

                        # Create a flat array for the data
                        # Note we make it float64 although it can contain integers.
                        nDataArr = numpy.prod(dataDimensionVector)
                        dataArr = numpy.zeros(nDataArr, dtype=numpy.float64, order="F")

                        # Get the data. Note the dataArr is populated when the routine exits
                        libcgns_utils.utils.getbcdataarray(
                            inFile, iBlock, iBoco, iBocoDataSet, iDir, flagDirNeu, dataArr, nDataArr
                        )

                        # Create a BocoDataSetArray object and return
                        return BocoDataSetArray(dataArrayName, dataType, nDimensions, dataDimensionVector, dataArr)

                    if nDirichletArrays > 0:
                        # Loop over Dirichlet data and get the actual data
                        for iDir in range(1, nDirichletArrays + 1):

                            # Get the data set
                            bcDSetArr = getBocoDataSetArray(BCDATATYPE["Dirichlet"], iDir)

                            # Append a BocoDataSetArray to the datasets
                            bcDSet.addDirichletDataSet(bcDSetArr)

                        # Append the Dirichlet BC dataset to the BC
                        bc.addBocoDataSet(bcDSet)

                    if nNeumannArrays > 0:
                        # Loop over Neumann data sets
                        for iDir in range(1, nNeumannArrays + 1):

                            # Get the data set
                            bcDSetArr = getBocoDataSetArray(BCDATATYPE["Neumann"], iDir)

                            # Append a BocoDataSetArray to the datasets
                            bcDSet.addNeumannDataSet(bcDSetArr)

                        # Append the Neumann BC dataset to the BC
                        bc.addBocoDataSet(bcDSet)

            blk.addBoco(bc)

        for iB2B in range(1, nB2B + 1):
            connectName, donorName, ptRange, donorRange, transform = libcgns_utils.utils.getb2binfo(
                inFile, iBlock, iB2B
            )
            blk.addB2B(B2B(connectName, donorName, ptRange, donorRange, transform))

        newGrid.addBlock(blk)

    # Read convergence history if available
    if nIterations > 0:
        for arrayID in range(nArrays):
            # Read array
            arrayName, arrayData = libcgns_utils.utils.getconvarray(inFile, nIterations, arrayID + 1)

            # Remove blank spaces
            arrayName = arrayName.strip()

            # Store results in the newGrid.convArray dictionary
            newGrid.addConvArray(arrayName, arrayData)

    libcgns_utils.utils.closefile(inFile)

    # Store grid dimension
    newGrid.cellDim = cellDim

    return newGrid


def convertPlot3d(plot3dFile, cgnsFile):
    """Read a multiblock, fortran big endiend unformatted plot3d. This
    routine is necessary becuase the supplied plot3d_to_cgns converter
    from the cgnslib doesn't always work properly.
    """
    # Full conversion is done in fortran.
    libcgns_utils.utils.convertplot3d(plot3dFile, cgnsFile)


def mirrorGrid(grid, axis, tol):
    """Method that takes a grid and mirrors about the axis. Boundary
    condition information is retained if possible"""

    # First make sure the grid is face matched:
    grid.split([])

    # Now copy original blocks
    newGrid = Grid()
    for blk in grid.blocks:
        blk.removeSymBCs()
        blk.B2Bs = []
        newGrid.addBlock(blk)

        mirrorBlk = copy.deepcopy(blk)
        mirrorBlk.flip(axis)
        newGrid.addBlock(mirrorBlk)

    # Now rename the blocks and redo-connectivity
    newGrid.renameBlocks()
    newGrid.renameBCs()
    newGrid.connect(tol)

    return newGrid


def divideGrid(grid):
    """Method that takes a grid and generates a new grid with 8 times
    as many blocks"""
    newGrid = Grid()
    for blk in grid.blocks:
        newBlks = blk.divide()
        for nblk in newBlks:
            newGrid.addBlock(nblk)

    # Now rename the blocks and redo-connectivity
    newGrid.renameBlocks()
    newGrid.renameBCs()
    newGrid.connect()

    return newGrid


def splitGrid(grid, splitFile):
    """Method that takes a grid and propagates any splits using
    connectivity information. This is a rewrite of the original
    Fortran implementation that is quite a bit simpler due to Python"""
    # Split the current grid
    extraSplits = []
    if splitFile is not None:
        f = open(splitFile, "r")
        for line in f:
            aux = line.split()
            extraSplits.append([int(aux[0]), int(aux[1]), int(aux[2])])
        f.close()
    grid.split(extraSplits=extraSplits)

    # New grid
    newGrid = Grid()
    for blk in grid.blocks:
        newBlks = blk.getSplitBlocks()
        for nblk in newBlks:
            newGrid.addBlock(nblk)

    # # Now rename the blocks, bcs and redo-connectivity
    newGrid.renameBlocks()
    newGrid.renameBCs()
    newGrid.connect()

    return newGrid


def mergeGrid(grid):
    """Method that that takes a grid with block to block connections and
    merges as many blocks as possible, reducing the total number of
    blocks in the mesh"""

    def fullFace(blk, ptRange):

        # Face size of the patch:
        fSize = abs(ptRange[:, 1] - ptRange[:, 0])

        fullFace = True
        for i in range(3):
            if fSize[i] != 0:
                if fSize[i] != blk.dims[i] - 1:
                    fullFace = False

        return fullFace

    def faceID(ptRange, blk):
        if ptRange[0, 0] == ptRange[0, 1] == 1:
            faceID = -1
        elif ptRange[0, 0] == ptRange[0, 1] == blk.dims[0]:
            faceID = 1

        elif ptRange[1, 0] == ptRange[1, 1] == 1:
            faceID = -2
        elif ptRange[1, 0] == ptRange[1, 1] == blk.dims[1]:
            faceID = 2

        elif ptRange[2, 0] == ptRange[2, 1] == 1:
            faceID = -3
        elif ptRange[2, 0] == ptRange[2, 1] == blk.dims[2]:
            faceID = 3

        return faceID

    # Outer Iterative Loop
    cont = True
    iteration = 0
    while cont:

        # First create a mapping of the blocks from the name to the index
        zoneMap = {}

        for i in range(len(grid.blocks)):
            blk = grid.blocks[i]
            zoneMap[blk.name] = i

        blockUsed = numpy.zeros(len(grid.blocks), "intc")
        newBlocks = []
        # Loop over each block:
        for i in range(len(grid.blocks)):
            blk = grid.blocks[i]

            # We haven't used this block yet:
            if blockUsed[i] == 0:

                # Loop over the B2B of this block:
                for b2b in blk.B2Bs:
                    otherIndex = zoneMap[b2b.donorName]
                    otherBlk = grid.blocks[otherIndex]

                    # Determine if this B2B is a full patch on my block
                    # *and* the other block
                    if (
                        fullFace(blk, b2b.ptRange)
                        and fullFace(otherBlk, b2b.donorRange)
                        and blockUsed[otherIndex] == 0
                        and i != otherIndex
                    ):

                        print("Merging:", i + 1, otherIndex + 1)

                        # Great! These block match. Let's make the new
                        # block

                        # Transform describes how the current block,
                        # blk, is related to the other block, otherBlk
                        transform = b2b.transform

                        # We need to determine which face we are
                        # dealing with (iLow, iHigh, etc) on each block

                        face = faceID(b2b.ptRange, blk)

                        dims = blk.dims.copy()
                        dims[abs(face) - 1] += otherBlk.dims[abs(transform[abs(face) - 1]) - 1] - 1
                        newCoords = numpy.zeros((dims[0], dims[1], dims[2], 3))

                        # Now transform the other coordinates to make
                        # them conform with the existing
                        # block. Essentially what we need to do
                        # perform a series of operations to convert
                        # the transfrom matrix back to [1, 2, 3],
                        # which means the block completely line up.
                        otherCoords = otherBlk.coords

                        tmp = abs(transform)
                        # Now swap the axes until we get 1,2,3. There
                        # are only 6 cases, so just do them individually
                        if tmp[0] == 1 and tmp[1] == 2 and tmp[2] == 3:
                            # Nothing to do:
                            pass
                        elif tmp[0] == 1 and tmp[1] == 3 and tmp[2] == 2:
                            otherCoords = numpy.swapaxes(otherCoords, 2, 1)

                        elif tmp[0] == 2 and tmp[1] == 1 and tmp[2] == 3:
                            otherCoords = numpy.swapaxes(otherCoords, 0, 1)

                        elif tmp[0] == 2 and tmp[1] == 3 and tmp[2] == 1:
                            otherCoords = numpy.swapaxes(otherCoords, 0, 1)
                            otherCoords = numpy.swapaxes(otherCoords, 1, 2)

                        elif tmp[0] == 3 and tmp[1] == 1 and tmp[2] == 2:
                            otherCoords = numpy.swapaxes(otherCoords, 0, 2)
                            otherCoords = numpy.swapaxes(otherCoords, 1, 2)

                        elif tmp[0] == 3 and tmp[1] == 2 and tmp[2] == 1:
                            otherCoords = numpy.swapaxes(otherCoords, 2, 0)

                        # This flips any axis not in the right order
                        if transform[0] < 0:
                            otherCoords = otherCoords[::-1, :, :, :]
                        if transform[1] < 0:
                            otherCoords = otherCoords[:, ::-1, :, :]
                        if transform[2] < 0:
                            otherCoords = otherCoords[:, :, ::-1, :]

                        if face > 0:
                            # blk then otherBlk
                            if face == 1:
                                newCoords[0 : blk.dims[0], :, :, :] = blk.coords
                                newCoords[blk.dims[0] - 1 :, :, :, :] = otherCoords

                            elif face == 2:
                                newCoords[:, 0 : blk.dims[1], :, :] = blk.coords
                                newCoords[:, blk.dims[1] - 1 :, :, :] = otherCoords

                            elif face == 3:
                                newCoords[:, :, 0 : blk.dims[2], :] = blk.coords
                                newCoords[:, :, blk.dims[2] - 1 :, :] = otherCoords

                        else:
                            # otherBlk then blk
                            if face == -1:
                                newCoords[0 : dims[0] - blk.dims[0] + 1, :, :, :] = otherCoords
                                newCoords[dims[0] - blk.dims[0] :, :, :] = blk.coords

                            elif face == -2:
                                newCoords[:, 0 : dims[1] - blk.dims[1] + 1, :, :] = otherCoords
                                newCoords[:, dims[1] - blk.dims[1] :, :, :] = blk.coords

                            elif face == -3:
                                newCoords[:, :, 0 : dims[2] - blk.dims[2] + 1, :] = otherCoords
                                newCoords[:, :, dims[2] - blk.dims[2] :, :] = blk.coords

                        # Create the new block
                        newBlk = Block("doesNotMatter", dims, newCoords)

                        # Now deal with the boundary conditions. These
                        # need to be corrected depending on how the
                        # blocks get added.
                        offset = numpy.zeros(3, "intc")
                        if face == 1:
                            offset[0] = blk.dims[0] - 1
                        elif face == 2:
                            offset[1] = blk.dims[1] - 1
                        elif face == 3:
                            offset[2] = blk.dims[2] - 1
                        elif face == -1:
                            offset[0] = dims[0] - blk.dims[0]
                        elif face == -2:
                            offset[1] = dims[1] - blk.dims[1]
                        elif face == -3:
                            offset[2] = dims[2] - blk.dims[2]

                        # Add all the bocos from the first block:
                        for boco in blk.bocos:
                            if face > 0:
                                # blk then otherBlk. BCs go in as is:
                                pass
                            else:
                                # BCs have to offset:
                                boco.ptRange[:, 0] += offset
                                boco.ptRange[:, 1] += offset

                            newBlk.addBoco(boco)

                        # Add all the bocos from the second
                        # block. This is tricky since we need to
                        # offset and potentially reorient them .
                        for boco in otherBlk.bocos:

                            tmp = boco.ptRange.copy()
                            if face > 0:
                                # blk then otherBlk. BCs need to be increemented by offset.
                                for idim in range(3):
                                    jdim = transform[idim]

                                    if jdim > 0:
                                        # Other block dim +ve, (ie in the same direction)
                                        boco.ptRange[idim, 0] = offset[idim] + tmp[jdim - 1, 0]
                                        boco.ptRange[idim, 1] = offset[idim] + tmp[jdim - 1, 1]
                                    else:
                                        # Other block dim -ve, (ie in the opposite direction)
                                        jdim = -jdim
                                        boco.ptRange[idim, 0] = (
                                            offset[idim] + otherBlk.dims[jdim - 1] - tmp[jdim - 1, 1] + 1
                                        )
                                        boco.ptRange[idim, 1] = (
                                            offset[idim] + otherBlk.dims[jdim - 1] - tmp[jdim - 1, 0] + 1
                                        )

                            else:
                                # otherBlk then blk. BCs need to be transformed only
                                for idim in range(3):
                                    jdim = transform[idim]

                                    if jdim > 0:
                                        # Other block dim +ve, (ie in the same direction)
                                        boco.ptRange[idim, 0] = tmp[jdim - 1, 0]
                                        boco.ptRange[idim, 1] = tmp[jdim - 1, 1]
                                    else:
                                        # Other block dim -ve, (ie in the opposite direction)
                                        jdim = -jdim
                                        boco.ptRange[idim, 0] = otherBlk.dims[jdim - 1] - tmp[jdim - 1, 1] + 1
                                        boco.ptRange[idim, 1] = otherBlk.dims[jdim - 1] - tmp[jdim - 1, 0] + 1

                            newBlk.addBoco(boco)

                        # Add the new block to the list
                        newBlocks.append(newBlk)

                        # Flag the two existing blocks as deleted
                        blockUsed[i] = 1
                        blockUsed[otherIndex] = 1

                        # We can't do anything else on this block so skip the rest of the b2b loop
                        break

                    # end if (matching)
                # end for (b2b loop)
            # end if (used check)
        # end for (block loop)

        if len(newBlocks) == 0:
            cont = False

        # Now loop back through the grids appending the blocks we
        # haven't used to "newBlocks"
        for i in range(len(grid.blocks)):
            if blockUsed[i] == 0:
                newBlocks.append(grid.blocks[i])

        # Set the new blocks
        grid.blocks = newBlocks
        print("New number of blocks:", len(grid.blocks))

        # Rename the blocks and remove any B2B info since it will all
        # be recomputed:
        grid.renameBlocks()
        grid.B2Bs = []
        grid.connect()
        iteration += 1
    # end for (outer loop)

    return grid


def combineGrids(grids, useOldNames=False):

    """Method that takes in a list of grids and returns a new grid object
    containing all zones from each grid. The blocks are renamed as
    there could (almost most certainly) be conflicts between the zone
    names in the different grids. This also means that the block to
    block connectivities need to be updated based on the new zone names.

    If useOldNames=True we will preserve the domain names after merging
    the blocks, otherwise, we will replace all names by the filenames.
    """

    # Create a dictionary to contain grid objects with their names
    # as the corresponding keys
    gridDict = {}
    for grid in grids:

        # Create a dictionary of copies so the original grids are not modified
        gridDict[grid.name] = copy.deepcopy(grid)

    # Alphabetically sort the list of grid object names
    nameList = list(sorted(gridDict.keys()))

    # First determine the total number of blocks
    newGrid = Grid()

    # Loop through each name for the grid objects and add their blocks
    # to the newGrid object

    for name in nameList:

        # Get the grid object corresponding to this name
        grid = gridDict[name]

        # Mapping of the old names to the new names
        zoneMap = {}

        nBlock = 0

        # Loop over the blocks and obtain the name mapping
        for blk in grid.blocks:
            nBlock += 1
            if not useOldNames:
                blockName = name
            else:
                try:
                    blockName = blk.name.split(".")[0]
                except TypeError:
                    blockName = blk.name.decode().split(".")[0]
            newName = blockName + f".{nBlock:05}"
            zoneMap[blk.name] = newName
            blk.name = newName

        # Now loop back over the blocks, replacing the donorName using
        # the map we defined above
        for blk in grid.blocks:
            for b2b in blk.B2Bs:
                b2b.donorName = zoneMap[b2b.donorName]

        # Now add the new processed blocks to our newgrid
        newGrid.blocks.extend(grid.blocks)

    return newGrid


def explodeGrid(grid, kMin=False):
    """Method that takes one multiblock grid and returns a list of grids, each
    one containing a single block of the original multiblock grid
    """

    # Reduce block size if just Kmin face is needed
    if kMin:
        for blk in grid.blocks:
            blk.dims[2] = 1
            blk.coords = blk.coords[:, :, 0:1, :]

    # Initialize list of single block grids
    gridList = []

    # Add one block to each grid
    for blk in grid.blocks:
        # Initialize new grid
        newGrid = Grid()
        # Add a single block to this grid
        newGrid.addBlock(blk)
        # Now rename the blocks, bcs and redo-connectivity, only if we have full mesh
        if kMin is False:
            newGrid.renameBlocks()
            newGrid.renameBCs()
            newGrid.connect()
        # Append this new grid to the grids list
        gridList.append(newGrid)

    # return list of grids
    return gridList


def explodeByZoneName(grid):
    """Method that takes one multiblock grid and returns a list of grids, each
    containing all zones that have similar naming.
    """

    # Initialize list of grids
    gridList = []
    nameList = []

    # Loop over each block in the input grid and obtain all zone names
    for blk in grid.blocks:
        name = blk.name.split(".")[:-1]
        name = ".".join(name)
        nameList.append(name)

    # Get only the unique zone names
    nameList = list(sorted(set(nameList)))

    gridDict = {}

    # Add keys corresponding to each name to the dict
    for name in nameList:
        gridDict[name] = Grid()
        gridDict[name].name = "_" + name

    # Add the blocks to the corresponding grid
    for blk in grid.blocks:
        name = blk.name.split(".")[:-1]
        name = ".".join(name)
        gridDict[name].addBlock(blk)

    # Loop over all keys and add the grid to the output list
    for name in nameList:

        # Now rename the blocks, bcs and redo-connectivity, only if we have full mesh
        gridDict[name].renameBlocks(actualName=True)
        gridDict[name].renameBCs()
        gridDict[name].connect()

        # Append this new grid to the grids list
        gridList.append(gridDict[name])

    # return list of grids
    return gridList, nameList


def write_tecplot_file(filename, title, variable_names, data_points):

    """
    Auxiliary function that writes tecplot files
    """

    # Open the data file
    fid = open(filename, "w")

    # Write the title
    fid.write("title = " + title + "\n")

    # Write tha variable names
    varnames_commas = ",".join(variable_names)  # Merge names in a single string separated by commas
    fid.write("variables = " + varnames_commas + ",\n")  # Write to the file

    # Write data points
    if type(data_points) is list:  # Check if user provided a list
        for point in data_points:
            str_points = [str(x) for x in point]  # Covert each entry to string
            str_points = " ".join(str_points)  # Merge all entries in a single string separated by whitespace
            fid.write(str_points + "\n")  # Write to file
    else:  # The user probably provided a numpy array
        for index in range(data_points.shape[0]):
            str_points = [str(x) for x in data_points[index, :]]
            str_points = " ".join(str_points)  # Merge all entries in a single string separated by whitespace
            fid.write(str_points + "\n")  # Write to file

    # Close file
    fid.close()
