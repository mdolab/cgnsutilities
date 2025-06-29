"""
This is the new gateway program to all of the cgns_utils.

Run cgns_utils -help to get a list of all available options. The basic
idea is as follows::

                                                  | write new file
    read cngs file -> Do some operations on it -> |     .or.
                                                  | write modified file

Developed by Dr. Gaetan K. W. Kenway

"""
import argparse
import os
import pickle
import shutil
import sys
import tempfile

import numpy as np

from cgnsutilities.cgnsutilities import (
    Block,
    Boco,
    Grid,
    combineGrids,
    convertPlot3d,
    divideGrid,
    explodeByZoneName,
    explodeGrid,
    libcgns_utils,
    mergeGrid,
    mirrorGrid,
    readGrid,
    simpleCart,
    splitGrid,
    write_tecplot_file,
)

# set width of printing for line wrap
os.environ["COLUMNS"] = "120"


def get_parser():
    # List out all of the possible options here.
    parser = argparse.ArgumentParser(prog="cgns_utils")

    subparsers = parser.add_subparsers(help="Choose one of the listed operations to perform", dest="mode")

    # ------------- Options for 'scale' mode --------------------
    p_scale = subparsers.add_parser("scale", help="Scale a grid by a constant factor")
    p_scale.add_argument("gridFile", help="Name of input CGNS file")
    p_scale.add_argument("scale", help="scale factor", type=float)
    p_scale.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------- Options for 'flip' mode --------------------
    p_flip = subparsers.add_parser("flip", help="Flip a grid about a plane defined by an axis")
    p_flip.add_argument("gridFile", help="Name of input CGNS file")
    p_flip.add_argument("axis", help="Flip the mesh about plane defined by axis: 'x', 'y', 'z'")
    p_flip.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------- Options for 'coarsen' mode --------------------
    p_coarsen = subparsers.add_parser("coarsen", help="Coarsen a grid uniformly")
    p_coarsen.add_argument("gridFile", help="Name of input CGNS file")
    p_coarsen.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------- Options for 'refine' mode --------------------
    p_refine = subparsers.add_parser("refine", help="Refine a grid uniformly")
    p_refine.add_argument("gridFile", help="Name of input CGNS file")
    p_refine.add_argument("outFile", nargs="?", default=None, help="Optional output file")
    p_refine.add_argument(
        "--axes",
        nargs="+",
        help="Refine mesh only along specified axis or axes (default: %(default)s)",
        default=["i", "j", "k"],
    )

    # ------------- Options for 'extractSurface' mode --------------------
    p_extract = subparsers.add_parser("extractSurface", help="Extract a wall surface from file")
    p_extract.add_argument("gridFile", help="Name of input CGNS file")
    p_extract.add_argument("surfFile", help="Name of plot3d surface file")

    # ------------- Options for 'extractSpecifiedSurface' mode --------------------
    p_extract_spec = subparsers.add_parser(
        "extractSpecifiedSurface", help="Extract a surface from a specified set of layers in a cgns block"
    )
    p_extract_spec.add_argument("gridFile", help="Name of input CGNS file")
    p_extract_spec.add_argument("surfFile", help="Name of plot3d surface file")
    p_extract_spec.add_argument("blockID", help="cgns block number to extract from")
    p_extract_spec.add_argument("imin", help="lower i bound,use 0-based numbering")
    p_extract_spec.add_argument("imax", help="upper i bound,use 0-based numbering")
    p_extract_spec.add_argument("jmin", help="lower j bound,use 0-based numbering")
    p_extract_spec.add_argument("jmax", help="upper j bound,use 0-based numbering")
    p_extract_spec.add_argument("kmin", help="lower k bound,use 0-based numbering")
    p_extract_spec.add_argument("kmax", help="upper k bound,use 0-based numbering")

    # ------------- Options for 'mirror' mode --------------------
    p_mirror = subparsers.add_parser(
        "mirror", help="Mirror a grid about a plane defined by an axis. This doubles the grid size"
    )
    p_mirror.add_argument("gridFile", help="Name of input CGNS file")
    p_mirror.add_argument("axis", help="Mirror about plane defined by axis: 'x', 'y', 'z'")
    p_mirror.add_argument("tol", nargs="?", default=1e-12, help="Tolerance for node merge")
    p_mirror.add_argument("useOldNames", nargs="?", default=False, help="Whether to use old block names")
    p_mirror.add_argument("surface", nargs="?", default=False, help="Whether this is a surface or volume grid")
    p_mirror.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------- Options for 'split' mode --------------------
    p_split = subparsers.add_parser(
        "split", help="Face-match a grid. If the grid is already faced matched, this will have no effect"
    )
    p_split.add_argument("gridFile", help="Name of input CGNS file")
    p_split.add_argument("outFile", nargs="?", default=None, help="Optional output file")
    p_split.add_argument(
        "--splitFile",
        nargs="?",
        default=None,
        help="""Add additional splits specified in split file. Each
                line must contain a block index (1 based), idim (1, 2, or 3),
                and a 1-based index of the block to split at""",
    )

    # ------------- Options for 'merge' mode --------------------
    p_merge = subparsers.add_parser(
        "merge",
        help="Automatically merge as many blocks as possible. Boundary conditions and family information is kept.",
    )
    p_merge.add_argument("gridFile", help="Name of input CGNS file")
    p_merge.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------- Options for 'connect' mode --------------------
    p_connect = subparsers.add_parser(
        "connect", help="Determine the block-to-block connectivity information for a point-matched grid"
    )
    p_connect.add_argument("gridFile", help="Name of input CGNS file")
    p_connect.add_argument("tol", nargs="?", default=1e-12, help="Tolerance for node merge")
    p_connect.add_argument(
        "--connectSelf",
        help="only check for connection on-block (periodic type)",
        action="store_true",
        dest="connectSelf",
        default=False,
    )
    p_connect.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------- Options for 'divide' mode --------------------
    p_divide = subparsers.add_parser("divide", help="Divide all blocks in the grid into 8 sub-blocks")
    p_divide.add_argument("gridFile", help="Name of input CGNS file")
    p_divide.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------- Options for 'autoBC' mode --------------------
    p_bc = subparsers.add_parser(
        "autoBC", help="Try to determine boundary conditions for blocks. Only suitable for external flow applications."
    )
    p_bc.add_argument("gridFile", help="Name of input CGNS file")
    p_bc.add_argument("sym", help="Normal for possible symmetry plane.", choices=["x", "y", "z"])
    p_bc.add_argument("radius", help="Radius of sphere containing bodies", type=float, default=10.0)
    p_bc.add_argument("outFile", nargs="?", default=None, help="Optional output file")
    p_bc.add_argument("--xOffset", nargs="?", default=0.0, type=float, help="x-coordinate of sphere origin")
    p_bc.add_argument("--yOffset", nargs="?", default=0.0, type=float, help="y-coordinate of sphere origin")
    p_bc.add_argument("--zOffset", nargs="?", default=0.0, type=float, help="z-coordinate of sphere origin")

    # ------------ Options for 'overwriteFamilies' mode --------------------
    p_fam = subparsers.add_parser(
        "overwriteFamilies", help="Overwrite family information", formatter_class=argparse.RawTextHelpFormatter
    )
    p_fam.add_argument("gridFile", help="Name of input CGNS file")
    p_fam.add_argument(
        "familyFile",
        help="""File containing additional family information. The file must consist of one or more lines contaning the following data:
<blockID> <faceID> <family>
where:
    blockID - is the block index *IN 1 BASED NUMBERING*
    faceID  - one of iLow, iHigh jLow, jHigh, kLow, or kHigh
    family  - the family name.

To find blockID of any mesh using Tecplot,
1. Load the mesh with Advanced options > One Tecplot zone per non-poly CGNS zone/solution
2. Use the Zone Number for the blockID

Examples:
    7 kLow wing
    4 jHigh sym
""",
    )

    p_fam.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'writeSubfaceFamily' mode --------------------
    p_fam = subparsers.add_parser("writeSubfaceFamily", help="""Overwrite the family information on a subface.""")
    p_fam.add_argument("gridFile", help="Name of inputCGNS file")
    p_fam.add_argument(
        "familyFile",
        help="""file containing data for the new families and face division.
                Format is 1st line: 1-based blockID, 2nd line: {ilow, ihigh, etc}, subsequent lines,
                one per line: p_extractnge (as 6 ints seperated by commas, not spaces), newFamilyName""",
    )
    p_fam.add_argument(
        "outFile",
        nargs="?",
        default=None,
        help="""Optional output file. ***NOTE*** It is highly recommended that an output file
                is specified as this method will overwrite existing boundary conditions on a face, and it is
                up to the user to supply subfaces which sufficiently replace it.""",
    )

    # ------------ Options for 'copyFamilyInfo' mode --------------------
    p_fam = subparsers.add_parser("copyFamilyInfo", help="Copy family information from two otherwise identical grids")
    p_fam.add_argument("gridFile", help="Name of CGNS file to which family information is to be copied")
    p_fam.add_argument("sourceFile", help="Name of output CGNS file which contains family information")
    p_fam.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'removeBC' mode --------------------
    p_rem = subparsers.add_parser("removeBC", help="Remove all BC")
    p_rem.add_argument("gridFile", help="Name of input CGNS file")
    p_rem.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'overwriteBCFamilyWithBC' mode --------------------
    p_sub = subparsers.add_parser(
        "overwriteBCFamilyWithBC",
        help="Overwrite boundary conditions based on family name",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_sub.add_argument("gridFile", help="Name of input CGNS file")
    p_sub.add_argument("familyName", help="The BC family to overwrite")
    p_sub.add_argument("newBCType", help="The new boundary condition to apply")
    p_sub.add_argument("outFile", nargs="?", default=None, help="Optional output file")
    p_sub.add_argument(
        "--blockIDs",
        type=int,
        nargs="+",
        default=None,
        help="The 1-based indices of the blocks to overwrite. By default, BCs are overwritten on all blocks.",
    )

    # ------------ Options for 'overwriteBC' mode --------------------
    p_sub = subparsers.add_parser(
        "overwriteBC", help="Overwrite boundary condition information", formatter_class=argparse.RawTextHelpFormatter
    )
    p_sub.add_argument("gridFile", help="Name of input CGNS file")
    bcFile_txt = """
The file must consist of one or more lines contaning the following data:
<blockID> <faceID> <BCType> <family> [dataset]

where:
    blockID - is the block index *IN 1 BASED NUMBERING*
    faceID  - one of iLow, iHigh jLow, jHigh, kLow or kHigh
    BCType  - one of the supported CGNS boundary conditions. See below for supported types
    family  - the family name.

Supported BC types are : bcfarfield, bcsymmetryplane bcwall, bcwallinviscid, bcwallviscous
bcwallviscousheatflux, bcwallviscousisothermal, bcoutflow, bcoutflowsubsonic
bcinflow, bcinflowsubsonic, bcinflowsupersonic

Optionally, additional datasets may be specified. These
can be used to set additional boundary condition data.
The format of the dataset line is as follows:
<BCSetName> <BCSetType> <DirNeuArr> <DataArrName1> <dataArr1>, ..., <DataArrNameN> <dataArrN>

where:
    BCSetName     - bc dataset name
    BCSetType     - bc dataset type. This is in most cases the same type as the BCType specified
    DirNeuArr     - can have one of two options: Dirichlet or Neumann
    DataArrNameN  - name of first property specified. This can be a range of things. Refer to ICEM or ADflow for supported BC properties
    dataArrN      - the actual data for the property. Either a scalar or a flattened nodal array. If an array is passed, the solver will convert the 1D array to the (possibly) 2D BC face.

To find blockID of any mesh using Tecplot,
1. Load the mesh with Advanced options > One Tecplot zone per non-poly CGNS zone/solution
2. Use the Zone Number for the blockID

Examples:

    7 kLow bcwallviscous wing
    4 jHigh bcsymmetryplane sym
    5 khigh bcoutflowsubsonic turb_inlet BCDataSet_1 BCInFlowSubsonic Dirichlet PressureStagnation 1234.0 TemperatureStagnation 4556.0
"""
    p_sub.add_argument(
        "bcFile",
        help="File containing additional bc info." + bcFile_txt,
    )
    p_sub.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'writebcinfo' mode --------------------
    p_sub = subparsers.add_parser(
        "writebcinfo",
        help="Writes boundary condition information to a file.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_sub.add_argument("gridFile", help="Name of input CGNS file")
    p_sub.add_argument(
        "bcOutFile",
        default=None,
        help="A file containing bc info." + bcFile_txt,
    )

    # ------------ Options for 'rebunch' mode --------------------
    p_bunch = subparsers.add_parser("rebunch", help="Rebunch offwall spacing (experimental)")
    p_bunch.add_argument("gridFile", help="Name of input CGNS file")
    p_bunch.add_argument("spacing", help="The desired off-wall spacing", type=float)
    p_bunch.add_argument("outFile", nargs="?", default=None, help="Optional output file")
    p_bunch.add_argument(
        "--extraCells",
        help="Number of additional cells to use in re-bunching. *SHOULD BE A MG NUMBER*.",
        type=int,
        default=0,
    )
    p_bunch.add_argument("--nodes", help="Only rebunch the first 'nodes' in the offwall direction", type=int, default=1)

    # ------------ Options for 'cgns2plot3d' mode --------------------
    p3d = subparsers.add_parser("cgns2plot3d", help="Convert a cgns file to a plain plot3d file")
    p3d.add_argument("gridFile", help="Name of input CGNS file")
    p3d.add_argument("plot3dFile", help="Name of output plot3d file")

    # ------------ Options for 'plot3dtocgns' mode --------------------
    p3dtoc = subparsers.add_parser(
        "plot3d2cgns",
        help="""Convert a multiblock, unformatted fortran, big-endian, multiblock plot3d file to a plain
                cgns file. This specific format is widely used at NASA and Boeing.""",
    )
    p3dtoc.add_argument("plot3dFile", help="Name of input plot3d file")
    p3dtoc.add_argument("gridFile", help="Name of output CGNS file")

    # ------------ Options for 'randomize' mode --------------------
    p_ran = subparsers.add_parser("randomize", help="Randomize the block orientation and order. Useful for testing.")
    p_ran.add_argument("gridFile", help="Name of input CGNS file")
    p_ran.add_argument(
        "seed",
        type=int,
        default=0,
        help="Seed for random generator. Specifying a seed will make process deterministic.",
    )
    p_ran.add_argument("--keepRHS", help="Keep right hand coordinate system", action="store_true", default=False)
    p_ran.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'reorder' mode --------------------
    p_reorder = subparsers.add_parser(
        "reorder",
        help="""Sort blocks in an alpha-numerical order. It can also add extra digits
                to the integers at the end of the block names to facilitate ordering.""",
    )
    p_reorder.add_argument("gridFile", help="Name of input CGNS file")
    p_reorder.add_argument(
        "intDigits",
        type=int,
        default=5,
        help="""Number of digits used for the integers. When CGNSlib generates a CGNS file
                (when converting from a plot3d file, for instance), it does not add extra digits to the integers
                when naming zones. This becomes a problem when you have more than 10 zones because the ordering will be:
                Zone1, Zone11, Zone12, ..., Zone19, Zone2, Zone21, ...
                This method will add extra digits to the zone names to give the correct ordering.""",
    )
    p_reorder.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'symmZero' mode --------------------
    p_sym = subparsers.add_parser("symmZero", help="Hard-zero any nodes on symmetry plane BCs.")
    p_sym.add_argument("gridFile", help="Name of input CGNS file")
    p_sym.add_argument("sym", help="Normal for possible symmetry plane.", choices=["x", "y", "z"])
    p_sym.add_argument("family", help="Family of the boundary condition.", default=None)
    p_sym.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'symmZeroNoBC' mode --------------------
    p_symnobc = subparsers.add_parser(
        "symmZeroNoBC",
        help="Hard-zero any nodes within a given tolerance of the symmetry plane. BCs are not taken into account.",
    )
    p_symnobc.add_argument("gridFile", help="Name of input CGNS file")
    p_symnobc.add_argument("sym", help="Normal for possible symmetry plane.", choices=["x", "y", "z"])
    p_symnobc.add_argument(
        "--tol",
        help="Distance tolerance to bring nodes to symmetry plane (default: %(default)s)",
        type=float,
        default=1.0e-5,
    )
    p_symnobc.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'timeCombine' mode  --------------------
    p_tc = subparsers.add_parser(
        "timeCombine", help="Combine cgns files from time accurate simulation into unsteady tecplot file."
    )
    p_tc.add_argument("baseName", help="baseName of the files. Use %%d to denote the counter.")
    p_tc.add_argument("outFile", nargs="?", default=None, help="Output file name. If not given, unsteady.plt is used")

    # ------------ Options for 'double2D' mode  --------------------
    p_dd = subparsers.add_parser("double2D", help="Take a 2d mesh one cell wide and make it two cells wide.")
    p_dd.add_argument("gridFile", help="Name of input CGNS file")
    p_dd.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'combine' mode  --------------------
    p_dd = subparsers.add_parser("combine", help="Take 2 or more cgns files and combine them into a single file.")
    p_dd.add_argument("gridFiles", metavar="files", type=str, nargs="+", help="Name of CGNS files to combine")
    p_dd.add_argument("outFile", type=str, help="Output CGNS file name")

    # ------------ Options for 'removeBlocks' mode  --------------------
    p_rm = subparsers.add_parser(
        "removeBlocks",
        help="""Remove blocks from a cgns file. The user should ensure that the final mesh
                is still valid in terms of boundary conditions and connectivities.""",
    )
    p_rm.add_argument("gridFile", help="Name of input CGNS file")
    p_rm.add_argument(
        "blockIDs",
        metavar="files",
        type=int,
        nargs="+",
        help="IDs (integers) of the blocks that will be removed. The integers should be 1-indexed",
    )
    p_rm.add_argument("outFile", type=str, help="Output CGNS file name")

    # ------------ Options for 'explode' mode  --------------------
    p_exp = subparsers.add_parser(
        "explode", help="Take one multiblock cgns file and explodes it into single-block cgns files."
    )
    p_exp.add_argument("gridFile", type=str, help="Name of input multiblock CGNS file")
    p_exp.add_argument(
        "outFile",
        nargs="?",
        default=None,
        help="""Optional reference to name output files. An integer will be added to the end.
                If none is given, the input filename will be used as reference.
                All connectivity information between different blocks is lost in this step, only
                internal connectivity remains.""",
    )

    # ------------ Options for 'explodeKmin' mode  --------------------
    p_expkmin = subparsers.add_parser(
        "explodeKmin",
        help="Take one multiblock cgns file and explodes it into single-block plot3d files that contains only the K=1 faces.",
    )
    p_expkmin.add_argument("gridFile", type=str, help="Name of input multiblock CGNS file")
    p_expkmin.add_argument(
        "outFile",
        nargs="?",
        default=None,
        help="""Optional reference to name output files. An integer will be added to the end.
                if none is given, the input filename will be used as reference.""",
    )

    # ------------ Options for 'explodeByZoneName' mode  --------------------
    p_expkmin = subparsers.add_parser(
        "explodeByZoneName",
        help="""Take one multiblock cgns file and explode it into multiple multiblock
                cgns files based on the zone name from the blocks.""",
    )
    p_expkmin.add_argument("gridFile", type=str, help="Name of input multiblock CGNS file")

    # ------------ Options for 'cartesian' mode --------------------
    p_cart = subparsers.add_parser(
        "cartesian", help="Generates a background cartesian mesh", formatter_class=argparse.RawTextHelpFormatter
    )
    p_cart.add_argument("gridFile", help="Name of input CGNS file")
    p_cart.add_argument(
        "cartFile",
        help="""File containing background mesh info. The file must consist of
4 lines contaning the following data:

<extensionXneg> <extensionYneg> <extensionZneg>
<extensionXpos> <extensionYpos> <extensionZpos>
<numNodesX> <numNodesY> <numNodesZ>
<weightGRX> <weightGRY> <weightGRZ>

where:
    extension is the distance of the cartesian box
    face to the corresponding bounding box face divided by the
    bounding box length. We need 2 values of extension per
    direction as we have two parallel faces for each one of them.

    numNodes is the number of nodes that should be used along the
    edges of the cartesian mesh. If you want one symmetry plane
    in the z direction, for instance, you need to set one of the
    extensionZ values to 0. If you want two symmetry planes in
    the z direction, (e.g. to run a 2D case) you need to set both
    extensionZ values to 0.

    weightGR are values between 0.0 and 1.0 used to balance edge
    growth ratio and cell volume resolution mismatch during the
    optimization. If weightGR = 0, the optimizer will not care
    about the growth ratios at the farfield and will just try
    to match the bounding box resolution. If weightGR = 1, the
    optimizer will not care about the bounding box resolution
    and will just try to get an uniform growth ratio. This
    results in an uniform mesh.

example:
    10 10 0
    10 10 10
    65 65 65
    0.1 0.1 0.1
""",
    )
    p_cart.add_argument(
        "outFile",
        help="""Name of output CGNS file. The output file contains only one cartesian block.
The input mesh is not included and BCs are applied.""",
    )

    # ------------ Options for 'simpleCart' mode --------------------
    p_sub = subparsers.add_parser(
        "simpleCart",
        help="""Generates a background cartesian mesh containing two
regions: one internal region with uniform spacing dh
and an extended region where mesh sizes grow following
a geometric progression. The internal region is
defined by the bounding box of the mesh present in gridFile.""",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_sub.add_argument("gridFile", help="Name of input CGNS file")
    p_sub.add_argument(
        "dh",
        help="""Uniform spacing size. If a single number is provided,
the same spacing is applied to all three dimensions (x,y,z).
But the user can also provide a set of three spacings
separated by commas, for example: 0.1,0.2,0.3.""",
        type=str,
    )
    p_sub.add_argument("hExtra", help="Length of the extended region", type=float)
    p_sub.add_argument("nExtra", help="Number of nodes of the extended region", type=int)
    p_sub.add_argument(
        "sym",
        help="""Normal for possible sym planes.
Possible options are xmin, ymin, zmin, xmax, ymax, zmax.
The user can also define a combination of them using commas,
for example: zmin,zmax. This is useful to generate 1-cell wide
meshes for 2D cases (remember to set dh to the span length
(mesh width) and mgcycle to 0).
""",
        type=str,
    )
    p_sub.add_argument("mgcycle", help="Minimum MG cycle to enforce", type=int)
    p_sub.add_argument("outFile", help="Name of output CGNS file")

    # ------------ Options for 'explicitCart' mode --------------------
    p_sub = subparsers.add_parser(
        "explicitCart",
        help="""Generates a background cartesian mesh containing two
regions: one internal region with uniform spacing dh
and an extended region where mesh sizes grow following
a geometric progression. The internal region is
defined by the x,y,z coordinates given as explicit arguments.""",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p_sub.add_argument("xmin", type=float, help="min x coordinate of the internal region")
    p_sub.add_argument("ymin", type=float, help="min y coordinate of the internal region")
    p_sub.add_argument("zmin", type=float, help="min z coordinate of the internal region")
    p_sub.add_argument("xmax", type=float, help="max x coordinate of the internal region")
    p_sub.add_argument("ymax", type=float, help="max y coordinate of the internal region")
    p_sub.add_argument("zmax", type=float, help="max z coordinate of the internal region")
    p_sub.add_argument(
        "dh",
        help="""Uniform spacing size. If a single number is provided,
the same spacing is applied to all three dimensions (x,y,z).
But the user can also provide a set of three spacings
separated by commas, for example: 0.1,0.2,0.3.""",
        type=str,
    )
    p_sub.add_argument("hExtra", help="Length of the extended region", type=float)
    p_sub.add_argument("nExtra", help="Number of nodes of the extended region", type=int)
    p_sub.add_argument(
        "sym",
        help="""Normal for possible sym planes.
Possible options are xmin, ymin, zmin, xmax, ymax, zmax.
The user can also define a combination of them using commas,
for example: zmin,zmax. This is useful to generate 1-cell wide
meshes for 2D cases (remember to set dh to the span length
(mesh width) and mgcycle to 0).
""",
        type=str,
    )
    p_sub.add_argument("mgcycle", help="Minimum MG cycle to enforce", type=int)
    p_sub.add_argument("outFile", help="Name of output CGNS file")

    # ------------ Options for 'translate' mode  --------------------
    p_t = subparsers.add_parser("translate", help="Translate a grid.")
    p_t.add_argument("gridFile", help="Name of input CGNS file")
    p_t.add_argument("dx", help="x-displacement", type=float)
    p_t.add_argument("dy", help="y-displacement", type=float)
    p_t.add_argument("dz", help="z-displacement", type=float)
    p_t.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'rotate' mode  --------------------
    p_rot = subparsers.add_parser("rotate", help="Rotate a grid around a given direction.")
    p_rot.add_argument("gridFile", help="Name of input CGNS file")
    p_rot.add_argument("vx", help="x-component of the rotation axis", type=float)
    p_rot.add_argument("vy", help="y-component of the rotation axis", type=float)
    p_rot.add_argument("vz", help="z-component of the rotation axis", type=float)
    p_rot.add_argument("theta", help="rotation angle [deg]", type=float)
    p_rot.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'autoOversetBC' mode
    p_obc = subparsers.add_parser(
        "autoOversetBC",
        help="Automatically generate connectivity and boundary conditions"
        "for an overset near field mesh generated by pyHyp. It assumes the surface is a "
        "BCWallViscous and the outer boundary is a BCOverset condition."
        "Only used with pyHyp hyperbolically generated meshes.",
    )
    p_obc.add_argument("gridFile", help="Name of input CGNS file")
    p_obc.add_argument("sym", help="Normal for possible symmetry plane.", choices=["x", "y", "z", "none"])
    p_obc.add_argument(
        "--connectSelf",
        help="only check for connection on-block (periodic type)",
        action="store_true",
        dest="connectSelf",
        default=False,
    )
    p_obc.add_argument("--tol", type=float, default=1e-12, help="Tolerance for connect")
    p_obc.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'autoNearfieldBC' mode
    p_anf = subparsers.add_parser(
        "autoNearfieldBC",
        help="Automatically generate connectivity and boundary conditions"
        "for an overset near field mesh with possible symmetry plane.",
    )
    p_anf.add_argument("gridFile", help="Name of input CGNS file")
    p_anf.add_argument("sym", help="Normal for possible symmetry plane.", choices=["x", "y", "z", "none"])
    p_anf.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'autoFarfieldBC' mode
    p_anf = subparsers.add_parser(
        "autoFarfieldBC",
        help="Automatically generate connectivity and boundary conditions"
        "for an overset farfield mesh with possible symmetry plane.",
    )
    p_anf.add_argument("gridFile", help="Name of input CGNS file")
    p_anf.add_argument("sym", help="Normal for possible symmetry plane.", choices=["x", "y", "z"])
    p_anf.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'fillOpenBCs' mode
    p_fbc = subparsers.add_parser(
        "fillOpenBCs",
        help="Adds a given BC to the faces that are not face-matched and also that do not have any previously-assigned BCs.",
    )
    p_fbc.add_argument("gridFile", help="Name of input CGNS file")
    p_fbc.add_argument(
        "bocoType",
        help="""Boundary condition type. Supported types are:
                bcfarfield, bcsymmetryplane, bcwall, bcwallinviscid, bcwallviscous, bcwallviscousheatflux,
                bcwallviscousisothermal, bcoutflow, bcoutflowsubsonic, bcoutflowsupersonic, bcinflow, bcinflowsubsonic,
                bcinflowsupersonic and bcoverset""",
    )
    p_fbc.add_argument("famName", help="Family name for the new BCs.")
    p_fbc.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'extractConv' mode  --------------------
    p_conv = subparsers.add_parser(
        "extractConv",
        help="Reads the convergence history node of a CGNS file and saves the data in a pickle or tecplot file.",
    )
    p_conv.add_argument("gridFile", type=str, help="Name of input CGNS file.")
    p_conv.add_argument("outType", choices=["pickle", "tecplot"], help="The type of convergence data output file.")
    p_conv.add_argument(
        "outFile",
        nargs="?",
        default=None,
        help="The convergence data will be saved to this filename. If none is given, the grid filename will be used as reference.",
    )

    # ------------ Options for 'include' mode
    p_inc = subparsers.add_parser(
        "include", help="Write a new file including only the zones specified the given numbers/ranges."
    )
    p_inc.add_argument("gridFile", help="Name of input CGNS file")
    p_inc.add_argument(
        "rangeSpec",
        help="""Range to extract. Comma separated list. Ranges can given like 6-8. Must be 1 based.
                Example: rangeSpec=4,5,9-16,19""",
    )
    p_inc.add_argument("outFile", help="Output file")

    # ------------ Options for 'section' mode
    p_sec = subparsers.add_parser(
        "section",
        help="For cgns files with 1 domain ONLY, write a subsection of the zone. Boundary conditions/B2Bs are deleted.",
    )
    p_sec.add_argument("gridFile", help="Name of input CGNS file")
    p_sec.add_argument("iStart", type=int)
    p_sec.add_argument("iEnd", type=int)
    p_sec.add_argument("jStart", type=int)
    p_sec.add_argument("jEnd", type=int)
    p_sec.add_argument("kStart", type=int)
    p_sec.add_argument("kEnd", type=int)
    p_sec.add_argument("outFile", nargs="?", default=None, help="Optional output file")

    # ------------ Options for 'info' mode
    p_info = subparsers.add_parser("info", help="Print some metrics for the CGNS file.")
    p_info.add_argument("gridFile", help="Name of input CGNS file")

    # ------------ Options for 'extrude' mode
    p_p2D = subparsers.add_parser(
        "extrude",
        help="""Takes a true 2D mesh (planar) and extrude it in one direction to make
                it a 3D mesh, one cell wide. This method assumes that BCs are already set in the CGNS file.
                BCs are retained and symmetry BCs are applied on planar surfaces.""",
    )
    p_p2D.add_argument(
        "gridFile",
        help="Name of input CGNS file. Note that the planar grid should not have symmetry BCs set on the plane.",
    )
    p_p2D.add_argument("direction", help="Direction which to extrude the grid", choices=["x", "y", "z"])
    p_p2D.add_argument("outFile", nargs="?", default=None, help="Optional output CGNS file")

    # ------------ Options for 'revolve' mode
    p_p2R = subparsers.add_parser(
        "revolve",
        help="""Takes a true 2D mesh (planar) and reloves it about specified axis to make
                a 3D axisymmetric mesh, one cell wide. This method assumes that BCs are already set in the CGNS file.
                BCs are retained and symmetry BCs are applied on planar surfaces. Output should be a wedge shape.""",
    )
    p_p2R.add_argument(
        "gridFile",
        help="Name of input CGNS file. Note that the planar grid should not have symmetry BCs set on the plane.",
    )
    p_p2R.add_argument(
        "normalDirection",
        help="""This is the direction in which the plane normal points in.
                Example: If supplied data is in xz plane, the normal points in y""",
        choices=["x", "y", "z"],
    )
    p_p2R.add_argument(
        "axis",
        help="Axis which to rotate about Example: If supplied data is in xz plane, you would give either x or z",
        choices=["x", "y", "z"],
    )
    p_p2R.add_argument(
        "startAngle", type=float, help="Rotation starting angle given in degrees. Typically this is a small quantity."
    )
    p_p2R.add_argument(
        "endAngle", type=float, help="Rotation ending angle given in degrees. Typically this is a small quantity."
    )
    p_p2R.add_argument("nThetas", type=int, help="number of angular locations to put points", default=2)
    p_p2R.add_argument("outFile", nargs="?", default=None, help="Optional output CGNS file")

    # ------------ Options for 'testBlock' mode
    p_test = subparsers.add_parser(
        "testBlock", help="Creates a single block mesh with specified dimensions. Used to quicky generate a test grid."
    )
    p_test.add_argument("nx", help="Number of nodes in x", type=int)
    p_test.add_argument("ny", help="Number of nodes in y", type=int)
    p_test.add_argument("nz", help="Number of nodes in z", type=int)
    p_test.add_argument("outFile", help="Name of output file")

    # ------------- Options for 'blockSizes' mode --------------------
    p_blockSizes = subparsers.add_parser("blockSizes", help="Print the sizes of each block in the mesh")
    p_blockSizes.add_argument("gridFile", help="Name of input CGNS file")

    # return the parser
    return parser


def main():
    parser = get_parser()
    # Get the arguments we need!
    args = parser.parse_args()

    # -------------------------------------------
    #         Selection of the task
    # -------------------------------------------

    # The time combine is special. First we generate the list of files we
    # need to deal with.
    if args.mode == "timeCombine":
        # Get the directory name where the baseName is:
        path = os.path.dirname(os.path.abspath(args.baseName))

        # Get the full list of files in this directory:
        allFiles = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
        files = []

        parts = args.baseName.split("%d")
        maxLength = 0
        for f in allFiles:
            if (parts[0] == "" or parts[0] in f) and (parts[1] == "" or parts[1] in f):
                # Make sure there is a .cgns in there somwhere
                if ".cgns" in f:
                    files.append(f)
                    maxLength = max(maxLength, len(f))
                    files = sorted(files)

        if args.outFile is None:
            outFile = "unsteady.plt"
        else:
            outFile = args.outFile

        # Now we make a character array of the file names, and hand if off to
        # fortran for all the actual reading/writing.
        fileNames = np.zeros((len(files), 256), "c")
        for i in range(len(files)):
            fileNames[i, 0 : len(files[i])] = files[i]

        libcgns_utils.utils.time_combine(fileNames, outFile)
        sys.exit(0)

    if args.mode == "testBlock":
        nx = args.nx
        ny = args.ny
        nz = args.nz
        X = np.zeros((nx, ny, nz, 3))
        Xcart = []
        Xcart.append(np.linspace(0, 1, nx))
        Xcart.append(np.linspace(0, 1, ny))
        Xcart.append(np.linspace(0, 1, nz))
        Xx, Xy, Xz = np.meshgrid(Xcart[0], Xcart[1], Xcart[2], indexing="ij")
        X[:, :, :, 0] = Xx
        X[:, :, :, 1] = Xy
        X[:, :, :, 2] = Xz
        b = Block("domain.00001", [nx, ny, nz], X)
        # Add bocos so we can run it:
        b.addBoco(Boco("iMin", "bcfarfield", [[1, 1], [1, ny], [1, nz]], "far"))
        b.addBoco(Boco("iMax", "bcfarfield", [[nx, nx], [1, ny], [1, nz]], "far"))
        b.addBoco(Boco("jMin", "bcsymmetryplane", [[1, nx], [1, 1], [1, nz]], "sym"))
        b.addBoco(Boco("jMax", "bcsymmetryplane", [[1, nx], [ny, ny], [1, nz]], "sym"))
        b.addBoco(Boco("kMin", "bcwallviscous", [[1, nx], [1, ny], [1, 1]], "wall"))
        b.addBoco(Boco("kMax", "bcfarfield", [[1, nx], [1, ny], [nz, nz]], "far"))
        g = Grid()
        g.addBlock(b)
        g.writeToCGNS(args.outFile)
        sys.exit(0)

    # The 'combine' function is done first sicne it is the only function
    # that reads multiple cgns files.
    if args.mode == "combine":
        grids = []
        for fName in args.gridFiles:
            grid = readGrid(fName)
            grids.append(grid)

        combinedGrid = combineGrids(grids)
        combinedGrid.writeToCGNS(args.outFile)

        # This task is now finished
        sys.exit(0)

    if args.mode == "explicitCart":
        # This task doesn't have args.gridFile so do it first
        xMin = [args.xmin, args.ymin, args.zmin]
        xMax = [args.xmax, args.ymax, args.zmax]

        # Change dh to a list of floats, since we had
        # to define it as a string in arg_parser to avoid
        # issues with variable number of inputs.
        # (The user could give a single dh, or a set of three values)
        dh = list(map(float, args.dh.split(",")))
        if len(dh) == 1:
            dh = dh[0]

        # A similar procedure is necessary for the symmetry planes argument
        sym = args.sym.split(",")

        simpleCart(xMin, xMax, dh, args.hExtra, args.nExtra, sym, args.mgcycle, args.outFile)

        sys.exit(0)

    if args.mode == "plot3d2cgns":
        # This doesn't read a cgns file so do this first too.
        convertPlot3d(args.plot3dFile, args.gridFile)
        sys.exit(0)

    # Get the current working grid 'curGrid' by reading the input
    curGrid = readGrid(args.gridFile)

    # The following are "special" and done first since they do not
    # have a CGNS output.

    if args.mode == "extractSurface":
        curGrid.extractSurface(args.surfFile)
        sys.exit(0)

    if args.mode == "extractSpecifiedSurface":
        curGrid.extractSpecifiedSurface(
            args.surfFile, args.blockID, args.imin, args.imax, args.jmin, args.jmax, args.kmin, args.kmax
        )
        sys.exit(0)

    if args.mode == "cgns2plot3d":
        curGrid.writePlot3d(args.plot3dFile)
        sys.exit(0)

    if args.mode == "blockSizes":
        curGrid.printBlockInfo()
        sys.exit(0)

    # Determine if we have an output file:
    try:
        if args.outFile is None:
            # Determine where to put a file:
            dirpath = tempfile.mkdtemp()

            # Define a temp output file
            outFileName = os.path.join(dirpath, "tmp.cgns")
        else:
            outFileName = args.outFile
    except Exception:
        outFile = None

    # Perform one of the following actions:
    if args.mode == "flip":
        curGrid.flip(args.axis)

    elif args.mode == "scale":
        curGrid.scale(args.scale)

    elif args.mode == "mirror":
        curGrid = mirrorGrid(curGrid, args.axis, args.tol, useOldNames=args.useOldNames, surface=args.surface)

    elif args.mode == "coarsen":
        curGrid.coarsen()

    elif args.mode == "refine":
        curGrid.refine(args.axes)

    elif args.mode == "split":
        curGrid = splitGrid(curGrid, args.splitFile)

    elif args.mode == "merge":
        curGrid = mergeGrid(curGrid)

    elif args.mode == "connect":
        if args.connectSelf:
            curGrid.connectSelfOnly(args.tol)
        else:
            curGrid.connect(args.tol)

    elif args.mode == "divide":
        curGrid = divideGrid(curGrid)

    elif args.mode == "autoBC":
        curGrid.autoBC(args.radius, args.sym, [args.xOffset, args.yOffset, args.zOffset])

    elif args.mode == "overwriteFamilies":
        curGrid.overwriteFamilies(args.familyFile)

    elif args.mode == "writeSubfaceFamily":
        curGrid.writeSubfaceFamily(args.familyFile)

    elif args.mode == "copyFamilyInfo":
        sourceGrid = readGrid(args.sourceFile)
        curGrid.copyFamilyInfo(sourceGrid)

    elif args.mode == "overwriteBCFamilyWithBC":
        curGrid.overwriteBCFamilyWithBC(args.familyName, args.newBCType, args.blockIDs)

    elif args.mode == "overwriteBC":
        print("overwriting bc with file", args.bcFile)
        curGrid.overwriteBCs(args.bcFile)

    elif args.mode == "writebcinfo":
        curGrid.writeBCs(args.bcOutFile)
        sys.exit(0)

    elif args.mode == "removebc":
        curGrid.removeBCs()

    elif args.mode == "rebunch":
        curGrid.rebunch(args.spacing, args.extraCells, args.nodes)

    elif args.mode == "randomize":
        curGrid.randomize(args.seed, args.keepRHS)

    elif args.mode == "reorder":
        curGrid.reorder(args.intDigits)

    elif args.mode == "symmZero":
        curGrid.symmZero(args.sym, args.family)

    elif args.mode == "symmZeroNoBC":
        curGrid.symmZeroNoBC(args.sym, args.tol)

    elif args.mode == "double2D":
        curGrid.double2D()

    elif args.mode == "removeBlocks":
        curGrid.removeBlocks(args.blockIDs)

    elif args.mode == "cartesian":
        found_overset = False
        for block in curGrid.blocks:
            for boco in block.bocos:
                if boco.internalType == "bcoverset":
                    found_overset = True
        if found_overset:
            curGrid.cartesian(args.cartFile, args.outFile)
        else:
            print("The CGNS file has no overset boundary conditions")
        sys.exit(0)

    elif args.mode == "simpleCart":
        # Change dh to a list of floats, since we had
        # to define it as a string in arg_parser to avoid
        # issues with variable number of inputs.
        # (The user could give a single dh, or a set of three values)
        dh = list(map(float, args.dh.split(",")))
        if len(dh) == 1:
            dh = dh[0]

        # A similar procedure is necessary for the symmetry planes argument
        sym = args.sym.split(",")

        curGrid.simpleCart(dh, args.hExtra, args.nExtra, sym, args.mgcycle, args.outFile)
        sys.exit(0)

    elif args.mode == "translate":
        curGrid.translate(args.dx, args.dy, args.dz)

    elif args.mode == "rotate":
        curGrid.rotate(args.vx, args.vy, args.vz, args.theta)

    elif args.mode == "autoOversetBC":
        curGrid.autoOversetBC(args.sym, args.connectSelf, args.tol)

    elif args.mode == "autoNearfieldBC":
        curGrid.autoNearfieldBC(args.sym)

    elif args.mode == "autoFarfieldBC":
        curGrid.autoFarfieldBC(args.sym)

    elif args.mode == "fillOpenBCs":
        curGrid.fillOpenBCs(args.bocoType, args.famName)

    elif args.mode == "include":
        toWrite = []
        for spec in args.rangeSpec.split(","):
            if "-" in spec:
                tmp = spec.split("-")
                start = int(tmp[0])
                end = int(tmp[1])
            else:
                start = int(spec)
                end = int(spec)
            for i in range(start, end + 1):
                toWrite.append(i)
        toWrite = np.unique(toWrite)
        toWrite.sort()
        curGrid.writeToCGNSSelected(args.outFile, toWrite)
        sys.exit(0)

    elif args.mode == "section":
        if len(curGrid.blocks) != 1:
            print("section command works only on grids with 1 block")
            sys.exit(0)
        curGrid.blocks[0].section(args.iStart, args.iEnd, args.jStart, args.jEnd, args.kStart, args.kEnd)

    elif args.mode == "explode":
        # Split original multiblock grid in a list of single-block grids
        gridList = explodeGrid(curGrid)

        # Check if the user gave a reference name. Otherwise, use the input name as reference
        if args.outFile is None:
            # Get the base name
            outFile = os.path.splitext(os.path.basename(args.gridFile))[0]
        else:
            outFile = args.outFile

        # Generate a list of names for the files by adding integers to the reference name
        fileNames = [outFile + "_%03d" % index + ".cgns" for index in range(1, len(gridList) + 1)]

        # Save each grid
        for index in range(len(gridList)):
            gridList[index].writeToCGNS(fileNames[index])

        # Finish execution
        sys.exit(0)

    elif args.mode == "explodeKmin":
        # Split original multiblock grid in a list of single-block grids
        # that contains just the K = 1 face
        gridList = explodeGrid(curGrid, kMin=True)

        # Check if the user gave a reference name. Otherwise, use the input name as reference
        if args.outFile is None:
            # Get the base name
            outFile = os.path.splitext(os.path.basename(args.gridFile))[0]
        else:
            outFile = args.outFile

        # Generate a list of names for the files by adding integers to the reference name
        fileNames = [outFile + "_%03d" % index + ".xyz" for index in range(1, len(gridList) + 1)]

        # Save each grid
        for index in range(len(gridList)):
            gridList[index].writePlot3d(fileNames[index])

        # Finish execution
        sys.exit(0)

    elif args.mode == "explodeByZoneName":
        # Split original multiblock grid into a list of multiblock grids
        # that correspond to each component based on zone names
        gridList, nameList = explodeByZoneName(curGrid)

        # Save each grid
        for grid in gridList:
            fileName = grid.name
            grid.writeToCGNS(fileName + ".cgns")

        # Finish execution
        sys.exit(0)

    elif args.mode == "info":
        curGrid.printInfo()
        sys.exit(0)

    elif args.mode == "extrude":
        curGrid.extrude(args.direction)

    elif args.mode == "revolve":
        if args.normalDirection == args.axis:
            print("ERROR: Normal direction and revolve axis cannot be the same. Exiting...")
            sys.exit(0)
        curGrid.revolve(args.normalDirection, args.axis, args.startAngle, args.endAngle, args.nThetas)

    elif args.mode == "extractConv":
        # extracts the convergence history contained in the CGNS file and saves it in a pickle file

        # Check if the user gave a reference name. Otherwise, use the input name as reference
        if args.outFile is None:
            # Get the base name
            outFile = os.path.splitext(os.path.basename(args.gridFile))[0]

            # Add extension based on output type
            if args.outType == "pickle":
                outFile = outFile + ".pickle"
            elif args.outType == "tecplot":
                outFile = outFile + ".dat"

        else:
            outFile = args.outFile

        # The function readGrid already read all the convergence history arrays.
        # Now we just need to save them in a file!
        if args.outType == "pickle":
            with open(outFile, "w") as fid:
                pickle.dump(curGrid.convArray, fid)

        elif args.outType == "tecplot":
            # Create a single 2D array that will contain all data
            data = []

            # Get the number of iterations
            numIter = len(curGrid.convArray[curGrid.convArray.keys()[0]])

            # Append iteration counter
            data.append(range(1, numIter + 1))

            for entry in curGrid.convArray.keys():
                data.append(curGrid.convArray[entry])

            # Convert data to array
            data = np.array(data).T

            # Write tecplot results
            write_tecplot_file(outFile, "Convergence", ["Iteration"] + curGrid.convArray.keys(), data)

        # Print log
        print("Saved convergence history into:")
        print(outFile)

        # Finish execution
        sys.exit(0)

    # Write out the grid.
    curGrid.writeToCGNS(outFileName)

    # Possibly copy back to the original:
    if args.outFile is None:
        shutil.copyfile(outFileName, args.gridFile)
        shutil.rmtree(dirpath)
