import os
import subprocess
from tkinter import Grid
import unittest
from matplotlib.pyplot import grid
from parameterized import parameterized
import numpy as np
from baseclasses import BaseRegTest
from cgnsutilities.cgnsutilities import readGrid, BC, combineGrids, mirrorGrid
import copy

baseDir = os.path.dirname(os.path.abspath(__file__))


class TestGrid(unittest.TestCase):
    def setUp(self):
        self.grid = readGrid(os.path.abspath(os.path.join(baseDir, "../examples/717_wl_L2.cgns")))

    def test_getTotalCellsNodes(self, train=False):
        totalCells, totalNodes = self.grid.getTotalCellsNodes()
        refFile = os.path.join(baseDir, "ref", "totalCellsNodes.ref")
        with BaseRegTest(refFile, train=train) as handler:
            handler.root_add_val("Total cells", totalCells, tol=0)
            handler.root_add_val("Total nodes", totalNodes, tol=0)

    def train_getTotalCellsNodes(self):
        self.test_getTotalCellsNodes(train=True)

    def test_getWallCellsNodes(self, train=False):
        nWallCells, nWallNodes = self.grid.getWallCellsNodes()
        refFile = os.path.join(baseDir, "ref", "wallCellsNodes.ref")
        with BaseRegTest(refFile, train=train) as handler:
            handler.root_add_val("Wall cells", nWallCells, tol=0)
            handler.root_add_val("Wall nodes", nWallNodes, tol=0)

    def train_getWallCellsNodes(self):
        self.test_getWallCellsNodes(train=True)

    def test_getBlockInfo(self, train=False):
        blockInfo = self.grid.getBlockInfo()
        refFile = os.path.join(baseDir, "ref", "blockInfo.ref")
        with BaseRegTest(refFile, train=train) as handler:
            handler.root_add_dict("Block info", blockInfo, tol=0)

    def train_getBlockInfo(self):
        self.test_getBlockInfo(train=True)

    def test_overwriteFamilies(self):
        # Find a specific BC and overwrite the family
        famFile = os.path.abspath(os.path.join(baseDir, "../examples/family_famFile"))
        # Check the family before overwriting.
        self.assertEqual(self.grid.blocks[0].bocos[0].family, "wall")
        self.grid.overwriteFamilies(famFile)
        self.assertEqual(self.grid.blocks[0].bocos[0].family, "wing1")

    def test_overwriteBCs(self):
        # Find a specific BC and overwrite the type and family
        bcFile = os.path.abspath(os.path.join(baseDir, "../examples/overwriteBCs_bcFile"))
        # Check the BC before overwriting. Note that the "updated" BC is first deleted and new appended
        self.assertEqual(self.grid.blocks[0].bocos[0].family, "wall")
        self.assertEqual(self.grid.blocks[0].bocos[0].type, BC["bcwallviscous"])
        self.grid.overwriteBCs(bcFile)
        self.assertEqual(self.grid.blocks[0].bocos[-1].family, "wall_inviscid")
        self.assertEqual(self.grid.blocks[0].bocos[-1].type, BC["bcwallinviscid"])

    def test_overwriteBCfamily(self):
        # Find a specific family and overwrite the BCs for the entire family
        # Check the BC before overwriting
        self.assertEqual(self.grid.blocks[0].bocos[1].family, "Far")
        self.assertEqual(self.grid.blocks[0].bocos[1].type, BC["bcfarfield"])
        self.assertEqual(self.grid.blocks[1].bocos[1].family, "Far")
        self.assertEqual(self.grid.blocks[1].bocos[1].type, BC["bcfarfield"])
        self.assertEqual(self.grid.blocks[1].bocos[0].family, "wall")
        self.assertEqual(self.grid.blocks[1].bocos[0].type, BC["bcwallviscous"])
        self.grid.overwriteBCFamilyWithBC("Far", "bcoverset", blockIDs=[2])

        # block 0 should be unchanged even though family matches
        self.assertEqual(self.grid.blocks[0].bocos[1].family, "Far")
        self.assertEqual(self.grid.blocks[0].bocos[1].type, BC["bcfarfield"])

        # block 1 wall family should be unchanged because family doesn't match
        self.assertEqual(self.grid.blocks[1].bocos[0].family, "wall")
        self.assertEqual(self.grid.blocks[1].bocos[0].type, BC["bcwallviscous"])

        # block 1 Far family should be overwritten with bcoverset
        self.assertEqual(self.grid.blocks[1].bocos[1].family, "Far")
        self.assertEqual(self.grid.blocks[1].bocos[1].type, BC["bcoverset"])

        # Check that using a non-existent blockID gives an error
        with self.assertRaises(IndexError):
            self.grid.overwriteBCFamilyWithBC("Far", "bcoverset", blockIDs=[0, 2])

    def test_overwriteBCs_array(self):
        self.grid.removeBCs()
        self.grid.overwriteBCs(os.path.abspath(os.path.join(baseDir, "../examples/hotwall_boco.info")))

        self.assertEqual(self.grid.blocks[4].bocos[0].type, BC["bcwallviscousisothermal"])
        np.testing.assert_array_equal(
            self.grid.blocks[4].bocos[0].dataSets[0].dirichletArrays[0].dataArr, np.array(range(300, 300 + 19 * 3))
        )

    def test_coarsen(self):
        coarse_bc_pt_range = {}

        for block in self.grid.blocks:
            coarse_bc_pt_range[block.name] = {}
            for boco in block.bocos:
                coarse_bc_pt_range[block.name][boco.name] = copy.deepcopy(boco.ptRange)

                # make new variable for convenience.
                # Base variable will be modified too because it is a mutable type
                c_range = coarse_bc_pt_range[block.name][boco.name]
                for idim in range(3):
                    c_range[idim, 0] = int(np.floor((c_range[idim, 0]) / 2)) + 1
                    c_range[idim, 1] = int(np.ceil((c_range[idim, 1]) / 2))

        self.grid.coarsen()
        totalCells = self.grid.getTotalCellsNodes()[0]
        self.assertEqual(15120 // 8, totalCells)

        # test that the point range was coarsened correctly.

        for block in self.grid.blocks:
            for boco in block.bocos:
                np.testing.assert_array_equal(coarse_bc_pt_range[block.name][boco.name], boco.ptRange)

        # we should be able to coarsen to a single cell in each block
        # and if a block already has a single cell, it should be skipped
        self.grid.coarsen()
        self.grid.coarsen()
        self.grid.coarsen()
        self.grid.coarsen()

        # 5 because there is one cell in each block
        totalCells = self.grid.getTotalCellsNodes()[0]
        self.assertEqual(5, totalCells)

    def test_refine(self):
        self.grid.refine("ijk")
        totalCells = self.grid.getTotalCellsNodes()[0]
        self.assertEqual(15120 * 8, totalCells)

    def test_mirror(self):
        newGrid = Grid()
        newGrid.name = []
        for blk in self.grid.blocks:
            blk.removeSymBCs()
            blk.B2Bs = []
        
            newGrid.name.append(blk.name)

            mirrorBlk = copy.deepcopy(blk)
            mirrorBlk.flip('z')
  
            mirrorBlk.name = blk.name + "_mirror"
            newGrid.name.append(mirrorBlk.name)
            

        newMirrorGrid = mirrorGrid(self.grid, "z", 1e-12, actualName=True)
        newNames = [blk.name for blk in newMirrorGrid.blocks]
        self.assertEqual(newNames, newGrid.name)

    def test_refine_axes(self):
        self.grid.refine("i")
        totalCells = self.grid.getTotalCellsNodes()[0]
        self.assertEqual(15120 * 2, totalCells)

        self.grid.refine("k")
        totalCells = self.grid.getTotalCellsNodes()[0]
        self.assertEqual(15120 * 4, totalCells)

        self.grid.refine("j")
        totalCells = self.grid.getTotalCellsNodes()[0]
        self.assertEqual(15120 * 8, totalCells)


class TestSimpleGrid(unittest.TestCase):
    def setUp(self):
        self.grid = readGrid(os.path.abspath(os.path.join(baseDir, "../examples/block_4x2x3.cgns")))

    def test_coarsen(self):

        self.grid.coarsen()
        for iDim in range(3):
            self.assertEqual(2, self.grid.blocks[0].dims[iDim])

        self.grid.coarsen()
        for iDim in range(3):
            self.assertEqual(2, self.grid.blocks[0].dims[iDim])


class TestCLI(unittest.TestCase):
    def setUp(self):
        self.grid = os.path.abspath(os.path.join(baseDir, "../examples/717_wl_L2.cgns"))

    def test_overwriteBCs_CLI(self):
        if os.path.isfile("717_wl_L2_overwriteBCs.cgns"):
            os.remove("717_wl_L2_overwriteBCs.cgns")

        cmd = "cgns_utils overwriteBC "
        cmd += self.grid + " "
        cmd += os.path.abspath(os.path.join(baseDir, "../examples/overwriteBCs_bcFile")) + " "
        cmd += "717_wl_L2_overwriteBCs.cgns"

        out = subprocess.run(cmd, shell=True)
        self.assertFalse(out.returncode)
        self.assertTrue(os.path.isfile("717_wl_L2_overwriteBCs.cgns"))
        os.remove("717_wl_L2_overwriteBCs.cgns")

    def test_overwriteFamilies(self):
        if os.path.isfile("717_wl_L2_overwriteFamilies.cgns"):
            os.remove("717_wl_L2_overwriteFamilies.cgns")

        cmd = "cgns_utils overwriteFamilies "
        cmd += self.grid + " "
        cmd += os.path.abspath(os.path.join(baseDir, "../examples/family_famFile")) + " "
        cmd += "717_wl_L2_overwriteFamilies.cgns"

        out = subprocess.run(cmd, shell=True)
        self.assertFalse(out.returncode)
        self.assertTrue(os.path.isfile("717_wl_L2_overwriteFamilies.cgns"))
        os.remove("717_wl_L2_overwriteFamilies.cgns")

    def test_overwriteBCFamilyWithBC(self):
        if os.path.isfile("717_wl_L2_overwriteBCFamilyWithBC.cgns"):
            os.remove("717_wl_L2_overwriteBCFamilyWithBC.cgns")

        cmd = "cgns_utils overwriteBCFamilyWithBC "
        cmd += self.grid + " Far bcoverset "
        cmd += "717_wl_L2_overwriteBCFamilyWithBC.cgns "
        cmd += "--blockIDs 2 4"

        out = subprocess.run(cmd, shell=True)
        self.assertFalse(out.returncode)
        self.assertTrue(os.path.isfile("717_wl_L2_overwriteBCFamilyWithBC.cgns"))
        os.remove("717_wl_L2_overwriteBCFamilyWithBC.cgns")


class TestBlock(unittest.TestCase):
    def setUp(self):
        self.grid = readGrid(os.path.abspath(os.path.join(baseDir, "../examples/717_wl_L2.cgns")))

    def test_getNumCells(self, train=False):
        refFile = os.path.join(baseDir, "ref", "block_getNumCells.ref")
        with BaseRegTest(refFile, train=train) as handler:
            # Just pick the first block from the grid to test
            numCells = self.grid.blocks[0].getNumCells()
            handler.root_add_val("Number of cells", numCells, tol=0)

    def train_getNumCells(self):
        self.test_getNumCells(train=True)

    def test_getNumNodes(self, train=False):
        refFile = os.path.join(baseDir, "ref", "block_getNumNodes.ref")
        with BaseRegTest(refFile, train=train) as handler:
            # Just pick the first block from the grid to test
            numCells = self.grid.blocks[0].getNumNodes()
            handler.root_add_val("Number of nodes in the first block", numCells, tol=0)

    def train_getNumNodes(self):
        self.test_getNumNodes(train=True)


class TestReturnFuncs(unittest.TestCase):
    def setUp(self):
        # Use the same grid file with different names
        self.grid1 = readGrid(os.path.abspath(os.path.join(baseDir, "../examples/717_wl_L2.cgns")))
        self.grid1.name = "grid1"
        self.grid2 = readGrid(os.path.abspath(os.path.join(baseDir, "../examples/717_wl_L2.cgns")))
        self.grid2.name = "grid2"
        self.grids = [self.grid1, self.grid2]

    def test_combineGrids(self, train=False):

        combinedNewNames = combineGrids(self.grids)
        combinedOldNames = combineGrids(self.grids, useOldNames=True)

        newNames = [blk.name for blk in combinedNewNames.blocks]
        oldNames = [blk.name for blk in combinedOldNames.blocks]

        # Blocks are named after the grid name
        self.assertEqual(
            newNames,
            [
                "grid1.00001",
                "grid1.00002",
                "grid1.00003",
                "grid1.00004",
                "grid1.00005",
                "grid2.00001",
                "grid2.00002",
                "grid2.00003",
                "grid2.00004",
                "grid2.00005",
            ],
        )

        # Block names from original grids are preserved
        self.assertEqual(
            oldNames,
            [
                "domain.00001",
                "domain.00002",
                "domain.00003",
                "domain.00004",
                "domain.00005",
                "domain.00001",
                "domain.00002",
                "domain.00003",
                "domain.00004",
                "domain.00005",
            ],
        )


class TestExamples(unittest.TestCase):

    # Get all example scripts in the example folder and its subfolders
    exampleDir = os.path.abspath(os.path.join(baseDir, "../examples"))
    examples = []
    for root, _, files in os.walk(exampleDir):
        for file in files:
            if file.endswith(".sh"):
                # Note: we cd into dir as each script assumes to be run in current folder
                cmd = f"cd {root} && bash {file}"
                examples.append([file, cmd])

    # Generate a custom test function name for each script that will be run
    def generateFuncName(testcase_func, _, param):
        return "{}_{}".format(
            testcase_func.__name__,
            parameterized.to_safe_name(param.args[0]),
        )

    @parameterized.expand(examples, name_func=generateFuncName)
    def test_example(self, _, cmd):
        """
        Extract and run all examples from the examples folder.
        We assume that all have a .sh extension. Note that this test
        does not guarantee that they make sense, only that they run.
        """
        out = subprocess.run(cmd, shell=True)
        self.assertFalse(out.returncode)


if __name__ == "__main__":
    unittest.main()
