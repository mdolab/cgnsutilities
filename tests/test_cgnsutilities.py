import unittest
from cgnsutilities.cgnsutilities import readGrid, BC


class TestGrid(unittest.TestCase):
    def setUp(self):
        self.grid = readGrid("../examples/717_wl_L2.cgns")
        return

    def test_getTotalCellsNodes(self):
        totalCells, totalNodes = self.grid.getTotalCellsNodes()
        self.assertEqual(totalCells, 15120)
        self.assertEqual(totalNodes, 18753)

    def test_getWallCellsNodes(self):
        nWallCells, nWallNodes = self.grid.getWallCellsNodes()
        self.assertEqual(nWallCells, 756)
        self.assertEqual(nWallNodes, 893)

    def test_overwriteFamilies(self):
        # Find a specific BC and overwrite the family
        famFile = "../examples/overwriteFamily_famFile"
        # Check the family before overwriting.
        self.assertEqual(self.grid.blocks[0].bocos[0].family.strip().decode("utf-8"), "wall")
        self.grid.overwriteFamilies(famFile)
        self.assertEqual(self.grid.blocks[0].bocos[0].family, "wing1")
        return

    def test_overwriteBCs(self):
        # Find a specific BC and overwrite the type and family
        bcFile = "../examples/overwriteBCs_bcFile"
        # Check the BC before overwriting. Note that the "updated" BC is first deleted and new appended
        self.assertEqual(self.grid.blocks[0].bocos[0].family.strip().decode("utf-8"), "wall")
        self.assertEqual(self.grid.blocks[0].bocos[0].type, BC["bcwallviscous"])
        self.grid.overwriteBCs(bcFile)
        self.assertEqual(self.grid.blocks[0].bocos[-1].family, "wall_invicid")
        self.assertEqual(self.grid.blocks[0].bocos[-1].type, BC["bcwallinviscid"])
        return
