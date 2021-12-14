# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
This sample loads model red blood cells and simulates its motion
in a complex geometry.
"""

import espressomd

required_features = ["LB_BOUNDARIES", "EXTERNAL_FORCES", "SOFT_SPHERE",
                     "MEMBRANE_COLLISION", "MASS"]
espressomd.assert_features(required_features)

from espressomd import lb
from espressomd import lbboundaries
from espressomd import shapes

import numpy as np
import os
import sys
import warnings

import object_in_fluid as oif
from object_in_fluid.oif_utils import output_vtk_rhomboid, output_vtk_cylinder

import time

t = time.localtime()
simNo = str(t.tm_year) + "-" + str(t.tm_mon) + "-" + str(t.tm_mday)  + "_" + str(t.tm_hour) + "-" + str(t.tm_min)
if not os.path.isdir("output/sim" + simNo):
    os.makedirs("output/sim" + simNo)
else:
    warnings.warn("Folder {} already exists, files will be overwritten"
                  .format("output/sim" + simNo))


boxX = 300.0
boxY = 1000.0
boxZ = 12.0
slrNum = 1
slrSX = 100
slrX = 40
slrYX = 87.5
slrYY = boxY
slrZ = 18
time_step = 0.05

system = espressomd.System(box_l=(boxX, boxY, boxZ + slrZ))
system.time_step = time_step
system.cell_system.skin = 0.2

# creating the template for RBCs
cell_type = oif.OifCellType(
    nodes_file="input/rbc374nodes.dat", triangles_file="input/rbc374triangles.dat",
    system=system, ks=0.04, kb=0.016, kal=0.02, kag=0.9, kv=1.0, check_orientation=False, resize=(2.0, 2.0, 2.0))

# creating the RBCs
cell0 = oif.OifCell(cell_type=cell_type,
                    particle_type=0, origin=[10.0, 25.0, 5.0])
cell1 = oif.OifCell(cell_type=cell_type,
                    particle_type=1, origin=[10.0, 25.0, 8.0])

cell2 = oif.OifCell(cell_type=cell_type,
                    particle_type=2, origin=[10.0, 50.0, 5.0])
cell3 = oif.OifCell(cell_type=cell_type,
                    particle_type=3, origin=[10.0, 50.0, 8.0])

cell4 = oif.OifCell(cell_type=cell_type,
                    particle_type=4, origin=[10.0, 75.0, 5.0])
cell5 = oif.OifCell(cell_type=cell_type,
                    particle_type=5, origin=[10.0, 75.0, 8.0])

cell6 = oif.OifCell(cell_type=cell_type,
                    particle_type=6, origin=[10.0, 100.0, 5.0])
cell7 = oif.OifCell(cell_type=cell_type,
                    particle_type=7, origin=[10.0, 100.0, 8.0])

cell8 = oif.OifCell(cell_type=cell_type,
                    particle_type=8, origin=[10.0, 125.0, 5.0])
cell9 = oif.OifCell(cell_type=cell_type,
                    particle_type=9, origin=[10.0, 125.0, 8.0])

cell10 = oif.OifCell(cell_type=cell_type,
                    particle_type=10, origin=[10.0, 150.0, 5.0])
cell11 = oif.OifCell(cell_type=cell_type,
                    particle_type=11, origin=[10.0, 150.0, 8.0])

cell12 = oif.OifCell(cell_type=cell_type,
                    particle_type=12, origin=[10.0, 175.0, 5.0])
cell13 = oif.OifCell(cell_type=cell_type,
                    particle_type=13, origin=[10.0, 175.0, 8.0])

cell14 = oif.OifCell(cell_type=cell_type,
                    particle_type=14, origin=[10.0, 200.0, 5.0])
cell15 = oif.OifCell(cell_type=cell_type,
                    particle_type=15, origin=[10.0, 200.0, 8.0])

cell16 = oif.OifCell(cell_type=cell_type,
                    particle_type=16, origin=[10.0, 225.0, 5.0])
cell17 = oif.OifCell(cell_type=cell_type,
                    particle_type=17, origin=[10.0, 225.0, 8.0])

cell18 = oif.OifCell(cell_type=cell_type,
                    particle_type=18, origin=[10.0, 250.0, 5.0])
cell19 = oif.OifCell(cell_type=cell_type,
                    particle_type=19, origin=[10.0, 250.0, 8.0])

cell20 = oif.OifCell(cell_type=cell_type,
                    particle_type=20, origin=[10.0, 275.0, 5.0])
cell21 = oif.OifCell(cell_type=cell_type,
                    particle_type=21, origin=[10.0, 275.0, 8.0])

cell22 = oif.OifCell(cell_type=cell_type,
                    particle_type=22, origin=[10.0, 300.0, 5.0])
cell23 = oif.OifCell(cell_type=cell_type,
                    particle_type=23, origin=[10.0, 300.0, 8.0])

cell24 = oif.OifCell(cell_type=cell_type,
                    particle_type=24, origin=[10.0, 325.0, 5.0])
cell25 = oif.OifCell(cell_type=cell_type,
                    particle_type=25, origin=[10.0, 325.0, 8.0])

cell26 = oif.OifCell(cell_type=cell_type,
                    particle_type=26, origin=[10.0, 350.0, 5.0])
cell27 = oif.OifCell(cell_type=cell_type,
                    particle_type=27, origin=[10.0, 350.0, 8.0])

cell28 = oif.OifCell(cell_type=cell_type,
                    particle_type=28, origin=[10.0, 375.0, 5.0])
cell29 = oif.OifCell(cell_type=cell_type,
                    particle_type=29, origin=[10.0, 375.0, 8.0])

cell30 = oif.OifCell(cell_type=cell_type,
                    particle_type=30, origin=[10.0, 400.0, 5.0])
cell31 = oif.OifCell(cell_type=cell_type,
                    particle_type=31, origin=[10.0, 400.0, 8.0])

cell32 = oif.OifCell(cell_type=cell_type,
                    particle_type=32, origin=[10.0, 425.0, 5.0])
cell33 = oif.OifCell(cell_type=cell_type,
                    particle_type=33, origin=[10.0, 425.0, 8.0])

cell34 = oif.OifCell(cell_type=cell_type,
                    particle_type=34, origin=[10.0, 450.0, 5.0])
cell35 = oif.OifCell(cell_type=cell_type,
                    particle_type=35, origin=[10.0, 450.0, 8.0])

cell36 = oif.OifCell(cell_type=cell_type,
                    particle_type=36, origin=[10.0, 475.0, 5.0])
cell37 = oif.OifCell(cell_type=cell_type,
                    particle_type=37, origin=[10.0, 475.0, 8.0])

cell38 = oif.OifCell(cell_type=cell_type,
                    particle_type=38, origin=[10.0, 500.0, 5.0])
cell39 = oif.OifCell(cell_type=cell_type,
                    particle_type=39, origin=[10.0, 500.0, 8.0])

cell40 = oif.OifCell(cell_type=cell_type,
                    particle_type=40, origin=[10.0, 525.0, 5.0])
cell41 = oif.OifCell(cell_type=cell_type,
                    particle_type=41, origin=[10.0, 525.0, 8.0])

cell42 = oif.OifCell(cell_type=cell_type,
                    particle_type=42, origin=[10.0, 550.0, 5.0])
cell43 = oif.OifCell(cell_type=cell_type,
                    particle_type=43, origin=[10.0, 550.0, 8.0])

cell44 = oif.OifCell(cell_type=cell_type,
                    particle_type=44, origin=[10.0, 575.0, 5.0])
cell45 = oif.OifCell(cell_type=cell_type,
                    particle_type=45, origin=[10.0, 575.0, 8.0])

cell46 = oif.OifCell(cell_type=cell_type,
                    particle_type=46, origin=[10.0, 600.0, 5.0])
cell47 = oif.OifCell(cell_type=cell_type,
                    particle_type=47, origin=[10.0, 600.0, 8.0])

cell48 = oif.OifCell(cell_type=cell_type,
                    particle_type=48, origin=[10.0, 625.0, 5.0])
cell49 = oif.OifCell(cell_type=cell_type,
                    particle_type=49, origin=[10.0, 625.0, 8.0])

cell50 = oif.OifCell(cell_type=cell_type,
                    particle_type=50, origin=[10.0, 650.0, 5.0])
cell51 = oif.OifCell(cell_type=cell_type,
                    particle_type=51, origin=[10.0, 650.0, 8.0])

cell52 = oif.OifCell(cell_type=cell_type,
                    particle_type=52, origin=[10.0, 675.0, 5.0])
cell53 = oif.OifCell(cell_type=cell_type,
                    particle_type=53, origin=[10.0, 675.0, 8.0])

cell54 = oif.OifCell(cell_type=cell_type,
                    particle_type=54, origin=[10.0, 700.0, 5.0])
cell55 = oif.OifCell(cell_type=cell_type,
                    particle_type=55, origin=[10.0, 700.0, 8.0])

cell56 = oif.OifCell(cell_type=cell_type,
                    particle_type=56, origin=[10.0, 725.0, 5.0])
cell57 = oif.OifCell(cell_type=cell_type,
                    particle_type=57, origin=[10.0, 725.0, 8.0])

cell58 = oif.OifCell(cell_type=cell_type,
                    particle_type=58, origin=[10.0, 750.0, 5.0])
cell59 = oif.OifCell(cell_type=cell_type,
                    particle_type=59, origin=[10.0, 750.0, 8.0])

cell60 = oif.OifCell(cell_type=cell_type,
                    particle_type=60, origin=[10.0, 775.0, 5.0])
cell61 = oif.OifCell(cell_type=cell_type,
                    particle_type=61, origin=[10.0, 775.0, 8.0])

cell62 = oif.OifCell(cell_type=cell_type,
                    particle_type=62, origin=[10.0, 800.0, 5.0])
cell63 = oif.OifCell(cell_type=cell_type,
                    particle_type=63, origin=[10.0, 800.0, 8.0])

cell64 = oif.OifCell(cell_type=cell_type,
                    particle_type=64, origin=[10.0, 825.0, 5.0])
cell65 = oif.OifCell(cell_type=cell_type,
                    particle_type=65, origin=[10.0, 825.0, 8.0])

cell66 = oif.OifCell(cell_type=cell_type,
                    particle_type=66, origin=[10.0, 850.0, 5.0])
cell67 = oif.OifCell(cell_type=cell_type,
                    particle_type=67, origin=[10.0, 850.0, 8.0])

cell68 = oif.OifCell(cell_type=cell_type,
                    particle_type=68, origin=[10.0, 875.0, 5.0])
cell69 = oif.OifCell(cell_type=cell_type,
                    particle_type=69, origin=[10.0, 875.0, 8.0])

cell70 = oif.OifCell(cell_type=cell_type,
                    particle_type=70, origin=[10.0, 900.0, 5.0])
cell71 = oif.OifCell(cell_type=cell_type,
                    particle_type=71, origin=[10.0, 900.0, 8.0])

cell72 = oif.OifCell(cell_type=cell_type,
                    particle_type=72, origin=[10.0, 925.0, 5.0])
cell73 = oif.OifCell(cell_type=cell_type,
                    particle_type=73, origin=[10.0, 925.0, 8.0])

cell74 = oif.OifCell(cell_type=cell_type,
                    particle_type=74, origin=[10.0, 950.0, 5.0])
cell75 = oif.OifCell(cell_type=cell_type,
                    particle_type=75, origin=[10.0, 950.0, 8.0])

cell76 = oif.OifCell(cell_type=cell_type,
                    particle_type=76, origin=[10.0, 950.0, 5.0])
cell77 = oif.OifCell(cell_type=cell_type,
                    particle_type=77, origin=[10.0, 950.0, 8.0])


# cell-wall interactions
system.non_bonded_inter[0, 10].soft_sphere.set_params(
    a=0.0001, n=1.2, cutoff=0.1, offset=0.0)
system.non_bonded_inter[1, 10].soft_sphere.set_params(
    a=0.0001, n=1.2, cutoff=0.1, offset=0.0)

# cell-cell interactions
system.non_bonded_inter[0, 1].membrane_collision.set_params(
    a=0.0001, n=1.2, cutoff=0.1, offset=0.0)


# fluid
lbf = espressomd.lb.LBFluid(agrid=1, dens=1.0, visc=1.5,
                            tau=0.1, ext_force_density=[0.002, 0.0, 0.0])
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, gamma=1.5)

# creating boundaries and obstacles in the channel
# OutputVtk writes a file
# lbboundaries created boundaries for fluid
# constraints created boundaries for the cells

boundaries = []

############################################################################################################
                                                    # MAIN BODY #
############################################################################################################
# bottom of the channel
bottom_shape = shapes.Rhomboid(corner=[0.0, 0.0, 0.0], a=[boxX, 0.0, 0.0], 
                               b=[0.0, boxY, 0.0], c=[0.0, 0.0, 1.0], direction=1)
boundaries.append(bottom_shape)
output_vtk_rhomboid(bottom_shape, out_file="output/sim" + str(simNo) + "/wallBottom.vtk")

# top of the channel
top_shape = shapes.Rhomboid(corner=[0.0, 0.0, boxZ - 1.0], a=[slrSX, 0.0, 0.0], 
                            b=[0.0, boxY, 0.0], c=[0.0, 0.0, 1.0], direction=1)
boundaries.append(top_shape)
output_vtk_rhomboid(top_shape, out_file="output/sim" + str(simNo) + "/wallTop.vtk")

# front wall of the channel
front_shape = shapes.Rhomboid(corner=[0.0, 0.0, 0.0], a=[boxX, 0.0, 0.0], 
                              b=[0.0, 1.0, 0.0], c=[0.0, 0.0, boxZ], direction=1)
boundaries.append(front_shape)
output_vtk_rhomboid(front_shape, out_file="output/sim" + str(simNo) + "/wallFront.vtk")

# back wall of the channel
back_shape = shapes.Rhomboid(corner=[0.0, boxY - 1.0, 0.0], a=[boxX, 0.0, 0.0], 
                             b=[0.0, 1.0, 0.0], c=[0.0, 0.0, boxZ], direction=1)
boundaries.append(back_shape)
output_vtk_rhomboid(back_shape, out_file="output/sim" + str(simNo) + "/wallBack.vtk")


############################################################################################################
                                            # Slant Ridges Structure #
############################################################################################################

new_shape_1 = shapes.Rhomboid(corner=[slrSX, 0.0, boxZ - 1.0], a=[slrYX, slrYY, 0.0], 
                              b=[1.0, 0.0, 0.0], c=[0.0, 0.0, slrZ], direction=1)
boundaries.append(new_shape_1)
output_vtk_rhomboid(new_shape_1, out_file="output/sim" + str(simNo) + "/wallnewa.vtk")

new_shape_2 = shapes.Rhomboid(corner=[slrSX + slrX, 0.0, boxZ - 1.0], a=[slrYX, slrYY, 0.0], 
                              b=[1.0, 0.0, 0.0], c=[0.0, 0.0, slrZ], direction=1)
boundaries.append(new_shape_2)
output_vtk_rhomboid(new_shape_2, out_file="output/sim" + str(simNo) + "/wallnewb.vtk")

new_shape_3 = shapes.Rhomboid(corner=[slrSX, 0.0, boxZ + slrZ - 1.0], a=[slrX + 1.0, 0.0, 0], 
                              b=[slrYX, slrYY, 0.0], c=[0.0, 0.0, 1.0], direction=1)
boundaries.append(new_shape_3)
output_vtk_rhomboid(new_shape_3, out_file="output/sim" + str(simNo) + "/wallnewc.vtk")

############################################################################################################
                                    # Slant Ridges Structure Supplement #
############################################################################################################

new_shape_4 = shapes.Rhomboid(corner=[slrSX - slrYX, 0.0, boxZ - 1.0], a=[slrYX + 1.0, 0.0, 0.0], 
                              b=[slrYX, slrYY, 0.0], c=[0.0, 0.0, 1], direction=1)
boundaries.append(new_shape_4)
output_vtk_rhomboid(new_shape_4, out_file="output/sim" + str(simNo) + "/wallnewd.vtk")

new_shape_5 = shapes.Rhomboid(corner=[slrSX + slrX, 0.0, boxZ - 1.0], a=[slrYX + 1.0, 0.0, 0.0], 
                              b=[slrYX, slrYY, 0.0], c=[0.0, 0.0, 1.0], direction=1)
boundaries.append(new_shape_5)
output_vtk_rhomboid(new_shape_5, out_file="output/sim" + str(simNo) + "/wallnewe.vtk")

new_shape_6 = shapes.Rhomboid(corner=[slrSX + slrX + slrYX, 0.0, boxZ - 1.0], 
                              a=[boxX - slrSX - slrX - slrYX, 0.0, 0.0], b=[0.0, boxY, 0.0], 
                              c=[0.0, 0.0, 1.0], direction=1)
boundaries.append(new_shape_6)
output_vtk_rhomboid(new_shape_6, out_file="output/sim" + str(simNo) + "/wallnewf.vtk")

new_shape_7 = shapes.Rhomboid(corner=[slrSX, 0.0, boxZ - 1.0], a=[slrX + 1.0, 0.0, 0.0], 
                              b=[0.0, 1.0, 0.0], c=[0.0, 0.0, slrZ + 1.0], direction=1)
boundaries.append(new_shape_7)
output_vtk_rhomboid(new_shape_7, out_file="output/sim" + str(simNo) + "/wallnewg.vtk")

new_shape_8 = shapes.Rhomboid(corner=[slrSX + slrYX, slrYY - 1.0, boxZ - 1.0], a=[slrX + 1.0, 0.0, 0.0], 
                              b=[0.0, 1.0, 0.0], c=[0.0, 0.0, slrZ + 1.0], direction=1)
boundaries.append(new_shape_8)
output_vtk_rhomboid(new_shape_8, out_file="output/sim" + str(simNo) + "/wallnewh.vtk")


for boundary in boundaries:
    system.lbboundaries.add(lbboundaries.LBBoundary(shape=boundary))
    system.constraints.add(shape=boundary, particle_type=10)

# simulation time is (maxCycle-1)*time_step
# e.g. maxCycle = 100 ===> simulation time = 9.9s
maxCycle = 200
totalCell = 78
# main integration loop
names = locals()
for j in range(0, totalCell):
    vtkName = "/cell" + str(j) + "_"
    names['cell' + str(j)].output_vtk_pos_folded(
        file_name="output/sim" + str(simNo) + vtkName + "0.vtk")


for i in range(1, maxCycle):
    system.integrator.run(steps=500)
    for j in range(0, totalCell):
        vtkName = "/cell" + str(j) + "_"
        names['cell' + str(j)].output_vtk_pos_folded(
            file_name="output/sim" + str(simNo) + vtkName + str(i) + ".vtk")
    print("time: {:.1f}".format(i * time_step))
print("Simulation completed.")
