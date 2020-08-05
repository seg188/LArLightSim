###############################
###############################
#THIS FILE WILL NOT RUN BY ITSELF, IT SERVES AS AN EXAMPLE FOR CREATING A NEW GEOMETRY


#
#make sure you import the geometry file
#

import geometry as geo 

####
#### define the center points of the modules you would like to create
#### 
#### we create 2 modules that will be seperated by 10 cm, and are 1m x 3m x 1m. The origion will be defined as the point exactly between them

module1_center = geo.point_3d(55.0*cm, 0.0, 0.0)
module2_center = geo.point_3d(-55.0*cm, 0.0, 0.0)

whole_geometry_center = geo.point_3d(0.0, 0.0, 0.0)

module1_box = geo.box(100.0*cm, 300.0*cm, 100.0*cm, module1_center, "module_1") #MAKE SURE YOU NAME EVERY SHAPE IN THE GEOMETRY
module2_box = geo.box(100.0*cm, 300.0*cm, 100.0*cm, module2_center, "module_2") #MAKE SURE YOU NAME EVERY SHAPE IN THE GEOMETRY

new_geometry = geo.geometry( [module1_box, module2_box], whole_geometry_center )

#the new_geometry object can now be passed to the simulator initializer as the geometry for the simulaton
