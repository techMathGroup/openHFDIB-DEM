# openHFDIB-DEM
Version with:
Velocity relaxation
Optional interpolation
Dynamic mesh refinment
Possibility to start from lastTime
Wall detection without stl findNearest
Valgrind optimization
	- lambda gradient computed only once when IBs moved
	- convex createIB calcInside minimized
Add models works with 2D non empty meshes by defining geometricD
	- this is also important for proper contact evaluation
	- default is set to mesh_.geometricD()
Rotation was repaired (I guess) 
Bodies can be switched to static after several timeSteps in contact with
wall or other static body
3D contact repaired. Get rid off SVD. Quicker and more stable. Precision not
so bad.