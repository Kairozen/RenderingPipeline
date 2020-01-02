# Small rendering engine

To compile : ```make```

To launch the engine : ```bin/tdsv data/file1.geo data/file2.geo```

Some .geo files are provided in the data folder of the project.

To move the camera in the program :

W, A, S, D : Move up, left, down, right respectively.

Z, X : Move forward, move backward.

Q and E : Unzoom and zoom.

Arrow keys : Rotate the objects.

Multiple .geo file can be displayed, objects will be put next to each other

Speed of the rotations, translations and zoom can be changed in ```include/libgraphic.hpp```

Known problems : 

- When loading multiple objects, the hidden faces are not properly applied on the objects and it seems all objects world position are counted as (0,0,0) when it's not the case