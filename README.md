Finite element (FE) test case representing a truss structure. The latter consists of 13 bars of square sections, connected by 8 nodes.

Two load cases are coded:
1) The extreme nodes are clamped and a vertical force is applied on node 3
2) Only node 1 is clamped and a horizontal force is applied on node 5.
In both cases, the value of the force applied is F = 1e4N.

Optimization problems are defined using three objective functions: the weight, the cost and the compliance of the structure.

The file 'TestClampedBeam.py' enables to test the FE code and visualize modifications of the mesh.

The file 'InitTest.py' is a configuration file to modify before running optimization with different methods.

