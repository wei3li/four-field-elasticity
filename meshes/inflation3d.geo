SetFactory("OpenCASCADE");
lc = 400 * 0.001;
//+
Sphere(1) = {0, 0, 0, 1.0, 0, Pi/2, Pi/2};
//+
Sphere(2) = {0, 0, 0, 0.5, 0, Pi/2, Pi/2};
//+
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };
//+
Mesh.CharacteristicLengthMax = lc;
Mesh.CharacteristicLengthMin = lc;
//+
Physical Surface("neum_surfaces", 12) = {1};
//+
Physical Surface("diri_surfaces_xyz", 13) = {5};
//+
Physical Surface("diri_surfaces_y", 14) = {4};
//+
Physical Surface("diri_surfaces_z", 15) = {3};
//+
Physical Surface("diri_surfaces_x", 16) = {2};
//+
Physical Curve("diri_curves_xz", 17) = {6};
//+
Physical Curve("diri_curves_yz", 18) = {8};
//+
Physical Curve("diri_curves_xy", 19) = {5};
//+
Physical Point("diri_points_xyz", 20) = {5, 6, 4};
//+
Physical Point("diri_points_yz", 21) = {3};
//+
Physical Point("diri_points_xz", 22) = {2};
//+
Physical Point("diri_points_xy", 23) = {1};
//+
Physical Curve("diri_curves_xyz", 24) = {9, 7, 10};
//+
Physical Curve("diri_curves_z", 25) = {3};
//+
Physical Curve("diri_curves_y", 26) = {4};
//+
Mesh.Algorithm = 1;    // 2D: Delaunay
Mesh.Algorithm3D = 1;  // 3D: Delaunay

Mesh 3;

For i In {1:5}
  OptimizeMesh "Netgen";
EndFor

// gmsh partial_ball.geo -3 -save_all -o hollow_ball_delaunay_400.msh
