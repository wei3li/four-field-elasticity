SetFactory("OpenCASCADE");
//+
lc = 57.28 * 0.00125;  // 41.848, 57.28
//+
Box(1) = {0, 0, 0, 0.5, 0.5, 0.5};
Sphere(2) = {0, 0, 0, 0.3};
v() = BooleanDifference { Volume{1}; Delete; }{ Volume{2}; Delete; };
//+
Physical Volume("domain", 17) = {1};
Physical Surface("diri_surfaces_xyz", 18) = {3};
Physical Surface("diri_surfaces_x", 19) = {1};
Physical Surface("diri_surfaces_y", 20) = {2};
Physical Surface("diri_surfaces_z", 21) = {5};
Physical Surface("neum_surfaces", 22) = {4, 6, 7};
Physical Curve("neum_curves", 23) = {13};
Physical Curve("diri_curves_xyz", 24) = {10, 11, 2, 8};
Physical Curve("diri_curves_xy", 25) = {1};
Physical Curve("diri_curves_xz", 26) = {4};
Physical Curve("diri_curves_yz", 27) = {6};
Physical Curve("diri_curves_x", 28) = {3, 5};
Physical Curve("diri_curves_y", 29) = {7, 9};
Physical Curve("diri_curves_z", 30) = {15, 12, 14};
Physical Point("diri_points_xyz", 31) = {2, 3, 8, 9};
Physical Point("diri_points_xy", 32) = {1};
Physical Point("diri_points_yz", 33) = {6, 7};
Physical Point("diri_points_xz", 34) = {5, 4};
Physical Point("neum_points", 35) = {10};
//+
Mesh.Smoothing = 10;
Mesh.CharacteristicLengthMax = lc;
Mesh.CharacteristicLengthMin = lc;
// MeshSize { PointsOf{ Volume{v()}; }} = lc;
Mesh.Algorithm = 5;    // 5 -> Delaunay,  8 -> Frontal-Q
Mesh.Algorithm3D = 10; // 1 -> Delaunay, 10 -> HXT
