SetFactory("OpenCASCADE");
lc = 250 * 0.01;
//+
Point(1) = {0, 0, 0, lc};
//+
Point(2) = {48, 44, 0, lc};
//+
Point(3) = {48, 60, 0, lc};
//+
Point(4) = {0, 44, 0, lc};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 5} {
  Surface{1};
}
//+
MeshSize {5, 6, 7, 8} = lc;
Mesh.Smoothing = 5;
//+
Physical Point("diri_points_xyz", 13) = {1, 4, 5, 8};
Physical Point("diri_points_z", 14) = {6, 7};
Physical Curve("diri_lines_xyz", 15) = {4, 5, 10, 12};
Physical Curve("diri_lines_z", 16) = {7, 9, 11};
Physical Surface("diri_surfaces_xyz", 17) = {5};
Physical Surface("diri_surfaces_z", 18) = {6};
Physical Surface("neum_surfaces", 19) = {1, 2, 3, 4};
