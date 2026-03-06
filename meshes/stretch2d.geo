SetFactory("OpenCASCADE");
lc = 78 * 0.0008;  // 31.3, 34.8, 48.2
r1 = 0.05; r2 = 0.1;

// Define the geometry
Rectangle(1) = {0, 0, 0, 1, 1, 0};
Disk(2) = {0.2, 0.15, 0, r2, r2};
Disk(3) = {0.15, 0.85, 0, r1, r1};
Disk(4) = {0.45, 0.8, 0, r1, r1};
Disk(5) = {0.5, 0.25, 0, r1, r1};
Disk(6) = {0.6, 0.5, 0, r1, r1};
Disk(7) = {0.25, 0.55, 0, r2, r2};
Disk(8) = {0.8, 0.75, 0, r2, r2};
Disk(9) = {0.8, 0.2, 0, r2, r2};
BooleanDifference(10) = { Surface{1}; Delete; }{ Surface{2:9}; Delete; };

// Physical groups
Physical Point("diri_points", 13) = {1, 2, 3, 4};
Physical Point("neum_points", 14) = {5, 6, 7, 8, 9, 10, 11, 12};
Physical Curve("diri_lines", 15) = {2, 3};
Physical Curve("neum_lines", 16) = {1, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Physical Surface("domain", 17) = {10};

// Mesh controls
MeshSize {1, 2, 3, 4, 5, 7, 8, 12} = lc;
MeshSize {6, 9, 10, 11} = lc/1.6;
Mesh.Smoothing = 10;
Mesh.Algorithm = 5;  // 5 -> Delaunay, 8 -> Frontal-Q
