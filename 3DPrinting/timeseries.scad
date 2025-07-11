difference() {scale([0.1, 0.1, 100]) difference () {surface(file="C:/Users/yolko/Desktop/turbulence-at-fronts/3DPrinting/images/b-timeseries.dat", convexity=3);
    translate([-1000, -1000, -2000]) cube(2000);
};
translate([40, 0, -1]) cylinder(r=1, h=8);
translate([40-4.3, 0, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
translate([40-2*4.3, 0, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
translate([40-3*4.3, 0, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
translate([40-4*4.3, 0, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
translate([40-5*4.3, 0, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
translate([40-6*4.3, 0, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
translate([40-7*4.3, 0, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
translate([40-8*4.3, 0, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
translate([40-9*4.3, 0, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
translate([40-10*4.3, 0, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
translate([40, 8, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
translate([40, 16, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
translate([40, 24, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
translate([40, 32, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
translate([40, 40, -1]) translate([-0.5, -0.5, 0]) cube([1, 1, 8]);
}


