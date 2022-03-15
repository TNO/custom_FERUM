h=%h%;
ech=7.0;

Point(1) = {0.0,0.0,0.0,h*ech};
Point(2) = {1.0*ech,0.0,0.0,h*ech};
Point(3) = {1.0*ech,1.0*ech,0.0,h*ech};
Point(4) = {0.0,1.0*ech,0.0,h*ech};
Point(5) = {3.0*ech,1.0*ech,-0.15*ech,h*ech};

Point(6) = {0.0,0.0,-0.4*ech,h*ech};
Point(7) = {1.0*ech,0.0,-0.4*ech,h*ech};
Point(8) = {1.0*ech,1.0*ech,-0.4*ech,h*ech};
Point(9) = {0.0,1.0*ech,-0.4*ech,h*ech};
Point(10) = {3.0*ech,1.0*ech,-0.25*ech,h*ech};

Point(99) = {3.0*ech,1.0*ech,-0.2*ech,h*ech};

Line(11)={1,2};
Line(12)={2,3};
Line(13)={3,4};
Line(14)={4,1};
Line(15)={1,3};
Line(16)={2,5};
Line(17)={3,5};
Line(18)={3,7};
Line(19)={4,6};

Line(51)={4,8};
Line(52)={3,10};
Line(53)={1,7};
Line(54)={2,10};
Line(55)={4,9};
Line(56)={3,8};
Line(57)={5,99};
Line(58)={1,6};
Line(59)={2,7};
Line(60)={99,10};

Line(81)={6,7};
Line(82)={7,8};
Line(83)={8,9};
Line(84)={9,6};
Line(85)={6,8};
Line(86)={7,10};
Line(87)={8,10};

Physical Point(21) = {1};
Physical Point(24) = {2};
Physical Point(25) = {5};

Physical Point(26) = {6};
Physical Point(29) = {7};
Physical Point(30) = {10};

Physical Point(90) = {99};

Physical Line(31) = {11};
Physical Line(32) = {12};
Physical Line(33) = {13};
Physical Line(34) = {14};
Physical Line(35) = {15};
Physical Line(36) = {16};
Physical Line(37) = {17};
Physical Line(38) = {18};
Physical Line(39) = {19};

Physical Line(61)={51};
Physical Line(62)={52};
Physical Line(63)={53};
Physical Line(64)={54};
Physical Line(65)={55};
Physical Line(66)={56};
Physical Line(67)={57};
Physical Line(68)={58};
Physical Line(69)={59};
Physical Line(70)={60};

Physical Line(41) = {81};
Physical Line(42) = {82};
Physical Line(43) = {83};
Physical Line(44) = {84};
Physical Line(45) = {85};
Physical Line(46) = {86};
Physical Line(47) = {87};

