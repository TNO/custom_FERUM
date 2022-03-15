function g = gfun_rigid_frame_system(M1, M2, M3, S1, S2, als)

switch als
    case 1
        g               = 2*M1 + 2*M2 + 0*M3 - 15*S1 - 0*S2;
    case 2
        g               = 1*M1 + 3*M2 + 2*M3 - 15*S1 - 10*S2;
    case 3
        g               = 2*M1 + 1*M2 + 1*M3 - 15*S1 - 0*S2;
    case 4
        g               = 1*M1 + 2*M2 + 1*M3 - 15*S1 - 0*S2;
    case 5
        g               = 1*M1 + 1*M2 + 2*M3 - 15*S1 - 0*S2;
    case 6
        g               = 1*M1 + 1*M2 + 4*M3 - 15*S1 - 10*S2;    
    case 999
        gg(1)           = 2*M1 + 2*M2 + 0*M3 - 15*S1 - 0*S2;
        gg(2)           = 1*M1 + 3*M2 + 2*M3 - 15*S1 - 10*S2;
        gg(3)           = 2*M1 + 1*M2 + 1*M3 - 15*S1 - 0*S2;
        gg(4)           = 1*M1 + 2*M2 + 1*M3 - 15*S1 - 0*S2;
        gg(5)           = 1*M1 + 1*M2 + 2*M3 - 15*S1 - 0*S2;
        gg(6)           = 1*M1 + 1*M2 + 4*M3 - 15*S1 - 10*S2;
        
        g               = min(gg);
end

end