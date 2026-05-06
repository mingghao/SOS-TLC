classdef Parameter
    properties (Constant)
        delta = [1,3,1,3]
        dt = 0.1
        color = ["c","b","m","r","g"]
        exps = ["HOCBF","ZOH-TLC-1","ZOH-TLC-3","rTLC-1","rTLC-3"] 
        g = 9.81
        m = 1650
        f_0 = 0.1
        f_1 = 5
        f_2 = 0.25
        v_d = 24
        v_max = 30
        v_min = 0
        v_0 = 13.89 
        sigma = 10
        eps = 10
        c_a = 0.4
        c_d = 0.7
        options = optimoptions('quadprog', 'Algorithm',...
            'interior-point-convex','Display','off')
    end
end

