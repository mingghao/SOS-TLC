classdef Parameter
    properties (Constant)
        delta = [1.5,3,1.5,3,1.5,3]
        dt = 0.1
        color = ["c","b","m","y","k","r","g"]
        exps = ["HOCBF","ZOH-TLC-1.5","ZOH-TLC-3","SOS-TLC-1.5(QP)","SOS-TLC-3(QP)","SOS-TLC-1.5(SDP)","SOS-TLC-3(SDP)"] 
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
        % options = optimoptions('quadprog', 'Algorithm',...
        %     'interior-point-convex','Display','off')
        options = optimset('Display','off')
    end
end

