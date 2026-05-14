classdef Parameter
    properties (Constant)
        m = 1
        I = 0.25
        g = 9.81
        x_0 = 0
        z_0 = 1.0
        theta_0 = 0.0
        x_d = 10
        z_d = 1.0
        theta_d = 0
        theta_max = pi/3;
        x_o = [5,8];
        z_o = [1.5,2.5];
        r = [1,1];
        delta = 0.1;
        F_max = 40
        M_max = 20
        uF_max = 70
        dt=0.01
        T = 6
        Delta_t = 0.1

        L = 0.6; 
        H = 0.2; 
        marker_color = [0.5, 0.5, 0];
        color = ['g','m','r'];

        exps = ["HOCBF","ZOH-TLC","SOS-TLC"]
        options = optimset('Display','off')
    end
end

