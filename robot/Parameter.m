classdef Parameter
    properties (Constant)
        vmax = 1.2
        x_0 = 0
        y_0 = 3.0
        vel_0 = 0.2
        acc_0 = 0;
        theta_0 = 0.0;
        dtheta_0 = 0.0;

        x_d = 10
        y_d = 1.0


        x_o = [5,10];
        y_o = [1.5,2.5];
        r = [1,0.9];
        delta = 0.1;

        jerk_max = 10
        M_max = 8
        dt = 0.01
        T = 20
        Delta_t = 0.1

        cbf_k = [1,2,3];
        tlc_t = [0.4,0.5,0.6];

        L = 0.6;
        H = 0.2;
        marker_color = [0.5, 0.5, 0];
        color = ['b','m','r'];
        linestyle = ["--","-", ":"]

        exps = ["HOCBF","ZOH-TLC","SOS-TLC"]
        options = optimset('Display','off')
    end
end
