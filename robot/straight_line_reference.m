function ref = straight_line_reference(t,p_start,p_goal)
    p_start=p_start(:);
    p_goal=p_goal(:);
    d = p_goal-p_start;
    L = norm(d,2);
    ref = struct('pos',p_start, ...
                 'vel',zeros(numel(d),1), ...
                 'acc',zeros(numel(d),1), ...
                 'jerk',zeros(numel(d),1));
    peak_factor = 315/128;
    T_move = peak_factor*(L/Parameter.vmax);

    if t>=T_move
        ref.pos=p_goal;
    end

    if t>0 && t<T_move
        s=t/T_move;
        q=1-s;
        a=min(s,q);
        S=a.^5.*(126+a.*(-420+a.*(540+a.*(-315+70*a))));
        right=s>0.5;
        S(right)=1-S(right);
        S1=630*s.^4.*q.^4;
        S2=2520*s.^3.*q.^3.*(1-2*s);
        S3=2520*s.^2.*q.^2.*(3-14*s+14*s.^2);
        ref.pos=p_start+d*S;
        ref.vel=d*(S1/T_move);
        ref.acc=d*(S2/T_move^2);
        ref.jerk=d*(S3/T_move^3);
    end
end
