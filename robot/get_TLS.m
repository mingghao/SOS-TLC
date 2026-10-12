function [A_tls,b_tls,u_nom]=get_TLS(state,ref)
    wn = 1.2;
    K = [wn^3, 3*wn^2, 3*wn];

    persistent cached_wn cached_P
    if isempty(cached_wn) || cached_wn~=wn
        A = [0, 1, 0;
             0, 0, 1;
             -K];
        rhs=-eye(3);
        cached_P=reshape((kron(eye(3),A')+kron(A',eye(3)))\rhs(:),3,3);
        cached_P=(cached_P+cached_P')/2;
        cached_wn=wn;
    end

    P = cached_P;
    v = state(3);
    a = state(4);
    theta = state(5);
    dtheta = state(6);
    n = [cos(theta); sin(theta)];
    np = [-sin(theta); cos(theta)];
    vel = v*n;
    acc = a*n+v*dtheta*np;
    jerk0 = 2*a*dtheta*np-v*dtheta^2*n;
    B = [n, v*np];
    E = [(state(1:2)-ref.pos(:))';
         (vel-ref.vel(:))';
         (acc-ref.acc(:))'];
    Edot0 = [E(2:3,:); (jerk0-ref.jerk(:))'];
    PE=P*E;
    V=sum(E.*PE,1)';
    LfV=2*sum(PE.*Edot0,1)';
    LgV=2*diag(PE(3,:))*B;
    D=Parameter.delta;
    A_tls=[D*LgV,-eye(2)];
    b_tls=-V-D*LfV;

    desired_jerk = ref.jerk(:) - (K*E)' - jerk0;
    v_reg = 0.1;
    u_nom = [n'*desired_jerk; ...
             v/(v^2+v_reg^2)*(np'*desired_jerk)];

end
