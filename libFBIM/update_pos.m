function x_new = update_pos(body_ref, new_q, new_theta)

    q0     = body_ref.q;
    theta0 = body_ref.theta;
    x0     = body_ref.x_pos; 
    np     = body_ref.np;
    nb     = length(new_theta);

    nstart = 0;
    for nbod = 1 : nb
        
        dtheta = new_theta(nbod) - theta0;
        Rot = [cos(dtheta), -sin(dtheta); ...
               sin(dtheta),  cos(dtheta)];
        
        for i = 1:np
            
            
            v = x0(i,:) - q0;
            
            v = Rot * v';
            
            x_new(i+nstart,:) = v' + new_q(nbod,:);
            
        end
        
        nstart = nstart + np;
        
    end
end
