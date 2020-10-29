function k = isFormClosure(contact)
    
    contact_coord = zeros(size(contact));
    contact_coord(:, 1:2) = contact(:, 1:2);
    angles = contact(:, 3);
    
    % find normal force
    normal = horzcat(cosd(angles),sind(angles));
    normal = horzcat(normal,zeros(size(normal,1),1));
    
    % find moment = r x F
    moment = cross(contact_coord, normal);
    moment = moment(:,3);
    
    wrench = horzcat(moment, normal(:, 1:2)).';    
    
    % check rank and calculate the linprog
    if rank(wrench) == 3
        f = ones(1, size(contact, 1));
        A = eye(size(contact, 1)).* -1;
        b = ones(1, size(contact, 1)) .* -1;
        Aeq = wrench;
        beq = [0 0 0];
        [k, ~, exitflag] = linprog(f, A, b, Aeq, beq);
        
        if exitflag == -2      
            error("The object is not in form closure.");
        else
            disp("The object is in form closure.");
        end
    end
end