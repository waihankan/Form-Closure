function k = isFormClosure(contact)
    % HELP SESSION :
    % _______________________________
    % for object not in form closure.    
    % contact = [1.5 1 225; -1.5 1 315; -1.5 -1 45; 1.5 -1 154];
    % k = isFormClosure(contact);
    % OUTPUT: The object is not in form closure; Error using
    % isFormClosure(line 38).
    % _______________________________
    % for object in form closure.
    % contact = [1.5 1 225; -1.5 1 315; -1.5 -1 45; 1.5 -1 135];
    % k = isFormClosure(contact);
    % OUTPUT: The object is in form closure; corresponding k values.
    
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
    if rank(wrench) >= 3
        f = ones(1, size(contact, 1));
        A = eye(size(contact, 1)).* -1;
        b = ones(1, size(contact, 1)) .* -1;
        Aeq = wrench;
        beq = [0 0 0];
        [k, ~, exitflag] = linprog(f, A, b, Aeq, beq);
        
        if exitflag == 1 
            disp("The object is in form closure.");
        else
            error("The object is not in form closure.");
        end
    else
        fprintf("Rank of wrench matrix = %d\n", rank(wrench));
        error("The object is not in form closure");  
    end
end

