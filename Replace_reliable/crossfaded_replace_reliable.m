function [output] = crossfaded_replace_reliable(x, xc, Mr, flag, w, transition_type, treat_short)  
% CROSSFADED_REPLACE_RELIABLE contains several methods for replacing 
% the reliable samples of signals obtained from declipping algorithms 
% producing solutions inconsistent in the reliable part.

% Input parameters:
%       x                   declipped signal using an inconsistent declipping method
%       xc                  clipped signal
%       Mr                  mask of reliable samples
%       flag                position of the crossfading transition
%       w                   length of the crossfading transition
%       transition_type     type of the crossfade
%       treat_short         way of treating segments shorter than w

% Pavel Záviška, Brno University of Technology, 2021

    if nargin < 4
        flag = 3; % 1 ... old traditional way of replacing all reliable samples (RR)
                  % 2 ... transition in the clipped part
                  % 3 ... transition in the reliable part
                  % 4 ... transition in both parts
                  
        w = 8; % length of the transition
        
        transition_type = 'soft'; % options: 'linear', 'soft'
        
        treat_short = 2; % options: 0 ... do not replace short reliable  
                         %          1 ... replace short reliable the traditional way 
                         %          2 ... shorten w accordingly
    end
    
     
       
    switch transition_type
        case 'linear'
            transition_up = linspace(0, 1, w+2).';
            transition_up = transition_up(2:end-1);
            transition_down = 1 - transition_up;
        case 'soft'
            transition_up = sin(linspace(0, pi/2, w+2).').^2;
            transition_up = transition_up(2:end-1);
            transition_down = 1 - transition_up;
        otherwise
            error('Invalid transition type!')
    end

    
    switch flag
        
        case 1  % traditional way
            output = x;
            output(Mr) = xc(Mr);
            
            
        case 2  % transition in clipped part          
            c_start_idxs = find(diff(Mr)==-1) + 1;
            c_stop_idxs = find(diff(Mr)==1);
            
            if length(c_stop_idxs) < length(c_start_idxs)
                c_stop_idxs = [c_stop_idxs; length(Mr)];
            elseif length(c_start_idxs) < length(c_stop_idxs)
                c_start_idxs = [1 c_start_idxs];
            end
            
            c_lengths = c_stop_idxs - c_start_idxs +1;
            
            yc = double(Mr);
            yr = 1 - yc;
            
            for n = 1:length(c_lengths) 
                if c_lengths(n) >= 2*w
                   yc(c_start_idxs(n):c_start_idxs(n)+w-1) = transition_down;
                   yr(c_start_idxs(n):c_start_idxs(n)+w-1) = transition_up;
                   
                   yc(c_stop_idxs(n)+1-w:c_stop_idxs(n)) = transition_up;
                   yr(c_stop_idxs(n)+1-w:c_stop_idxs(n)) = transition_down;
                   
                else
                    switch treat_short
                        case 0 % Does not make much sense, it throws away the restoration peaks 
                            yc(c_start_idxs(n):c_start_idxs(n) + c_lengths(n) - 1) = ones(c_lengths(n), 1);
                            yr(c_start_idxs(n):c_start_idxs(n) + c_lengths(n) - 1) = zeros(c_lengths(n), 1);
                        case 1
                            % No need to do anything
                        case 2                            
                            transition_r = linspace(0,1, ceil(c_lengths(n)/2)+2).';
                            transition_r = transition_r(2:end-1);
                            if c_lengths(n) > 1
                                if mod(c_lengths(n),2) == 0
                                    transition_r = [transition_r; transition_r(end:-1:1)];
                                else
                                    transition_r = [transition_r; transition_r(end-1:-1:1)];
                                end
                            end
                            transition_c = 1 - transition_r;

                            yc(c_start_idxs(n):c_start_idxs(n) + c_lengths(n) - 1) =  transition_c; 
                            yr(c_start_idxs(n):c_start_idxs(n) + c_lengths(n) - 1) =  transition_r;
                    end       
                end
            end
            
            output = xc.*yc + x.*yr;
            
            
        case 3  % transition in reliable part (probably the most reasonable way)          
            difference = diff(Mr);
            
            c_start_idxs = find(diff(Mr)==-1) + 1;
            c_stop_idxs = find(diff(Mr)==1);
            
            if length(c_stop_idxs) < length(c_start_idxs)
                c_stop_idxs = [c_stop_idxs; length(Mr)];
            elseif length(c_start_idxs) < length(c_stop_idxs)
                c_start_idxs = [1 c_start_idxs];
            end
            
            r_lengths = [c_start_idxs; 0] - [0; c_stop_idxs] -1;
            r_lengths(end) = length(x) - c_stop_idxs(end);
            
            yc = double(Mr);
            yr = 1 - yc;
            
            N = 1;
            for n = 1:length(x)-1
                if difference(n) == -1  % transition from reliable to clipped
                    if r_lengths(N) >= 2*w  % reliable part is long enough
                        yc(n-w+1 : n) = transition_down;
                        yr(n-w+1 : n) = transition_up;
                    % No need to solve if reliable part is not long enough, it will be solved for difference(n) == 1
                    end
                elseif difference(n) == 1  % transition from clipped to reliable
                    N = N+1;
                    if r_lengths(N) >= 2*w % reliable part is long enough
                        yc(n+1 : n+w) = transition_up;
                        yr(n+1 : n+w) = transition_down;
                    else
                        switch treat_short
                            case 0  % do not replace 
                                yc(n+1 : n+r_lengths(N)) = zeros(r_lengths(N), 1);
                                yr(n+1 : n+r_lengths(N)) = ones(r_lengths(N), 1);
                            case 1  % replace all
                                % don't do anything, reliable samples will be replaced
                            case 2  % shorten w
                                yc(n+1 : n+r_lengths(N)) = 1 - comp_short_transition(r_lengths(N), transition_type);
                                yr(n+1 : n+r_lengths(N)) = comp_short_transition(r_lengths(N), transition_type);
                            otherwise
                                error('Invalid ''treat_rel_short'' option!')
                        end  
                    end
                end
            end
            
            output = xc.*yc + x.*yr;
        
        
        case 4  % transition in the middle of the clipped and reliable part
            difference = diff(Mr);
            
            c_start_idxs = find(diff(Mr)==-1) + 1;
            c_stop_idxs = find(diff(Mr)==1);
            
            if length(c_stop_idxs) < length(c_start_idxs)
                c_stop_idxs = [c_stop_idxs; length(Mr)];
            elseif length(c_start_idxs) < length(c_stop_idxs)
                c_start_idxs = [1 c_start_idxs];
            end
            
            c_lengths = c_stop_idxs - c_start_idxs +1;
            r_lengths = [c_start_idxs; 0] - [0; c_stop_idxs] -1;
            r_lengths(end) = length(x) - c_stop_idxs(end);
            
            yc = double(Mr);
            yr = 1 - yc;
            
            N = 1;
            for n = 1:length(x)-1 

                if difference(n) == -1  % transition from reliable to clipped
                    limit = min(r_lengths(N), c_lengths(N)) + (r_lengths(N) - c_lengths(N) > 1);

                    if limit >= w  % reliable part is long enough
                        yc(n - ceil(w/2) + 1 : n + floor(w/2)) = transition_down;
                        yr(n - ceil(w/2) + 1 : n + floor(w/2)) = transition_up;
                    else
                        switch treat_short
                            case 0
                                % don't do anything, it will be solved in the transition from clipped to reliable
                            case 1
                                % don't do anything, reliable samples will be replaced
                            case 2
                                yc(n - ceil(limit/2) + 1 : n + floor(limit/2)) = 1 - comp_short_transition_up(limit, transition_type);
                                yr(n - ceil(limit/2) + 1 : n + floor(limit/2)) = comp_short_transition_up(limit, transition_type);
                        end
                    end
                elseif difference(n) == 1  % transition from clipped to reliable
                    N = N+1;
                    limit = min(r_lengths(N), c_lengths(N-1)) + (r_lengths(N) - c_lengths(N-1) > 1);
                    
                    if limit >= w % reliable part is long enough
                        yc(n - floor(w/2) + 1 : n + ceil(w/2)) = transition_up;
                        yr(n - floor(w/2) + 1 : n + ceil(w/2)) = transition_down;
                    else
                        switch treat_short
                            case 0
                                if r_lengths(N) < c_lengths(N-1)
                                    yc(n + 1 : n + r_lengths(N)) = zeros(r_lengths(N), 1);
                                    yr(n + 1 : n + r_lengths(N)) = ones(r_lengths(N), 1);
                                end
                            case 1
                                % don't do anything, reliable samples will be replaced
                            case 2
                                yc(n - floor(limit/2) + 1 : n + ceil(limit/2)) = comp_short_transition_up(limit, transition_type);
                                yr(n - floor(limit/2) + 1 : n + ceil(limit/2)) = 1 - comp_short_transition_up(limit, transition_type);
                        end
                    end
                end
            end
            
            output = xc.*yc + x.*yr;
    end
   
end


function transition = comp_short_transition(w, transition_type)

    if w == 1
        transition = 0.5;
    elseif w == 2
        transition = [0.5; 0.5];
    else
        if strcmp(transition_type, 'linear')
            transition_up = linspace(0, 1, floor(w/2)+2).';
            transition_up = transition_up(2:end-1);
            transition_down = 1 - transition_up;

        elseif strcmp(transition_type, 'soft')
            transition_up = sin(linspace(0, pi/2, floor(w/2)+2).').^2;
            transition_up = transition_up(2:end-1);
            transition_down = 1 - transition_up;
        end

        if mod(w,2) == 0
            transition = [transition_down; transition_up];
        else
            transition = [transition_down; 0; transition_up];
        end
    end

end


function transition = comp_short_transition_up(w, transition_type)
    
    if w == 1
        transition = 0.5;
    else
        if strcmp(transition_type, 'linear')
            transition = linspace(0, 1, w+2).';
            transition = transition(2:end-1);

        elseif strcmp(transition_type, 'soft')
            transition = sin(linspace(0, pi/2, w+2).').^2;
            transition = transition(2:end-1);
        end
    end

end
