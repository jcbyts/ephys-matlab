function B = concomptest(A)
%// Step #1
visited = false(size(A));
[rows,cols] = size(A);
B = zeros(rows,cols);
ID_counter = 1;

%// Step 2
%// For each location in your matrix...
for row = 1 : rows
    for col = 1 : cols
        %// Step 2a
        %// If this location is not 1, mark as visited and continue
        if A(row,col) == 0
            visited(row,col) = true;

        %// Step 2b
        %// If we have visited, then continue
        elseif visited(row,col)
            continue;

        %// Step 2c
        %// Else...
        else
            %// Step 3
            %// Initialize your stack with this location
            stack = [row col];

            %// Step 3a
            %// While your stack isn't empty...
            while ~isempty(stack)
                %// Step 3b
                %// Pop off the stack
                loc = stack(1,:);
                stack(1,:) = [];

                %// Step 3c
                %// If we have visited this location, continue
                if visited(loc(1),loc(2))
                    continue;
                end

                %// Step 3d
                %// Mark location as true and mark this location to be
                %// its unique ID
                visited(loc(1),loc(2)) = true;
                B(loc(1),loc(2)) = ID_counter;

                %// Step 3e
                %// Look at the 8 neighbouring locations
                [locs_y, locs_x] = meshgrid(loc(2)-1:loc(2)+1, loc(1)-1:loc(1)+1);
                locs_y = locs_y(:);
                locs_x = locs_x(:);

                 %%%% USE BELOW IF YOU WANT 4-CONNECTEDNESS
                 % See bottom of answer for explanation
                 %// Look at the 4 neighbouring locations
                 % locs_y = [loc(2)-1; loc(2)+1; loc(2); loc(2)];
                 % locs_x = [loc(1); loc(1); loc(1)-1; loc(1)+1];

                %// Get rid of those locations out of bounds
                out_of_bounds = locs_x < 1 | locs_x > rows | locs_y < 1 | locs_y > cols;

                locs_y(out_of_bounds) = [];
                locs_x(out_of_bounds) = [];

                %// Step 3f
                %// Get rid of those locations already visited
                is_visited = visited(sub2ind([rows cols], locs_x, locs_y));

                locs_y(is_visited) = [];
                locs_x(is_visited) = [];

                %// Get rid of those locations that are zero.
                is_1 = A(sub2ind([rows cols], locs_x, locs_y));
                locs_y(~is_1) = [];
                locs_x(~is_1) = [];

                %// Step 3g
                %// Add remaining locations to the stack
                stack = [stack; [locs_x locs_y]];
            end

            %// Step 4
            %// Increment counter once complete region has been examined
            ID_counter = ID_counter + 1;
        end
    end %// Step 5
 end   