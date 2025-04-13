clc;
clear;
close all;

%%  %%%%%%%%%%%%%%%%%%%%%%%%% Introduction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This file reorders element connectivity based on the three
% dimensions of the brick element to ensure that the third natural
% coordinate is always in the out-of-plane direction

% Option 1: reorder based on three dimensinon of element, taking zeta
% direction as direction with minimum dimension, useful if the thickness is
% smaller than in-plane dimension of individual element

% Option 2: reorder based on direction of local coordinates, useful if the
% in-plane dimension is smaller than out-of-plane direction, but requires
% knowledge on original mesh ordering

% Option 3: reorder based on refernece point, such point could be centre of
% circle or anything else, and the zeta diection is defined along element 
% centre to reference point, useful if in-plane dimension is smaller than out-of-plane direction
% and shell is curved or flat

% Dimension
% a: direction along node pair 1,4
% b: direction along node pair 1,2
% c: direction along node pair 1,5

%%  %%%%%%%%%%%%%%%%%%%%%%%%% Reading stage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method specification
option = 1;
direction = 'a'; % Only used for option 2
refcoord = [0,0,-1000];  % Only used for option 3

%%  %%%%%%%%%%%%%%%%%%%%%%%%% Reading stage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part reads nodal coordinate and element connectivities from target
% data file

% Open the file for reading
fid = fopen('Example.dat', 'r');
if fid == -1
    error('Cannot open the file');
end
i = 5;
nod_rdlflag = 0;           % Flag indicating reading node coordinate
elm_rdlflag = 0;           % Flag indicating reading element connectivity
node_coord = zeros(1,4);
elm_conn = zeros(1,9);

% Loop through the file to find the nodal coordiante and element
% connectivity

while ~feof(fid)
    line = fgetl(fid);

    if contains(line,'#')
         nod_rdlflag = 0;
         elm_rdlflag = 0;
    end

    if nod_rdlflag == 1
          cellArr = strsplit(line,'   '); % Split the string by commas
          node_coord(end+1,:) = str2double(cellArr); % Convert cell array to matrix
    end

    if elm_rdlflag == 1
          cellArr = strsplit(line,' '); % Split the string by commas
          cellArr = str2double(cellArr); % Convert cell array to matrix
          cellArr(isnan(cellArr)) = [];
          elm_conn(end+1,:) = cellArr;
    end

     if contains(line, 'nod.name     x   y   z')
       nod_rdlflag = 1;
     end

     if contains(line, 'elm.name')
       elm_rdlflag = 1;
     end
end

node_coord = node_coord(2:end,:);
node_coord(node_coord(:,4)>1e-8) = 0.0001;
elm_conn = elm_conn(2:end,:);


%%  %%%%%%%%%%%%%%%%%%%%%%%% Processing stage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reorder node order in element connectivity based on three average dimensions of
% the brick element

% Initialise
elm_conn_mod = zeros(size(elm_conn));

% 1st dimension
index_a = [2 3;
           1 4;
           6 7;
           5 8];
a1v = zeros(3,1);

% 2nd dimension
index_b = [1 2;
           4 3;
           5 6;
           8 7];
b1v = zeros(3,1);

% 3rd dimension
index_c = [2 6;
           1 5;
           3 7;
           4 8];
c1v = zeros(3,1);


for i = 1:size(elm_conn,1)

    % Compute three dimensions of the brick element
    a1 = 0; b1 = 0; c1 = 0;
    for j = 1:4
        a = elm_conn(i,index_a(j,1)+1);
        b = elm_conn(i,index_a(j,2)+1);
        a1v = (node_coord(a,2:4) - node_coord(b,2:4))';
        a1 = a1 + norm(a1v) / 4;

        a = elm_conn(i,index_b(j,1)+1);
        b = elm_conn(i,index_b(j,2)+1);
        b1v = (node_coord(a,2:4) - node_coord(b,2:4))';
        b1 = b1 + norm(b1v) / 4;

        a = elm_conn(i,index_c(j,1)+1);
        b = elm_conn(i,index_c(j,2)+1);
        c1v = (node_coord(a,2:4) - node_coord(b,2:4))';
        c1 = c1 + norm(c1v) / 4;
    end

    switch option
        case 1
            % Reorder depending on smallest dimension of the brick element
            if (b1<a1) && (b1<c1)
                index_seq = [1,4,8,5,2,3,7,6];
            elseif (a1<c1) && (a1<b1)
                index_seq = [4,3,7,8,1,2,6,5];
            else
                index_seq = 1:8;
            end

        case 2
            % Reorder based on specified direction, e.g. b denotes
            % direction along fibre b is zeta direction
            switch direction
                case 'b'
                    index_seq = [1,4,8,5,2,3,7,6];
                case 'a'
                    index_seq = [4,3,7,8,1,2,6,5];
                case 'c'
                    index_seq = 1:8;
                otherwise
                    error('Specify valid direction')
            end

        case 3
            % Reorder based on reference point
            % Compute centre of element
            a = elm_conn(i,2:9);
            d1 = [sum(node_coord(a,2)); sum(node_coord(a,3)); sum(node_coord(a,4))]/8;

            % Compute vector connecting element centre and reference point
            d1v = d1 - refcoord;

            % Compute dot product with the three dimension length connected
            % to point 1 and reorder based on which dot product is max
            a1v = (node_coord(elm_conn(i,2),2:4) - node_coord(elm_conn(i,5),2:4))';
            b1v = (node_coord(elm_conn(i,2),2:4) - node_coord(elm_conn(i,3),2:4))';
            c1v = (node_coord(elm_conn(i,2),2:4) - node_coord(elm_conn(i,6),2:4))';
            a1 = abs(d1v' * a1v); b1 = abs(d1v' * b1v); c1 = abs(d1v' * c1v);
            if a1 == max([a1,b1,c1])
                index_seq = [4,3,7,8,1,2,6,5];
            elseif b1 == max([a1,b1,c1])
                index_seq = [1,4,8,5,2,3,7,6]; 
            else
                index_seq = 1:8;
            end

    end 
    
    % Stores modified element connectivity
    elm_conn_mod(i,1) = i;
    for j = 1:8
        a = index_seq(j);
        elm_conn_mod(i,j+1) = elm_conn(i,a+1);
    end

end

% Close the file
fclose(fid);
%%  %%%%%%%%%%%%%%%%%%%%%%%%% Writing stage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rewrite element connectivity of the data file
fid2 = fopen('new element conn.txt', 'w');
writematrix(elm_conn_mod,'new element conn.txt','Delimiter',',');
