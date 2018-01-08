% This script is to generate the SPD for the propsoed technique


clear all;
clc;
%% Initialize camera parameter
img_height = 1024;
img_width = 1024;
pixelsize= 13.3/1024;
FOV =15;
angle = 45;
test = 0;
bin_size = 512;

file_path='spr_ham_spear/simulate/SKY2000_Magnitude6_doublestars_0.12.txt';
[SKYMAP_No,star_RA,star_DEC,star_MAG]= textread(file_path,'%d %f %f %f');

% Si is coordiante of star in Earth reference frame,
% the 3 column are X, Y,and Z
no_stars=length(star_RA);
Si = zeros(no_stars, 3);
catalog = zeros(no_stars, 4);

for i=1: no_stars
    % Convert RA, and DEC into ECI unit vector
    ECI_vector =[cosd(star_DEC(i))* cosd(star_RA(i)) cosd(star_DEC(i))* sind(star_RA(i))  sind(star_DEC(i))];
    Si(i,:)= ECI_vector;
end

% SPD stored
SPD_matrix = [];

for i=1: no_stars
    
    RA= star_RA(i);
    DEC= star_DEC(i);
    angle_hor_list = [];
    star_ids = [];
    star_ids_sorted= [];
    angle_final_list = [];
    sorted_angle_final_list = [];
    angle_rad_list = [];
    dist_list= [];
    dist_list_sorted= [];
    dist_list_sorted_final = [];
    a = [];
    
    
    C= Convert_Axis_2_AttitudeMatrix(RA,DEC,angle);
    [R_camera_to_earth,star_matrix]= Find_neighbor_star_FOV(C,FOV, img_height, img_width, pixelsize);
    
    
    for j = 1 : size(star_matrix,1)
        
        
        if(star_matrix(j,1) ~= i)
            
            x_cord = star_matrix(j,11);
            y_cord = star_matrix(j,12);
            
            % x_cord corresponds to the row number of the star in the image.
            % y_cord corresponds to the column number of the star in the image.
            
            dim_x = abs(x_cord - (img_height/2));
            dim_y = abs(y_cord - (img_width/2));
            dist = sqrt((dim_x).^ 2 + (dim_y).^ 2);
            angle_hor = atand(dim_x/dim_y);
            
            % See the documentation for quadrant information.
            % Star lies in the first quadrant.
            if( x_cord > 0 && x_cord < (img_height/2))
                if( y_cord > (img_height/2) && y_cord < (img_height))
                    angle_hor_list = [angle_hor_list angle_hor];
                end
            end
            
            % Star lies in the second quadrant.
            if(x_cord > 0 && x_cord < (img_height/2))
                if(y_cord > 0 && y_cord < (img_height/2))
                    angle_hor = 180 - angle_hor;
                    angle_hor_list = [angle_hor_list angle_hor];
                end
            end
            
            % Star lies in the third quadrant.
            if( x_cord > (img_height/2) && x_cord < img_height)
                if( y_cord > 0 && y_cord < (img_height/2))
                    angle_hor = 180 + angle_hor;
                    angle_hor_list = [angle_hor_list angle_hor];
                end
            end
            
            % Star lies in the fourth quadrant.
            if( x_cord > (img_height/2) && x_cord < img_height)
                if( y_cord > (img_height/2) && y_cord < (img_height))
                    angle_hor = 360 - angle_hor;
                    angle_hor_list = [angle_hor_list angle_hor];
                end
            end
            
            % Covering the extreme conditions: 0,90,180,270.
            
            if(x_cord == (img_height/2) && y_cord > (img_height/2))
                angle_hor_list = [angle_hor_list 0];
            end
            
            if(x_cord == (img_height/2) && y_cord < (img_height/2))
                angle_hor_list = [angle_hor_list 180];
            end
            
            if(y_cord == (img_height/2) && x_cord < (img_height/2))
                angle_hor_list = [angle_hor_list 90];
            end
            
            if(y_cord == (img_height/2) && x_cord > (img_height/2))
                angle_hor_list = [angle_hor_list 270];
            end
            star_ids = [star_ids; star_matrix(j,1) angle_hor];
            dist_list = [dist_list; dist angle_hor];
        end
    end
    
    sorted_angle_hor_list = sort(angle_hor_list);
    
    for l = 1 : length(sorted_angle_hor_list)
        
        if(l ~= length(sorted_angle_hor_list))
            angle_between = sorted_angle_hor_list(l+1) - sorted_angle_hor_list(l);
        end
        
        if(l == length(sorted_angle_hor_list))
            angle_between = 360 - sorted_angle_hor_list(l) + sorted_angle_hor_list(1);
        end
        
        angle_final_list = [angle_final_list angle_between];
    end
    
    % Preparing the star_ids SPD.
    count = 1;
    sorted_angle_hor_list = sort(angle_hor_list);
    
    for jj = 1 : length(sorted_angle_hor_list)
        for kk = 1 : (size(star_matrix,1) - 1)
            if(star_ids(kk,1) ~= i)
                if(star_ids(kk,2) == sorted_angle_hor_list(jj))
                    star_ids_sorted = [star_ids_sorted star_ids(kk,1)];
                    count = count + 1;
                end
            end
        end
    end
    
    %Preparing the radial distance SPD.
    count_dist = 1;
    for ll = 1 : length(sorted_angle_hor_list)
        for mm = 1 : (size(star_matrix,1) - 1)
            if(dist_list(mm,2) == sorted_angle_hor_list(ll))
                dist_list_sorted = [dist_list_sorted dist_list(mm,1)];
                count_dist = count_dist + 1;
            end
        end
    end
    
    % Arranging all the radial distances in an increasing manner.
    
    for ll = 1 : length(angle_final_list)
        angle_rad_list = [angle_rad_list ; angle_final_list(ll) dist_list_sorted(ll)];
    end
    
    dist_list_sorted_final = sort(dist_list_sorted);
    
    for ll = 1 : length(dist_list_sorted_final)
        for mm = 1 : length(dist_list_sorted_final)
            if(angle_rad_list(mm,2) == dist_list_sorted_final(ll))
                sorted_angle_final_list = [sorted_angle_final_list angle_rad_list(mm,1)];
            end
        end
    end
    
    bin_size = 512; % Number of bins
    for ll = 1 : bin_size
        count = 0;
        a(ll) = 0;
        for mm = 1 : length(dist_list_sorted_final)
            if((512/bin_size)*(ll-1) < dist_list_sorted_final(mm) && dist_list_sorted_final(mm) < (512/bin_size)*ll)
                a(ll) = a(ll) + sorted_angle_final_list(mm);
                count = count + 1;
            end
        end
        if(a(ll) ~= 0)
            a(ll) = a(ll) / count;
        end
    end    
    SPD_matrix = [SPD_matrix ; a];
end
