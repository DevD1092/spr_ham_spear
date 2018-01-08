% This script is to test the performance of the proposed technique in the scenario of false stars present in the image.


clear all;

clc;

close all;

clear classes;


% Initialize camera parameter

% SI STS parameters

FOV = 15;

img_height = 1024;

img_width = 1024;

pixel_size = 13.3/1024;

f = (img_height)*pixel_size /2/ tand(FOV/2);

% g=4;
%
% pr=ceil(384/g);

img_dimension=[img_width img_height];

%Load SC

catalog_path='spr_ham_spear/simulate/SKY2000_Magnitude6_doublestars_0.12.txt';

% Initialize camera attitude

int=2; %2.5: 10368, 3.75: 4608, 5:2592, 7.5: 1152

RAini=1;%from [0, 360];

RAend=359;

RA=[RAini:int:RAend]'; %row-wise

DECini=-89; %from(-90,90);

DECend=89;

DEC=[DECini:int:DECend]; %Column-wise

angle=45;

var_test = 0;

%feature_vector tolerance

% tol=1;

% fn_threshold=0.7;


%Simulator parameters

cent_variance=0;

no_ran_star=0;

SNR=10; % SNR = 10 : ideal.

background_noise=0.0; % background_noise = 0.0 : ideal

PSF_set=0;          % Change this setting to (1,2 or 3) if you want to calculate the centroid. PSF_set = 0 means exact co-ordinate of the star given with no magnitude consideration.

no_miss = 0;

%centroider parameters

thresh=0.3;


% Read star coordinates in Earth reference frame from star catalog

[SKYMAP_No,star_RA,star_DEC,star_MAG]= textread(catalog_path,'%d %f %f %f');


% Si is coordiante of star in Earth reference frame,

% the 3 column are X, Y,and Z

Si = [cosd(star_DEC).*cosd(star_RA) cosd(star_DEC).*sind(star_RA)  sind(star_DEC)];

% for i=1: length(star_RA)

%     % Convert RA, and DEC into ECI unit vector

%     ECI_vector =[cosd(star_DEC(i))* cosd(star_RA(i)) cosd(star_DEC(i))* sind(star_RA(i))  sind(star_DEC(i))];

%     Si(i,:)= ECI_vector;

% end


catalog=struct('SKYMAP_No',SKYMAP_No,'star_RA',star_RA,'star_DEC',star_DEC,'star_MAG',star_MAG,'Si',Si);

% Testing   variables
acc = 0;
match_acc=0;
count_acc = 0;
time_elapsed = [];
Star_ref_empty_cases = 0;
Star_ref_empty_cases_list = [];
not_successful_RA_empty = [];
not_successful_DEC_empty = [];
multiple_star_id = 0;
failed_images_number_of_stars = [];
false_match_RA = [];
false_match_DEC = [];
avg_number_star_return = [];
missing_stars_list = [];
x = 1;
bin_size = 512;
sp_pho_top_50_2_count = [];
sp_pho_top_50_max = [];
sp_pho_id = [];
match_check = [];

% starID=zeros(size(RA,1),size(DEC,2));

% star_q_fov=zeros(size(RA,1),size(DEC,2));

% star_q_ad=zeros(size(RA,1),size(DEC,2));

% star_num_returned_match=zeros(size(RA,1),size(DEC,2));

% starnum_ad_frm_sra=zeros(size(RA,1),size(DEC,2));

% starnum_fov_frm_ad=zeros(size(RA,1),size(DEC,2));

%% THE LOOP FOR GENERATING IMAGES

for ii = 1:size(RA,1)
    
  for jj = 1:size(DEC,2)
        
        id=[];
        
        fn=[];
        
        temp=[];
        
                
        % Generate sky image at predetermined attitude
        
        Reci2body= Convert_Axis_2_AttitudeMatrix(RA,DEC,angle);
        
        [star_matrix, I]= Plot_sky_images(Reci2body, FOV, img_height, img_width, pixel_size,cent_variance, no_ran_star, SNR, background_noise, PSF_set,catalog);
        
        [row_star_matrix,col_star_matrix] = size(star_matrix);
        
        for x = 1:row_star_matrix
            id = [id star_matrix(x,1)];
        end
        
        id_centroid  =[];
        
        for x = 1 : size(star_matrix,1)
            id_centroid = [id_centroid ; star_matrix(x,1) star_matrix(x,11) star_matrix(x,12)];
        end
        
        len_id = length(id);
        
        if (len_id > 3)
            
            tic
            
%             imshow(I);
            
            % --- This is for exact co-ordinates of the star with PSF_setting = 0 ----
            
            id = [];
            centroid = [];
            Sc = [];
            count_valid_stars = 0;
            valid_stars_list = [];
            for i = 1 : len_id
                if(star_matrix(i,4) <= 6.0)
                    count_valid_stars = count_valid_stars + 1;
                    valid_stars_list = [valid_stars_list i];
                    id = [id  star_matrix(i,1)];
                end
            end
            
            id
            
            missing_stars_list = [missing_stars_list (len_id - count_valid_stars)];
            
            for i = 1 : count_valid_stars
                centroid(i,1) = star_matrix(valid_stars_list(i),11);
                centroid(i,2) = star_matrix(valid_stars_list(i),12);
                Sc(i,1) = star_matrix(valid_stars_list(i),8);
                Sc(i,2) = star_matrix(valid_stars_list(i),9);
                Sc(i,3) = star_matrix(valid_stars_list(i),10);
            end

% --- This is for exact co-ordinates of the star with PSF_setting = 0 ----
            
%             if(length(centroid) <= 2)
%                 continue;
%             end
%             
            % choose a reference star that is nearest to the center
            d=zeros(size(centroid,1),1);
            for i=1:length(centroid)
                temp= centroid(i,:)-[img_height/2 img_width/2];
                d(i)= sum(temp.*temp);
            end
            k= find(d==min(d));
            k=k(1);            
            
            % Construct the new id            
            
            % Co-ordinates of the star nearest to the center of the image.
            % Shifting the co-oridnates
            Star_ref = [];
            centroid_new = [];
            centroid_newer = [];
            dist_new = [];
            center_star_x_cord = (centroid(k,1));
            center_star_y_cord = (centroid(k,2));
            shift_x_cord = abs((img_height/2) - center_star_x_cord);
            shift_y_cord = abs((img_height/2) - center_star_y_cord);
            j = 1;
            for i = 1 : length(centroid)
                if(i ~= k)
                    if(center_star_x_cord < 512)
                        centroid_new(i,1) = centroid(i,1) + shift_x_cord;
                    end
                    if (center_star_x_cord > 512)
                        centroid_new(i,1) = centroid(i,1) - shift_x_cord;
                    end
                    if(center_star_y_cord < 512)
                        centroid_new(i,2) = centroid(i,2) + shift_y_cord;
                    end
                    if(center_star_y_cord > 512)
                        centroid_new(i,2) = centroid(i,2) - shift_y_cord;
                    end
                    if(center_star_x_cord == 512)
                        centroid_new(i,1) = centroid(i,1);
                    end
                    if(center_star_y_cord == 512)
                        centroid_new(i,2) = centroid(i,2);
                    end
                end
                
                if(i == k)
                    centroid_new(i,1) = 512;
                    centroid_new(i,2) = 512;
                end
                
                dist_new(i) = sqrt((centroid_new(i,1) - 512) ^ 2 + (centroid_new(i,2) - 512) ^ 2);
                
                if(centroid_new(i,1) >= 0 && centroid_new(i,1) <= 1024)
                    if(centroid_new(i,2) >= 0 && centroid_new(i,2) <= 1024)
                        if(dist_new(i) < 510)     % Just to be on the safe side.
                            centroid_newer(j,1) = centroid_new(i,1);
                            centroid_newer(j,2) = centroid_new(i,2);
                            j = j + 1;
                        end
                    end
                end
            end
            
            
            if(length(centroid_newer) > 4)
                % Listing the angles of the stars in the image the same anticlockwise order
                % These angles are of the stars in the image.
                angle_hor_list = [];
                sorted_angle_hor_list = [];
                angle_final_list = [];
                angle_rad_list = [];
                sorted_angle_final_list = [];
                dist_list = [];
                dist_list_sorted = [];
                dist_list_sorted_final = [];
                Image_vector = [];
                
                for i=1: length(centroid_newer)
                    
                    if(centroid_newer(i,1) ~= (img_height/2) && centroid_newer(i,2) ~= (img_height/2))
                        
                        x_cord = centroid_newer(i,1);
                        y_cord = centroid_newer(i,2);
                        
                        % x_cord corresponds to the row number of the star in the image.
                        % y_cord corresponds to the column number of the star in the image.
                        
                        dim_x = abs(x_cord - (img_height/2));
                        dim_y = abs(y_cord - (img_height/2));
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
                            if(y_cord > 0 && y_cord <  (img_height/2))
                                angle_hor = 180 - angle_hor;
                                angle_hor_list = [angle_hor_list angle_hor];
                            end
                        end
                        
                        % Star lies in the third quadrant.
                        if( x_cord >  (img_height/2) && x_cord < img_height)
                            if( y_cord > 0 && y_cord <  (img_height/2))
                                angle_hor = 180 + angle_hor;
                                angle_hor_list = [angle_hor_list angle_hor];
                            end
                        end
                        
                        % Star lies in the fourth quadrant.
                        if( x_cord >  (img_height/2) && x_cord < img_height)
                            if( y_cord >  (img_height/2) && y_cord < (img_height))
                                angle_hor = 360 - angle_hor;
                                angle_hor_list = [angle_hor_list angle_hor];
                            end
                        end
                        
                        % Covering the extreme conditions: 0,90,180,270.
                        
                        if(x_cord ==  (img_height/2) && y_cord >  (img_height/2))
                            angle_hor_list = [angle_hor_list 0];
                        end
                        
                        if(x_cord ==  (img_height/2) && y_cord <  (img_height/2))
                            angle_hor_list = [angle_hor_list 180];
                        end
                        
                        if(y_cord ==  (img_height/2) && x_cord <  (img_height/2))
                            angle_hor_list = [angle_hor_list 90];
                        end
                        
                        if(y_cord ==  (img_height/2) && x_cord >  (img_height/2))
                            angle_hor_list = [angle_hor_list 270];
                        end
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
                    angle_final_list(l) =  angle_between;
                end
                
                count_dist = 1;
                for ll = 1 : length(sorted_angle_hor_list)
                    for mm = 1 : length(sorted_angle_hor_list)
                        if(dist_list(mm,2) == sorted_angle_hor_list(ll))
                            dist_list_sorted = [dist_list_sorted dist_list(mm,1)];
                            count_dist = count_dist + 1;
                        end
                    end
                end
                
                % Arranging all the angles in an increasing manner.
                
                for ll = 1 : length(angle_final_list)
                    angle_rad_list = [angle_rad_list ; angle_final_list(ll) dist_list_sorted(ll)];
                end
                
                sorted_angle_final_list = sort(angle_final_list);
                
                for ll = 1 : length(sorted_angle_final_list)
                    for mm = 1 : length(sorted_angle_final_list)
                        if(angle_rad_list(mm,1) == sorted_angle_final_list(ll))
                            dist_list_sorted_final = [dist_list_sorted_final angle_rad_list(mm,2)];
                        end
                    end
                end
                
                % Matching
                
                for ll = 1 : bin_size
                    count = 0;
                    a(ll) = 0;
                    for mm = 1 : length(sorted_angle_final_list)
                        if((180/bin_size)*(ll-1) < sorted_angle_final_list(mm) && sorted_angle_final_list(mm) < (180/bin_size)*ll)
                            a(ll) = a(ll) + dist_list_sorted_final(mm);
                            count = count + 1;
                        end
                    end
                    if(a(ll) ~= 0)
                        a(ll) = a(ll) / count;
                    end
                end
                Image_vector = [Image_vector ; a];               
                
                
                aa = Image_vector;
                std_im = std(aa,1);
                cc = [ zeros(1,length(aa)-1) aa];
                ee = fft(cc);
                
                SPD = dlmread('spr_ham_spear/SPD/SPD.txt');                
                
                cross_corr = [];
                sorted_cross_corr = [];
                cross_corr_top_10 = [];
                diff = [];
                sp_pho = [];
                sp_pval =[];
                sp_pho_top_50 = [];
                sp_pho_pval = [];
                sp_pho_top_50_2 = [];
                sp_pho_pval_top_50 = [];
                
                for comp_ind = 1 : length(SPD)                    
                    SPD_vector = SPD(comp_ind,2: length(SPD(comp_ind,:)));
                    bb = SPD_vector;
                    [sp_pho_1 sp_pval_1] = corr(aa',bb','Type','Spearman');
                    sp_pho = [sp_pho sp_pho_1];
                    sp_pval = [sp_pval sp_pval_1];
                    sp_pho_pval = [sp_pho_pval ; sp_pho_1 sp_pval_1];
                end
            end            
            
            sp_pho_sorted = sort(sp_pho,'descend');
            
            for ind_pho = 1 : 10
                sp_pho_top_50 = [sp_pho_top_50 sp_pho_sorted(ind_pho)];
            end
            
            for ind_pho = 1 : length(sp_pho_top_50)
                for ind_pval = 1 : length(sp_pho_pval)
                    if(sp_pho_pval(ind_pval,1) == sp_pho_top_50(ind_pho))
                        sp_pho_pval_top_50 = [sp_pho_pval_top_50; sp_pval(ind_pval) ind_pval];
                    end
                end
            end
            
            sp_pho_top_50_2 = sp_pho_top_50;
            
            sp_pho_top_50_max = [sp_pho_top_50_max ; max(sp_pho_top_50)];
            sp_pho_id = [sp_pho_id ; sp_pho(id(k))];
            
            toc
            
            if(ismember(sp_pho(id(k)),sp_pho_top_50_2) == 1)
                match = 1;
                match_acc = match_acc + 1;
            else
                false_match_RA = [false_match_RA RA(ii)];
                false_match_DEC = [false_match_DEC DEC(jj)];
            end        
            
            match_check = [match_check ; match];      
            
           count_acc = count_acc + 1;            
        end
        if(count_acc == 16200)
            return;
        end
        pause(0.2);
        
     end   
end



























