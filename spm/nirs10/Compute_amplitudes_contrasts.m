%Stats for SD - Patient 1
TOPO = [];
load('TOPO.mat');
v = 6; %left view
%TOPO.v{v}.s


%h=1;
%c=1;
%number of sessions
ns = length(TOPO.v{v}.s);
%number of contrasts divided by 2 (for positive and negative contrasts)
nc = length(TOPO.v{v}.group.hb{1}.c)/2;
beta_max = zeros(nc,3,ns+1);
beta_min = zeros(nc,3,ns+1);
beta_max1 = zeros(3,ns+1);
beta_min1 = zeros(3,ns+1);
beta_max2 = zeros(3,ns+1);
beta_min2 = zeros(3,ns+1);
c_max = zeros(nc,3,3); %triplets row, col, val
c_min = zeros(nc,3,3);

%TOPO.v{v}.group.hb{h}.c{c}.Tmap(row,col)
%loop over chromophore
for h=1:3
    %loop over contrasts 
    %for c=1:nc
    c=1;
    h_choice_max = 2; %
    h_choice_min = 2; %look at HbR
    %before:
    %h_choice_max = h;
    %h_choice_min = h;
    
        %careful, in group, there are twice as many contrasts alternating
        %between positive and negative ones
        [val1 col_max] = max(max(TOPO.v{v}.group.hb{h_choice_max}.c{2*c-1}.Tmap,[],1));
        [val2 row_max] = max(max(TOPO.v{v}.group.hb{h_choice_max}.c{2*c-1}.Tmap,[],2));
        if val1==val2
            val_max = val1;
        else
            disp('Problem with val_max');
        end
        [val1 col_min] = min(min(TOPO.v{v}.group.hb{h_choice_min}.c{2*c-1}.Tmap,[],1));
        [val2 row_min] = min(min(TOPO.v{v}.group.hb{h_choice_min}.c{2*c-1}.Tmap,[],2));
        if val1==val2
            val_min = val1;
        else
            disp('Problem with val_min');
        end
        %store
        c_max(1,h,1) = row_max;
        c_max(1,h,2) = col_max;
        c_max(1,h,3) = val_max;
        c_min(1,h,1) = row_min;
        c_min(1,h,2) = col_min;
        c_min(1,h,3) = val_min;
        c_max(2,h,1) = row_max;
        c_max(2,h,2) = col_max;
        c_max(2,h,3) = TOPO.v{v}.group.hb{h}.c{2*c+1}.Tmap(row_max,col_max);
        c_min(2,h,1) = row_min;
        c_min(2,h,2) = col_min;
        c_min(2,h,3) = TOPO.v{v}.group.hb{h}.c{2*c+1}.Tmap(row_min,col_min);
        %loop over sessions
        for s=1:ns
            beta_max(1,h,s) = TOPO.v{v}.s{s}.hb{h}.c_interp_beta(1,row_max,col_max);
            %at the site of max activation, give 
            %the second volterra at the same location
            beta_max(2,h,s) = TOPO.v{v}.s{s}.hb{h}.c_interp_beta(2,row_max,col_max);
            
            beta_min(1,h,s) = TOPO.v{v}.s{s}.hb{h}.c_interp_beta(1,row_min,col_min);
            beta_min(2,h,s) = TOPO.v{v}.s{s}.hb{h}.c_interp_beta(2,row_min,col_min);
            
            beta_max1(h,s) = TOPO.v{v}.s{s}.hb{h}.c_interp_beta(1,row_max,col_max);
            beta_max2(h,s) = TOPO.v{v}.s{s}.hb{h}.c_interp_beta(2,row_max,col_max);            
            beta_min1(h,s) = TOPO.v{v}.s{s}.hb{h}.c_interp_beta(1,row_min,col_min);
            beta_min2(h,s) = TOPO.v{v}.s{s}.hb{h}.c_interp_beta(2,row_min,col_min);
        end
        beta_max(1,h,ns+1) = TOPO.v{v}.group.hb{h}.c{1}.beta_group(row_max,col_max);
        beta_min(1,h,ns+1) = TOPO.v{v}.group.hb{h}.c{1}.beta_group(row_min,col_min);
        beta_max(2,h,ns+1) = TOPO.v{v}.group.hb{h}.c{3}.beta_group(row_max,col_max);
        beta_min(2,h,ns+1) = TOPO.v{v}.group.hb{h}.c{3}.beta_group(row_min,col_min);
        beta_max1(h,ns+1) = TOPO.v{v}.group.hb{h}.c{1}.beta_group(row_max,col_max);
        beta_min1(h,ns+1) = TOPO.v{v}.group.hb{h}.c{1}.beta_group(row_min,col_min);
        beta_max2(h,ns+1) = TOPO.v{v}.group.hb{h}.c{3}.beta_group(row_max,col_max);
        beta_min2(h,ns+1) = TOPO.v{v}.group.hb{h}.c{3}.beta_group(row_min,col_min);
    %end
end