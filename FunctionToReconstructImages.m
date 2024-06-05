

function [final_image, final_recon, num_total_patterns,N] = FunctionToReconstructImages (N, imagemeasurements) 


%% SPECIFY THE NUMBER OF TOTAL MEASUREMENTS THAT WERE TAKEN (POSITIVE + NEGATIVE PATTERN MEASUREMENTS)

num_total_patterns =2*N^2; % N^2 patterns times 2, because there is a positive and negative measurement for each pattern (differential method)


%% CONVERT THE TABLE WITH ALL THE MEASUREMENTS TO AN ARRAY

array_total_measurements = table2array(readtable(imagemeasurements, 'Sheet', 'Result')); 

%% STORE BACKGROUND, WHITE LIGTH, POSITIVE, AND NEGATIVE MEASUREMENTS IN SEPARATE ARRAYS

background=mean(array_total_measurements(2:101, 2:end)); % Rows 2 to 100 correspond to background measurements
whitelight=mean(array_total_measurements(102:201, 2:end)); % Rows 101:200 correspond to white light measurements
pos_measurements=zeros(num_total_patterns/2, 512); % 512 are the number of intensities at different wavelengths that TOKUSPEC captures for each pattern
neg_measurements=zeros(num_total_patterns/2, 512);

pos_index=0; % To store positive measurements in order
neg_index=0; % To store negative measurements in order

for i=202:(2*N^2+201)
    if mod(i, 2) == 0 % Measurements in even rows correspond to the positive patterns
        pos_index=pos_index + 1;
        pos_measurements(pos_index,:)=array_total_measurements(i, 2:end);

    else % Measurements in odd rows correspond to the negative patterns
        neg_index=neg_index + 1;
        neg_measurements(neg_index,:)=array_total_measurements(i, 2:end);
    end
end


%% SUBSTRACT BACKGROUND TO ALL POSITIVE AND NEGATIVE DATA AND DIVIDE BY THE WHITE LIGHT CONTRIBUTION OF THE PROJECTOR 

for i=1:num_total_patterns/2
    pos_measurements(i,:) = (pos_measurements(i,:) - background)./(whitelight + 1) ;
    neg_measurements(i,:) = (neg_measurements(i,:) - background)./(whitelight + 1);
end

%% MAKE NEGATIVE MEASUREMENTS EQUAL TO 0 SINCE THEY DO NOT CONTAIN ANY INFORMATION

for i=1:num_total_patterns/2
    for j=1:512 % Go through all the measurements and, if they are negative, make them equal to 0
        if pos_measurements(i,j) < 0
            pos_measurements(i,j) = 0;
        end
        if neg_measurements(i,j) < 0
            neg_measurements(i,j) = 0;
        end
    end
end

%% INTERPOLATE DATA FROM WAVELENGTHS 430-829 nm IN INTERVALS OF 1 nm

range=430:1:829;
measurements=length(range); % Number of measurements that will result from the interpolation

interpolated_pos_measurements=zeros(num_total_patterns/2, measurements);
interpolated_neg_measurements=zeros(num_total_patterns/2, measurements);

for i = 1:num_total_patterns/2
    interpolated_pos_measurements(i,:)=interp1(array_total_measurements(1,2:end), pos_measurements(i,:), range, 'spline');
    interpolated_neg_measurements(i,:)=interp1(array_total_measurements(1,2:end), neg_measurements(i,:), range, 'spline');
end


%% SUBSTRACT THE POSITIVE MEASUREMENTS MINUS THEIR CORRESPONDING NEGATIVE MEASUREMENTS

substracted_measurements=zeros(num_total_patterns/2, measurements);

for i=1:num_total_patterns/2
    substracted_measurements(i,:) = interpolated_pos_measurements(i,:) - interpolated_neg_measurements(i,:);
end


%% CREATE A HYPERCUBE FOR EACH PATTERN AND STORE THEM IN A CELL ARRAY

% Hadamard matrix from which we obtain each pattern
had=hadamard(N^2); 

% Cell array to store the hypercubes belonging to each pattern
hypercubes=cell(N*N, measurements); 

for i=1:N*N
    W=had(i,:);
    pattern=reshape(W,N,N);

    hypercube=zeros(N, N, measurements); % To go through all the measurements stored in final_added_measurements. Row by Row (N^2 in total)
       
    for k=1:N
        for t=1:N
            for z=1:measurements
                hypercube(k, t, z) = pattern(k,t).*substracted_measurements(i, z); %Fill the hypercube with measurements
            end
        end
    end


    for z=1:measurements
        hypercubes{i,z}=hypercube(:, :, z);% Store the current cube in the cell array that contains all the hypercubes
    end

end

%% FINAL IMAGE (FINAL HYPERSPECTRAL IMAGE) OBTAINED BY ADDING ALL THE STORED HYPERCUBES COMPONENT BY COMPONENT

final_image = zeros(N, N, measurements);

% Iterate over each cell
for i = 1:N*N
    for z=1:measurements
        final_image(:,:,z) = final_image(:,:,z) + hypercubes{i,z}; % Sum the data in each cell to reconstruct the final image
    end
end


%% FINAL IMAGE (FINAL HYPERSPECTRAL IMAGE) OBTAINED BY DOING HADAMARD INVERSE

final_recon=cell(1,measurements);
for z=1:measurements
    final_recon_img=had\substracted_measurements(:,z); % Divide the final measurements by the original Hadamard matrix to reconstruct the final image
    final_recon{z}=reshape(final_recon_img,N,N); % Reshape to obtain a squared image
end



end