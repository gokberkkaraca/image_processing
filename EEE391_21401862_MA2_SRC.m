% Read the image from the file
image_rgb = imread('fruits.png');
image = im2double(rgb2gray(image_rgb));
figure('Name','Gray Scale Image');
imshow(image, []);

% Apply Fourier Transform to the image
FT_image = FourierTransform(image);

% Plot magnitude and phase of the image
FT_magnitude = log(1+abs(FT_image));
FT_angle = angle(FT_image);

figure('Name', 'Magnitude');
imshow(FT_magnitude, []);

figure('Name', 'Angle');
imshow(FT_angle, []);

% --- QUESTION 1 --- 

% Creating and appliyng low pass filter
LP = zeros(size(image));
bandwidth = pi/4;
for i = 1:size(image,1)
    for j = 1:size(image,2)
        wx = -pi + 2*pi*(j-1)/512;
        wy = -pi + 2*pi*(i-1)/512;
        wx_condition = -bandwidth <= wx && wx <= bandwidth;
        wy_condition = -bandwidth <= wy && wy <= bandwidth;
        LP(i,j) = wx_condition && wy_condition;
    end
end

Filtered_Image = InverseFourierTransform(LP.*FT_image);
figure('Name','Low Pass Filter');
imshow(real(Filtered_Image),[]);

% Creating and appliyng high pass filter
HP = ones(512,512) - LP;
Filtered_Image = InverseFourierTransform(HP.*FT_image);
figure('Name','High Pass Filter');
imshow(real(Filtered_Image),[])

% --- QUESTION 2 ---
A = image;

% Calculate A*hx
hx = [0, 0, 0; 1, 0, -1; 0, 0, 0];
res = my_conv2(A, hx);
figure('Name','A*hx');
imshow(res, []);

% Calculate A*hy
hy = [0, 1, 0; 0, 0, 0; 0, -1, 0];
res = my_conv2(A, hy);
figure('Name','A*hy');
imshow(res, []);

% Calculate A*hx*hy
hxhy = my_conv2(hx, hy);
res = my_conv2(A, hxhy);
figure('Name','A*hx*hy');
imshow(res, []);

% Calculate hxy analytically then calculate A*hxy
hxy = [1, 0, -1; 0, 0, 0; -1, 0, 1];
res = my_conv2(A, hxy);
figure('Name','A*hxy');
imshow(res, []);

% --- QUESTION 3 ---
hx = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
gx = my_conv2(A, hx);

hy = [-1, -2, -1; 0, 0, 0; 1, 2, 1];
gy = my_conv2(A, hy);

g = sqrt(gx.^2 + gy.^2);
figure('Name','G');
imshow(g, []);

function result = my_conv2(A, B)
    % Assume that A is an axb array
    row_a = size(A,1);
    column_a = size(A,2);
    
    % Assume that B is an cxd array
    row_b = size(B,1);
    column_b = size(B,2);
    
    % Rotate B 180 degrees
    B = rot180(B);
    
    % Add zero paddings to the original matrix
    A = add_zero_padding(A, row_b-1, column_b-1);
    
    % Create and fill the result matrix
    result = zeros(row_a , column_a);
    for i = 1:row_a
        for j = 1:column_a
            filtered_area_of_A = A(i:i+row_b-1, j:j+column_b-1);
            result(i, j) = result(i, j) + sum(filtered_area_of_A.*B, 'all');
        end
    end
end

function result = add_zero_padding(A, row_padding, column_padding)
    % Calculate new dimensions and create the matrix
    new_row_size = size(A,1) + row_padding;
    new_column_size = size(A,2) + column_padding;
    result = zeros(new_row_size, new_column_size);
    % Calculate left and top paddings
    top_padding = floor(row_padding/2);
    left_padding = floor(column_padding/2);
    % Copy the matrix into its field in new one
    result(top_padding + 1:top_padding + size(A,1), left_padding + 1:left_padding + size(A,2)) = A;
end

% Improves readability
function result = rot180(A)
    result = rot90(A,2);
end