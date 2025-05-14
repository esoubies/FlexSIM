function dark_channel = get_dark_channel(image, win_size)

[m, n, ~] = size(image);

pad_size = floor(win_size/2);

%padded_image = padarray(image, [pad_size pad_size], Inf);

dark_channel = zeros(m, n); 

% for j = 1 : m
%     for i = 1 : n
%         patch = padded_image(j : j + (win_size-1), i : i + (win_size-1), :);
%         dark_channel(j,i) = min(patch(:));
%      end
% end
% 
% parfor k =  1:m*n       
%         i = mod(k+m-1,m)+1;
%         j = floor((k+m-1)/m);
%         patch = padded_image(j : j + (win_size-1), i : i + (win_size-1), :);
%         dark_channel_temp(k) = min(patch(:));
% end

dark_channel_temp = movmin(movmin(image,win_size,1),win_size,2);
dark_channel_temp = dark_channel_temp';

for k = 1:m*n
    dark_channel(floor((k+m-1)/m),mod(k+m-1,m)+1) = dark_channel_temp(k);
end


end