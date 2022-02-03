clc;
clear;
close all;

%utils
my_utils = utils;

%% Choose input video
inp_yuv = 'foreman_qcif.yuv';
%inp_yuv = 'mother-daughter_qcif.yuv';

%% Import YUV file
frames = yuv_import_y(inp_yuv,[176 144],50);

%% Variables
step_sizes = [2^3, 2^4, 2^5, 2^6];
bitrates_arr = zeros(size(step_sizes,1),1) ;
psnr_arr = zeros(size(step_sizes,1),1) ;
recon_frames = cell(size(step_sizes, 1),1);

%block counter
b=1;

%% Load statistics of intra frames from part 1
load('IntraStats.mat')

%% run for all step sizes
for i=1:numel(step_sizes)
    
    %% Uniform Quantizer
    step = step_sizes(i) ;
    
    %% Load statistics for range
    block_VLC_stat = intraStats{b};
       
    %% Rate Calculation
    fps = 30 ;
    lambda = 0.2   ; % Lambda = c*Q^2
    [recon_frames{b},bits_tot] = cond_rep_process(my_utils, frames,step,block_VLC_stat,lambda);
    rate = bits_tot/size(frames,1)*fps/1000 ; %bitrate given 30 fps
    bitrates_arr(b) = rate ;
    
    %% PSNR Calculation
    psnr_arr(b) = my_utils.get_vid_PSNR(frames,recon_frames{b}) ;
    
    b=b+1 ;
end

%% Plotting
figure ;
plot(bitrates_arr,psnr_arr,'-o')
xlabel('Rate[kbit/s]') ;
ylabel('PSNR') ;
rates3 = bitrates_arr;
psnr3 = psnr_arr;
vidReconstructed3 = recon_frames;
figure;
subplot(1,3,1); imshow(uint8(frames{1})); title('Original Frame');
subplot(1,3,2); imshow(uint8(vidReconstructed3{2}{1})); title('Frame with quantization step 16');
subplot(1,3,3); imshow(uint8(vidReconstructed3{4}{1})); title('Frame with quantization step 64');

save('ex3.mat','rates3','psnr3','vidReconstructed3');
