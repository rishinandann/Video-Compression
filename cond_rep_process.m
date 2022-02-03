function [rec_vid,bits_tot] = cond_rep_process(my_utils, vid,step,block_VLC_stat,c)
    % Extrabits indicate mode for each block in one frame
    extraBits_frame = 144*176/(16*16) ;
    rec_vid = cell(size(vid));
    
    % Take first frame of the original video
    first_frame = vid{1};
    
    % Intra code the first frame
    first_frame_dct_coeffs = my_utils.quantize_frame(my_utils.dct_frame_process(first_frame),step);
    row = size(first_frame,1) ;
    col = size(first_frame,2) ;
    bits_tot = 0;
    
    % get the length of every coefficient in the first frame
    for i_1=1:row/16
        for j_1 = 1:col/16
            % Get the corresponding block in the frame
            block = first_frame_dct_coeffs((i_1-1)*16+1:i_1*16,(j_1-1)*16+1:j_1*16);
            % Get the bitrate/block of this specific Block
            R_block = my_utils.br_per_block(block,block_VLC_stat);
            bits_tot = bits_tot + R_block;
        end
    end
    
    rec_vid{1} = my_utils.idct_frame_process(first_frame_dct_coeffs);
    bits_tot = bits_tot + extraBits_frame;
    
    % Frames after first frame
    for i = 2:size(vid,1)
        [rec_vid{i} ,R_frame] = rep_frame_process(my_utils,vid{i},rec_vid{i-1},step,block_VLC_stat,c) ;
        bits_tot = bits_tot+R_frame+extraBits_frame ;
    end
end

function [rep_frame, R_frame] = rep_frame_process(my_utils,actual_frame,prev_frame,step,block_VLC_stat,c)
    
    row = size(actual_frame,1);
    col = size(actual_frame,2);
    
    % Calculate the Intra frame 
    intra_frame_dct_coeffs = my_utils.quantize_frame(my_utils.dct_frame_process(actual_frame),step);
    intra_frame = my_utils.idct_frame_process(intra_frame_dct_coeffs);
    copy_frame = prev_frame;
    
    % Calculate squared errors for whole image
    intra_sq_err = (actual_frame-intra_frame).^2 ;
    copy_sq_err = (actual_frame-copy_frame).^2 ;
    
    ind_frame = zeros(row,col) ;
    lambda = c*step^2;
    R_frame = 0;
    % Mode Bits
    mode_bit = 1;
    for i=1:row/16
        for j = 1:col/16
            % Get intra blocks DCT
            intra_block_dct_coeffs = intra_frame_dct_coeffs((i-1)*16+1:i*16,(j-1)*16+1:j*16);
            % D_intra
            intra_err = intra_sq_err((i-1)*16+1:i*16,(j-1)*16+1:j*16);
            D_intra = mean(intra_err(:));
            %D_copy
            copy_err = copy_sq_err((i-1)*16+1:i*16,(j-1)*16+1:j*16); % per pixel
            D_copy = mean(copy_err(:));
            
            %R_intra
            R_block = my_utils.br_per_block(intra_block_dct_coeffs,block_VLC_stat);
            R_intra = (R_block+mode_bit)/256;
            % R_copy
            R_copy = mode_bit/256; % per pixel
            
            % Calculate Laplacian
            j_intra = D_intra + lambda*R_intra ;
            j_copy = D_copy + lambda*R_copy ;
            if(j_intra<j_copy)
                ind_frame((i-1)*16+1:i*16,(j-1)*16+1:j*16) = ones(16,16) ;
                R_frame = R_frame+R_intra*256;
            else
                R_frame = R_frame+R_copy*256;
            end
        end
    end
    rep_frame = intra_frame.*ind_frame +prev_frame.*(~ind_frame) ;
end


