classdef utils
   methods
        function [out] = quantize_frame(obj,inp,stepsize)
            %scalar uniform quantizer
            out = stepsize*round(inp/stepsize);
        end
        function [frame_dct] = dct_frame_process(obj,frame)
            % Get dimensions of 1 frame
            row = size(frame,1) ;
            col = size(frame,2) ;
            % Divide the cell/frame into matrices of size 8x8
            frame_cell = mat2cell(frame,8*ones(row/8,1),8*ones(col/8,1)) ;
            % Create new array for DCT coefficients 
            frame_dct_cell = cell(size(frame_cell)) ;
            % Iterate over every created block and save coefficient
            for i = 1:row/8 
                for j  = 1:col/8
                    frame_dct_cell{i,j} = dct2(frame_cell{i,j}) ;
                end
            end
            % Create a final matrix
            frame_dct = cell2mat(frame_dct_cell) ;
        end
        function [frame] = idct_frame_process(obj,frame_dct_cell)
            row = size(frame_dct_cell,1) ;
            col = size(frame_dct_cell,2) ;

            frame_dct_cell = mat2cell(frame_dct_cell,8*ones(row/8,1),8*ones(col/8,1)) ;
            frame_cell = cell(size(frame_dct_cell)) ;

            for i = 1:row/8 
                for j  = 1:col/8
                    frame_cell{i,j} = idct2(frame_dct_cell{i,j}) ;
                end
            end

            frame = cell2mat(frame_cell) ;
        end
        function R_block = br_per_block(obj,block,block_VLC_stat)
            vec = block(:);
            R_block = 0;
            for k = 1:length(vec)
                coeffStats = block_VLC_stat{k};
                coeff = vec(k);
                [symIdx,~] = find(coeffStats(:,1)==coeff);
                if isempty(symIdx)
                    symLength = 8; %bit
                    disp('Error')
                else
                    symLength = coeffStats(symIdx,3);
                end
                R_block = R_block + symLength;
            end
        end 
        function [psnr] = get_vid_PSNR(obj, ref_vid, rec_vid)
            psnr_arr = zeros(size(ref_vid,1),1) ;
            for i = 1:size(ref_vid,1)
                psnr_arr(i) = get_frame_PSNR(ref_vid{i},rec_vid{i});
            end
            psnr = mean(psnr_arr) ;
        end
   end
end

function [psnr] = get_frame_PSNR(ref_frame, rec_frame)
    dist = mean(((ref_frame(:)-rec_frame(:)).^2)) ;
    psnr = 10*log10((255^2)/dist) ;
end