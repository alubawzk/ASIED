function []  = LCS_ellipse_multi()
%% parameters illustration
%1) Tac: 
%The threshold of elliptic angular coverage which ranges from 0~360. 
%The higher Tac, the more complete the detected ellipse should be.
%2) Tr:-0.109400
%The ratio of support inliers to ellipse which ranges from 0~1.
%The higher Tr, the more sufficient the support inliers are.
%3) specified_polarity: 
%1 means detecting the ellipses with positive polarity;
%-1 means detecting the ellipses with negative polarity; 
%0 means detecting all ellipses from image
warning('off');
close all;
clc; clear;
%% parameters
tau_p = 0.095;
tau_h = 0.95;

%% 批量处理
dataset_name = [{'Prasad Images - Dataset Prasad'},...%1
                {'Random Images - Dataset #1'},...%2
                {'Smartphone Images - Dataset #2'},...%3
                {'Satellite Images - Dataset Meng #1'},...%4
                {'Satellite Images - Dataset Meng #2'},...%5
                {'Prasad_plus'}, ...%6
                {'Tableware'},...%7
                {'Synthetic Images - Occluded Ellipses'},...%8
                {'Synthetic Images - Overlap Ellipses'},... %9
                {'Concentric Ellipses - Dataset Synthetic'},...%10
                {'Concurrent Ellipses - Dataset Synthetic'},...%11
                {'Industrial PCB Image Dataset'} ...%12
                ];
            
dataset_type = [{'.jpg'},{'.jpg'},{'.bmp'},{'.jpg'},{'.jpg'},{'.jpg'}, ...
                {'.jpg'},{'.jpg'},{'.jpg'},{'.png'},{'.png'},{'.bmp'}];

% dataset_num = 10;
for data_n = 1  % [1,2,3,4,5,6,7,8,9,10,11]
    Dataset = dataset_name{data_n};
    imgType = dataset_type{data_n};
    % imgPath=['F:\Program\Datasets\', Dataset, '\images\'];
    imgPath=['D:\Program\Elp\Datasets\', Dataset, '\images\'];
    imgData_root = [imgPath,'*',imgType];
    imgDir = dir(imgData_root);
    img_num = 0;
    time_total = 0;
    group_total = 0; clust_total = 0; verfi_total = 0;
    for Mahala_tolerance=0.076 %0.076
        for i=1:length(imgDir)
            I = imread([imgPath imgDir(i).name]);
            scale = 1; %Prasad:1.3
            I = imresize(I, scale);
            img_num = img_num + 1;
            disp([Dataset, ': ', num2str(img_num)]);
%             if img_num==49
%                 beep;
%             end

            % add noise
            % I = imnoise(I,'gaussian',0,0.05);
            
            % detecting ellipses from real-world images
            [ellipses, ~, ~, ~, time_PI, time_combine, time_cluster, time_verifiction] = ellipseDetectionByArcSupportLSs(I, tau_p, tau_h);
            ellipses(:,1:4) = ellipses(:,1:4) / scale;
            time_total = time_total+time_PI;
            group_total = group_total+time_combine; 
            clust_total = clust_total+time_cluster; 
            verfi_total = verfi_total+time_verifiction;
            %% 显性融合
            % ellipses_label = true(1,size(ellipses,1));
            % for k=1:size(ellipses,1)-1
            %     for j=k+1:size(ellipses,1)
            %         if sqrt((ellipses(k,1)-ellipses(j,1))^2+(ellipses(k,2)-ellipses(j,2))^2) < 2 &&...
            %                 (ellipses(k,3)-ellipses(j,3)) < 2 && (ellipses(k,4)-ellipses(j,4)) < 2 &&...
            %                 (ellipses(k,5)-ellipses(j,5)) < 0.2
            %             ellipses_label(1,j)=false;
            %         end
            %     end
            % end
            % ellipses = ellipses(ellipses_label,:);
            %% write the result ellipses txt
            resultname=strrep(imgDir(i).name,imgType,'.txt');
            new_folder = ['D:\Program\Elp\Datasets\',Dataset,'\ASIED_GitHub'];
            if exist(new_folder,'dir')==0
                mkdir(new_folder);% 或者用 mkdir data,在当前目录下，生成一个data文件夹
            end
            fid=fopen([new_folder,'\',resultname],'wt');
            for m=1:size(ellipses,1)+1
                if m==1
                    fprintf(fid,'%d\n',size(ellipses,1));
                    continue
                end
                for n=1:5
                    fprintf(fid,'%f',ellipses(m-1,n));
                    fprintf(fid,'\t');
                end
                fprintf(fid,'\n');
            end
            fclose(fid);
        end
    end
    disp('-----------------------------------------------------------');
    disp(['average total time:',num2str(time_total/img_num),'s']);
    disp(['average group time:',num2str(group_total/img_num),'s']);
    disp(['average clust time:',num2str(clust_total/img_num),'s']);
    disp(['average verfi time:',num2str(verfi_total/img_num),'s']);

end

    
% for Mahala_tolerance=0.076
%     for i=1:length(imgDir)
%         I = imread([imgPath imgDir(i).name]);
%         img_num = img_num + 1;
%         disp(img_num);
%         % add noise
% %         I = imnoise(I,'gaussian',0,0.05);
%         %% detecting ellipses from real-world images
%         [ellipses, ~, ~, ~, ~, ~, ~, ~] = ellipseDetectionByArcSupportLSs(I, Tac, Tr, specified_polarity, Mahala_tolerance);        
% %        %% 显性融合
% %         ellipses_label = true(1,size(ellipses,1));
% %         for k=1:size(ellipses,1)-1
% %             for j=k+1:size(ellipses,1)
% %                 if sqrt((ellipses(k,1)-ellipses(j,1))^2+(ellipses(k,2)-ellipses(j,2))^2) < 2 &&...
% %                         (ellipses(k,3)-ellipses(j,3)) < 2 && (ellipses(k,4)-ellipses(j,4)) < 2 &&...
% %                         (ellipses(k,5)-ellipses(j,5)) < 0.2
% %                     ellipses_label(1,j)=false;
% %                 end
% %             end
% %         end
% %         ellipses = ellipses(ellipses_label,:);
%         %% write the result ellipses txt
%         resultname=strrep(imgDir(i).name,imgType,'.txt');
% %         new_folder = ['D:\360MoveData\Users\ylg\Desktop\ProgramLab\Datasets\',Dataset,'\HBCED_Ma\',num2str(Mahala_tolerance)];
% %         new_folder = ['D:\360MoveData\Users\ylg\Desktop\ProgramLab\Datasets\',Dataset,'\HBCED_Ma'];
% %         new_folder = ['D:\360MoveData\Users\ylg\Desktop\ProgramLab\Datasets\',Dataset,'\HBCED_Noise','_',num2str(dev)];
%         new_folder = ['D:\360MoveData\Users\ylg\Desktop\ProgramLab\Datasets\',Dataset,'\ASIED'];
%         if exist(new_folder,'dir')==0
%             mkdir(new_folder);% 或者用 mkdir data,在当前目录下，生成一个data文件夹
%         end
%         fid=fopen([new_folder,'\',resultname],'wt');
%         for m=1:size(ellipses,1)+1
%             if m==1
%                 fprintf(fid,'%d\n',size(ellipses,1));
%                 continue
%             end
%             for n=1:5
%                 fprintf(fid,'%f',ellipses(m-1,n));
%                 fprintf(fid,'\t');
%             end
%             fprintf(fid,'\n');
%         end
%         fclose(fid);
%     end
%     
%     
% % end
% 
% end

%%
% beep;
% load gong
% sound(y,Fs)
end