function [] = Produce_fig_from3dbin(path_in)

% Produce_fig_from3dbin('/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4')

% clearvars;

path_in = '/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4';
path_out = '/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs';

%vars
NSIDE = 4;
filename_in = '_1_2dproj_z3_data_slAll';
filename_out = '2dproj_z3_ts32_sl';
slices = 32;
nc = 512;

specs_path_list_in =dir(strcat(path_in,'/*'));

specs_path_list_nowake = strcat(path_in,'/',specs_path_list_in(3).name);
sample_list_nowake=dir(strcat(specs_path_list_nowake,'/sample5*'));
sample_list_nowake={sample_list_nowake.name};

specs_path_list_wake = strcat(path_in,'/',specs_path_list_in(4).name);
sample_list_wake=dir(strcat(specs_path_list_wake,'/sample5*'));
sample_list_wake_out={sample_list_wake.name};
sample_list_wake_in=strcat(sample_list_wake_out,'/half_lin_cutoff_half_tot_pert_nvpw_v0p6');

sample_list_nowake_out_path=strcat(path_out,'/',specs_path_list_in(3).name);
sample_list_wake_out_path=strcat(path_out,'/',specs_path_list_in(4).name);

for w_nw=1:2
% for w_nw=1   
    if w_nw==1
        specs_path_list=specs_path_list_nowake;
        sample_list_in=sample_list_nowake;
        sample_list_out=sample_list_nowake;
        sample_list_out_path = sample_list_nowake_out_path;
    end    
    if w_nw==2
        specs_path_list=specs_path_list_wake;
        sample_list_in=sample_list_wake_in;
        sample_list_out=sample_list_wake_out;
        sample_list_out_path = sample_list_wake_out_path;
    end
    mkdir(sample_list_out_path);
        
    for sample = 1:length(sample_list_in)
%     for sample = 1
        list_of_angle_paths=dir(char(strcat(specs_path_list,'/',string(sample_list_in(sample)),'/data','/1lf_0.5rf','/NSIDE_',num2str(NSIDE),'/anglid_*')));       
        list_of_angle_paths={list_of_angle_paths.name};
%         mkdir(strcat(sample_list_out,'/',string(sample_list(sample))));

        for angle_id=1:length(list_of_angle_paths)
%       for angle_id=1
            path_prefix_name = strcat(sample_list_out_path,'/',string(sample_list_out(sample)),'-',list_of_angle_paths(angle_id));
%             mkdir(prefix_name);

            path1=strcat(specs_path_list,'/',string(sample_list_in(sample)),'/data','/1lf_0.5rf','/NSIDE_',num2str(NSIDE),'/',list_of_angle_paths(angle_id),'/');               
            path2=dir(char(strcat(path1,"*pv*")));
            path2={path2.name};
            path2=path2(1);
            file_in = strcat(specs_path_list,'/',string(sample_list_in(sample)),'/data','/1lf_0.5rf','/NSIDE_',num2str(NSIDE),'/',list_of_angle_paths(angle_id),'/',string(path2),'/2dproj/dm/',filename_in,'.bin');
            fid = fopen(file_in);

            for slice_id=1:slices
%             for slice_id=1
                file_out = strcat(path_prefix_name,'-',filename_out,num2str(slice_id),'.png');
                skip = 4*(nc*nc)*(slice_id-1);     %each number has 4 `bytes to skip
                fseek(fid,skip,'bof');
                slice_2d = fread(fid,[512 512], 'float32','l') ; 
                slice_2d=(slice_2d-mean(slice_2d(:)))/mean(slice_2d(:));
                slice_2d=(atan((slice_2d+1)*16));
                % do figure
                fig=figure('Visible', 'off');
%                 fig=figure;
                set(gcf, 'Position', [0 0 nc-1 nc-1]);
                hold on;
                axes('Position',[0 0 1 1],'Visible','off');
                imagesc([0 nc-1],[0 nc-1],slice_2d);
                saveas(fig,file_out);
                close(fig);
            end
            fclose(fid);
        end
    end
end






end

