function [] = CopyRelevantFigure_eachslice_all_groups()

% this function should be run after obtaining the signal of each figure. It
% will generate the fist of relevant figures and copy then to a diferent
% directory

% info_ = info(:)';

% slice_names = repmat([1:3]',1,5);
% slice_names = repmat([5001:5100]',1,3072);

% "info" matriz (96x32)

%  Z = X+Y'    %same as plus(X,Y')

threshold = 30;


simulations = "sample"+string([5001:5100]);

angles = "-anglid_"+string([1:96]);
slices = "-2dproj_z3_ts32_sl"+string([1:32]);

info = angles'+slices;
info = info(:)';

slice_names = simulations'+info;


select_wake = signal_sample_w>threshold;

%list of figure with wakes higher than threshould
selected_slice_names_wake = slice_names(select_wake(:));  




sorted_signal_sample_nw = sort(signal_sample_nw(:));


thresh_no_wake = sorted_signal_sample_nw(end-sum(select_wake(:))-1);
select_no_wake = signal_sample_nw>thresh_no_wake;
%list of figure without wakes with higher threshould, maching the number of
%figures selected to the wake case (may diuffer by a few, since there could be more than one figure with a given threshold)
selected_slice_names_nowake = slice_names(select_no_wake(:));  

size(selected_slice_names_nowake(:))
size(selected_slice_names_wake(:))


% dlmwrite('/home/asus/Dropbox/extras/storage/graham/ht/labels_eachSlice_nw_thr50.txt',selected_slice_names_nowake,'delimiter','\t');
% dlmwrite('/home/asus/Dropbox/extras/storage/graham/ht/labels_eachSlice_w_thr50.txt',selected_slice_names_wake,'delimiter','\t');


writematrix(selected_slice_names_nowake,strcat('/home/asus/Dropbox/extras/storage/graham/ht/labels_eachSlice_nw_thr',num2str(threshold),'.txt'));
writematrix(selected_slice_names_wake,strcat('/home/asus/Dropbox/extras/storage/graham/ht/labels_eachSlice_w_thr',num2str(threshold),'.txt'));


% in the second part we will move the relevant figures to a folder



end