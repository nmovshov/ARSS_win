%% Make anaglyph
clear all
close all
clc

%% Get and sort file names
rfiles=dir('*.R.jpg');
lfiles=dir('*.L.jpg');
numbers=nan(size(rfiles));
for k=1:length(rfiles)
    ix=strfind(rfiles(k).name,'.');
    numbers(k)=str2double(rfiles(k).name(ix(1)+1:ix(2)-1));
end
[~,IDX]=sort(numbers);
rnames={rfiles(IDX).name}';
lnames={lfiles(IDX).name}';

%% Read the frames, remove red from right frame and cyan from left frame and combine
for k=1:length(lnames)
    lframe=imread(lnames{k});
    rframe=imread(rnames{k});
    lframe(:,:,2:3)=0;
    rframe(:,:,1)=0;
    frame=lframe+rframe;
    fname=['frame' num2str(k) '.jpg'];
    imwrite(frame,fname);
end