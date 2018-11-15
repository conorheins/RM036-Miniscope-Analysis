function [imData,sizx,sizy]=sequence_bigread(directory,sframe,num2read)
% reads sequence of tiff files into matlab array, allows reading from sframe to sframe+num2read-1 frames of the tiff 
% assumes individual tifs are labeled in increasing order of capture time (e.g. has
% the form '...X.tif', where X is the frame number.
% Assumes uncompressed, non-negative (i.e. unsigned) data.
% imData will be your [M,N,frames] array.
% Conor Heins 2017

imgs = dir([directory, filesep,'*.tif']);
frames = zeros(1,length(imgs));
for file = 1:length(imgs)
    fnam = imgs(file).name;
    IDs = regexp(fnam,'\d*','match'); % use a regular expression to find frame index within each file name
    frames(file) = str2num(IDs{1});
end

% sort frames in ascending order before reading
[~,srt] = sort(frames,'ascend');
imgs = imgs(srt);

%get image info from first image of sequence
path_to_file= fullfile(directory,imgs(1).name);

info = imfinfo(path_to_file);
sizx = info.Height;
sizy = info.Width;
        
num_tot_frames=length(imgs);

if (num2read+sframe<= num_tot_frames+1)
    lastframe=num2read+sframe-1;
else
    num2read=numFrames-sframe+1;
    lastframe=num_tot_frames;
    display('Hmmm...just reading from starting frame until the end');
end
        
imData = zeros(sizx,sizx,num2read);
  
sframemsg = ['Reading from frame ',num2str(sframe),' to frame ',num2str(lastframe),' of ',num2str(num_tot_frames), ' total frames'];
disp(sframemsg)

for t = 1:num2read
    imData(:,:,t)=imread(fullfile(directory,imgs(t+sframe-1).name));
end
       
display('Finished reading images')