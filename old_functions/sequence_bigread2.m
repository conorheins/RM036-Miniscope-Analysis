function [imData]=sequence_bigread2(directory, images,sframe,num2read)
%reads sequence of tiff files in Matlab bigger than 4GB, allows reading from sframe to sframe+num2read-1 frames of the tiff 
%Conor Heins 2017
%Assumes uncompressed, non-negative (i.e. unsigned) data.
% imData will be your [M,N,frames] array.

%get image info from first image of sequence
path_to_file= fullfile(directory,images(1).name);

info = imfinfo(path_to_file);
x = info.Height;
y = info.Width;
        
num_tot_frames=length(images);

if (num2read+sframe<= num_tot_frames+1)
    lastframe=num2read+sframe-1;
else
    num2read=numFrames-sframe+1;
    lastframe=num_tot_frames;
    display('Hmmm...just reading from starting frame until the end');
end
        
imData = zeros(x,y,num2read,'uint16');
  
sframemsg = ['Reading from frame ',num2str(sframe),' to frame ',num2str(lastframe),' of ',num2str(num_tot_frames), ' total frames'];
disp(sframemsg)

for t = 1:num2read
    imData(:,:,t)=imread(fullfile(directory,images(t+sframe-1).name));
end
       
display('Finished reading images')