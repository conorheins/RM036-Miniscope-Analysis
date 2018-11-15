%% script to pre-process data (get spectral estimate of noise and set signals below (sig * noise) to 0, save as Ythresh in data object

tic

Yfilt_temp = M1.Yfilt;
Yfilt_temp = reshape(Yfilt_temp,d1*d2,numFrames);
Yfilt_temp = bsxfun(@minus,Yfilt_temp,median(Yfilt_temp,2));
Ynoise = get_noise_fft(Yfilt_temp,options_preprocess);
Yfilt_temp(bsxfun(@lt,Yfilt_temp,Ynoise*sig)) = 0;
data.Ythresh = Yfilt_temp;

clear Yfilt_temp Ynoise

fprintf('Time taken to denoise & threshold data: %.2f minutes\n',(toc/60))