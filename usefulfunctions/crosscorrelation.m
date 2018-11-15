%
%crosscorrelation.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%                       Network Identification                        %
%                by  Erik Smedler MSc, Phd student                    %
%                                                                     %
%                        Version 3.0 Publication                      %
%                                                                     %
%           Copyright 2014 by Erik Smedler All Rights Reserved        %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%  Program file: NetworkIdentification.m                              %
%  Related needed files: crosscorrelation.m; pickcells.m;             %
%  Required file: oscillationdata.dat                                 %  
%  Creates file: correlationdata.mat                                  %
%  1st version: 2008-12-12, Karolinska Institutet, Stockholm, Sweden  %
%  2nd version: 2011-09-24, Karolinska Institutet, Stockholm, Sweden  %
%  3nd version: 2014-01-23, Karolinska Institutet, Stockholm, Sweden  %
%                                                                     %
%                                                                     %
%  Contact: erik.smedler@ki.se                                        %
%                                                                     %
%  Last updated: 2014-01-23                                           %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - Main file for the NetworkIdentification program - - - %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%--------------------------------------------------------------------------
%
%Maximum correlation coefficient
%Function that calculates the correlation coefficient between two signals
%by using a predefined lag in the cross covariance function.
%
%--------------------------------------------------------------------------
%Input:
%t: the time points of the signals.
%y1: the fluorescence intensities for cell i at every time point.
%y2: the fluorescence intensities for cell j at every time point.
%maxlags: the maximum lag used by xcov.m.
%i1: index for the first cell.
%j1: index for the second cell.
%
%Output:
%z: maximum correlation coefficient between two signals and the lag z=[Rmax lag].
%
%--------------------------------------------------------------------------
function z=crosscorrelation(t,y1,y2,maxlags,i1,j1)
Y2=y2;
y2old=y2;
y1old=y1;
told=t;
Maxlag=0;
if maxlags~=0
    if (maxlags/length(y1))>=0.25 
        maxlags=floor(0.25*length(y1));
    end
    %Cross covariance between the two signals.
    [Ry1y2,lags]=xcov(y1,y2,maxlags,'coeff');

    %Finds for which lag Maxlag the maximum correlation coefficient is.
    N=length(Ry1y2);
    maxlag=0;
    minlag=0;
    maxlagindex=1;
    minlagindex=1;
    for i=1:N
        MAX=max(Ry1y2);
        MIN=min(Ry1y2);
        if Ry1y2(i)==MAX
            maxlag=lags(i);
            maxlagindex=i;
        end
        if Ry1y2(i)==MIN
            minlag=lags(i);
            minlagindex=i;
        end
    end

    %Finds the lag with the maximum absolute value.
    if abs(Ry1y2(maxlagindex))>=abs(Ry1y2(minlagindex))
        Maxlag=maxlag;
    else
        Maxlag=minlag;
    end
    %Moves the signal y2 in order to get the maximum correlation with y1.
    if Maxlag~=0
        N1=length(y2);
        Y2=pi*ones(1,N1);
        for i=1:N1
            if Maxlag>0 %If the lag is to the right.
                if i-Maxlag>0
                    Y2(i)=y2(i-Maxlag);
                end
            end
            if Maxlag<0 %If the lag is to the left.
                if i-Maxlag<=N1
                    Y2(i)=y2(i-Maxlag);
                end
            end
        end
    %Removes those part of the signals that do not overlap.
    %Finds the concerned indices.
        indexvector=[];
        for i=1:N1
            if Y2(i)==pi   
                indexvector=[indexvector i];
            end
        end
        removeindex=indexvector;
        Y2(removeindex)=[];
        y1(removeindex)=[];
        t(removeindex)=[];
    else
        Y2=y2;
    end
end

%The maximum correlation coefficient.
M1=corrcoef(y1,Y2);
maxcorrcoef=M1(1,2);

if isnan(maxcorrcoef)
    disp(['Warning! correlation coefficient is NaN for cells ',num2str(i1),' and ',num2str(j1),'!'])
end

%Output from the function file.
z=[maxcorrcoef; Maxlag];