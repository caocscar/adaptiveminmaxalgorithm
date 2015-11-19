function adaptminmax040210
% This program pre-processes the Raman data (single & mapping files) by denoising the spectra,
% performing background subtraction and normalizes the data
% It also pauses every ~100 figures to allow for examination
% Weights the overfitted points so that overfitting is minimized or eliminated
% altogether
% Calculates residual of median filter, wavelet smoothing to original data
warning off MATLAB:polyfit:RepeatedPointsOrRescale
close all

graph = 1;

% Open & read txt file with names of files to work with
filelist = 'filelist.txt';
filenames = textread(filelist,'%s');
nfiles = length(filenames);
disp([filelist, ' has ', num2str(nfiles), ' files.']);
figcount = 0;

% Processing txt files
for i = 1:nfiles

    % checking what kind of file (single or mapping) it is
    filename = char(filenames(i));
    fprintf('%s: single file\n', filename);
    [x,y] = textread(filename,'%f%f','headerlines',1);     % Reads in txt file
    data = [x y];
    FS = process(data, graph, filename);
    figcount = figcount + 1;
    % pausing program to look at figures
    if graph == 1 && figcount > 50
        fprintf('PAUSED after displaying %d figures\n\n', figcount);
        pause
        figcount = 0;
    end

end

end

%------------------------------------------------------
% Pre-processing data file
function FS = process(data, graph, filename)

data = sortrows(data,1);            % sort rows in ascending order
rawdata = data;
medfiltdata(:,2) = medfilt1(data(:,2),5);  % applies median filter to data
data(:,2) = wden(medfiltdata(:,2), 'sqtwolog', 's', 'sln', 5, 'sym8'); %denoises spectra

residual = data(:,2) - rawdata(:,2);
meanres = mean(residual);
stdres = std(residual);
cutoff = 3;
badindex = find(abs(residual) > cutoff*stdres);
goodindex = find(abs(residual) <= cutoff*stdres);
figure(101);
plot(data(:,1),residual(:),'b',data(badindex,1),residual(badindex),'ro'); hold on
plot(data(:,1),-cutoff*stdres,'k',data(:,1),cutoff*stdres,'k')
axis tight
ylabel('Residual')
title(['Cutoff of ',num2str(cutoff),' standard deviations']);
close

[yfit_final, chosen_order, FS] = minmax(data);
unique(chosen_order);
order = FS(2);

ysub_final = data(:,2) - yfit_final;
ysub_final(ysub_final < 0) = 0;
y_normalized = normalize(ysub_final);
BGS_data = [data(:,1) y_normalized];
datastore(BGS_data, filename);                % saves to txt file

[yfit_final_weighted, chosen_order, FS] = minmax_weighted(data);
ysub_final_weighted = data(:,2) - yfit_final_weighted;
ysub_final_weighted(ysub_final_weighted < 0) = 0;
y_normalized_weighted = normalize(ysub_final_weighted);
BGS_data_weighted = [data(:,1) y_normalized_weighted];

if graph

    figure; handle = gcf;
    set(handle,'Position',[100 80 1100 650])       % maximizes window

    h(1) = subplot(211);
    plot(data(:,1),data(:,2),'linewidth',2.5); hold on      % plots denoised spectra
    title(filename);
    plot(data(:,1),yfit_final,'k','linewidth',2);
    plot(data(:,1),yfit_final_weighted,'m','linewidth',1.5);    
    axis tight;

    h(2) = subplot(212);
    plot(data(:,1),y_normalized,'k','linewidth',2); hold on
    plot(data(:,1),y_normalized_weighted,'m','linewidth',1.5);    
    title('Subtracted Normalized Spectrum');
    axis tight;
    
    linkaxes(h,'x');

end

end


%------------------------------------------------------
% Minmax process
function [yfit_final, chosen_order, FS] = minmax_weighted(data)

initial_order = 4;
initial_constraint = 0;
[ysub_init yfit_init] = backsub_weighted(data,initial_order,initial_constraint);        % modified polynomial fit
FSratio = (max(yfit_init)-min(yfit_init))/(max(ysub_init)-min(ysub_init));
order = findorder(FSratio);
FS = [FSratio order];
YFITS = zeros(size(yfit_init,1),4);
delta = 1;
col = 0;
for constraint = 0:1
    for z = order:delta:order+delta
        col = col + 1;
        if z == initial_order && constraint == initial_constraint
            YFITS(:,col) = yfit_init;
        else
            [ysub1 yfit1] = backsub_weighted(data,z,constraint);
            YFITS(:,col) = yfit1;
        end
    end
end
[yfit_final, chosen_fit] = max(YFITS,[],2);
chosen_order = order + delta - delta*rem(chosen_fit,size(YFITS,2)/2);

end
%------------------------------------------------------
% Background subtraction with weights for overfitted regions
function [ysub,yfit] = backsub_weighted(data,n,constraint)

% compute n-th order polynomial 
x_weighted = data(:,1);
y_weighted = data(:,2);
count = 0;
maxcount = 200;
% pts4convergence = 10;
percentage_increased = 0.03;
index_overfittedpts_old = 0;
while(1)
    count = count + 1;
    P = polyfit(x_weighted,y_weighted,n);  % shape-preserving interpolant
    yfit = polyval(P,data(:,1));
    ymin = min(data(:,2),yfit);
    if constraint == 1
        endpts = 1;
        ymin(1:endpts) = data(1:endpts,2);
        ymin(end-endpts+1:end) = data(end-endpts+1:end,2);
    end
    ysub = data(:,2) - yfit;                                % original - bestfit curve
    index_overfittedpts = find(ysub < 0);
    if ~isempty(index_overfittedpts) && count <= 200
        if isequal(index_overfittedpts,index_overfittedpts_old)
            if percentage_increased < 0.5
                percentage_increased = 2*percentage_increased;
            else
                [count -1]
                break;
            end
        else
            percentage_increased = 0.03;
        end
        data_overfittedpts = data(index_overfittedpts,:);
        if count < 50
            percentage = max(0.01,percentage_increased);
        elseif count < 100
            percentage = max(0.1,percentage_increased);
        elseif count < 150
            percentage = max(0.2,percentage_increased);
        else
            percentage = max(0.4,percentage_increased);
        end
        weightsize = [round(percentage*size(data,1)) 1];
        x_weighted = [data(:,1); repmat(data_overfittedpts(:,1),weightsize)];
        y_weighted = [ymin; repmat(data_overfittedpts(:,2),weightsize)];
        index_overfittedpts_old = index_overfittedpts;       
    else
        count
        break;
    end
end
if count >= maxcount; disp('MAX number of iterations reached'); end

end

%------------------------------------------------------
% Minmax process
function [yfit_final, chosen_order, FS] = minmax(data)

initial_order = 4;
initial_constraint = 0;
[ysub_init yfit_init] = backsub(data,initial_order,initial_constraint);        % modified polynomial fit
FSratio = (max(yfit_init)-min(yfit_init))/(max(ysub_init)-min(ysub_init));
order = findorder(FSratio);
FS = [FSratio order];
YFITS = zeros(size(yfit_init,1),4);
delta = 1;
col = 0;
for constraint = 0:1
    for z = order:delta:order+delta
        col = col + 1;
        if z == initial_order && constraint == initial_constraint
            YFITS(:,col) = yfit_init;
        else
            [ysub1 yfit1] = backsub(data,z,constraint);
            YFITS(:,col) = yfit1;
        end
    end
end
[yfit_final, chosen_fit] = max(YFITS,[],2);
chosen_order = order + delta - delta*rem(chosen_fit,size(YFITS,2)/2);

end
%------------------------------------------------------
% Background subtraction
function [ysub,yfit] = backsub(data,z,constraint)

% compute z-th order polynomial 
baseline = data(:,2);
count = 0;
maxcount = 200;
pts4convergence = 10;
while(1)
    P = polyfit(data(:,1),baseline,z);  % shape-preserving interpolant
    yfit = polyval(P,data(:,1));
    count = count + 1;
    baseline = min(data(:,2),yfit);
    if constraint == 1
        endpts = 1;
        baseline(1:endpts) = data(1:endpts,2);
        baseline(end-endpts+1:end) = data(end-endpts+1:end,2);
    end
    ysub = data(:,2) - yfit;                                % original - bestfit curve       
    R(count) = length(find(ysub < 0));
    if count > pts4convergence
        S = R(end-pts4convergence+1:end);
        if (length(unique(S))==1 || count >= maxcount);
            count
            break; end
    end
end
if count >= maxcount; disp('MAX number of iterations reached'); end

end

%------------------------------------------------------
% Choose order based on fluorescence-to-signal ratio
function order = findorder(FSratio)

if FSratio < 0.2
    order = 1;
elseif FSratio < 0.75
    order = 2;
elseif FSratio < 8.5
    order = 3;
elseif FSratio < 55 
    order = 4;
elseif FSratio < 240
    order = 5;
else
    order = 6;
end

end
%------------------------------------------------------
% Normalize data from [0,1]
function y = normalize(x)

z = x - min(x);
y = z/max(z);

end


%------------------------------------------------------
% Writes data to txt file
function datastore(data,filename)

if size(data,2) == 2
    Nfile = strcat('AMM-',filename);               
    fout = fopen(char(Nfile),'w');
    fprintf(fout,'\t%s\n','x           y_normal');
else
    return
end
for i = 1:size(data,1)
    for j = 1:size(data,2)
        fprintf(fout,'\t%f',data(i,j));
    end
    fprintf(fout,'\n');
end
fclose(fout);

end
