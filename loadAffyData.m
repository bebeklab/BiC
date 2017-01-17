% This method loads mRNA data associated with Affy chips
% with ALL BLANKS IN COLUMN 2 OF IDFILE REPLACED WITH '-' dashes
% It median centers each array

function mRNA = loadAffyData(gsefile, gplfile )

%%
% clear
% gsefile=[pwd '/Data/Apc RMA.txt'];
% gplfile=[pwd '/Data/Apc p21 - ids.txt']
%%
raw_data = importdata(gsefile);
raw_data.textdata = raw_data.textdata(2:end,1);

fid = fopen(gplfile);
fgetl(fid);
ids = textscan(fid, '%s %s %*[^\n]', 'delimiter','\t');
fclose(fid);

% % Rescale mRNA values to be median centered
% for i=1:size(raw_data.data,2)
%     meds(i) = median(raw_data.data(~isnan(raw_data.data(:,i)),i));
%     raw_data.data(:,i) = raw_data.data(:,i) - meds(i);
% end

% Map gene names from array id file onto HUGO/MGI names
for i = 1:length(ids{1})
 
    wholestring = ids{2}{i};
    bars = strfind( wholestring,' /// '  );
    if isempty(bars)
        ids{3}{i} = wholestring;
        ids{4}{i} = wholestring;
    else
        ids{3}{i} = upper(wholestring(1:bars(1)-1));
        ids{4}{i} = upper(wholestring(bars(1)+5:end));
    end
end
%%

mRNA.data = raw_data.data; 
mRNA.probes = ids{1};
mRNA.genest = ids{3};
mRNA.genes = mRNA.genest';
mRNA.syn = ids{4}';

