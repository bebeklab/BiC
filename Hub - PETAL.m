% Calculate coexpression linkages from a group (Apcp21net) to another group
% Calculate significance based on connections to randomly chosen stuff

%
% Need files: 'dChip Mouse Intestine Matrix.txt', 'Apc p21 dChip.txt'
% Calls: CoexpNeighborhood, 

%%
%Load the ppi
clear

%load([pwd '/Data/' 'Interaction Structure.mat']);

% mRNA's gonna have some NaNs (b/c dChip was made to leave blanks)
% Make sure the dChip file has blanks replaced with NaN
% mRNA = loadAffyData('dChip Mouse Intestine Matrix.txt', [pwd '/Data/Apc p21 - ids.txt'], ppi);
%mRNA = loadAffyData([pwd '/Data/Apc dChip.txt'], [pwd '/Data/Apc p21 - ids.txt'], ppi);
mRNA_t = loadAffyData([pwd '/Data/Apc RMA.txt'], [pwd '/Data/Apc p21 - ids.txt'] );
mRNA=mRNA_t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load up the proteomic targets and the Apc-p21 network
fid1 = fopen([pwd '/Data/Apc prots.txt']);
ApcProts_raw = textscan(fid1, '%s');
fclose(fid1);

fid2 = fopen([pwd '/Data/petal_nodes_p_0.05.txt']);
NumCols = 37; % Replace blank cells with '-'
i=1;
while ~feof(fid2)
    a_net = textscan(fid2, ['%*s\t %*s\t' repmat('%s\t',1, NumCols-2) '\n']);
    if ~isempty([a_net{:}])
        cut = min(strmatch('-',[a_net{:}]));
        if isempty(cut)
            CANnets{i} = [a_net{:}];
        else
            CANnets{i} = [a_net{1:cut-1}];
        end
        i=i+1;
    end
end
fclose(fid2);
%%
[junk junk ApcProts]= intersect(upper(ApcProts_raw{:}), upper(mRNA_t.genes));

[junk leftovers]    = setdiff(  upper(ApcProts_raw{:}), upper(mRNA_t.genes));

[junk junk hits] = intersect(upper(ApcProts_raw{1}(leftovers)), upper(mRNA_t.syn));

%ApcProts = [ApcProts mRNA_t.syn(hits)];

%[junk junk ApcLocs] = intersect(ApcProts, mRNA.genes);

ApcLocs=ApcProts;


%%
clear NetPPI NetLocs;
for i=1:length(CANnets)
    [junk junk NetLocs{i}] = intersect(upper(CANnets{i}), upper(mRNA_t.genes));

end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For starters, calculate coexpression from (A) each of the CANnets to (B) prot targets
% Crypts
% tscore = (mean(mRNA_t.data(:,5:8)')-mean(mRNA_t.data(:,13:16)'))./sqrt(var(mRNA_t.data(:,5:8)')/4+var(mRNA_t.data(:,13:16)')/4);
% All
tscore = (mean(mRNA_t.data(:,1:8)')-mean(mRNA_t.data(:,9:16)'))./sqrt(var(mRNA_t.data(:,1:8)')/8+var(mRNA_t.data(:,9:16)')/8);
ts = abs(tscore)./max(abs(tscore));
ts_bin = zeros(length(ts),1);
[x y] = sort(ts); cut = round(.75*length(ts));
ts_bin(find(ts>=x(cut))) = 1;
p=1;
%%
for p=1:length(CANnets)
    rho.all{p} = (corr(mRNA.data(NetLocs{p},:)', mRNA.data', 'type','Pearson','rows','pairwise'));
    rho.all{p}(isnan(rho.all{p})) = 0;
    rho.t{p} = bsxfun(@times, rho.all{p},(ts(NetLocs{p})'));
%     rho.t{p} = rho.all{p}(find(ts_bin(NetLocs{p})),:);
    
%     rho.apc{p} = rho.t{p}(:,ApcLocs);

    torque.net(p) = CalculateTorque(rho.t{p}, ApcLocs, 0);
    maxperm = 1000;
    D = zeros(maxperm,1);
    for m=1:maxperm
        r_group= round(1 + (size(rho.t{p},2)-1).*rand(length(ApcLocs),1));
        D(m) = CalculateTorque(rho.t{p}, r_group, 0);
    end
    torque.pval(p) = length(find(D<torque.net(p)))/length(find(D));
    p
end
save([pwd '/Results/torque.mat'],'torque','rho');
%%
fid =fopen([pwd '/Results/bimodality.txt'],'w');
fprintf(fid,'%s\t','Bimodality, B');
fprintf(fid,'%s\t','p-value');
fprintf(fid,'%s\n','Network');
for p=1:length(torque.pval)
    fprintf(fid,'%f\t',torque.net(p));
    fprintf(fid,'%f\t',torque.pval(p));
    for z=1:length(CANnets{p})
        fprintf(fid, '%s\t', CANnets{p}{z});
    end
    fprintf(fid, '\n');
end
fclose(fid);

%%
p=6; ecdf(rho.apc{p}(:)); hold on; ecdf(rho.t{p}(:)); hold off

fid =fopen([pwd '/Results/bimodality - edge list.sif'],'w');
% fprintf(fid,'%s\t','Node 1');
% fprintf(fid,'%s\t','Type');
% fprintf(fid,'%s\n','Node 2');
hit = 6;
for p=1:length(CANnets{hit})
    for a=1:length(ApcLocs)
        fprintf(fid,'%s\t',CANnets{hit}{p});
        fprintf(fid,'%s\t','coexp');
        fprintf(fid,'%s\n',ppi.genes{mRNA.pos(ApcLocs(a))});
    end
end
fclose(fid);

fid =fopen([pwd '/Results/bimodality - edge attributes.eda'],'w');
fprintf(fid,'%s\n','Interaction Strength');
hit = 6;
for p=1:length(CANnets{hit})
    for a=1:length(ApcLocs)
        fprintf(fid,'%s',CANnets{hit}{p});
        fprintf(fid,'%s',[' (coexp) ' ppi.genes{mRNA.pos(ApcLocs(a))} ' = ']);
        fprintf(fid,'%f\n',rho.apc{hit}(p,a));
    end
end
fclose(fid);


fid =fopen([pwd '/Results/bimodality - node attributes.noa'],'w');
fprintf(fid,'%s\n','Expression');
hit = 6;
for p=1:length(CANnets{hit})
    fprintf(fid,'%s ',[CANnets{hit}{p} '= ']);
    fprintf(fid,'%f\n',ts(NetLocs{hit}(p)));
end
fclose(fid);

fid =fopen([pwd '/Results/bimodality - node degree.noa'],'w');
fprintf(fid,'%s\n','Degrees');
hit = 6;
for a=1:length(ApcLocs)
    fprintf(fid,'%s',[ppi.genes{mRNA.pos(ApcLocs(a))} ' = ']);
    fprintf(fid,'%f\n', round(sum(abs(rho.apc{hit}(:,a)))*100));
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now examine physical distance on ppi same way:
% Construct distance matrix 20 x n
% ApcProts, p21Prots, Apcp21net
% There are a few connections, but not many, so check 2-step distance, too
physical = zeros(length(Apcp21net),length(ppi.interax));
for i=1:length(Apcp21net)
%     physical(i,:) = ppi.interax(Apcp21net(i),:); % direct links
%     physical(i,:) = ppi.pfam(Apcp21net(i),:); % inferred links
    ind = ppi.interax*ppi.interax(:,Apcp21net(i)); % indirect links
    physical(i,setdiff(find(ind),find(ppi.interax(:,Apcp21net(i))))) = 1;
end

% What is the chance of grabbing a random handful of genes with this much
% connectivity?
this_group = ApcProts;
tot_conn = sum(physical(:,this_group)');
has_vals = find(sum(physical)~=0);
D = zeros(size(physical,1), 1000);
for p=1:10000
    r_group= round(1 + (length(has_vals)-1).*rand(length(this_group),1));
    D(:,p) = sum(physical(:,has_vals(r_group))');
end
for g=1:size(physical,1)
    pvals(g) = length(find(D(g,:)>tot_conn(g)))/length(find(D(g,:)));
    if pvals(g)<threshold
        apcphys.prots{g} = ApcProts(find(physical(g, ApcProts)));
    end
end

apcphys.p = pvals;
apcphys.votes = zeros(length(ApcLocs),1);
for v=1:length(apcphys.prots)
    [junk x y] = intersect(apcphys.prots{v},ApcProts);
    if ~isempty(junk)
        apcphys.votes(y) = apcphys.votes(y) + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine p-values
% Either Fisher's (1932) or Bailey+Gribskov (1998)
p=prod(p21coexp.p); [1-chi2cdf(-2*log(p),2*20) sum(((-log(p)).^a)./factorial(a))*p]
% or 
this_group = p21coexp.p; 1-binocdf(length(find(this_group<threshold)), length(this_group), threshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT CYTOSCAPE SIF FILES
% Apcp21net in the center
% Proteomic targets on the outside
% New edges between significant nodes
% Edges = 
fid = fopen([pwd '/Results/Apcp21net-prots.sif'],'w');
hits = find(apcphys.p<threshold);
for p=1:length(hits)
    output{1} = ppi.genes{Apcp21net(hits(p))};
    output{2} = 'ApcIndirect';
    output{3} = ppi.genes(apcphys.prots{hits(p)});
    fprintf(fid, '%s\t', output{1});
    fprintf(fid, '%s\t', output{2});
    for z=1:length(output{3})
        fprintf(fid, '%s\t', output{3}{z});
    end
    fprintf(fid, '\n');
end

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid =fopen([pwd '/Results/Just direct links.sif'],'w');
for p=1:length(Apcp21net)
    if sum(physical(p,ApcProts)>0)
        fprintf(fid,'%s ',ppi.genes{Apcp21net(p)});
        fprintf(fid,'%s ','ApcDirect');
        output = ppi.genes(ApcProts(find(physical(p,ApcProts))));
        for z=1:length(output)
            fprintf(fid, '%s ', output{z});
        end
        fprintf(fid, '\n');
    end
end

fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen([pwd '/Results/Apcp21net-coexp.sif'],'w');
hits = find(apccoexp.p<threshold);
for p=1:length(hits)
    output{1} = ppi.genes{mRNA.pos(NetLocs(hits(p)))};
    output{2} = 'ApcCoexp';
    output{3} = ppi.genes(mRNA.pos(apccoexp.prots{hits(p)}));
    fprintf(fid, '%s ', output{1});
    fprintf(fid, '%s ', output{2});
    for z=1:length(output{3})
        fprintf(fid, '%s ', output{3}{z});
    end
    fprintf(fid, '\n');
end

fclose(fid);

