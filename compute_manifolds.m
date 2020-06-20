function [manifolds,partitions]=compute_manifolds(G,couplingtype,varargin);

% manifolds: all possible partial synchronization manifolds
% partitions: all possible partitions of the systems

% use: 
% 1) adjacency matrix consisting of integer elments:
%
% manifolds=compute_manifolds(G,couplingtype);
% G: adjacency matrix 
% coupling type: 1 (invasive coupling) or 2 (non-invasive coupling)
%
% 
% 2) real adjacency matrix not consisting of integer elments:
%
% manifolds=compute_manifolds(G,couplingtype,tolerance);
% G: adjacency matrix 
% coupling type: 1 (invasive coupling) or 2 (non-invasive coupling)
% tolerance: tolerance in comparison of row-sums (finite precision
% arithmetic)
% 
%
%
% 3) adjacency matrix that can be decomposed as
%   sum_{i=1}^m G_i r_i   with {r1,... r_m} rationally independent
%    and G_i consisting of positive integer elments, i=1,...,m
%
% manifolds=compute_manifolds(G,couplingtype);
% G: structured array of length m; G{i}=G_i 
% coupling type: 1 (invasive coupling) or 2 (non-invasive coupling)
% tolerance: tolerance in comparison of row-sums
% 
%
% Based on the the input type of G, the function automatically
% decides between case 3) and cases 1)-2)





%% Initialization
manifolds=[];

% Dimensions
if iscell(G)==1, 
    N=length(G{1});
else N=length(G);
end

%% part 1: create all combinations
% G is a cell 
if couplingtype ==1
    if nargin ==3
        tolerance = varargin{1};
    [A, partitions]=generate_partition_invasive(G,tolerance);
    else 
    [A, partitions]=generate_partition_invasive(G);
    end
else
    [A, partitions]=generate_partition_universal(G);
end

%disp('click to continue')
%pause

%% part 2: check manifolds - conditions with row sums

if ~iscell(G),
N=length(G);

if norm(mod(G,1))==0,
disp('integer adjacency matrix')

manifolds=[];

if couplingtype==1;
    disp('invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N;
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    Gs=EE*G*EE';
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs);
        if rijs(m)~=rijs(m-1),
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1);
        for tel2=1:(length(indices)-1);
          %if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
            block= Gs(positionsh,positionsv); 
            blocksum=rowsum(block);
            if blocksum==0,
               isman=0;
            end         
           if isman==0, break, end
        end
      if isman==0, break, end
    end
    if isman==1,
        manifolds=[manifolds;rij]; 
    end
end

elseif couplingtype==2;
        disp('non-invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N;
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    Gs=EE*G*EE';
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs);
        if rijs(m)~=rijs(m-1),
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1);
        for tel2=1:(length(indices)-1);
          if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
            block= Gs(positionsh,positionsv); 
            blocksum=rowsum(block);
            if blocksum==0,
               isman=0;
            end         
           if isman==0, break, end
          end
        end
      if isman==0, break, end
    end
    if isman==1,
        manifolds=[manifolds;rij]; 
    end
end



else disp('invalid coupling type number')
return
end
    



else % non-integer row sum

if nargin~=3
disp('incorrected number of inputs (tolerance for non-integer adjacency matrix to be specified)');
return
else
epsilon=varargin{1};



manifolds=[];

if couplingtype==1;
    disp('invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N;
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    Gs=EE*G*EE';
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs);
        if rijs(m)~=rijs(m-1),
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1);
        for tel2=1:(length(indices)-1);
          %if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
            block= Gs(positionsh,positionsv); 
            blocksum=rowsumnint(block);
            if blocksum==0,
               isman=0;
            end         
           if isman==0, break, end
        end
      if isman==0, break, end
    end
    if isman==1,
        manifolds=[manifolds;rij]; 
    end
end

elseif couplingtype==2;
        disp('non-invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N;
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    Gs=EE*G*EE';
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs);
        if rijs(m)~=rijs(m-1),
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1);
        for tel2=1:(length(indices)-1);
          if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
            block= Gs(positionsh,positionsv); 
            blocksum=rowsumnint(block);
            if blocksum==0,
               isman=0;
            end         
           if isman==0, break, end
          end
        end
      if isman==0, break, end
    end
    if isman==1,
        manifolds=[manifolds;rij]; 
    end
end



else disp('invalid coupling type number')
end
    


end
end




%
% iscell loop
else

mm=length(G);


manifolds=[];

if couplingtype==1;
    disp('invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N;
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrix
    %%
    for kk=1:mm,
          Gs{kk}=EE*G{kk}*EE';
    end
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs);
        if rijs(m)~=rijs(m-1),
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1);
        for tel2=1:(length(indices)-1);
          %if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
		 %%	            
		  for kk=1:mm,	
		        block= Gs{kk}(positionsh,positionsv); 
                blocksum=rowsum(block);
                 if blocksum==0,
                 isman=0;
                end         
            end
            %%
            %%
            if isman==0, break, end
        end
      if isman==0, break, end
    end
    if isman==1,
        manifolds=[manifolds;rij]; 
    end
end

elseif couplingtype==2;
        disp('non-invasive coupling')

for k=1:size(A,1)
    % create permutation matrix
    rij=A(k,:);
    [rijs,E]=sort(rij);
    EE=zeros(N,N);
    for l=1:N;
        EE(l,E(l))=1;
    end
    % reordened adjacancy matrices
    %%
    for kk=1:mm,
         Gs{kk}=EE*G{kk}*EE';
    end
    %%
    % indices= positions where new block begins
    indices=1;
    for m=2:length(rijs);
        if rijs(m)~=rijs(m-1),
            indices=[indices m];
        end
    end
    indices=[indices N+1];
    isman=1;
    for tel1=1:(length(indices)-1);
        for tel2=1:(length(indices)-1);
          if tel1~=tel2,
            positionsh=indices(tel1):(indices(tel1+1)-1);
            positionsv=indices(tel2):(indices(tel2+1)-1);
		 %%
            for kk=1:mm,
                 block= Gs{kk}(positionsh,positionsv); 
                 blocksum=rowsum(block);
                 if blocksum==0,
                    isman=0;
                 end
            end
            %%         
           if isman==0, break, end
          end
        end
      if isman==0, break, end
    end
    if isman==1,
        manifolds=[manifolds;rij]; 
    end
end



else disp('invalid coupling type number')
return
end



%
% iscell loop
%
end


    
%% inline functions    
    
function y=rowsum(bl);
    [nn]=size(bl);
    bl2=bl*ones(nn(2),1);
    bl2=abs(bl2-bl2(1)*ones(nn(1),1));
   if max(bl2)==0,
       y=1;
   else
       y=0;
   end 
end
    

function y=rowsumnint(bl)
    [nn]=size(bl);
    bl2=bl*ones(nn(2),1);
    bl2=abs(bl2-bl2(1)*ones(nn(1),1));
   if max(bl2)<epsilon,
       y=1;
   else
       y=0;
   end 
end

function [A,partitions]=generate_partition_universal(G)
% Dimension
    if iscell(G)==1, 
    N=length(G{1});
    else N=length(G);
    end
    
    A=[0];
    for k=2:N;
    B=[];
    for l=1:size(A,1)
       for m=0:max(A(l,:))+1 
       B=[B; A(l,:) m];
       end
    end
    A=B;

    end
    partitions=A;
    A=A(1:end-1,:);
end

function [A,partitions]=generate_partition_invasive(G, varargin)
    % Dimension
    if iscell(G)==1, 
    N=length(G{1});
    else N=length(G);
    end
    % G is a matrix
    if ~iscell(G)
    Grow = G*ones(size(G,2),1);
        if norm(mod(G,1))==0
        %disp('integer adjacency matrix')
        [Grows,IA,IC] = unique(Grow);
        else % non-integer row sum
            if nargin~=2
            disp('incorrected number of inputs (tolerance for non-integer adjacency matrix to be specified)');
            return
            else
        tolerance=varargin{1};
        tol=tolerance/max(abs(Grow(:)));
        [Grows,IA,IC] = uniquetol(Grow,tol);
            end
        end

        if length(IA)==N
            disp('No partial synchronization manifolds exist')
            return
        end

        A=[0];
        for k=2:N
            B=[];
            for l=1:size(A,1)
                index=(IC(k)==IC(1:k-1));  % check the system k have the same row sum with systems 1,2,...,k-1 
                Al=A(l,:)';                        %  the values of system k can be chosen from by the row sums
                P=unique(Al(index));        %  remove duplicated entries
                P=[P;max(Al)+1];             %  add the situation where the k system is in a new group 
                for m=1:length(P)
                    B=[B;A(l,:) P(m)];
                end
            end
            A=B;
        end
        partitions=A;
        A=A(1:end-1,:);
  % G is a cell 
    else 
        r=size(G,2);  % G={G1, G2, ... Gr} 
        Grow=[]; % matrix to store sum of row in G1, G2, ..., Gr
        for ii=1:r
        if norm(mod(G{1,ii},1))==0
        Gp=G{1,ii}*ones(N,1);
        else 
        if nargin~=2
            disp('incorrected number of inputs (tolerance for non-integer adjacency matrix to be specified)');
            return
        else
        Gp=G{1,ii}*ones(N,1);
        tolerance=varargin{1};
        tol=tolerance/max(abs(Gp(:)));
        [Gu,IA,Gp] = uniquetol(Gp,tol) % Gp indicates which sums of row are equal w.r.t. the given tolerance
            if length(IA)==N
            disp('No partial synchronization manifolds exist')
            return
            end
        end
        end
        Grow=[Grow Gp];
        end
        IC=[];
        Skip=[];
        for ii=1:N
            if any(ii==Skip)
                continue;
            end
            index=[];
            for jj=1:r
            temp=(Grow(ii,jj)==Grow(:,jj));
            index=[index temp];
            end
            index_com=prod(index,2); % find common indicator of constant row sum of all matrices 
            Skip=[Skip;find(index_com)]; % store which nodes can be skipped
            IC=[IC index_com];
        end
            if size(IC,2)==N
            disp('No partial synchronization manifolds exist')
            return
            end
        IC=IC*[1:size(IC,2)]';

        
        A=[0];
        for k=2:N
            B=[];
            for l=1:size(A,1)
                index=(IC(k)==IC(1:k-1));
                Al=A(l,:)';
                P=unique(Al(index));
                P=[P;max(Al)+1];
                for m=1:length(P)
                    B=[B;A(l,:) P(m)];
                end
            end
            A=B;
        end
        partitions=A;
        A=A(1:end-1,:);
    end 
end 
% end main function
end