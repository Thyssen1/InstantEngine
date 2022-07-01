function f=dsd(nf,ncat,varargin)
%   DSD calculates definitive screening design conditions given an number  
%   of continuous (nf) and categorical (ncat) factors, based on:
%
%   Jones, B and Nachtsheim, C. (2011) 
%       "A Class of Three-Level Designs for Definitive Screening in the 
%       Presence of Second-Order Effects" 
%       Journal of Quality Technology, 43, 1-15
%
%   Jones, B and Nachtsheim, C. (2013) 
%       "Definitive Screening Designs with Added Two-Level Categorical 
%           Factors"
%       Journal of Quality Technology, 45, 121-129

%   Usage:  f=dsd(nf,ncat,designChoice)
%
%   Inputs: nf: Number of continuous factors
%           ncat: number of categorical factors
%           designChoice: 
%                   1 or 'DSD': De-alias all two-factor interactions with 
%                       categorical factors. (Default)
%                   2 or 'ORTH': Make orthogonal main-effects plan.
%
%   Outputs: f: design matrix (-1,0,or 1) with a column for each of the 
%               three level continuous variables, 1 or 2 for two level 
%               categorical variables.
%
%   Validated against the equivalent JMP10 addin for all 1736 permutations   
%   of (nf=3:30, ncat=0:30,designChoice=1:2)
%
%   See also ROWEXCH, DAUGMENT, DCOVARY, X2FX, CORDEXCH
%
%   Ported from JMP script 04-Mar-2015
%   Jacob Albrecht, BMS

if nargin==3
   if varargin{1}==1 || strcmpi('dsd',varargin{1})
       designChoice=1;
   elseif varargin{1}==2 || strcmpi('orth',varargin{1})
       designChoice=2;
   else
      error('Design Choice must be ''DSD'' (or 1) , ''ORTH'' (or 2)')
   end
elseif nargin>3
    error('Too many input arguments');
else
    designChoice=1;
end

f10 = [0 1 1 1 1 1 1 1 1 1
    1 0 -1 -1 -1 -1 1 1 1 1
    1 -1 0 -1 1 1 -1 -1 1 1
    1 -1 -1 0 1 1 1 1 -1 -1
    1 -1 1 1 0 -1 -1 1 -1 1
    1 -1 1 1 -1 0 1 -1 1 -1
    1 1 -1 1 -1 1 0 -1 -1 1
    1 1 -1 1 1 -1 -1 0 1 -1
    1 1 1 -1 -1 1 -1 1 0 -1
    1 1 1 -1 1 -1 1 -1 -1 0
    0 -1 -1 -1 -1 -1 -1 -1 -1 -1
    -1 0 1 1 1 1 -1 -1 -1 -1
    -1 1 0 1 -1 -1 1 1 -1 -1
    -1 1 1 0 -1 -1 -1 -1 1 1
    -1 1 -1 -1 0 1 1 -1 1 -1
    -1 1 -1 -1 1 0 -1 1 -1 1
    -1 -1 1 -1 1 -1 0 1 1 -1
    -1 -1 1 -1 -1 1 1 0 -1 1
    -1 -1 -1 1 1 -1 1 -1 0 1
    -1 -1 -1 1 -1 1 -1 1 1 0];

f16 = [-0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    -1 0 1 1 -1 1 -1 -1 1 -1 1 1 -1 1 -1 -1
    -1 -1 0 1 1 -1 1 -1 1 -1 -1 1 1 -1 1 -1
    -1 -1 -1 0 1 1 -1 1 1 -1 -1 -1 1 1 -1 1
    -1 1 -1 -1 0 1 1 -1 1 1 -1 -1 -1 1 1 -1
    -1 -1 1 -1 -1 0 1 1 1 -1 1 -1 -1 -1 1 1
    -1 1 -1 1 -1 -1 0 1 1 1 -1 1 -1 -1 -1 1
    -1 1 1 -1 1 -1 -1 0 1 1 1 -1 1 -1 -1 -1
    -1 -1 -1 -1 -1 -1 -1 -1 0 1 1 1 1 1 1 1
    -1 1 1 1 -1 1 -1 -1 -1 0 -1 -1 1 -1 1 1
    -1 -1 1 1 1 -1 1 -1 -1 1 0 -1 -1 1 -1 1
    -1 -1 -1 1 1 1 -1 1 -1 1 1 0 -1 -1 1 -1
    -1 1 -1 -1 1 1 1 -1 -1 -1 1 1 0 -1 -1 1
    -1 -1 1 -1 -1 1 1 1 -1 1 -1 1 1 0 -1 -1
    -1 1 -1 1 -1 -1 1 1 -1 -1 1 -1 1 1 0 -1
    -1 1 1 -1 1 -1 -1 1 -1 -1 -1 1 -1 1 1 0];
f16 = [f16 ; -1 .* f16];

minList2=ones(ncat,1);  % define levels for categorical variables
maxList2=2*ones(ncat,1);

% utility functions
% generate legendre symbol given i,j and fld
    function legendre = Function_legendre( i, j, q )
        m = j - i;
        modq = mod( q .^ 2,size( q,2 ) );
        if m == 0
            legendre= 0;
        elseif any( modq == m )
            legendre= -1;
        else
            legendre=1; 
        end
    end

% construct paley matrix given prime number
    function paley_matrix = Function_paley_matrix(q)
     
        m = zeros( q, q);
        fld = 0 : (q - 1);
        for jj = 1: q
            for j = jj: q
                m(jj, j) = Function_legendre( jj - 1, j - 1, fld );
            end
        end
        mt = m';
        if mod( q, 4 ) == mod( 3, 4 )
            factor= -1;
        else
            factor= 1;
        end
        m = m+(mt .* factor);
        paley_matrix=m; % paley construction of hadamard matrix
    end

nctn = nf;
nf = nctn + ncat; % total number of factors

if( 2 * floor( nf / 2 ) == nf)
    p = nf - 1;
else
    p = nf;
end

done = 0;
while( ~done)
    if( isprime( p ))
        c = [[0 ; ones( p, 1)] , [ones( 1, p);Function_paley_matrix( p )]];
        f = [c ; -c];
        done = 1;
    else
        p =p+ 2;
    end
end
if( nf == 10)
    f = f10;
end
if( nf == 9)
    f = f10;
    f(:, 10) = [];
end
if( nf == 16)
    f = f16 ;
end
if( nf == 15)
    f = f16;
    f(:, 16) = [];
end
if( nf == 26 || nf == 25)
    a =Function_paley_matrix( 13 );
    %% starter vector for B
    strt = [-1, -1, 1, -1, 1, 1, 1, 1, 1, -1, 1, 1, 1];
    b = [];
    %% construct B
    
    for ii = 1: 13
        b = [b'; strt ]';
        strt = circshift( strt ,[0,-1]);
    end
    c = [[a , b,] ; [b' , -1 * a]];
    if( nf == 26)
        f = [c ; -c];
    else
        f = [c ; -c];
    end
end
nc = size( f ,2);
nr = size( f ,1);
if( nc > nf)
    f(:, nf + 1 : nc) = [];
end
if( ncat == 0)
    f = [f ; zeros( 1, nf)];
end
if( ncat == 1)
    f = [f ; zeros( 2, nf)];
end
if( ncat > 1 && designChoice==2)
    f = [f ; zeros( 4, nf)];
end
if( ncat > 1&&designChoice==1)
    f =[ f ; zeros( 2, nf)];
end
tmpf = f;
for rowidx=1:nr/2
    tmpf(2*rowidx-1,:) = f(rowidx,:);
    tmpf(2*rowidx,:) = f(rowidx+nr/2,:);
end
f = tmpf;
if (designChoice == 2)
    B = [-1 -1 -1 1
        -1 -1 1 -1
        -1 1 -1 -1
        1 -1 -1 -1];
    for fidx = (nctn+1):nf
        if (ncat>1)
            colidx = mod(fidx-nctn-1,4)+1;
            f((nr+1):(nr+4),fidx)=B(:,colidx);
        end
    end
else
    B = [-1 -1 -1
        1 1 1];
    for fidx = (nctn+1):nf
        if (ncat>1)
            colidx = mod(fidx-nctn-1,3)+1;
            f((nr+1):(nr+2),fidx)=B(:,colidx);
        end
    end
end
nr = size( f ,1);
if( ncat > 0)
    for fidx = (nctn + 1): nf
        for rowidx = 1: nr
            if( f(rowidx, fidx) == 1)
                 f(rowidx,fidx) = maxList2(fidx - nctn);
            end
            if( f(rowidx, fidx) == -1)
                 f(rowidx,fidx ) = minList2(fidx - nctn);
            end
            if( f(rowidx, fidx) == 0)
                if (designChoice == 2)
                    f(rowidx,fidx ) = maxList2(fidx - nctn);
                else
                    % even rows are mins, odd rows are maxes
                    if( mod( rowidx, 2 ) == 0) 
                        f(rowidx,fidx ) = minList2(fidx - nctn);
                    else
                        f(rowidx,fidx ) = maxList2(fidx - nctn);
                    end
                end
            end
        end
    end
end


end