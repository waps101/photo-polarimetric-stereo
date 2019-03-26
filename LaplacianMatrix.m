%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Output param:
%   @m => The sparse matrix of laplacian matrix
%
%   Amos.zhu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m=LaplacianMatrix(z,mask)

rows=size(z,1);
colums=size(z,2);

offsetm=maskoffset(mask);
noofPixelsUsed=sum(mask(:));

count=1;
jidx_v=zeros(noofPixelsUsed*5,1); %jidx_v_count=1; % Jacobian matrix variable colume idx
jidx_f=zeros(noofPixelsUsed*5,1); %jidx_f_count=1; % Jacobian martrix function row idx
jidx_value=ones(noofPixelsUsed*5,1);

% jidx_v=[]; % Jacobian matrix variable colume idx
% jidx_f=[]; % Jacobian martrix function row idx

for m=1:rows
    for n=1:colums
        if mask(m,n)
%             fidx=(n-1)*rows+m-offsetm(m,n); % jacobian function idx;
%             vidx=(n-1)*rows+m-offsetm(m,n); % jacobian variable idx;
%             jidx_v(count)=vidx;
%             jidx_f(count)=fidx;
%             jidx_value(count)=-4;
%             count=count+1;
            if n-1>0&&n+1<=colums&&m-1>0&& m+1<=rows
                if mask(m,n-1)&&mask(m,n+1)&&mask(m-1,n)&&mask(m+1,n)
                    
                    fidx=(n-1)*rows+m-offsetm(m,n); % jacobian function idx;
                    vidx=(n-1)*rows+m-offsetm(m,n); % jacobian variable idx;
                    jidx_v(count)=vidx; jidx_f(count)=fidx;jidx_value(count)=-4;count=count+1;
                    
                    vidx=(n-2)*rows+m-offsetm(m,n-1); % jacobian variable idx;
                    jidx_v(count)=vidx; jidx_f(count)=fidx;jidx_value(count)=1;count=count+1;
                    vidx=n*rows+m-offsetm(m,n+1); % jacobian variable idx;
                    jidx_v(count)=vidx; jidx_f(count)=fidx;jidx_value(count)=1;count=count+1;
                    vidx=(n-1)*rows+m-1-offsetm(m-1,n); % jacobian variable idx;
                    jidx_v(count)=vidx; jidx_f(count)=fidx;jidx_value(count)=1;count=count+1;
                    vidx=(n-1)*rows+m+1-offsetm(m+1,n); % jacobian variable idx;
                    jidx_v(count)=vidx; jidx_f(count)=fidx;jidx_value(count)=1;count=count+1;
                end
            end
            
%             if n+1<=colums
%                 if mask(m,n+1)
%                     vidx=n*rows+m-offsetm(m,n+1); % jacobian variable idx;
%                     jidx_v(count)=vidx;
%                     jidx_f(count)=fidx;
%                     count=count+1;
%                 end
%             end
%             
%             if m-1>0
%                 if mask(m-1,n)
%                     vidx=(n-1)*rows+m-1-offsetm(m-1,n); % jacobian variable idx;
%                     jidx_v(count)=vidx;
%                     jidx_f(count)=fidx;
%                     count=count+1;
%                 end
%             end
%             
%             if m+1<=rows
%                 if mask(m+1,n)
%                     vidx=(n-1)*rows+m+1-offsetm(m+1,n); % jacobian variable idx;
%                     jidx_v(count)=vidx;
%                     jidx_f(count)=fidx;
%                     count=count+1;
%                 end
%             end
        end
    end
end


jidx_f=jidx_f(1:count-1);
jidx_v=jidx_v(1:count-1);
jidx_value=jidx_value(1:count-1);
m=sparse(jidx_f,jidx_v,jidx_value,noofPixelsUsed,noofPixelsUsed);

end