%support altering employed "global" method and exchange one index pair
function out_x=sa_cwmax(in_x,R)%in_x is a PC-loading vector. R is the covariance matrix
%find the support of in_x;
supp_init=find(in_x);
ze_init=find(in_x==0);
card_init=length(supp_init);%the cardinality of in_x;
if card_init~=1
    %enumerate all non-zero elements
    in_x_supp_init=full(in_x(supp_init));
    sum2_init=in_x_supp_init'*in_x_supp_init;
    sum3_init=in_x_supp_init'*R(supp_init,supp_init)*in_x_supp_init;
    DIAG_R=diag(R);
    RQQ=DIAG_R(ze_init);
    RQQ_ALL=repmat(RQQ,[1,card_init]);
    ze_num=length(ze_init);
    SUM1_PART1=R(ze_init,supp_init)*in_x_supp_init;
    SUM1_MINUS=R(ze_init,supp_init).*repmat(in_x_supp_init',[ze_num,1]);
    SUM1_ALL=repmat(SUM1_PART1,[1,card_init])-SUM1_MINUS;
    
    SUM2_INIT=repmat(sum2_init,[ze_num,card_init]);
    in_x_supp_init_2=in_x_supp_init.^2;
    SUM2_ALL=SUM2_INIT-repmat(in_x_supp_init_2',[ze_num,1]);
    
    SUM3_INIT=repmat(sum3_init,[ze_num,card_init]);
    SUPP=repmat(supp_init,[1,card_init]);
    SUPP_LOWER=tril(SUPP,-1);
    SUPP_UPPER=triu(SUPP,1);
    SUPP=SUPP_LOWER(2:end,:)+SUPP_UPPER(1:end-1,:);
    R_SUPP_INIT_SUPP=zeros(card_init,card_init-1);
    for k=1:card_init
        R_SUPP_INIT_SUPP(k,:)=R(supp_init(k),SUPP(:,k));
    end
    IN_X_SUPP=repmat(in_x_supp_init,[1,card_init]);
    IN_X_SUPP_LOWER=tril(IN_X_SUPP,-1);
    IN_X_SUPP_UPPER=triu(IN_X_SUPP,1);
    IN_X_SUPP=IN_X_SUPP_LOWER(2:end,:)+IN_X_SUPP_UPPER(1:end-1,:);
    SUM3_ALL_PART2=zeros(ze_num,card_init);
    for k=1:card_init
        SUM3_ALL_PART2(:,k)=repmat(2*in_x_supp_init(k)*(R_SUPP_INIT_SUPP(k,:)*IN_X_SUPP(:,k)),[ze_num,1]);
    end
    SUM3_ALL=SUM3_INIT-repmat((DIAG_R(supp_init).*in_x_supp_init_2)',[ze_num,1])-SUM3_ALL_PART2;
    
    B_ALL=SUM3_ALL-RQQ_ALL.*SUM2_ALL;
    C_ALL=-(SUM1_ALL.*SUM2_ALL);
    DELTA_ALL=B_ALL.^2-4*SUM1_ALL.*C_ALL;
    ROOT_DELTA_ALL=sqrt(DELTA_ALL);
    ROOT2_ALL=(-B_ALL+ROOT_DELTA_ALL)./(2*SUM1_ALL);
    SUM1eq0=SUM1_ALL==0;    %b=0
    Bge0=B_ALL>=0;  %ad-c<=0
    ROOT_ALL=~SUM1eq0.*ROOT2_ALL+SUM1eq0.*~Bge0*2^64;  
    VALUE_ALL=(SUM3_ALL+2*SUM1_ALL.*ROOT_ALL+RQQ_ALL.*ROOT_ALL.^2)./(SUM2_ALL+ROOT_ALL.^2);
    
    %find the desired exchange
    if length(ze_init)~=1
        [column_max,id_in]=max(VALUE_ALL);
        [~,id_out]=max(column_max);
        id_in=id_in(id_out);
    else
        id_in=1;
        [~,id_out]=max(VALUE_ALL);
    end
    %desired pc-loading
    out_x=in_x;%initialize
    out_x(supp_init(id_out))=0;
    out_x(ze_init(id_in))=ROOT_ALL(id_in,id_out);
    %normalize
    out_x=out_x/norm(out_x);
else
    [~,id]=max(diag(R));
    out_x=zeros(length(R),1);
    out_x(id)=1;
end
end