%support altering employed "partial" method
function [out_x,out_X]=sa_partial(in_x,R,change_num)%in_x is a PC-loading vector. R is the covariance matrix
%find the support of in_x;
supp_init=find(in_x);
ze_init=find(in_x==0);
card_init=length(supp_init);%the cardinality of in_x;
if card_init~=1
    %initialize change records
    in_x_supp_init=full(in_x(supp_init));
    sum2_init=in_x_supp_init'*in_x_supp_init;
    sum3_init=in_x_supp_init'*R(supp_init,supp_init)*in_x_supp_init;
    sum2=sum2_init;
    sum3=sum3_init;
    Q=ze_init;
    DIAG_R=diag(R);
    RQQ=DIAG_R(ze_init);
    out_x=in_x;%init value
    out_X=zeros(length(in_x),change_num);
    supp_now=supp_init;
    [~,smallest_id_vec]=mink(abs(out_x(supp_init)),change_num);
    SUM1=R(ze_init,supp_init)*out_x(supp_init);
    for i=1:change_num
        %set the i-th smallest element of in_x to zero
        %initialize
        supp=supp_now;
        smallest_id=smallest_id_vec(i);
        id1=supp_init(smallest_id)==supp;
        supp(id1)=[];
        %update sum2 & sum3
        sum2=sum2-out_x(supp_init(smallest_id))^2;
        sum3=sum3-(R(supp_init(smallest_id),supp_init(smallest_id))*out_x(supp_init(smallest_id))^2+2*out_x(supp_init(smallest_id))*(R(supp_init(smallest_id),supp)*out_x(supp)));
        %SUM1=R(Q,supp)*out_x(supp);
        SUM1=SUM1-out_x(supp_init(smallest_id))*R(Q,supp_init(smallest_id));
        out_x(supp_init(smallest_id))=0;
        %
        length_Q=length(Q);
        SUM3=repmat(sum3,[length_Q,1]);
        SUM2=repmat(sum2,[length_Q,1]);
        B=SUM3-RQQ.*SUM2;
        C=-(SUM1.*SUM2);
        DELTA=B.^2-4*SUM1.*C;
        ROOT_DELTA=sqrt(DELTA);
        ROOT2=(-B+ROOT_DELTA)./(2*SUM1);    %the larger root of the quadratic formula
        SUM1eq0=SUM1==0;    %b=0
        Bge0=B>=0;  %ad-c<=0
        ROOT=~SUM1eq0.*ROOT2+SUM1eq0.*~Bge0*2^64;
        VALUE=(SUM3+2*SUM1.*ROOT+RQQ.*ROOT.^2)./(SUM2+ROOT.^2);
        [~,id]=max(VALUE);
        %out_x
        out_x(Q(id))=ROOT(id);
        out_x=out_x/norm(out_x);
        out_X(:,i)=out_x;
        if i~=change_num
            %update supp_now
            supp_now=sort([supp;Q(id)]);
            %update sum2 sum3
            sum2=sum2+ROOT(id)^2;
            sum3=sum3+(R(Q(id),Q(id))*ROOT(id)^2+2*ROOT(id)*(R(Q(id),supp)*out_x(supp)));
            %add to SUM1
            SUM1=SUM1+ROOT(id)*R(Q,Q(id));
            %delete No. id row of SUM1
            SUM1(id,:)=[];
            %update Q
            Q(id)=[];
            %update RQQ
            RQQ(id)=[];
        end
    end
else
    [~,id]=max(diag(R));
    out_x=zeros(length(R),1);
    out_x(id)=1;
    out_X=out_x;
end
end