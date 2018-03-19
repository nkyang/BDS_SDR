function channel = calCN0( channel,trackResults,subFrameStart,settings )
%¼ÆËãÔØÔë±ÈCN0
bitToCal = floor((settings.msToProcess - max(subFrameStart))/5000);
NHcode = [-1,-1,-1,-1,-1,1,-1,-1,1,1,-1,1,-1,1,-1,-1,1,1,1,-1]';
for i=1:settings.numberOfChannels
    if ~isnan(subFrameStart(i))
        Q = trackResults(i).Q_P(subFrameStart(i) : subFrameStart(i) + (bitToCal * 5000) -1);
        I = trackResults(i).I_P(subFrameStart(i) : subFrameStart(i) + (bitToCal * 5000) -1);
        if trackResults(i).PRN > 5 && trackResults(i).PRN ~= 17
            M = 20;
            I = reshape(I,M,size(I,2)/M);
            I = bsxfun(@times,I,NHcode);
            bits = 2*(sum(I)>0)-1;
            I = bsxfun(@times,I,bits);
        else
            M = 2;
            I = reshape(I,M,size(I,2)/M);
            bits = 2*(sum(I)>0)-1;
            I = bsxfun(@times,I,bits);
        end
        I = reshape(I,5000,numel(I)/5000);
        Q = reshape(Q,5000,numel(Q)/5000);
        NBP = sum(I).^2 + sum(Q).^2;
        WBP = sum(I.^2) + sum(Q.^2);
%         NBP = reshape(NBP,10000/M,size(NBP,2)*M/10000);
%         WBP = reshape(WBP,10000/M,size(WBP,2)*M/10000);
        Z   = (NBP./WBP);
        if Z > 1
            channel(i).CN0 = 10*log10(1000*(Z-1)./(5000-Z));
        else
            channel(i).CN0 = 0;
        end
    end
end
% plot(CN0');

