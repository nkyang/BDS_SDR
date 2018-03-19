function cumSignal = cumData( longSignal, interTime, settings )
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
codeDopple = -4:0.005:4;
blockSize  = round(settings.codeFreqBasis ./ codeDopple);
blockSize((blockSize>interTime * samplesPerCode)&(~isinf(blockSize))) = [];
codeDopple = settings.codeFreqBasis ./ blockSize;
blockSize  = abs(blockSize);
m = ceil(interTime * samplesPerCode ./ blockSize);
cumSignal = inf(length(codeDopple),samplesPerCode);
for i = 1:length(codeDopple)
    if codeDopple(i) < 0
        corrSignal = reshape(longSignal(1:m(i)*(blockSize(i)+1)),blockSize(i)+1,m(i));
        corrSignal = reshape(corrSignal(1:blockSize(i),:),1,m(i)*blockSize(i));
        corrSignal = sum(reshape(corrSignal(1:interTime*samplesPerCode),...
            samplesPerCode,interTime),2)';        
    elseif codeDopple(i) > 0
        corrSignal = reshape(longSignal(1:m(i)*(blockSize(i)-1)),blockSize(i)-1,m(i));
        corrSignal = reshape([corrSignal; corrSignal(blockSize(i)-1,:)],1,m(i)*blockSize(i));
        corrSignal = sum(reshape(corrSignal(1:interTime*samplesPerCode),...
            samplesPerCode,interTime),2)';
    else
        corrSignal = sum(reshape(longSignal(1:interTime*samplesPerCode),...
            samplesPerCode,interTime),2)';
    end
    cumSignal(i,:) = corrSignal;
end
end 