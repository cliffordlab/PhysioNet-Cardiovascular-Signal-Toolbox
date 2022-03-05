function dlY = spatialDropout(dlX,dropoutFactor,fmt,doTraining)

if doTraining
    maskSize = size(dlX);
    maskSize(fmt=='S') = 1;
    
    dropoutScaleFactor = single(1 - dropoutFactor);
    dropoutMask = (rand(maskSize,'like',dlX) > dropoutFactor) / dropoutScaleFactor;
    
    dlY = dlX .* dropoutMask;
else
    dlY = dlX;
end

end