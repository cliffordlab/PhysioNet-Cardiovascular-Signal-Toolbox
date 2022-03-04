function dlY = residualBlock(dlX,dilationFactor,dropoutFactor,parametersBlock,doTraining)

% Convolution options.
filterSize = size(parametersBlock.Conv1.Weights,1);
paddingSize = (filterSize - 1) * dilationFactor;

% Convolution.
weights = parametersBlock.Conv1.Weights;
bias =  parametersBlock.Conv1.Bias;
dlY = dlconv(dlX,weights,bias, ...
    'DataFormat','CBS', ...
    'Stride', 1, ...
    'DilationFactor', dilationFactor, ...
    'Padding', [paddingSize; 0] );

% Instance normalization, ReLU, spatial dropout.
dlY = instanceNormalization(dlY,'CBS');
dlY = relu(dlY);
dlY = spatialDropout(dlY,dropoutFactor,'CBS',doTraining);

% Convolution.
weights = parametersBlock.Conv2.Weights;
bias = parametersBlock.Conv2.Bias;
dlY = dlconv(dlY,weights,bias, ...
    'DataFormat','CBS', ...
    'Stride', 1, ...
    'DilationFactor', dilationFactor, ...
    'Padding',[paddingSize; 0] );

% Instance normalization, ReLU, spatial dropout.
dlY = instanceNormalization(dlY,'CBS');
dlY = relu(dlY);
dlY = spatialDropout(dlY,dropoutFactor,'CBS',doTraining);

% Optional 1-by-1 convolution.
if ~isequal(size(dlX),size(dlY))
    weights = parametersBlock.Conv3.Weights;
    bias = parametersBlock.Conv3.Bias;
    dlX = dlconv(dlX,weights,bias,'DataFormat','CBS');
end

% Addition and ReLU
dlY = relu(dlX+dlY);

end