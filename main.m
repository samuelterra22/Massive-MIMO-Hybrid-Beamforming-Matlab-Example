s = rng(67);                        % Set RNG state for repeatability

%%%%%%%%%%%%%%% DEFINE SYSTEM PARAMETERS FOR THE EXAMPLE %%%%%%%%%%%%%%%%%%

% Multi-user system with single/multiple streams per user
prm.numUsers = 4;                   % Number of users
prm.numSTSVec = [3 2 1 2];          % Number of independent data streams per user
prm.numSTS = sum(prm.numSTSVec);    % Must be a power of 2
prm.numTx = prm.numSTS * 8;         % Number of BS transmit antennas (power of 2)
prm.numRx = prm.numSTSVec * 4;      % Number of receive antennas, per user (any >= numSTSVec)

% Each user has the same modulation
prm.bitsPerSubCarrier = 4;          % 2: QPSK, 4: 16QAM, 6: 64QAM, 8: 256QAM
prm.numDataSymbols = 10;            % Number of OFDM data symbols

% MS positions: assumes BS at origin
%   Angles specified as [azimuth;elevation] degrees
%   az in range [-180 180], el in range [-90 90], e.g. [45;0]
maxRange = 1000;                    % all MSs within 1000 meters of BS
prm.mobileRanges = randi([1 maxRange], 1, prm.numUsers);
prm.mobileAngles = [rand(1, prm.numUsers) * 360 - 180; ...
                    rand(1, prm.numUsers) * 180 - 90];

prm.fc = 28e9;                      % 28 GHz system
prm.chanSRate = 100e6;              % Channel sampling rate, 100 Msps
prm.ChanType = 'Scattering';        % Channel options: 'Scattering', 'MIMO'
prm.NFig = 8;                       % Noise figure (increase to worsen, 5-10 dB)
prm.nRays = 500;                    % Number of rays for Frf, Fbb partitioning

%%%%%%% DEFINE OFDM MODULATION PARAMETERS USEDS FOR THE SYSTEM   %%%%%%%%%%

prm.FFTLength = 256;
prm.CyclicPrefixLength = 64;
prm.numCarriers = 234;              % Number of carries
prm.NullCarrierIndices = [1:7 129 256-5:256]'; % Guards and DC
prm.PilotCarrierIndices = [26 54 90 118 140 168 204 232]';
nonDataIdx = [prm.NullCarrierIndices; prm.PilotCarrierIndices];
prm.CarriersLocations = setdiff((1:prm.FFTLength)', sort(nonDataIdx));

numSTS = prm.numSTS;
numTx = prm.numTx;
numRx = prm.numRx;
numSTSVec = prm.numSTSVec;
codeRate = 1/3;                     % same code rate per user
numTails = 6;                       % number of termination tail bits
prm.numFrmBits = numSTSVec.* (prm.numDataSymbols * prm.numCarriers * ...
                 prm.bitsPerSubCarrier * codeRate) - numTails;
prm.modMode = 2 ^ prm.bitsPerSubCarrier; % Modulation order
% Account for channel filter delay
numPadSym = 3;                      % number of symbols to zeropad
prm.numPadZeros = numPadSym* (prm.FFTLength + prm.CyclicPrefixLength);

%%%%%%%%%%%%%%% DEFINE TRANSMITE AND RECEIVE ARRAYS AND %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% POSITIONAL PARAMETERS FOR THE SYSTEM %%%%%%%%%%%%%%%%%%%%

prm.cLight = physconst('LightSpeed');
prm.lambda = prm.cLight / prm.fc;

% Get transmit and receive array information
[isTxURA, expFactorTx, isRxURA, expFactorRx] = helperArrayInfo(prm, true);

% Transmit antenna array definition
%   Array locations and angles
prm.posTx = [0; 0; 0];              % BS/Transmit array position, [x; y; z], meters
if isTxURA
    % Uniform Rectangular array
    txArray = phased.PartitionedArray(...
        'Array', phased.URA([expFactorTx numSTS], 0.5 * prm.lambda), ...
        'SubarraySelection', ones(numSTS, numTx), 'SubarraySteering', 'Custom');
else
    % Uniform Linear array
    txArray = phased.ULA(numTx, 'ElementSpacing', 0.5 * prm.lambda, ...
        'Element', phased.IsotropicAntennaElement('BackBaffled', false));
end
prm.posTxElem = getElementPosition(txArray) / prm.lambda;

spLoss = zeros(prm.numUsers, 1);
prm.posRx = zeros(3, prm.numUsers);
for uIdx = 1:prm.numUsers

    % Receive arrays
    if isRxURA(uIdx)
        % Uniform Rectangular array
        rxarray = phased.PartitionedArray(...
            'Array', phased.URA([expFactorRx(uIdx) numSTSVec(uIdx)], ...
            0.5*prm.lambda), 'SubarraySelection', ones(numSTSVec(uIdx), ...
            numRx(uIdx)), 'SubarraySteering', 'Custom');
        prm.posRxElem = getElementPosition(rxarray) / prm.lambda;
    else
        if numRx(uIdx) > 1
            % Uniform Linear array
            rxarray = phased.ULA(numRx(uIdx), ...
                'ElementSpacing', 0.5 * prm.lambda, ...
                'Element', phased.IsotropicAntennaElement);
            prm.posRxElem = getElementPosition(rxarray)/prm.lambda;
        else
            rxarray = phased.IsotropicAntennaElement;
            prm.posRxElem = [0; 0; 0]; % LCS
        end
    end

    % Mobile positions
    [xRx, yRx, zRx] = sph2cart(deg2rad(prm.mobileAngles(1, uIdx)), ...
                             deg2rad(prm.mobileAngles(2, uIdx)), ...
                             prm.mobileRanges(uIdx));
    prm.posRx(:,uIdx) = [xRx; yRx; zRx];
    [toRxRange, toRxAng] = rangeangle(prm.posTx, prm.posRx(:, uIdx));
    spLoss(uIdx) = fspl(toRxRange, prm.lambda);
end


%%%%%%%%%%%%%%%%%%%%%%% CHANNEL STATE INFORMATION %%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the preamble signal
prm.numSTS = numTx;                 % set to numTx to sound out all channels
preambleSig = helperGenPreamble(prm);

% Transmit preamble over channel
prm.numSTS = numSTS;                % keep same array config for channel
[rxPreSig, chanDelay] = helperApplyMUChannel(preambleSig, prm, spLoss);

% Channel state information feedback
hDp = cell(prm.numUsers, 1);
prm.numSTS = numTx;                 % set to numTx to estimate all links
for uIdx = 1:prm.numUsers

    % Front-end amplifier gain and thermal noise
    rxPreAmp = phased.ReceiverPreamp( ...
        'Gain', spLoss(uIdx), ...    % account for path loss
        'NoiseFigure', prm.NFig,'ReferenceTemperature',290, ...
        'SampleRate', prm.chanSRate);
    rxPreSigAmp = rxPreAmp(rxPreSig{uIdx});
    %   scale power for used sub-carriers
    rxPreSigAmp = rxPreSigAmp * (sqrt(prm.FFTLength - ...
        length(prm.NullCarrierIndices)) / prm.FFTLength);

    % OFDM demodulation
    rxOFDM = ofdmdemod(rxPreSigAmp(chanDelay(uIdx) + 1: ...
        end - (prm.numPadZeros-chanDelay(uIdx)), :), prm.FFTLength, ...
        prm.CyclicPrefixLength, prm.CyclicPrefixLength, ...
        prm.NullCarrierIndices, prm.PilotCarrierIndices);

    % Channel estimation from preamble
    %       numCarr, numTx, numRx
    hDp{uIdx} = helperMIMOChannelEstimate(rxOFDM(:, 1:numTx, :), prm);

end


%%%%%%%%%%%%%%%%%%%%%%%%%% HYBRID BEAMFORMING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate the hybrid weights on the transmit side
if prm.numUsers == 1
    % Single-user OMP
    %   Spread rays in [az;el]=[-180:180;-90:90] 3D space, equal spacing
    %   txang = [-180:360/prm.nRays:180; -90:180/prm.nRays:90];
    txang = [rand(1, prm.nRays) * 360 - 180; rand(1, prm.nRays) * 180 - 90]; % random
    At = steervec(prm.posTxElem, txang);

    Fbb = complex(zeros(prm.numCarriers, numSTS, numSTS));
    Frf = complex(zeros(prm.numCarriers, numSTS, numTx));
    for carrIdx = 1:prm.numCarriers
        [Fbb(carrIdx, :, :), Frf(carrIdx, :, :)] = helperOMPTransmitWeights( ...
            permute(hDp{1}(carrIdx, :, :), [2 3 1]), numSTS, numSTS, At);
    end
    v = Fbb;    % set the baseband precoder (Fbb)
    % Frf is same across subcarriers for flat channels
    mFrf = permute(mean(Frf, 1), [2 3 1]);

else
    % Multi-user Joint Spatial Division Multiplexing
    [Fbb, mFrf] = helperJSDMTransmitWeights(hDp, prm);

    % Multi-user baseband precoding
    %   Pack the per user CSI into a matrix (block diagonal)
    steeringMatrix = zeros(prm.numCarriers, sum(numSTSVec), sum(numSTSVec));
    for uIdx = 1:prm.numUsers
        stsIdx = sum(numSTSVec(1:uIdx-1)) + (1:numSTSVec(uIdx));
        steeringMatrix(:, stsIdx, stsIdx) = Fbb{uIdx};  % Nst-by-Nsts-by-Nsts
    end
    v = permute(steeringMatrix, [1 3 2]);

end

% Transmit array pattern plots
if isTxURA
    % URA element response for the first subcarrier
    pattern(txArray, prm.fc, -180:180, -90:90, 'Type', 'efield', ...
            'ElementWeights', mFrf.' * squeeze(v(1, :, :)), ...
            'PropagationSpeed', prm.cLight);
else % ULA
    % Array response for first subcarrier
    wts = mFrf.' * squeeze(v(1, :, :));
    pattern(txArray, prm.fc, -180:180, -90:90, 'Type', 'efield', ...
            'Weights', wts(:, 1), 'PropagationSpeed', prm.cLight);
end

prm.numSTS = numSTS;                 % revert back for data transmission
