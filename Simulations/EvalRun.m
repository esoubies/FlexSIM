function res = EvalRun(params, res)
%--------------------------------------------------------------------------
% Function eval = EvalRun(params)
% 
% Returns a structure that contains the results of the simulation. This
% inlcudes the accuracy of the parameter estimation, and the SNR of the 
% reconstruction as measured against the original image.
%
% Inputs : params -> Structure with fields  
%                         - or: Orientations of the SIM images [rad]
%                         - ph: Phases of the SIM
%                         - DataPath: Path of the original image
%          res    -> Structure with fields  
%                         - k: Estimated wavevectors 
%                         - phase: Estimated phases
%                         - rec: reconstruction 
%
% Outputs: res    -> Input structure with the additional fields
%                         - kGt: Groundtruth wavevector of the acquisition
%                         - kErr: Wavevector error expressed as norm
%                         - phGt: Groundtruth phase of the acqu
%                         - phErr: Phase error expressed as norm
%                         - MEP: Maximum expected number of photons during acq
%                         - snr: SNR of the reconstruction
%                         - ssi: SSI of the reconstruction
%
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

kGt = zeros(params.nbOr, 2);                          % Calculate the ground truth wavevector
for or = 1:params.nbOr                                % Pattern generation
    kGt(or,:) = 2*pi*params.ns/params.lamb*[cos(params.or(or)), sin(params.or(or))]*params.Na/params.nl;
end
res.kGt=kGt.*sign(kGt(:,1)).*sign(res.k(:,1));      % Store gt wavevector with sign convention on the gt as on the results       
res.kErr=vecnorm(res.kGt - res.k, 2, 2);                  % Store error as norm

res.phGt = params.ph;                               % Store phase GT
if params.method == 2                               % And error, depending on the chosen method
    res.phErr = abs(res.phase - params.ph(1));
else     
    res.phErr = abs(res.phase - repmat(params.ph, ([params.nbOr, 1])));
end

x0 = loadtiff('phantom.tif');
y = res.rec;
y = cast(y, 'double');
y = (y - min(y(:)))/(max(y(:)) - min(y(:)));
res.snr = 20*log10(norm(x0(:))/norm(x0(:)-y(:)));
res.MEP = params.MEP; 
if params.sav
    prefix = char(params.DataPath); 
    prefix = prefix(1:end-4);
    save(strcat(prefix,'_Results.m'),'res');
end
end