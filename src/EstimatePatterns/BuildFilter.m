function [att_filt, filt] = BuildFilter(k, params, displ)
    kPix = k * params.sz(1) * params.res / pi; 
    if isfield(params,'ringMaskLim')
        rr = min(0.5,1-params.ringMaskLim(1))* norm(kPix);
    else
        rr = 0.5 * norm(kPix);
    end
    OTF = params.OTF; I = params.I; J = params.J; sz = params.sz;
    OTF0=double(OTF.*ifftshift((sqrt((I-kPix(1)-floor(sz(1)/2)-1).^2+(J-kPix(2)-floor(sz(2)/2)-1).^2)<rr)+(sqrt((I+kPix(1)-floor(sz(1)/2)-1).^2+(J+kPix(2)-floor(sz(2)/2)-1).^2)<rr))>0);
    OTFshift=ifftshift(imtranslate(fftshift(OTF),kPix(:)')+imtranslate(fftshift(OTF),-kPix(:)'));
    att_filt = OTFshift.*OTF0;
    filt = OTF0.*OTF; 
    if displ == 3
        quickPlot(fftshift(filt))
    end
end