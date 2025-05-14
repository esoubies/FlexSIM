function im = ConvertDataType(im,datatype)
switch datatype
    case 'uint8'
        im = uint8(im/max(im(:))*(2^8-1));
    case 'uint16'
        im = uint16(im/max(im(:))*(2^16-1));
    case 'uint32'
        im = uint32(im/max(im(:))*(2^32-1));
    case 'single'
        im = single(im);
end
end