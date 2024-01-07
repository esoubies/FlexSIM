function finfo = nd2finfo(file, varargin)
%ND2FINFO analyses Nikon ND2 file structure and extracts metadata
% finfo = nd2finfo(file) returns file structure ( including header of data
% segment, data start position, and data length ) and basic image
% attributes including image width, image height, channel count, and stage
% positions. 
%
% Version 1.2
% Copyright, Chao-Yuan Yeh, 2017

fid = fopen(file, 'r');

% Each datasegment begins with these signature bytes
sigbyte = [218, 206, 190, 10];

count = 1;
fs = struct;
signature = fread(fid, 4, '*uint8')';
if ~isempty(signature) && strfind(signature, sigbyte)
    [fs, count, ~] = readHeader(fid, fs, count);
end

% Second segment always begins at 4096
fseek(fid, 4096, 'bof');

flag = 0;
while flag == 0
    signature = fread(fid, 4, '*uint8')';
    if ~isempty(signature) && sum(signature == sigbyte) == 4
        [fs, count, flag] = readHeader(fid, fs, count);
        if strfind(fs(count-1).nameAttribute, 'ImageDataSeq')
            break;
        end
    else
        break
    end
end

% In a large file, there is often a seemingly random number of zeros
% between the ImageDataSeq segment and its following segment. Other
% segments usually don't have these randome padding bytes. 
temp = fread(fid, 10000, '*uint8')';
idx_nxt = strfind(temp, sigbyte);
fseek(fid, fs(count-1).dataStartPos + fs(count-1).dataLength + ...
    idx_nxt(1)-1, 'bof');

flag = 0;
while flag == 0
  signature = fread(fid, 4, '*uint8')';
  if ~isempty(signature) && sum(signature == sigbyte) == 4
      [fs, count, flag] = readHeader(fid, fs, count);   
  else
      if isempty(signature)
          break
      end
  end
end

for ii = 1 : length(fs)
  if strfind(fs(ii).nameAttribute, 'ImageAttributesLV!')
    img_attrib_idx = ii;
  elseif strfind(fs(ii).nameAttribute, 'ImageDataSeq')
    img_data_idx = ii;
  elseif strfind(fs(ii).nameAttribute, 'ImageMetadataSeq')
    img_metadata_idx = ii;
  elseif strfind(fs(ii).nameAttribute, 'ImageCalibration')
    img_calib_idx = ii;
  elseif strfind(fs(ii).nameAttribute, 'ImageTextInfo')
    img_txt_idx = ii;  
  end
end

finfo = struct;
finfo.file_structure = fs;
fseek(fid, fs(img_attrib_idx).dataStartPos, 'bof');
img_attrib = fread(fid, fs(img_attrib_idx).dataLength, '*char')';


% fseek(fid, fs(img_attrib_idx).dataStartPos + ...
%     strfind(img_attrib, insert0('uiWidth')) + length(insert0('uiWidth')), 'bof');
strloc(fid,fs, img_attrib_idx, img_attrib, 'uiWidth');
finfo.img_width = fread(fid, 1, '*uint32');
strloc(fid, fs, img_attrib_idx, img_attrib, 'uiHeight');
finfo.img_height = fread(fid, 1, '*uint32');
strloc(fid, fs, img_attrib_idx, img_attrib, 'uiSequenceCount');
finfo.img_seq_count = fread(fid, 1, '*uint32');

fseek(fid, fs(img_metadata_idx).dataStartPos, 'bof');
img_metadata = fread(fid, fs(img_metadata_idx).dataLength, '*char')';
strloc(fid, fs, img_metadata_idx, img_metadata, 'XPos');
finfo.center_x = fread(fid, 1, 'float64');
strloc(fid, fs, img_metadata_idx, img_metadata, 'YPos');
finfo.center_y = fread(fid, 1, 'float64');
strloc(fid, fs, img_metadata_idx, img_metadata, 'dCalibration');
finfo.calib_factor = fread(fid, 1, 'float64');

pos = strfind(img_metadata, insert0('ChannelIsActive'));
finfo.ch_count = length(pos);

% fseek(fid, fs(img_calib_idx).dataStartPos, 'bof');
% img_calib = fread(fid, fs(img_calib_idx).dataLength, '*uint8')';
% fseek(fid, fs(img_calib_idx).dataStartPos + ...
%     strfind(img_calib, insert0('dCalibration')) + ...
%     length(insert0('dCalibration')), 'bof');
% finfo.calib_factor = fread(fid, 1, 'float64');

finfo.padding_bytes = fs(img_data_idx).dataLength - 8 - ...
    finfo.img_width * finfo.img_height *finfo.ch_count *2;

if finfo.padding_bytes == finfo.img_height * 2
  finfo.padding_style = 1;
elseif finfo.padding_bytes == 0
  finfo.padding_style = 2;
else
  finfo.padding_style = 3;
end

if ~isempty(varargin) && strcmpi(varargin{1}, 'getmeta')
  fseek(fid, fs(img_calib_idx).dataStartPos, 'bof');
  finfo.meta.img_calib = fread(fid, fs(img_calib_idx).dataLength, '*char')';
  
  fseek(fid, fs(img_txt_idx).dataStartPos, 'bof');
  finfo.meta.img_txt = fread(fid, fs(img_txt_idx).dataLength, '*char')';
  
  finfo.meta.img_meta = img_metadata;
end

function [attrib, count, flag] = readHeader(fid, attrib, count)
attrib(count).nameLength = fread(fid, 1, 'uint32');
attrib(count).dataLength = fread(fid, 1, 'uint64');
attrib(count).nameAttribute = fread(fid, attrib(count).nameLength, '*char')';
attrib(count).dataStartPos = ftell(fid);
flag = fseek(fid, attrib(count).dataLength, 'cof');
count = count + 1;
end

% ND2 file has text sequence with 0 between characters. 
function out = insert0(in)
num = in + 0;
out = char([reshape(cat(1, num, zeros(size(num))), [1, length(num)*2]), 0]);
end

function strloc( fid, fs, fsidx, text, str )
idx = strfind(text, insert0(str)) + length(insert0(str));
fseek(fid, fs(fsidx).dataStartPos + idx(1), 'bof');
end
end