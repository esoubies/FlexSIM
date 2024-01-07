function im = nd2read(filename, varargin)
%
%
if nargin==2
    ch=varargin{1};
else
    ch=[];
end

% Get info
finfo = nd2finfo(filename);

% Create image


if isempty(ch)
    im =  zeros([finfo.img_height, finfo.img_width, finfo.ch_count, finfo.img_seq_count]);
    ch=1: finfo.ch_count;
else
    im =  zeros([finfo.img_height, finfo.img_width, 1, finfo.img_seq_count]);
end


% Open file
tic
for tt=1:finfo.img_seq_count
    fid = fopen(filename, 'r');
    fseek(fid, finfo.file_structure(strncmp(['ImageDataSeq|',num2str(tt-1),'!'], ...
        {finfo.file_structure(:).nameAttribute}, 15)).dataStartPos+8, 'bof');

    % Image extracted from ND2 has image width defined by its first dimension.
    if finfo.padding_style == 1
        for ii = 1: finfo.img_height
            temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
                [finfo.ch_count finfo.img_width]);
            jid=1;
            for jj=ch
                im(ii,:,jid,tt) = temp(jj,:);
                jid=jid+1;
            end
            fseek(fid, 2, 'cof');
        end
    else
        for ii = 1: finfo.img_height
            temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
                [finfo.ch_count finfo.img_width]);
            jid=1;
            for jj=ch
                im(ii,:,jid,tt) = temp(jj,:);
                jid=jid+1;
            end
        end
    end
    fclose(fid);
end
im=squeeze(im);

end