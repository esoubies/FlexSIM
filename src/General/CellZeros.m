function c = CellZeros(c_ref,fac,dims)
[I,J]=size(c_ref);
c=cell([I,J]);
for i=1:I
    for j=1:J
        sz=size(c_ref{i,j});
        c{i,j}=zeros(sz(dims).*fac);
    end
end
end