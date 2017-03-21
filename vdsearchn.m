function I=vdsearchn(fullvec,subvec)
%vdsearchn(fullvec,subvec)
for i=1:length(subvec)
I(i)=dsearchn(fullvec',subvec(i));
end