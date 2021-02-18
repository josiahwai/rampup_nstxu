sys = load('nstxu_sys.mat').nstxu_sys;
struct_to_ws(sys);

mxx = Pxx'*mxx*Pxx;

mcc = mxx(1:ncx,1:ncx);
mcv = mxx(1:ncx,ncx+1:ncx+nvx);

figure
contourf([mcc mcv]); colorbar



figure
contourf(mcv); colorbar


mean(mcc(:))
mean(mcv(:))

median(mcc(:))
median(mcv(:))

max(mcc(:))
max(mcv(:))











