function [A, B, C, D] = SystemDynamics(Mxx, Rxx,fit_coils,...
    ncx, Rext_mOhm, Lext_mH, Ts,vargin)

nx = size(Rxx,1);
nvx = nx-ncx;

for i=1:length(fit_coils)
    Mxx(fit_coils(i),fit_coils(i)) = Mxx(fit_coils(i),fit_coils(i)) + Lext_mH(i)/1000;
    Rxx(fit_coils(i),fit_coils(i)) = Rxx(fit_coils(i),fit_coils(i)) + Rext_mOhm(i)/1000;
end

Ac = -inv(Mxx)*Rxx;

coil_mask = zeros(ncx,length(fit_coils));
for i=1:length(fit_coils)
    coil_mask(fit_coils(i),i) = 1;
end

Bc = Mxx\[eye(ncx); zeros(nvx,ncx)]*coil_mask;

nx = length(Rxx);

A = inv(eye(nx)-Ts*Ac);
%A = 

B = A*Bc*Ts;

C = coil_mask'*[eye(ncx) zeros(ncx,nvx)];

D = zeros(size(C,1),size(B,2));

end
