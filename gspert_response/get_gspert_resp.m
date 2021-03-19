shot = 204660;
time = 0.100;
tree = 'EFIT01';
tokamak = 'nstxu';
server = 'skylark.pppl.gov:8501';

tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
mpc = tok_data_struct.mpc;
mpv = tok_data_struct.mpv;

circ = nstxu2016_circ(tok_data_struct);

eq = read_eq(shot,time,tree,tokamak,server);

response = gspert(eq, tok_data_struct);

xmats = response.xmats;
xmatp = [mpc'; mpv'] * reshape(response.dcphidip,[], 1);
Ldum = 0;
xmatx = [xmats xmatp; xmatp' Ldum];
xmat = circ.Pxx' * xmatx * circ.Pxx;
xmat_cv = xmat(1:end-1, :);

xi = reshape(xmat_cv, [], 1);

