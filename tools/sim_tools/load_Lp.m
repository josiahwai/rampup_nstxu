modeldir = '/Users/jwai/Research/rampup_nstxu/buildmodel/built_models/std/';

load('nstxu_obj_config2016_6565.mat')
circ = nstxu2016_circ(tok_data_struct);

t = [30:10:960]/1e3;
shot = 204660;

Lp = zeros(length(t), 1);
for i = 1:length(t)
  time_ms = t(i)*1e3;
  sys = load([modeldir num2str(shot) '_' num2str(time_ms) '_sys.mat']).sys;      
  Lp(i) = sys.Lp;
  Mpc(i,:) = sys.lstar(circ.iipx, circ.iicx);
  Mpv(i,:) = sys.lstar(circ.iipx, circ.iivx);
end


Mpc = smoothdata(Mpc,1, 'movmedian', 11);
Mpv = smoothdata(Mpv,1, 'movmedian', 11);
Mpc = smoothdata(Mpc,1, 'movmean', 11);
Mpv = smoothdata(Mpv,1, 'movmean', 11);

figure
plot(t,Mpc)

figure
plot(t,Mpv)

% cftool(t,Lp)













