% make_mdsipmex.m
% make sure ../mdsplus points to current version of mdsplus distribution
% to make the debug version, debug=1;make_mdsipmex;
% to make the shared version, shared=1;make_mdsipmex;
% Note; if shared, install MdsIpShr library or point $LD_LIBRARY_PATH to shareables
% Basil P. DUVAL, May 2000

MDSPLUS=getenv('MDSPLUS_DIR');
if(length(MDSPLUS)==0)
  disp('shell variable MDSPLUS must point to /f/mdsplus/linux/mdsplus before compilation');end

if(exist('debug','var'));DEBUG = '-DDEBUG';else;DEBUG='';end

if(~strcmp(computer,'VMS'));PASSWD = '-DPASSWD';else;PASSWD='mdsipmex.opt';end

if(~exist('shared','var'))
comm = sprintf('mex -v %s %s mdsipmex.c %s/mdstcpip/mdsipshr.c %s/mdstcpip/mdsiputil.c %s/mdstcpip/mdsip_socket_io.o -DNOCOMPRESSION -I%s/include -I%s/mdstcpip',...
	       DEBUG,PASSWD,MDSPLUS,MDSPLUS,MDSPLUS,MDSPLUS);
else
comm = sprintf('mex -v %s %s mdsipmex.c -DNOCOMPRESSION -I%s/include -I%s/mdstcpip -L%s/lib -lMdsIpShr',...
	       DEBUG,PASSWD,MDSPLUS,MDSPLUS,MDSPLUS);
end
disp(comm);
comm
eval(comm);
