// this loader is not generated 
// we have added the global=%t option 
// because some other dynamic libraries 
// may want to use symbols defined internally here.

libcplex_path=file('join',['.','libcplex'+%shext]);
addinter(libcplex_path,'libcplex',global=%t);
clear libcplex_path;
