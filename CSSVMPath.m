function out=CSSVMPath(x,y)
global NumSingular
    nu=initial_nu(y);
    tic
    NumSingular=0;
    path=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%core code
        initial_solution=initial(x,y,nu);
        [local_path]=regularization_path(initial_solution,x,y);
        path=[path;local_path];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc
%     time=toc ;
    steps=length(path);
    out.Steps=steps;
    out.Path=path;
end