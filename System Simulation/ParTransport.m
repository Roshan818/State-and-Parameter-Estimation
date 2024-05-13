%% A function for parallel transporting a vector from point xold to xnew 
function PAR = ParTransport(xold,xnew)
PAR = (eye(3) - ((eye(3)-xold*xold')*xnew*xnew') ...
      /(1+xold'*xnew) - xold*xnew');
end