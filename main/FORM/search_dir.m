function  d = search_dir(G,grad_G,u)


alpha = -grad_G / norm(grad_G);

d = ( G / norm(grad_G) + alpha' * u ) * alpha - u;
