function G = gfunwithbulge(lsf,x,probdata,gfundata,G)

nx = size(x,2);

m_minus_1 = gfundata(lsf).m_minus_1;
si = gfundata(lsf).si;
ri = gfundata(lsf).ri;
dpstui = gfundata(lsf).dsptui;

u = x_to_u(x,probdata);

Bi = zeros(m_minus_1,nx);

for i = 1:m_minus_1
    
   disti = ( dot( u-dpstui(:,i)*ones(1,nx) , u-dpstui(:,i)*ones(1,nx) ) ).^0.5;
   
   jbulge = find( disti <= ri(i) );
   
   if ~isempty(jbulge)
      Bi(i,jbulge) = si(i) * ( ri(i)^2 - disti(jbulge).^2 ).^2;
   end
   
end

G = G + sum( Bi, 1 );