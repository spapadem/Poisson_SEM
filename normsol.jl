function  normsol(nodes,el,uh,hx,hy,r)


itomn = zeros((r+1)^2,2);

pref = [ -1 -1 0;
          1 -1 0;
          1  1 0;
         -1  1 0];
elref = [1 2 3 4];

glp,glw = GLpw(r);


nu  = zeros(size(el,1),1);

for l = 1 : size(el,1)
glpx = nodes[Int(el[l,1]),1] .+ .5*hx*(glp .+ 1);
glpy = nodes[Int(el[l,1]),2] .+ .5*hy*(glp .+ 1);
indp = 1;
	for p = 1 : r + 1
		for q = 1 : r + 1
			xp = findall(x -> x .< 1e-10, vec(abs.(nodes[Int(el[l,indp]),1] .- glpx)));
		    yq = findall(y -> y .< 1e-10, vec(abs.(nodes[Int(el[l,indp]),2] .- glpy)));
            nu[l]  = nu[l]  + .25*hx*hy*glw[p]*glw[q]*
			(ex_sol(glpx[xp[1]],glpy[yq[1]])).^2;
			indp = indp + 1;
        end
	end
end
nu   = sqrt(sum(nu));

return nu

end

