function  Fvec(r,hx,hy,t,c,node,osc_f,xl,yb)

#F = zeros((r+1)^2,1);
glp,glw = GLpw(r);




glpx = xl .+ .5*hx*(glp .+ 1);
glpy = yb .+ .5*hy*(glp .+ 1);

	itomnx = findall(abs.(node[1] .- glpx) .<1e-10);
	itomny = findall(abs.(node[2] .- glpy) .<1e-10);
F = 0;
	for p = 1 : r + 1
		for q = 1 : r + 1
			F = F + .5*hx*hy*glw[p]*glw[q]*
				f_ex(osc_f,t,c,glpx[p],glpy[q])*
				(itomny[1] == q)*(itomnx[1] == p);
		end

	end


return F
end



