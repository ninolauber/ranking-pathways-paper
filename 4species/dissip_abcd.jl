##----load packages----------------------------------------
using JuMP
using Ipopt
using DelimitedFiles
##---------------------------------------------------------


##----define functions-------------------------------------
"""
"arrhenius" type function to calculate
reaction-rate constant

Eb: reaction energy-barrier
E1: educt energy of formation
E2: product energy of formation
"""
function f(Eb,E1,E2)
	## calculate reaction free-energy
	dG = E2-E1;
	## check if reaction is exergonic (dG<=0)
	if dG <= 0.0
		## activation energy = energy-barrier
		return exp(-Eb)
	## check if reaction is endergonic (dG>0)
	elseif dG > 0.0
		## activation energy = energy-barrier + free-energy
		return exp(-(Eb+dG))
	end
end


function cost_full(kp,km,vext,q0)
	## init model
	model = Model(Ipopt.Optimizer);
	## disable terminal output
	set_silent(model);
	
	## init concentration variables
	@variable(model,1e-2 <= a);
	@variable(model,1e-2 <= b);
	@variable(model,1e-2 <= c);
	@variable(model,1e-2 <= d);
	## init slack variable
	@variable(model,slack==0.0);
	
	## calculate mass-action forward- and backward-fluxes
	j1p,j1m = kp[1]*a,km[1]*b;
	j2p,j2m = kp[2]*a,km[2]*c;
	j3p,j3m = kp[3]*b,km[3]*d;
	j4p,j4m = kp[4]*c,km[4]*d;

	## set flux-constraints
	@constraint(model,-j1p+j1m-j2p+j2m == vext[1]);
	@constraint(model,j1p-j1m-j3p+j3m == vext[2]);
	@constraint(model,j2p-j2m-j4p+j4m == vext[3]);
	@constraint(model,j3p-j3m+j4p-j4m == vext[4]);
	## set conservation-constraints
	@constraint(model,a+b+c+d+slack == q0);

	## calculate dissipation
	Sigma = log(j1m/j1p)*(j1m-j1p) +
		log(j2m/j2p)*(j2m-j2p) +
		log(j3m/j3p)*(j3m-j3p) +
		log(j4m/j4p)*(j4m-j4p);
	Delta = 0.0;
	Chi = Sigma + Delta;

	@objective(model,Min,Chi);
	## DEBUG: print model details
	# print(model)

	## solve optimization problem
	optimize!(model)
	## DEBUG: print termination statusus
	# println(termination_statusus(model))

	costs = [value(Sigma),value(Delta),value(Chi)];

	q = [
		value(a),
		value(b),
		value(c),
		value(d)
	];

	j = zeros(4,3);
	j[:,1] =[ 
		value(j1p),
		value(j2p),
		value(j3p),
		value(j4p)
	]
	j[:,2] =[ 
		value(j1m),
		value(j2m),
		value(j3m),
		value(j4m)
	]
	j[:,3] = j[:,1] .- j[:,2]
	
	## return optimal dissipation and corresponding concentrations
	return termination_status(model),costs,q,j
end


function cost_b(kp,km,vext,q0)
	## init model
	model = Model(Ipopt.Optimizer);
	## disable terminal output
	set_silent(model);
	
	## init concentration variables
	@variable(model,1e-5 <= a);
	@variable(model,1e-5 <= b);
	@variable(model,1e-5 <= c);
	@variable(model,1e-5 <= d);
	
	## calculate mass-action forward- and backward-fluxes
	j1p,j1m = kp[1]*a,km[1]*b;
	j2p,j2m = kp[2]*a,km[2]*c;
	j3p,j3m = kp[3]*b,km[3]*d;
	j4p,j4m = kp[4]*c,km[4]*d;
	
	## set flux-constraints
	@constraint(model,-j1p+j1m == vext[1]);
	@constraint(model,j1p-j1m-j3p+j3m == vext[2]);
	@constraint(model,j3p-j3m == vext[4]);
	## set conservation-constraints
	@constraint(model,a+b+c+d == q0);

	## calculate dissipation
	Sigma = log(j1m/j1p)*(j1m-j1p)+log(j3m/j3p)*(j3m-j3p);
	Delta = (sqrt(j2p)-sqrt(j2m))^2+(sqrt(j4p)-sqrt(j4m))^2;
	Chi = Sigma + Delta
	## set objective function: minimize dissipation
	@objective(model,Min,Chi);
	## DEBUG: print model details
	# print(model)

	## solve optimization problem
	optimize!(model)
	## DEBUG: print termination statusus
	# println(termination_statusus(model))

	costs = [value(Sigma),value(Delta),value(Chi)];

	q = [
		value(a),
		value(b),
		value(c),
		value(d)
	];

	j = zeros(4,3);
	j[:,1] =[ 
		value(j1p),
		value(j2p),
		value(j3p),
		value(j4p)
	]
	j[:,2] =[ 
		value(j1m),
		value(j2m),
		value(j3m),
		value(j4m)
	]
	j[:,3] = j[:,1] .- j[:,2]
	
	## return optimal dissipation and corresponding concentrations
	return termination_status(model),costs,q,j
end


function cost_c(kp,km,vext,q0)
	## init model
	model = Model(Ipopt.Optimizer);
	## disable terminal output
	set_silent(model);
	
	## init concentration variables
	@variable(model,1e-5 <= a);
	@variable(model,1e-5 <= b);
	@variable(model,1e-5 <= c);
	@variable(model,1e-5 <= d);
	
	## calculate mass-action forward- and backward-fluxes
	j1p,j1m = kp[1]*a,km[1]*b;
	j2p,j2m = kp[2]*a,km[2]*c;
	j3p,j3m = kp[3]*b,km[3]*d;
	j4p,j4m = kp[4]*c,km[4]*d;
	
	## set flux-constraints
	@constraint(model,-j2p+j2m == vext[1]);
	@constraint(model,j2p-j2m-j4p+j4m == vext[3]);
	@constraint(model,j4p-j4m == vext[4]);
	## set conservation-constraints
	@constraint(model,a+b+c+d == q0);

	## calculate dissipation
	Sigma = log(j2m/j2p)*(j2m-j2p)+log(j4m/j4p)*(j4m-j4p);
	Delta = (sqrt(j1p)-sqrt(j1m))^2+(sqrt(j3p)-sqrt(j3m))^2;
	Chi = Sigma + Delta
	## set objective function: minimize dissipation
	@objective(model,Min,Chi);
	## DEBUG: print model details
	# print(model)

	## solve optimization problem
	optimize!(model)
	## DEBUG: print termination statusus
	# println(termination_statusus(model))

	costs = [value(Sigma),value(Delta),value(Chi)];

	q = [
		value(a),
		value(b),
		value(c),
		value(d)
	];

	j = zeros(4,3);
	j[:,1] =[ 
		value(j1p),
		value(j2p),
		value(j3p),
		value(j4p)
	]
	j[:,2] =[ 
		value(j1m),
		value(j2m),
		value(j3m),
		value(j4m)
	]
	j[:,3] = j[:,1] .- j[:,2]
	
	## return optimal dissipation and corresponding concentrations
	return termination_status(model),costs,q,j
end


function spec4_costs(eb,ef,fname)
	## set constraints
	q0 = 50.0;
	Vext = [-1.0,0.0,0.0,1.0];
	Vclose = [0.0,0.0,0.0,0.0];

	## calculate the forward reaction-rates
	kp = [
		f(eb[1],ef[1],ef[2]),
		f(eb[2],ef[1],ef[3]),
		f(eb[3],ef[2],ef[4]),
		f(eb[4],ef[3],ef[4])
	];
	## calculate the backward reaction-rates
	km = [
		f(eb[1],ef[2],ef[1]),
		f(eb[2],ef[3],ef[1]),
		f(eb[3],ef[4],ef[2]),
		f(eb[4],ef[4],ef[3])
	];

	println("------------------------------------------")
	println("Full CRN Detailed Balance Eq.")
	println("------------------------------------------")
	status,costs_eq,q_eq,Jeq = cost_full(kp,km,Vclose,q0);

	println("Optimization: ",status)

	println("------------------------------------------")
	println("Full CRN NESS")
	println("------------------------------------------")
	status,costs_full,q_full,Jfull = cost_full(kp,km,Vext,q0);
	println("Optimization: ",status)

	println("------------------------------------------")
	println("ABD-Pathway")
	println("------------------------------------------")
	status,costs_b,q_b,Jb = cost_b(kp,km,Vext,q0);
	println("Optimization: ",status)

	println("------------------------------------------")
	println("ACD-Pathway")
	println("------------------------------------------")
	status,costs_c,q_c,Jc = cost_c(kp,km,Vext,q0);
	println("Optimization: ",status)

	k = zeros(2,4);
	k[1,:] = round.(kp,digits=4);
	k[2,:] = round.(km,digits=4);

	q = zeros(4,4);
	q[1,:] = round.(q_eq,digits=4);
	q[2,:] = round.(q_full,digits=4);
	q[3,:] = round.(q_b,digits=4);
	q[4,:] = round.(q_c,digits=4);

	Jm = zeros(3,4);
	Jm[1,:] = round.(Jfull[:,1]',digits=4);
	Jm[2,:] = round.(Jb[:,1]',digits=4);
	Jm[3,:] = round.(Jc[:,1]',digits=4);

	Jp = zeros(3,4);
	Jp[1,:] = round.(Jfull[:,2]',digits=4);
	Jp[2,:] = round.(Jb[:,2]',digits=4);
	Jp[3,:] = round.(Jc[:,2]',digits=4);

	J = zeros(3,4);
	J[1,:] = round.(Jfull[:,3]',digits=4);
	J[2,:] = round.(Jb[:,3]',digits=4);
	J[3,:] = round.(Jc[:,3]',digits=4);

	costs = zeros(3,3);
	costs[1,:] = round.(costs_full,digits=4);
	costs[2,:] = round.(costs_b,digits=4);
	costs[3,:] = round.(costs_c,digits=4);

	if !isdir("data")
		mkdir("data")
	end
	open("data/"*fname*"_Jm.dat","w") do f
		write(f,"r1,r2,r3,r4\n");
		writedlm(f,Jm,',');
	end
	open("data/"*fname*"_Jp.dat","w") do f
		write(f,"r1,r2,r3,r4\n");
		writedlm(f,Jp,',');
	end
	open("data/"*fname*"_J.dat","w") do f
		write(f,"r1,r2,r3,r4\n");
		writedlm(f,J,',');
	end
	open("data/"*fname*"_q.dat","w") do f
		write(f,"a,b,c,d\n");
		writedlm(f,q,',');
	end
	open("data/"*fname*"_k.dat","w") do f
		write(f,"r1,r2,r3,r4\n");
		writedlm(f,k,',');
	end
	open("data/"*fname*"_costs.dat","w") do f
		write(f,"Sigma,Delta,Chi\n");
		writedlm(f,costs,',');
	end

	return 0
end
##---------------------------------------------------------


##----wrapper functions for results------------------------
function sym_results()
	## set energy landscape
	eb = [2.0,2.0,2.0,2.0];
	ef = [-2.0,-4.0,-4.0,-6.0];
	spec4_costs(eb,ef,"sym");
	return 0
end

function asym11_results()
	## set energy landscape
	eb = [2.0,2.0,2.0,2.0];
	ef = [-2.0,-3.0,-4.0,-6.0];
	spec4_costs(eb,ef,"asym11");
	return 0
end

function asym12_results()
	## set energy landscape
	eb = [0.0,2.0,0.0,2.0];
	ef = [-2.0,-3.0,-4.0,-6.0];
	spec4_costs(eb,ef,"asym12");
	return 0
end

function asym21_results()
	## set energy landscape
	eb = [2.0,2.0,2.0,2.0];
	ef = [-2.0,-1.0,-4.0,-6.0];
	spec4_costs(eb,ef,"asym21");
	return 0
end

function asym22_results()
	## set energy landscape
	eb = [0.0,2.0,0.0,2.0];
	ef = [-2.0,-1.0,-4.0,-6.0];
	spec4_costs(eb,ef,"asym22");
	return 0
end
##---------------------------------------------------------
