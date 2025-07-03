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
	@variable(model,1e-2 <= e);
	## init slack variable
	@variable(model,slack==0.0);
	
	## calculate mass-action forward- and backward-fluxes
	j1p,j1m = kp[1]*a,km[1]*b;
	j2p,j2m = kp[2]*a,km[2]*c;
	j3p,j3m = kp[3]*a,km[3]*d;
	j4p,j4m = kp[4]*b,km[4]*e;
	j5p,j5m = kp[5]*c,km[5]*e;
	j6p,j6m = kp[6]*d,km[6]*e;

	## set flux-constraints
	@constraint(model,-j1p+j1m-j2p+j2m-j3p+j3m == vext[1]);
	@constraint(model,j1p-j1m-j4p+j4m == vext[2]);
	@constraint(model,j2p-j2m-j5p+j5m == vext[3]);
	@constraint(model,j3p-j3m-j6p+j6m == vext[4]);
	@constraint(model,j4p-j4m+j5p-j5m+j6p-j6m == vext[5]);
	## set conservation-constraints
	@constraint(model,a+b+c+d+e+slack == q0);

	## calculate dissipation
	Sigma = log(j1m/j1p)*(j1m-j1p) +
		log(j2m/j2p)*(j2m-j2p) +
		log(j3m/j3p)*(j3m-j3p) +
		log(j4m/j4p)*(j4m-j4p) +
		log(j5m/j5p)*(j5m-j5p) +
		log(j6m/j6p)*(j6m-j6p);
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
		value(d),
		value(e)
	];

	## return optimal dissipation and corresponding concentrations
	return termination_status(model),costs,q
end


function cost_bc(kp,km,vext,q0)
	## init model
	model = Model(Ipopt.Optimizer);
	## disable terminal output
	set_silent(model);
	
	## init concentration variables
	@variable(model,1e-2 <= a);
	@variable(model,1e-2 <= b);
	@variable(model,1e-2 <= c);
	@variable(model,1e-2 <= d);
	@variable(model,1e-2 <= e);
	
	## calculate mass-action forward- and backward-fluxes
	j1p,j1m = kp[1]*a,km[1]*b;
	j2p,j2m = kp[2]*a,km[2]*c;
	j3p,j3m = kp[3]*a,km[3]*d;
	j4p,j4m = kp[4]*b,km[4]*e;
	j5p,j5m = kp[5]*c,km[5]*e;
	j6p,j6m = kp[6]*d,km[6]*e;

	## set flux-constraints
	@constraint(model,-j1p+j1m-j2p+j2m == vext[1]);
	@constraint(model,j1p-j1m-j4p+j4m == vext[2]);
	@constraint(model,j2p-j2m-j5p+j5m == vext[3]);
	@constraint(model,j4p-j4m+j5p-j5m == vext[5]);
	## set conservation-constraints
	@constraint(model,a+b+c+d+e == q0);

	## calculate dissipation
	Sigma = log(j1m/j1p)*(j1m-j1p)+
		log(j2m/j2p)*(j2m-j2p) +
		log(j4m/j4p)*(j4m-j4p) +
		log(j5m/j5p)*(j5m-j5p);
	Delta = (sqrt(j3p)-sqrt(j3m))^2+(sqrt(j6p)-sqrt(j6m))^2;
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
		value(d),
		value(e)
	];
	
	## return optimal dissipation and corresponding concentrations
	return termination_status(model),costs,q
end


function cost_bd(kp,km,vext,q0)
	## init model
	model = Model(Ipopt.Optimizer);
	## disable terminal output
	set_silent(model);
	
	## init concentration variables
	@variable(model,1e-2 <= a);
	@variable(model,1e-2 <= b);
	@variable(model,1e-2 <= c);
	@variable(model,1e-2 <= d);
	@variable(model,1e-2 <= e);
	
	## calculate mass-action forward- and backward-fluxes
	j1p,j1m = kp[1]*a,km[1]*b;
	j2p,j2m = kp[2]*a,km[2]*c;
	j3p,j3m = kp[3]*a,km[3]*d;
	j4p,j4m = kp[4]*b,km[4]*e;
	j5p,j5m = kp[5]*c,km[5]*e;
	j6p,j6m = kp[6]*d,km[6]*e;

	## set flux-constraints
	@constraint(model,-j1p+j1m-j3p+j3m == vext[1]);
	@constraint(model,j1p-j1m-j4p+j4m == vext[2]);
	@constraint(model,j3p-j3m-j6p+j6m == vext[4]);
	@constraint(model,j4p-j4m+j6p-j6m == vext[5]);
	## set conservation-constraints
	@constraint(model,a+b+c+d+e == q0);

	## calculate dissipation
	Sigma = log(j1m/j1p)*(j1m-j1p)+
		log(j3m/j3p)*(j3m-j3p) +
		log(j4m/j4p)*(j4m-j4p) +
		log(j6m/j6p)*(j6m-j6p);
	Delta = (sqrt(j2p)-sqrt(j2m))^2+(sqrt(j5p)-sqrt(j5m))^2;
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
		value(d),
		value(e)
	];
	
	## return optimal dissipation and corresponding concentrations
	return termination_status(model),costs,q
end


function cost_cd(kp,km,vext,q0)
	## init model
	model = Model(Ipopt.Optimizer);
	## disable terminal output
	set_silent(model);
	
	## init concentration variables
	@variable(model,1e-2 <= a);
	@variable(model,1e-2 <= b);
	@variable(model,1e-2 <= c);
	@variable(model,1e-2 <= d);
	@variable(model,1e-2 <= e);
	
	## calculate mass-action forward- and backward-fluxes
	j1p,j1m = kp[1]*a,km[1]*b;
	j2p,j2m = kp[2]*a,km[2]*c;
	j3p,j3m = kp[3]*a,km[3]*d;
	j4p,j4m = kp[4]*b,km[4]*e;
	j5p,j5m = kp[5]*c,km[5]*e;
	j6p,j6m = kp[6]*d,km[6]*e;

	## set flux-constraints
	@constraint(model,-j2p+j2m-j3p+j3m == vext[1]);
	@constraint(model,j2p-j2m-j5p+j5m == vext[3]);
	@constraint(model,j3p-j3m-j6p+j6m == vext[4]);
	@constraint(model,j5p-j5m+j6p-j6m == vext[5]);
	## set conservation-constraints
	@constraint(model,a+b+c+d+e == q0);

	## calculate dissipation
	Sigma = log(j2m/j2p)*(j2m-j2p)+
		log(j3m/j3p)*(j3m-j3p) +
		log(j5m/j5p)*(j5m-j5p) +
		log(j6m/j6p)*(j6m-j6p);
	Delta = (sqrt(j1p)-sqrt(j1m))^2+(sqrt(j4p)-sqrt(j4m))^2;
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
		value(d),
		value(e)
	];
	
	## return optimal dissipation and corresponding concentrations
	return termination_status(model),costs,q
end


function cost_b(kp,km,vext,q0)
	## init model
	model = Model(Ipopt.Optimizer);
	## disable terminal output
	set_silent(model);
	
	## init concentration variables
	@variable(model,1e-2 <= a);
	@variable(model,1e-2 <= b);
	@variable(model,1e-2 <= c);
	@variable(model,1e-2 <= d);
	@variable(model,1e-2 <= e);
	
	## calculate mass-action forward- and backward-fluxes
	j1p,j1m = kp[1]*a,km[1]*b;
	j2p,j2m = kp[2]*a,km[2]*c;
	j3p,j3m = kp[3]*a,km[3]*d;
	j4p,j4m = kp[4]*b,km[4]*e;
	j5p,j5m = kp[5]*c,km[5]*e;
	j6p,j6m = kp[6]*d,km[6]*e;

	## set flux-constraints
	@constraint(model,-j1p+j1m == vext[1]);
	@constraint(model,j1p-j1m-j4p+j4m == vext[2]);
	@constraint(model,j4p-j4m == vext[5]);
	## set conservation-constraints
	@constraint(model,a+b+c+d+e == q0);

	## calculate dissipation
	Sigma = log(j1m/j1p)*(j1m-j1p)+log(j4m/j4p)*(j4m-j4p);
	Delta = (sqrt(j2p)-sqrt(j2m))^2 +
		(sqrt(j3p)-sqrt(j3m))^2 +
		(sqrt(j5p)-sqrt(j5m))^2 +
		(sqrt(j6p)-sqrt(j6m))^2;
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
		value(d),
		value(e)
	];
	
	## return optimal dissipation and corresponding concentrations
	return termination_status(model),costs,q
end


function cost_c(kp,km,vext,q0)
	## init model
	model = Model(Ipopt.Optimizer);
	## disable terminal output
	set_silent(model);
	
	## init concentration variables
	@variable(model,1e-2 <= a);
	@variable(model,1e-2 <= b);
	@variable(model,1e-2 <= c);
	@variable(model,1e-2 <= d);
	@variable(model,1e-2 <= e);
	
	## calculate mass-action forward- and backward-fluxes
	j1p,j1m = kp[1]*a,km[1]*b;
	j2p,j2m = kp[2]*a,km[2]*c;
	j3p,j3m = kp[3]*a,km[3]*d;
	j4p,j4m = kp[4]*b,km[4]*e;
	j5p,j5m = kp[5]*c,km[5]*e;
	j6p,j6m = kp[6]*d,km[6]*e;

	## set flux-constraints
	@constraint(model,-j2p+j2m == vext[1]);
	@constraint(model,j2p-j2m-j5p+j5m == vext[2]);
	@constraint(model,j5p-j5m == vext[5]);
	## set conservation-constraints
	@constraint(model,a+b+c+d+e == q0);

	## calculate dissipation
	Sigma = log(j2m/j2p)*(j2m-j2p)+log(j5m/j5p)*(j5m-j5p);
	Delta = (sqrt(j1p)-sqrt(j1m))^2 +
		(sqrt(j3p)-sqrt(j3m))^2 +
		(sqrt(j4p)-sqrt(j4m))^2 +
		(sqrt(j6p)-sqrt(j6m))^2;
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
		value(d),
		value(e)
	];
	
	## return optimal dissipation and corresponding concentrations
	return termination_status(model),costs,q
end


function cost_d(kp,km,vext,q0)
	## init model
	model = Model(Ipopt.Optimizer);
	## disable terminal output
	set_silent(model);
	
	## init concentration variables
	@variable(model,1e-2 <= a);
	@variable(model,1e-2 <= b);
	@variable(model,1e-2 <= c);
	@variable(model,1e-2 <= d);
	@variable(model,1e-2 <= e);
	
	## calculate mass-action forward- and backward-fluxes
	j1p,j1m = kp[1]*a,km[1]*b;
	j2p,j2m = kp[2]*a,km[2]*c;
	j3p,j3m = kp[3]*a,km[3]*d;
	j4p,j4m = kp[4]*b,km[4]*e;
	j5p,j5m = kp[5]*c,km[5]*e;
	j6p,j6m = kp[6]*d,km[6]*e;

	## set flux-constraints
	@constraint(model,-j3p+j3m == vext[1]);
	@constraint(model,j3p-j3m-j6p+j6m == vext[2]);
	@constraint(model,j6p-j6m == vext[5]);
	## set conservation-constraints
	@constraint(model,a+b+c+d+e == q0);

	## calculate dissipation
	Sigma = log(j3m/j3p)*(j3m-j3p)+log(j6m/j6p)*(j6m-j6p);
	Delta = (sqrt(j1p)-sqrt(j1m))^2 +
		(sqrt(j2p)-sqrt(j2m))^2 +
		(sqrt(j4p)-sqrt(j4m))^2 +
		(sqrt(j5p)-sqrt(j5m))^2;
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
		value(d),
		value(e)
	];
	
	## return optimal dissipation and corresponding concentrations
	return termination_status(model),costs,q
end


function spec5_costs(eb,ef,fname)
	## set constraints
	q0 = 20.1;
	Vext = [-1.0,0.0,0.0,0.0,1.0];
	Vclose = [0.0,0.0,0.0,0.0,0.0];

	## calculate the forward reaction-rates
	kp = [
		f(eb[1],ef[1],ef[2]),
		f(eb[2],ef[1],ef[3]),
		f(eb[3],ef[1],ef[4]),
		f(eb[4],ef[2],ef[5]),
		f(eb[5],ef[3],ef[5]),
		f(eb[6],ef[4],ef[5])
	];
	## calculate the backward reaction-rates
	km = [
		f(eb[1],ef[2],ef[1]),
		f(eb[2],ef[3],ef[1]),
		f(eb[3],ef[4],ef[1]),
		f(eb[4],ef[5],ef[2]),
		f(eb[5],ef[5],ef[3]),
		f(eb[6],ef[5],ef[4])
	];

	println("------------------------------------------")
	println("Full CRN Detailed Balance Eq.")
	println("------------------------------------------")
	status,costs_eq,q_eq = cost_full(kp,km,Vclose,q0);
	println("Optimization: ",status)

	println("------------------------------------------")
	println("Full CRN NESS")
	println("------------------------------------------")
	status,costs_full,q_full = cost_full(kp,km,Vext,q0);
	println("Optimization: ",status)

	println("------------------------------------------")
	println("ABCE-Pathway")
	println("------------------------------------------")
	status,costs_bc,q_bc = cost_bc(kp,km,Vext,q0);
	println("Optimization: ",status)

	println("------------------------------------------")
	println("ABDE-Pathway")
	println("------------------------------------------")
	status,costs_bd,q_bd = cost_bd(kp,km,Vext,q0);
	println("Optimization: ",status)

	println("------------------------------------------")
	println("ACDE-Pathway")
	println("------------------------------------------")
	status,costs_cd,q_cd = cost_cd(kp,km,Vext,q0);
	println("Optimization: ",status)

	println("------------------------------------------")
	println("ABE-Pathway")
	println("------------------------------------------")
	status,costs_b,q_b = cost_b(kp,km,Vext,q0);
	println("Optimization: ",status)

	println("------------------------------------------")
	println("ACE-Pathway")
	println("------------------------------------------")
	status,costs_c,q_c = cost_c(kp,km,Vext,q0);
	println("Optimization: ",status)

	println("------------------------------------------")
	println("ADE-Pathway")
	println("------------------------------------------")
	status,costs_d,q_d = cost_d(kp,km,Vext,q0);
	println("Optimization: ",status)

	q = zeros(8,5);
	q[1,:] = round.(q_eq,digits=4);
	q[2,:] = round.(q_full,digits=4);
	q[3,:] = round.(q_bc,digits=4);
	q[4,:] = round.(q_bd,digits=4);
	q[5,:] = round.(q_cd,digits=4);
	q[6,:] = round.(q_b,digits=4);
	q[7,:] = round.(q_c,digits=4);
	q[8,:] = round.(q_d,digits=4);

	costs = zeros(7,3);
	costs[1,:] = round.(costs_full,digits=4);
	costs[2,:] = round.(costs_bc,digits=4);
	costs[3,:] = round.(costs_bd,digits=4);
	costs[4,:] = round.(costs_cd,digits=4);
	costs[5,:] = round.(costs_b,digits=4);
	costs[6,:] = round.(costs_c,digits=4);
	costs[7,:] = round.(costs_d,digits=4);

	if !isdir("data")
		mkdir("data")
	end
	# writedlm("data/sym_k.dat",k,',');
	open("data/"*fname*"_q.dat","w") do f
		write(f,"a,b,c,d,e\n");
		writedlm(f,q,',');
	end
	open("data/"*fname*"_costs.dat","w") do f
		write(f,"Sigma,Delta,Chi\n");
		writedlm(f,costs,',');
	end

	return 0
end
##---------------------------------------------------------


##----wrapper functions for results------------------------
function asym_results()
	## set energy landscape
	eb = [2.0,2.0,2.0,2.0,2.0,2.0];
	ef = [-2.0,-3.0,-4.0,-5.0,-6.0];
	spec5_costs(eb,ef,"asym");
	return 0
end
##---------------------------------------------------------
