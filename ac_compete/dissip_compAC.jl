##----load packages----------------------------------------
using JuMP
using Ipopt
using DelimitedFiles
##---------------------------------------------------------


##----define functions-------------------------------------
function cost_full(kp,km,vext,q0,init)
	## init model
	model = Model(Ipopt.Optimizer);
	## disable terminal output
	set_silent(model);
	
	## init concentration variables
	@variable(model,0.0 <= f,start=init[1]);
	@variable(model,0.0 <= a,start=init[2]);
	@variable(model,0.0 <= b,start=init[3]);
	@variable(model,0.0 <= c1,start=init[4]);
	@variable(model,0.0 <= c2,start=init[5]);
	@variable(model,0.0 <= d,start=init[6]);
	## init slack variable
	@variable(model,slack==0.0);
	
	## calculate mass-action forward- and backward-fluxes
	j1p,j1m = kp[1]*a*f,km[1]*b;
	j2p,j2m = kp[2]*b,km[2]*c1;
	j3p,j3m = kp[3]*b,km[3]*c2;
	j4p,j4m = kp[4]*c1*f,km[4]*d;
	j5p,j5m = kp[5]*c2*f,km[5]*d;
	j6p,j6m = kp[6]*d,km[6]*a*a;

	## set flux-constraints
	@constraint(model,-j1p+j1m-j4p+j4m-j5p+j5m == vext[1]);
	@constraint(model,-j1p+j1m+2*j6p-2*j6m == vext[2]);
	@constraint(model,j1p-j1m-j2p+j2m-j3p+j3m == vext[3]);
	@constraint(model,j2p-j2m-j4p+j4m == vext[4]);
	@constraint(model,j3p-j3m-j5p+j5m == vext[5]);
	@constraint(model,j4p-j4m+j5p-j5m-j6p+j6m == vext[6]);
	## set conservation-constraints
	@constraint(model,f+2*a+3*b+3*c1+3*c2+4*d+slack == q0);

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
		value(f),
		value(a),
		value(b),
		value(c1),
		value(c2),
		value(d)
	];

	j = zeros(6,3);
	j[:,1] =[ 
		value(j1p),
		value(j2p),
		value(j3p),
		value(j4p),
		value(j5p),
		value(j6p)
	]
	j[:,2] =[ 
		value(j1m),
		value(j2m),
		value(j3m),
		value(j4m),
		value(j5m),
		value(j6m)
	]
	j[:,3] = j[:,1] .- j[:,2]
	
	## return optimal dissipation and corresponding concentrations
	return termination_status(model),costs,q,j
end


function cost_c1(kp,km,vext,q0,init)
	## init model
	model = Model(Ipopt.Optimizer);
	## disable terminal output
	set_silent(model);
	
	## init concentration variables
	@variable(model,0.0 <= f,start=init[1]);
	@variable(model,0.0 <= a,start=init[2]);
	@variable(model,0.0 <= b,start=init[3]);
	@variable(model,0.0 <= c1,start=init[4]);
	@variable(model,0.0 <= c2,start=init[5]);
	@variable(model,0.0 <= d,start=init[6]);

	## calculate mass-action forward- and backward-fluxes
	j1p,j1m = kp[1]*a*f,km[1]*b;
	j2p,j2m = kp[2]*b,km[2]*c1;
	j3p,j3m = kp[3]*b,km[3]*c2;
	j4p,j4m = kp[4]*c1*f,km[4]*d;
	j5p,j5m = kp[5]*c2*f,km[5]*d;
	j6p,j6m = kp[6]*d,km[6]*a*a;

	## set flux-constraints
	@constraint(model,-j1p+j1m-j4p+j4m == vext[1]);
	@constraint(model,-j1p+j1m+2*j6p-2*j6m == vext[2]);
	@constraint(model,j1p-j1m-j2p+j2m == vext[3]);
	@constraint(model,j2p-j2m-j4p+j4m == vext[4]);
	@constraint(model,j4p-j4m-j6p+j6m == vext[6]);
	## set conservation-constraints
	@constraint(model,f+2*a+3*b+3*c1+3*c2+4*d == q0);

	## calculate dissipation
	Sigma = log(j1m/j1p)*(j1m-j1p) +
		log(j2m/j2p)*(j2m-j2p) +
		log(j4m/j4p)*(j4m-j4p) +
		log(j6m/j6p)*(j6m-j6p);
	Delta = (sqrt(j3p)-sqrt(j3m))^2+(sqrt(j5p)-sqrt(j5m))^2;
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
		value(f),
		value(a),
		value(b),
		value(c1),
		value(c2),
		value(d)
	];

	j = zeros(6,3);
	j[:,1] =[ 
		value(j1p),
		value(j2p),
		value(j3p),
		value(j4p),
		value(j5p),
		value(j6p)
	]
	j[:,2] =[ 
		value(j1m),
		value(j2m),
		value(j3m),
		value(j4m),
		value(j5m),
		value(j6m)
	]
	j[:,3] = j[:,1] .- j[:,2]
	
	## return optimal dissipation and corresponding concentrations
	return termination_status(model),costs,q,j
end


function cost_c2(kp,km,vext,q0,init)
	## init model
	model = Model(Ipopt.Optimizer);
	## disable terminal output
	set_silent(model);
	
	## init concentration variables
	@variable(model,0.0 <= f,start=init[1]);
	@variable(model,0.0 <= a,start=init[2]);
	@variable(model,0.0 <= b,start=init[3]);
	@variable(model,0.0 <= c1,start=init[4]);
	@variable(model,0.0 <= c2,start=init[5]);
	@variable(model,0.0 <= d,start=init[6]);
	
	## calculate mass-action forward- and backward-fluxes
	j1p,j1m = kp[1]*a*f,km[1]*b;
	j2p,j2m = kp[2]*b,km[2]*c1;
	j3p,j3m = kp[3]*b,km[3]*c2;
	j4p,j4m = kp[4]*c1*f,km[4]*d;
	j5p,j5m = kp[5]*c2*f,km[5]*d;
	j6p,j6m = kp[6]*d,km[6]*a*a;

	## set flux-constraints
	@constraint(model,-j1p+j1m-j5p+j5m == vext[1]);
	@constraint(model,-j1p+j1m+2*j6p-2*j6m == vext[2]);
	@constraint(model,j1p-j1m-j3p+j3m == vext[3]);
	@constraint(model,j3p-j3m-j5p+j5m == vext[5]);
	@constraint(model,j5p-j5m-j6p+j6m == vext[6]);
	## set conservation-constraints
	@constraint(model,f+2*a+3*b+3*c1+3*c2+4*d == q0);

	## calculate dissipation
	Sigma = log(j1m/j1p)*(j1m-j1p) +
		log(j3m/j3p)*(j3m-j3p) +
		log(j5m/j5p)*(j5m-j5p) +
		log(j6m/j6p)*(j6m-j6p);
	Delta = (sqrt(j2p)-sqrt(j2m))^2+(sqrt(j4p)-sqrt(j4m))^2;
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
		value(f),
		value(a),
		value(b),
		value(c1),
		value(c2),
		value(d)
	];

	j = zeros(6,3);
	j[:,1] =[ 
		value(j1p),
		value(j2p),
		value(j3p),
		value(j4p),
		value(j5p),
		value(j6p)
	]
	j[:,2] =[ 
		value(j1m),
		value(j2m),
		value(j3m),
		value(j4m),
		value(j5m),
		value(j6m)
	]
	j[:,3] = j[:,1] .- j[:,2]
	
	## return optimal dissipation and corresponding concentrations
	return termination_status(model),costs,q,j
end


function compAC_costs(init,fname)
	## set constraints
	q0 = 20.0;
	Vext = [-2.0,1.0,0.0,0.0,0.0,0.0];
	Vclose = [0.0,0.0,0.0,0.0,0.0,0.0];

	## set forward- and backward reaction-rates
	kp = ones(6);
	km = ones(6);

	println("------------------------------------------")
	println("Full CRN Detailed Balance Eq.")
	println("------------------------------------------")
	status,costs_eq,q_eq,Jeq = cost_full(kp,km,Vclose,q0,init);

	println("Optimization: ",status)

	println("------------------------------------------")
	println("Full CRN NESS")
	println("------------------------------------------")
	status,costs_full,q_full,Jfull = cost_full(kp,km,Vext,q0,init);
	println("Optimization: ",status)

	println("------------------------------------------")
	println("C1-Pathway")
	println("------------------------------------------")
	status,costs_b,q_b,Jb = cost_c1(kp,km,Vext,q0,init);
	println("Optimization: ",status)

	println("------------------------------------------")
	println("C2-Pathway")
	println("------------------------------------------")
	status,costs_c,q_c,Jc = cost_c2(kp,km,Vext,q0,init);
	println("Optimization: ",status)

	q = zeros(4,6);
	q[1,:] = round.(q_eq,digits=4);
	q[2,:] = round.(q_full,digits=4);
	q[3,:] = round.(q_b,digits=4);
	q[4,:] = round.(q_c,digits=4);

	Jm = zeros(3,6);
	Jm[1,:] = round.(Jfull[:,1]',digits=4);
	Jm[2,:] = round.(Jb[:,1]',digits=4);
	Jm[3,:] = round.(Jc[:,1]',digits=4);

	Jp = zeros(3,6);
	Jp[1,:] = round.(Jfull[:,2]',digits=4);
	Jp[2,:] = round.(Jb[:,2]',digits=4);
	Jp[3,:] = round.(Jc[:,2]',digits=4);

	J = zeros(3,6);
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
		write(f,"r1,r2,r3,r4,r5,r6\n");
		writedlm(f,Jm,',');
	end
	open("data/"*fname*"_Jp.dat","w") do f
		write(f,"r1,r2,r3,r4,r5,r6\n");
		writedlm(f,Jp,',');
	end
	open("data/"*fname*"_J.dat","w") do f
		write(f,"r1,r2,r3,r4,r5,r6\n");
		writedlm(f,J,',');
	end
	open("data/"*fname*"_q.dat","w") do f
		write(f,"f,a,b,c1,c2,d\n");
		writedlm(f,q,',');
	end
	open("data/"*fname*"_costs.dat","w") do f
		write(f,"Sigma,Delta,Chi\n");
		writedlm(f,costs,',');
	end

	return 0
end
##---------------------------------------------------------


##---------------------------------------------------------
##----wrapper functions for results------------------------
function ac_results0()
	init = [20.0,0.0,0.0,0.0,0.0,0.0]
	compAC_costs(init,"compAC0");
	return 0
end

function ac_results1()
	init = [1.0766,1.1591,1.2479,1.2479,1.2479,1.3435]
	compAC_costs(init,"compAC1");
	return 0
end
##---------------------------------------------------------
